##___________________________________________________________________________##
##
## v3 - now works on deduplicated records from the cleaned wiSDM data
## Deduplication combination taxon keys, lat, lon, and uncertainty fields
##
## PARALLEL: the per-species work is now distributed across worker processes
## (outer loop parallelized with future.apply). Outputs are identical to the
## sequential version - same data structures, same row order - because
## future_lapply preserves input order and each task is deterministic.
##
##___________________________________________________________________________##


library(tidyverse)
library(qs)
library(terra)
library(future)
library(future.apply)

safe_ratio <- function(numerator, denominator) {
  if (denominator > 0) numerator / denominator else NA_real_
}

##___________________________________________________________________________##


paths <- read_csv("./data/external/file_paths_custom_data.csv", 
                  show_col_types = FALSE)
file_path <- paths$file_path[1]
r <- terra::rast(file_path)
# values(r) <- 0
# plot(r)


##___________________________________________________________________________##

## Collect info on the GBIF original data files and datasets!

fl <- list.files("D:/wisdm_v2/gbif_onestop_raw_data",
                 pattern="occurrence.txt",
                 recursive = TRUE,
                 full.names = TRUE)

tkeys <- c()
sp_names <- c()

cli::cli_progress_bar(total = length(fl),
                      name = "Progress files")

for(i in 1:length(fl)){
  
  dt <- read_delim(fl[i], n_max = 100, show_col_types = FALSE) |> 
    suppressWarnings() |> 
    suppressMessages()
  
  tkeys[i] <- dt$acceptedTaxonKey[1]
  sp_names[i] <- dt$species[1]
  
}

# List of original GBIF datafiles and paths to data
# Sorted by species name to allow simple joins (checked and it works)
gbif_raw_paths <- data.frame(tkey    = tkeys,
                             species = sp_names,
                             path    = fl) |> 
  arrange(species)



##___________________________________________________________________________##

sp_data_files <- list.files(
  "./data/projects",
  pattern = "_processed_occurrences.qs",
  recursive = TRUE,
  full.names = TRUE
) |>
  sort()

if (length(sp_data_files) == 0) {
  stop("No processed occurrence files were found in ./data/projects.")
}

# Paths to the shared spatial layers. Kept as variables so the worker
# processes can rebuild live terra objects from disk (external pointers cannot
# be shipped from the main session).
bounds_path <- "./data/external/gadm/europe_selected_countries_wgs84cea_v3-1.gpkg"
bbox_path   <- "./data/external/gadm/europe_selected_countries_bbox_wgs84cea_v2.gpkg"

eu_bounds_vec <- terra::vect(bounds_path)

eu_bbox_vec <- terra::vect(bbox_path)

# Spatial operations require a common CRS. The reference raster CRS is also the
# target CRS used to assign occurrence points to 1 km cells.
if (!terra::same.crs(eu_bounds_vec, r)) {
  eu_bounds_vec <- terra::project(eu_bounds_vec, terra::crs(r))
}

countries <- eu_bounds_vec$NAME_0

if (anyDuplicated(countries)) {
  stop("Country names in the European boundary file must be unique.")
}



##___________________________________________________________________________##

sp_names <- c()
sp_taxon_keys <- c()

for (i in seq_along(sp_data_files)) {
  
  sp_data_path <- sp_data_files[[i]]
  sp_data <- qread(sp_data_path)$cleaned_1km
  sp_tkey <- sp_data$acceptedTaxonKey[1]
  sp_name <- sp_data$species[1]
  
  sp_names[i] <- sp_name
  sp_taxon_keys[i] <- sp_tkey
  # 
  # if(!sp_tkey %in% gbif_raw_paths$tkey){
  #   cat("Key not matched for:", sp_name,"\n\n")
  # }
}

wisdm_sp_data <- data.frame(species_name = sp_names, 
                            acceptedTaxonKey = sp_taxon_keys) |> 
  arrange(species_name)

# Map the keys between sources: GBIF raw data and wiSDM
# Species names enforce link the two
comb_sp_data_tkeys <- cbind(gbif_raw_paths,
                            wisdm_sp_data)

##___________________________________________________________________________##
##
## PER-SPECIES WORKER FUNCTION
##
## This is the exact body of the original for-loop, wrapped in a function so it
## can run in a separate worker process. It returns the same three objects the
## loop used to store per species (sp_data, tmp_counts, tmp_by_country_counts),
## so the assembled outputs downstream are unchanged.
##
##___________________________________________________________________________##

process_species <- function(sp_data_path,
                            comb_sp_data_tkeys,
                            countries,
                            raster_path,
                            bounds_path,
                            bbox_path) {
  
  # Rebuild the shared spatial layers as LIVE terra objects inside this worker.
  # terra objects are external pointers and cannot be sent from the main
  # session. They are cached in the worker process (via options()) so this
  # load happens once per worker, not once per species.
  cache <- getOption("gbif_spatial_cache")
  if (is.null(cache)) {
    r <- terra::rast(raster_path)                 # lazy: header only
    eu_bounds_vec <- terra::vect(bounds_path)
    if (!terra::same.crs(eu_bounds_vec, r)) {
      eu_bounds_vec <- terra::project(eu_bounds_vec, terra::crs(r))
    }
    eu_bbox_vec <- terra::vect(bbox_path)
    cache <- list(r = r, eu_bounds_vec = eu_bounds_vec, eu_bbox_vec = eu_bbox_vec)
    options(gbif_spatial_cache = cache)
  }
  r             <- cache$r
  eu_bounds_vec <- cache$eu_bounds_vec
  eu_bbox_vec   <- cache$eu_bbox_vec
  
  sp_data <- qs::qread(sp_data_path)$cleaned_1km
  
  sp_tkey <- sp_data$acceptedTaxonKey[1]
  
  sp_data <- sp_data |>
    mutate(uid = paste(acceptedTaxonKey,
                       decimalLongitude,
                       decimalLatitude,
                       coordinateUncertaintyInMeters, sep = "|"))
  
  # Clean duplicates base on
  sp_data <- sp_data |>
    group_by(uid) |>
    slice(1) |>
    ungroup()
  
  
  # Get GBIF raw data by getting the original data
  # for all records with lat/lon coords
  # --------------------------------------------
  gbif_raw_dt_path <- comb_sp_data_tkeys |>
    filter(acceptedTaxonKey == sp_tkey) |>
    pull(path)
  
  gbif_dt <- read_delim(gbif_raw_dt_path,
                        show_col_types = FALSE) |>
    suppressWarnings() |>
    suppressMessages() |>
    select(decimalLongitude,
           decimalLatitude,
           scientificName,
           acceptedTaxonKey,
           species,
           coordinateUncertaintyInMeters)
  
  raw_ntotal_records <- nrow(gbif_dt)
  
  gbif_dt <- gbif_dt |>
    filter(!is.na(decimalLongitude),
           !is.na(decimalLatitude))
  
  required_columns <- c("decimalLongitude", "decimalLatitude", "species")
  missing_columns <- setdiff(required_columns, names(sp_data))
  
  if (length(missing_columns) > 0) {
    stop(
      "Missing required columns in ",
      sp_data_path,
      ": ",
      paste(missing_columns, collapse = ", ")
    )
  }
  
  if (nrow(sp_data) == 0) {
    stop("No cleaned 1 km occurrence records found in ", sp_data_path, ".")
  }
  
  sp_name <- sp_data$species[[1]]
  
  # Keep only a stable record identifier as point metadata. It is used to join
  # country memberships to raster cells without relying on row order.
  point_data <- sp_data |>
    transmute(
      .record_id = row_number(),
      decimalLongitude,
      decimalLatitude
    )
  
  sp_data_global_vect <- terra::vect(
    point_data,
    geom = c("decimalLongitude", "decimalLatitude"),
    crs = "EPSG:4326"
  ) |>
    terra::project(terra::crs(r))
  
  
  ##################
  ## GBIF RAW DATA
  ##################
  
  raw_point_data <- gbif_dt |>
    transmute(
      .record_id = row_number(),
      decimalLongitude,
      decimalLatitude
    )
  
  raw_sp_data_global_vect <- terra::vect(
    raw_point_data,
    geom = c("decimalLongitude", "decimalLatitude"),
    crs = "EPSG:4326"
  ) |>
    terra::project(terra::crs(r))
  
  ##______________________________________________________________##
  
  country_intersections <- terra::intersect(
    sp_data_global_vect,
    eu_bounds_vec
  )
  
  bbox_intersections <- terra::intersect(
    sp_data_global_vect,
    eu_bbox_vec
  )
  
  ##______________________________________________________________##
  
  
  # One country intersection supplies both European membership and national
  # totals. Points on a shared border retain one membership for each country.
  raw_country_intersections <- terra::intersect(
    raw_sp_data_global_vect,
    eu_bounds_vec
  )
  
  raw_bbox_intersections <- terra::intersect(
    raw_sp_data_global_vect,
    eu_bbox_vec
  )
  
  ##______________________________________________________________##
  
  
  country_memberships <- country_intersections |>
    as.data.frame() |>
    transmute(
      .record_id,
      country = NAME_0
    )
  
  # Extract cells once for all records. drop_na() preserves the previous
  # definition of a useful record: its reference-raster value must be present.
  extracted_global_cells <- terra::extract(
    r,
    sp_data_global_vect,
    cells = TRUE,
    ID = TRUE
  )
  
  names(extracted_global_cells)[names(extracted_global_cells) == "ID"] <-
    ".record_id"
  
  valid_global_cells <- extracted_global_cells |>
    drop_na() |>
    select(.record_id, cell)
  
  european_valid_cells <- country_memberships |>
    inner_join(valid_global_cells, by = ".record_id")
  
  n_total_global <- nrow(sp_data)
  #
  # n_total_global_dd <- length(unique(paste(sp_data$decimalLatitude,
  #                                          sp_data$decimalLongitude,
  #                                          sp_data$coordinateUncertaintyInMeters)))
  
  n_total_eur <- n_distinct(country_memberships$.record_id)
  n_total_eur_bbox <- nrow(bbox_intersections)
  
  n_unique_global <- n_distinct(valid_global_cells$cell)
  n_unique_eur <- n_distinct(european_valid_cells$cell)
  
  # Temporary data frame holding count data and ratios
  tmp_counts <- data.frame(
    
    # Species name
    sp_name = sp_name,
    
    raw_ntotal_global = raw_ntotal_records,
    raw_ntotal_eur = nrow(raw_country_intersections),
    raw_ntotal_eur_bbox = nrow(raw_bbox_intersections),
    
    # Nr of records in two versions: raw and deduplicated after catching
    # a bug related to having higher counts post filtering (duplications?) than
    # number of checked/retrieved records from GBIF
    n_total_global = n_total_global,
    #n_total_global_dd = n_total_global_dd, # deduplicated records
    
    # Dual count, countries in study area and in major Eur bbox
    n_total_eur = n_total_eur,
    n_total_eur_bbox = n_total_eur_bbox,
    
    # Unique records based on 1x1 km cells
    n_unique_global = n_unique_global,
    n_unique_eur = n_unique_eur,
    
    # Ratios
    r_uni_global = safe_ratio(n_unique_global, n_total_global),
    #r_uni_global_dd = safe_ratio(n_unique_global, n_total_global_dd),
    
    r_uni_eur = safe_ratio(n_unique_eur, n_total_eur),
    
    r_conc_eur = safe_ratio(n_unique_eur, n_unique_global)
  )
  
  country_record_counts <- country_memberships |>
    count(country, name = "n_records")
  
  country_unique_counts <- european_valid_cells |>
    group_by(country) |>
    summarize(
      n_records_unique = n_distinct(cell),
      .groups = "drop"
    )
  
  tmp_by_country_counts <- tibble(country = countries) |>
    left_join(country_record_counts, by = "country") |>
    left_join(country_unique_counts, by = "country") |>
    mutate(
      sp_name = sp_name,
      n_records = replace_na(n_records, 0L),
      n_records_unique = replace_na(n_records_unique, 0L)
    ) |>
    select(sp_name, country, n_records, n_records_unique)
  
  list(
    sp_data               = sp_data,
    tmp_counts            = tmp_counts,
    tmp_by_country_counts = tmp_by_country_counts,
    sp_name               = sp_name
  )
}

##___________________________________________________________________________##
##
## PARALLEL DISPATCH OVER SPECIES
##
##___________________________________________________________________________##

n_species <- length(sp_data_files)

# Number of parallel workers. Each worker holds its own copy of the spatial
# layers AND reads a full GBIF occurrence file into memory, so peak RAM is
# roughly n_workers x (largest occurrence file). Lower this if you hit memory
# limits; raise it toward your core count if you have RAM to spare.
n_workers <- 8L

# multisession = separate R processes (works on Windows; the D:/ paths above
# show this is Windows, where fork-based multicore is unavailable).
plan(multisession, workers = n_workers)

message(sprintf("Processing %d species across %d workers ...",
                n_species, n_workers))

# Optional live progress bar: uncomment the two `progressr` lines and wrap the
# future_lapply call in with_progress(), then call p() inside the worker.
# requires install.packages("progressr")
# progressr::handlers("cli")

results <- future_lapply(
  sp_data_files,
  process_species,
  comb_sp_data_tkeys = comb_sp_data_tkeys,
  countries          = countries,
  raster_path        = file_path,
  bounds_path        = bounds_path,
  bbox_path          = bbox_path,
  future.packages    = c("dplyr", "readr", "tidyr", "tibble", "terra", "qs"),
  future.seed        = TRUE,   # silences the parallel-RNG check (no RNG used)
  future.chunk.size  = 1L      # one species per task -> good load balancing
)

# Release the workers.
plan(sequential)

# Reassemble the same per-species lists the sequential loop produced. Order is
# preserved by future_lapply, so downstream rbind/bind_rows are identical.
sp_data_list              <- lapply(results, `[[`, "sp_data")
sp_counts_list            <- lapply(results, `[[`, "tmp_counts")
sp_counts_by_country_list <- lapply(results, `[[`, "tmp_by_country_counts")

for (res in results) {
  cat("Finished processing species:", res$sp_name, "\n")
}

# Some source files store acceptedTaxonKey with different atomic types. A
# single base rbind call preserves the coercion behavior of the original script
# without repeatedly growing the result inside the loop.
sp_data_all <- do.call(base::rbind, sp_data_list)
sp_counts <- bind_rows(sp_counts_list)
sp_counts_by_country <- bind_rows(sp_counts_by_country_list)

##___________________________________________________________________________##


output_dir <- "./data/post_model_outputs/gbif_records_stats"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

write_csv(
  sp_data_all,
  file.path(output_dir, "sp_data_all_v3.csv")
)
write_rds(
  sp_data_all,
  file.path(output_dir, "sp_data_all_v3.rds")
)

##___________________________________________________________________________##


write_csv(
  sp_counts,
  file.path(output_dir, "sp_counts_v3.csv")
)
write_rds(
  sp_counts,
  file.path(output_dir, "sp_counts_v3.rds")
)


##___________________________________________________________________________##


write_csv(
  sp_counts_by_country,
  file.path(output_dir, "sp_counts_by_country_v3.csv")
)
write_rds(
  sp_counts_by_country,
  file.path(output_dir, "sp_counts_by_country_v3.rds")
)