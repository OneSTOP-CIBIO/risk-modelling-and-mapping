##___________________________________________________________________________##
##
## v4 - deduplicated records from cleaned wiSDM data + PARALLEL over species
##      + verbose per-species timing to console (live) and a txt log file.
##
## Logging note: with multisession, workers are separate R processes. Printing
## from a worker does NOT reliably reach the console, and several workers
## appending to one file can interleave. So each worker MEASURES its own time
## and RETURNS it; all logging is done from the MAIN process (progressr's
## handlers also run in the main process, so the live bar is concurrency-safe).
##
## Outputs are identical to the sequential/v3 version - same data structures,
## same row order - because future_lapply preserves input order.
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

# Format a number of seconds as H:MM:SS for the log summary.
fmt_hms <- function(sec) {
  sec <- round(as.numeric(sec))
  sprintf("%d:%02d:%02d", sec %/% 3600, (sec %% 3600) %/% 60, sec %% 60)
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
## Exact body of the original loop, wrapped in a function so it can run in a
## worker process. Also measures its own wall-clock time and (optionally)
## signals a progressr progressor `p`. Returns the same three data objects the
## loop stored per species, plus timing/diagnostics used only for logging.
##
##___________________________________________________________________________##

process_species <- function(sp_data_path,
                            comb_sp_data_tkeys,
                            countries,
                            raster_path,
                            bounds_path,
                            bbox_path,
                            p = NULL) {
  
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
  
  # Time the per-species work only (spatial-layer load above is one-off per
  # worker and would otherwise skew the first species handled by each worker).
  work_start <- Sys.time()
  
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
  
  
  #####################################
  ## Save data
  
  
  global_filt_file_path <- paste0("D:/wisdm_v2/gbif_onestop_raw_geo_data/",
                                  sp_name,"_filt_wisdm.gpkg")
  
  writeVector(sp_data_global_vect, global_filt_file_path,overwrite=TRUE)
  
  raw_file_path <- paste0("D:/wisdm_v2/gbif_onestop_raw_geo_data/",
                          sp_name,"_raw_gbif.gpkg")
  
  writeVector(raw_sp_data_global_vect, raw_file_path,overwrite=TRUE)
  
  
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
  
  # --- timing + live progress signal ----------------------------------------
  work_end <- Sys.time()
  elapsed_sec <- as.numeric(difftime(work_end, work_start, units = "secs"))
  
  # Signal the (main-process) progress bar, carrying a per-species message.
  if (!is.null(p)) {
    p(message = sprintf("%s (%.1fs)", sp_name, elapsed_sec))
  }
  
  list(
    sp_data               = sp_data,
    tmp_counts            = tmp_counts,
    tmp_by_country_counts = tmp_by_country_counts,
    # --- diagnostics used only for logging ---
    sp_name               = sp_name,
    elapsed_sec           = elapsed_sec,
    n_1km                 = nrow(sp_data),
    raw_ntotal            = raw_ntotal_records,
    pid                   = Sys.getpid(),
    finished_at           = work_end
  )
}

##___________________________________________________________________________##
##
## PARALLEL DISPATCH OVER SPECIES  (with live progress + timing log)
##
##___________________________________________________________________________##

n_species <- length(sp_data_files)

# Number of parallel workers. Each worker holds its own copy of the spatial
# layers AND reads a full GBIF occurrence file into memory, so peak RAM is
# roughly n_workers x (largest occurrence file). Lower this if you hit memory
# limits; raise it toward your core count if you have RAM to spare.
n_workers <- 10L

# Logging goes to a timestamped txt file, written by the MAIN process only.
log_dir  <- "./data/post_model_outputs/gbif_records_stats/logs"
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)
log_path <- file.path(log_dir, format(Sys.time(), "processing_log_%Y%m%d_%H%M%S.txt"))

# multisession = separate R processes (works on Windows; the D:/ paths above
# show this is Windows, where fork-based multicore is unavailable).
plan(multisession, workers = n_workers)

# Live progress bar is provided by progressr if available; otherwise the run
# still works and the log file is still written at the end.
use_progress <- requireNamespace("progressr", quietly = TRUE)

message(sprintf("Processing %d species across %d workers ...",
                n_species, n_workers))
message("Timing log: ", log_path)
if (!use_progress) {
  message("Tip: install.packages('progressr') for a live progress bar.")
}

# Common dispatch call, reused by both the progress and no-progress branches.
dispatch <- function(p = NULL) {
  future_lapply(
    sp_data_files,
    process_species,
    p                  = p,
    comb_sp_data_tkeys = comb_sp_data_tkeys,
    countries          = countries,
    raster_path        = file_path,
    bounds_path        = bounds_path,
    bbox_path          = bbox_path,
    future.packages    = c("dplyr", "readr", "tidyr", "tibble", "terra", "qs"),
    future.seed        = TRUE,   # silences the parallel-RNG check (no RNG used)
    future.chunk.size  = 1L      # one species per task -> good load balancing
  )
}

wall_start <- Sys.time()

if (use_progress) {
  progressr::handlers("cli")            # live one-line bar with per-species msg
  results <- progressr::with_progress({
    p <- progressr::progressor(steps = n_species)
    dispatch(p)
  })
} else {
  results <- dispatch(NULL)
}

wall_end <- Sys.time()

# Release the workers.
plan(sequential)

##___________________________________________________________________________##
##
## WRITE THE TIMING LOG + CONSOLE SUMMARY (main process, ordered)
##
##___________________________________________________________________________##

sp_names_done <- vapply(results, `[[`, character(1), "sp_name")
elapsed       <- vapply(results, `[[`, numeric(1),   "elapsed_sec")
n1km          <- vapply(results, `[[`, integer(1),   "n_1km")
rawn          <- vapply(results, `[[`, integer(1),   "raw_ntotal")
pids          <- vapply(results, `[[`, integer(1),   "pid")
finished_at   <- lapply(results, `[[`, "finished_at")

wall_elapsed <- as.numeric(difftime(wall_end, wall_start, units = "secs"))
sum_elapsed  <- sum(elapsed)
slow_i       <- which.max(elapsed)
fast_i       <- which.min(elapsed)

log_lines <- c(
  "==============================================================================",
  sprintf(" GBIF record-count run - %s", format(wall_start, "%Y-%m-%d %H:%M:%S")),
  sprintf(" Species: %d    Workers: %d    Plan: multisession", n_species, n_workers),
  "==============================================================================",
  sprintf(" %-4s %-42s %9s %9s %9s %8s %-9s",
          "idx", "species", "raw_n", "n_1km", "secs", "pid", "finished"),
  "------------------------------------------------------------------------------"
)

for (i in seq_along(results)) {
  log_lines <- c(log_lines, sprintf(
    " %-4d %-42s %9d %9d %9.1f %8d %-9s",
    i, substr(sp_names_done[i], 1, 42), rawn[i], n1km[i], elapsed[i],
    pids[i], format(finished_at[[i]], "%H:%M:%S")
  ))
}

log_lines <- c(
  log_lines,
  "------------------------------------------------------------------------------",
  sprintf(" Wall-clock time (parallel) : %s  (%.1f s)", fmt_hms(wall_elapsed), wall_elapsed),
  sprintf(" Sum of per-species times   : %s  (%.1f s)", fmt_hms(sum_elapsed),  sum_elapsed),
  sprintf(" Approx. achieved speedup   : %.2fx  (on %d workers)",
          if (wall_elapsed > 0) sum_elapsed / wall_elapsed else NA_real_, n_workers),
  sprintf(" Mean / median per species  : %.1f s / %.1f s", mean(elapsed), median(elapsed)),
  sprintf(" Slowest : %-42s %.1f s", substr(sp_names_done[slow_i], 1, 42), elapsed[slow_i]),
  sprintf(" Fastest : %-42s %.1f s", substr(sp_names_done[fast_i], 1, 42), elapsed[fast_i]),
  "=============================================================================="
)

writeLines(log_lines, log_path)

# Echo the summary (not the full per-species table) to the console.
message("")
message(paste(tail(log_lines, 9), collapse = "\n"))
message("Full per-species timing written to: ", log_path)

##___________________________________________________________________________##
##
## ASSEMBLE OUTPUTS (identical structures/order to the sequential version)
##
##___________________________________________________________________________##

sp_data_list              <- lapply(results, `[[`, "sp_data")
sp_counts_list            <- lapply(results, `[[`, "tmp_counts")
sp_counts_by_country_list <- lapply(results, `[[`, "tmp_by_country_counts")

# Some source files store acceptedTaxonKey with different atomic types. A
# single base rbind call preserves the coercion behavior of the original script
# without repeatedly growing the result inside the loop.
sp_data_all <- do.call(base::rbind, sp_data_list)
sp_counts <- bind_rows(sp_counts_list)
sp_counts_by_country <- bind_rows(sp_counts_by_country_list)

##___________________________________________________________________________##

version_nr <- 4

output_dir <- "./data/post_model_outputs/gbif_records_stats"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

write_csv(
  sp_data_all,
  file.path(output_dir, paste0("sp_data_all_v",version_nr,".csv"))
)
write_rds(
  sp_data_all,
  file.path(output_dir, paste0("sp_data_all_v",version_nr,".rds"))
)

##___________________________________________________________________________##


write_csv(
  sp_counts,
  file.path(output_dir, paste0("sp_counts_v",version_nr,".csv"))
)
write_rds(
  sp_counts,
  file.path(output_dir, paste0("sp_counts_v",version_nr,".rds"))
)


##___________________________________________________________________________##


write_csv(
  sp_counts_by_country,
  file.path(output_dir, paste0("sp_counts_by_country_v",version_nr,".csv"))
)
write_rds(
  sp_counts_by_country,
  file.path(output_dir, paste0("sp_counts_by_country_v",version_nr,".rds"))
)
