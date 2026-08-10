
##___________________________________________________________________________##
##
## v3 - now works on deduplicated records from the cleaned wiSDM data
## Deduplication combination taxon keys, lat, lon, and uncertainty fields
##
##___________________________________________________________________________##


library(tidyverse)
library(qs)
library(terra)

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

eu_bounds_vec <- terra::vect(
  "./data/external/gadm/europe_selected_countries_wgs84cea_v3-1.gpkg"
)

eu_bbox_vec <- terra::vect(
  "./data/external/gadm/europe_selected_countries_bbox_wgs84cea_v2.gpkg"
)

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


n_species <- length(sp_data_files)
sp_data_list <- vector("list", n_species)
sp_counts_list <- vector("list", n_species)
sp_counts_by_country_list <- vector("list", n_species)

#pb <- txtProgressBar(min = 0, max = n_species, style = 3)

for (i in seq_along(sp_data_files)) {
  
  sp_data_path <- sp_data_files[[i]]
  sp_data <- qread(sp_data_path)$cleaned_1km

  sp_tkey <- sp_data$acceptedTaxonKey[1]
  
  sp_data <- sp_data |> 
    mutate(uid = paste(acceptedTaxonKey,
                       decimalLongitude,
                       decimalLatitude,
                       coordinateUncertaintyInMeters,sep="|"))
  
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

  sp_data_list[[i]] <- sp_data
  sp_counts_list[[i]] <- tmp_counts
  sp_counts_by_country_list[[i]] <- tmp_by_country_counts

  #setTxtProgressBar(pb, i)
  cat("[",i,"] Finished processing species:", sp_name,"\n\n")
}

#close(pb)

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
