library(tidyverse)
library(qs)
library(terra)

##___________________________________________________________________________##


paths <- read_csv("./data/external/file_paths_custom_data.csv")
file_path <- paths$file_path[1]
r <- terra::rast(file_path)
# values(r) <- 0
# plot(r)

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

# Spatial operations require a common CRS. The reference raster CRS is also the
# target CRS used to assign occurrence points to 1 km cells.
if (!terra::same.crs(eu_bounds_vec, r)) {
  eu_bounds_vec <- terra::project(eu_bounds_vec, terra::crs(r))
}

countries <- eu_bounds_vec$NAME_0

if (anyDuplicated(countries)) {
  stop("Country names in the European boundary file must be unique.")
}

safe_ratio <- function(numerator, denominator) {
  if (denominator > 0) numerator / denominator else NA_real_
}

##___________________________________________________________________________##


n_species <- length(sp_data_files)
sp_data_list <- vector("list", n_species)
sp_counts_list <- vector("list", n_species)
sp_counts_by_country_list <- vector("list", n_species)

pb <- txtProgressBar(min = 0, max = n_species, style = 3)

for (i in seq_along(sp_data_files)) {
  sp_data_path <- sp_data_files[[i]]
  sp_data <- qread(sp_data_path)$cleaned_1km

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

  # One country intersection supplies both European membership and national
  # totals. Points on a shared border retain one membership for each country.
  country_intersections <- terra::intersect(
    sp_data_global_vect,
    eu_bounds_vec
  )

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
  n_total_eur <- n_distinct(country_memberships$.record_id)
  n_unique_global <- n_distinct(valid_global_cells$cell)
  n_unique_eur <- n_distinct(european_valid_cells$cell)

  tmp_counts <- data.frame(
    sp_name = sp_name,
    n_total_global = n_total_global,
    n_total_eur = n_total_eur,
    n_unique_global = n_unique_global,
    n_unique_eur = n_unique_eur,
    r_uni_global = safe_ratio(n_unique_global, n_total_global),
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

  setTxtProgressBar(pb, i)
}

close(pb)

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
  file.path(output_dir, "sp_data_all.csv")
)
write_rds(
  sp_data_all,
  file.path(output_dir, "sp_data_all.rds")
)

##___________________________________________________________________________##


write_csv(
  sp_counts,
  file.path(output_dir, "sp_counts.csv")
)
write_rds(
  sp_counts,
  file.path(output_dir, "sp_counts.rds")
)


##___________________________________________________________________________##


write_csv(
  sp_counts_by_country,
  file.path(output_dir, "sp_counts_by_country.csv")
)
write_rds(
  sp_counts_by_country,
  file.path(output_dir, "sp_counts_by_country.rds")
)
