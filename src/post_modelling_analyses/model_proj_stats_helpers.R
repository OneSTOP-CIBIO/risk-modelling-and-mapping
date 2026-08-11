model_proj_mode_config <- function(calculation_mode = "binary_maps",
                                   output_dir = "./data/processed") {
  allowed_modes <- c("binary_maps", "favourability_maps")

  if (length(calculation_mode) != 1L ||
      is.na(calculation_mode) ||
      !calculation_mode %in% allowed_modes) {
    stop(
      "`calculation_mode` must be one of: ",
      paste(sprintf('"%s"', allowed_modes), collapse = ", "),
      ".",
      call. = FALSE
    )
  }

  path_list_field <- switch(
    calculation_mode,
    binary_maps = "paths_list_bin",
    favourability_maps = "paths_list_fav"
  )
  map_used <- switch(
    calculation_mode,
    binary_maps = "binary10pct",
    favourability_maps = "favourability"
  )

  list(
    calculation_mode = calculation_mode,
    path_list_field = path_list_field,
    map_used = map_used,
    out_file = file.path(
      output_dir,
      sprintf("suitability_by_country_%s.rds", calculation_mode)
    ),
    out_csv_file = file.path(
      output_dir,
      sprintf("suitability_by_country_%s.csv", calculation_mode)
    ),
    region_out_file = file.path(
      output_dir,
      sprintf("suitability_by_region_%s.rds", calculation_mode)
    ),
    region_out_csv_file = file.path(
      output_dir,
      sprintf("suitability_by_region_%s.csv", calculation_mode)
    ),
    err_file = file.path(
      output_dir,
      sprintf("suitability_errors_%s.rds", calculation_mode)
    ),
    perf_file = file.path(
      output_dir,
      sprintf("suitability_perf_%s.rds", calculation_mode)
    )
  )
}

model_proj_make_region_mask <- function(zone_raster) {
  region_mask <- terra::ifel(!is.na(zone_raster), 1L, NA)
  terra::trim(region_mask)
}

model_proj_crop_mask <- function(value_raster, region_mask) {
  cropped <- terra::crop(value_raster, region_mask, snap = "near")
  terra::mask(cropped, region_mask)
}

model_proj_get_path <- function(model_entry, period, scenario, path_list_field) {
  path <- tryCatch(
    model_entry[[path_list_field]][[period]][[scenario]],
    error = function(e) NULL
  )

  if (length(path) == 0L) {
    NA_character_
  } else {
    as.character(path)[[1L]]
  }
}

model_proj_summarise_values <- function(values, zones, n_zone,
                                        calculation_mode) {
  calculation_mode <- model_proj_mode_config(calculation_mode)$calculation_mode
  zones <- as.integer(zones)
  keep <- !is.na(values) &
    !is.na(zones) &
    zones >= 1L &
    zones <= n_zone
  values <- values[keep]
  zones <- zones[keep]

  if (identical(calculation_mode, "binary_maps")) {
    return(list(
      n0 = tabulate(zones[values == 0], nbins = n_zone),
      n1 = tabulate(zones[values == 1], nbins = n_zone),
      total = tabulate(zones, nbins = n_zone)
    ))
  }

  out <- list(
    total = integer(n_zone),
    mean_favourability = rep(NA_real_, n_zone),
    sd_favourability = rep(NA_real_, n_zone),
    median_favourability = rep(NA_real_, n_zone),
    mad_favourability = rep(NA_real_, n_zone)
  )
  if (!length(values)) {
    return(out)
  }

  values_by_zone <- split(values, zones)
  for (zone_name in names(values_by_zone)) {
    zone_id <- as.integer(zone_name)
    zone_values <- values_by_zone[[zone_name]]
    n_values <- length(zone_values)

    out$total[[zone_id]] <- n_values
    out$mean_favourability[[zone_id]] <- mean(zone_values)
    out$sd_favourability[[zone_id]] <- if (n_values > 1L) {
      stats::sd(zone_values)
    } else {
      NA_real_
    }
    out$median_favourability[[zone_id]] <- stats::median(zone_values)
    out$mad_favourability[[zone_id]] <- stats::mad(zone_values)
  }

  out
}

model_proj_zonal_vector <- function(zonal_result, n_zone,
                                    default = NA_real_) {
  out <- rep(default, n_zone)
  if (is.null(zonal_result) || !nrow(zonal_result)) {
    return(out)
  }

  zone_ids <- suppressWarnings(as.integer(zonal_result[[1L]]))
  zone_values <- suppressWarnings(as.numeric(zonal_result[[2L]]))
  keep <- !is.na(zone_ids) &
    zone_ids >= 1L &
    zone_ids <= n_zone
  out[zone_ids[keep]] <- zone_values[keep]
  out
}

model_proj_summarise_zonal <- function(value_raster, zone_raster, n_zone,
                                       calculation_mode) {
  calculation_mode <- model_proj_mode_config(calculation_mode)$calculation_mode

  if (identical(calculation_mode, "binary_maps")) {
    sums <- terra::zonal(
      value_raster,
      zone_raster,
      fun = "sum",
      na.rm = TRUE
    )
    valid_counts <- terra::zonal(
      value_raster,
      zone_raster,
      fun = "notNA"
    )
    total <- as.integer(
      model_proj_zonal_vector(valid_counts, n_zone, default = 0)
    )
    n1 <- as.integer(
      model_proj_zonal_vector(sums, n_zone, default = 0)
    )

    return(list(total = total, n1 = n1, n0 = total - n1))
  }

  valid_counts <- terra::zonal(
    value_raster,
    zone_raster,
    fun = "notNA"
  )
  means <- terra::zonal(
    value_raster,
    zone_raster,
    fun = "mean",
    na.rm = TRUE
  )
  standard_deviations <- terra::zonal(
    value_raster,
    zone_raster,
    fun = stats::sd,
    na.rm = TRUE
  )
  medians <- terra::zonal(
    value_raster,
    zone_raster,
    fun = stats::median,
    na.rm = TRUE
  )
  median_absolute_deviations <- terra::zonal(
    value_raster,
    zone_raster,
    fun = stats::mad,
    na.rm = TRUE
  )

  list(
    total = as.integer(
      model_proj_zonal_vector(valid_counts, n_zone, default = 0)
    ),
    mean_favourability = model_proj_zonal_vector(means, n_zone),
    sd_favourability = model_proj_zonal_vector(
      standard_deviations,
      n_zone
    ),
    median_favourability = model_proj_zonal_vector(medians, n_zone),
    mad_favourability = model_proj_zonal_vector(
      median_absolute_deviations,
      n_zone
    )
  )
}
