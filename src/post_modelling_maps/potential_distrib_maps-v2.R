# ==============================================================================
# wiSDM 2.0 potential-distribution maps from cached binary model projections
# ==============================================================================
#
# From RStudio, open the repository project and run the complete batch with:
#
#   source("src/post_modelling_maps/potential_distrib_maps-v2.R")
#
# Interactive source() calls autorun by default. To load only the functions:
#
#   options(wisdm.run_potential_maps_on_source = FALSE)
#   source("src/post_modelling_maps/potential_distrib_maps-v2.R")
#
# The batch can then be started explicitly with:
#
#   run_potential_distribution_maps()
#
# Direct Rscript execution also remains supported, but RStudio source() is the
# primary entry point.
#
# The v1 renderer is intentionally left unchanged.

# ==============================================================================
# Configuration
# ==============================================================================

potential_maps_config <- list(
  cache_path = file.path(
    "data",
    "processed",
    "best_mod_set_paths.rds"
  ),
  bounds_gpkg = file.path(
    "data",
    "external",
    "gadm",
    "europe_selected_countries_wgs84cea_v3-1.gpkg"
  ),
  plot_bounds_gpkg = file.path(
    "data",
    "external",
    "gadm",
    "europe_selected_countries_wgs84cea_v3-1_simplified_1km.gpkg"
  ),
  plot_bounds_tolerance_m = 1000,
  basemap_gpkg = file.path(
    "data",
    "external",
    "gadm",
    "eu_council_b100_bbox_gadm41_full_by_country_wgs84cea.gpkg"
  ),
  basemap_cache_gpkg = file.path(
    "data",
    "external",
    "gadm",
    paste0(
      "eu_council_b100_bbox_gadm41_full_by_country_wgs84cea",
      "_simplified_2km.gpkg"
    )
  ),
  basemap_tolerance_m = 2000,
  grid_path = file.path(
    "data",
    "external",
    "grids",
    "Grid_20km_Eur_OneSTOP_EPSG6933_v1.shp"
  ),
  study_grid_gpkg = file.path(
    "data",
    "external",
    "grids",
    "Grid_20km_Eur_OneSTOP_EPSG6933_v1_study_area.gpkg"
  ),
  grid_id_field = "id",
  occurrence_dir = "D:/wisdm_v2/gbif_onestop_raw_geo_data",
  output_dir = file.path(
    "data",
    "post_model_maps",
    "sp_maps_by_period_scenario"
  ),
  occurrence_suffix = "_filt_wisdm.gpkg",
  map_dpi = 300,
  map_width = 11,
  map_height = 11 * 0.536,
  display_aggregate_factor = 3L,
  display_aggregate_fun = "modal",
  disclaimer_text = paste0(
    "OneSTOP SDM Gallery | ⚠️Disclaimer: maps represent potential ",
    "distribution"
  ),
  basemap_disclaimer = "Basemap: GADM-v4.1 (EPSG:6933)",
  version = "0.2",
  ref_date = "05/08/2026",
  legend_position_inside = c(0.015, 0.60),
  color_suitable = "#2E7D32",
  color_unsuitable = "#EEEEEE",
  color_outside_area = "#FFFFFF",
  color_border = "#757575",
  color_basemap_border = "#B0B0B0",
  color_grid = "#3C64AD",
  color_grid_outside = "#BDBDBD",
  color_legend_border = "#777777"
)

potential_maps_condition_levels <- c(
  "Favourable conditions",
  "Unfavourable conditions",
  "Outside modelled area"
)

potential_maps_occurrence_levels <- c(
  "Within project area",
  "Outside project area"
)

# Keep the requested period/scenario names and ordering. Projection codes retain
# their leading and trailing underscores because the legacy filename builder
# adds another separator before the code.
per_scn <- data.frame(
  period = c(
    rep(c("mid_2041_2070", "late_2071_2100"), each = 3L),
    "hist"
  ),
  scenario = c(
    rep(c("ssp126", "ssp370", "ssp585"), times = 2L),
    "hist"
  ),
  stringsAsFactors = FALSE
)
per_scn$ps_name <- paste(per_scn$period, per_scn$scenario, sep = "_")
per_scn$proj_code <- c(
  "_2041_2070_ssp126_bin_",
  "_2041_2070_ssp370_bin_",
  "_2041_2070_ssp585_bin_",
  "_2071_2100_ssp126_bin_",
  "_2071_2100_ssp370_bin_",
  "_2071_2100_ssp585_bin_",
  "_hist_bin_"
)
per_scn$proj_name <- c(
  "SSP1/2.6 (2070)",
  "SSP3/7.0 (2070)",
  "SSP5/8.5 (2070)",
  "SSP1/2.6 (2100)",
  "SSP3/7.0 (2100)",
  "SSP5/8.5 (2100)",
  "Historical/baseline"
)

potential_maps_summary_columns <- c(
  "sp_name",
  "proj_code",
  "total_cells",
  "na_cells",
  "valid_cells",
  "suitable_cells",
  "unsuitable_cells"
)

# ==============================================================================
# Validation and naming helpers
# ==============================================================================

potential_maps_require_packages <- function() {
  required <- c(
    "terra",
    "sf",
    "ggplot2",
    "ggspatial",
    "cli",
    "ggtext"
  )
  missing <- required[!vapply(
    required,
    requireNamespace,
    logical(1),
    quietly = TRUE
  )]

  if (length(missing) > 0L) {
    stop(
      "Missing required R package(s): ",
      paste(missing, collapse = ", "),
      ".",
      call. = FALSE
    )
  }

  invisible(required)
}

potential_maps_progress_bar <- function(name, total,
                                        .envir = parent.frame()) {
  cli::cli_progress_bar(
    name = name,
    total = total,
    format = paste(
      "{cli::pb_spin}",
      "{cli::pb_name}",
      "{cli::pb_current}/{cli::pb_total}",
      "{cli::pb_bar}",
      "{cli::pb_percent}",
      "| {cli::pb_status}",
      "| ETA {cli::pb_eta}"
    ),
    current = FALSE,
    auto_terminate = FALSE,
    clear = FALSE,
    .envir = .envir
  )
}

potential_maps_progress_update <- function(id, set, status) {
  cli::cli_progress_update(
    id = id,
    set = set,
    inc = 0,
    status = status,
    force = TRUE
  )
  invisible(id)
}

potential_maps_assert_scalar_logical <- function(value, argument) {
  if (!is.logical(value) || length(value) != 1L || is.na(value)) {
    stop("`", argument, "` must be TRUE or FALSE.", call. = FALSE)
  }
  invisible(value)
}

potential_maps_assert_positive_number <- function(value, argument,
                                                  integer = FALSE) {
  if (!is.numeric(value) ||
      length(value) != 1L ||
      is.na(value) ||
      !is.finite(value) ||
      value <= 0 ||
      (integer && value != as.integer(value))) {
    stop(
      "`",
      argument,
      "` must be one positive ",
      if (integer) "integer" else "number",
      ".",
      call. = FALSE
    )
  }
  invisible(value)
}

potential_maps_assert_path <- function(path, argument, must_exist = FALSE) {
  if (!is.character(path) ||
      length(path) != 1L ||
      is.na(path) ||
      !nzchar(trimws(path))) {
    stop("`", argument, "` must be one non-empty path.", call. = FALSE)
  }
  if (must_exist && !file.exists(path)) {
    stop("Required path does not exist: ", path, call. = FALSE)
  }
  invisible(path)
}

potential_maps_species_tokens <- function(sp_name) {
  if (!is.character(sp_name) ||
      length(sp_name) != 1L ||
      is.na(sp_name) ||
      !nzchar(trimws(sp_name))) {
    stop("`sp_name` must be one non-empty species name.", call. = FALSE)
  }

  tokens <- strsplit(trimws(sp_name), "[[:space:]]+")[[1L]]
  tokens <- tokens[nzchar(tokens)]
  if (length(tokens) < 2L) {
    stop(
      "Cannot derive a binomial name from species name: ",
      sp_name,
      call. = FALSE
    )
  }

  epithet_index <- 2L
  hybrid_markers <- c("x", "X", "×")
  while (epithet_index <= length(tokens) &&
         tokens[[epithet_index]] %in% hybrid_markers) {
    epithet_index <- epithet_index + 1L
  }
  if (epithet_index > length(tokens)) {
    stop(
      "Cannot derive a species epithet from hybrid name: ",
      sp_name,
      call. = FALSE
    )
  }

  c(tokens[[1L]], tokens[[epithet_index]])
}

potential_maps_species_display_name <- function(sp_name) {
  paste(potential_maps_species_tokens(sp_name), collapse = " ")
}

potential_maps_species_slug <- function(sp_name) {
  slug <- paste(potential_maps_species_tokens(sp_name), collapse = "_")
  slug_ascii <- iconv(slug, from = "UTF-8", to = "ASCII//TRANSLIT")
  if (is.na(slug_ascii)) {
    slug_ascii <- slug
  }
  slug_ascii <- tolower(slug_ascii)
  slug_ascii <- gsub("[^a-z0-9]+", "_", slug_ascii)
  slug_ascii <- gsub("_+", "_", slug_ascii)
  gsub("^_|_$", "", slug_ascii)
}

potential_maps_output_filename <- function(sp_name, proj_code) {
  if (!is.character(proj_code) ||
      length(proj_code) != 1L ||
      is.na(proj_code) ||
      !proj_code %in% per_scn$proj_code) {
    stop("Unknown legacy projection code: ", proj_code, call. = FALSE)
  }

  paste0(
    "sdm_map_",
    potential_maps_species_slug(sp_name),
    "_",
    proj_code,
    ".png"
  )
}

potential_maps_occurrence_path <- function(
    sp_name,
    occurrence_dir = potential_maps_config$occurrence_dir,
    suffix = potential_maps_config$occurrence_suffix) {
  file.path(occurrence_dir, paste0(sp_name, suffix))
}

potential_maps_validate_catalog <- function(catalog) {
  if (!is.list(catalog) || length(catalog) == 0L) {
    stop("The best-model cache must be a non-empty list.", call. = FALSE)
  }

  required_fields <- c("sp_name", "tkey", "paths_list_bin")
  normalized_names <- character(length(catalog))

  for (i in seq_along(catalog)) {
    entry <- catalog[[i]]
    if (!is.list(entry)) {
      stop("Catalog entry ", i, " is not a list.", call. = FALSE)
    }

    missing_fields <- setdiff(required_fields, names(entry))
    if (length(missing_fields) > 0L) {
      stop(
        "Catalog entry ",
        i,
        " is missing required field(s): ",
        paste(missing_fields, collapse = ", "),
        ".",
        call. = FALSE
      )
    }

    sp_name <- entry$sp_name
    if (!is.character(sp_name) ||
        length(sp_name) != 1L ||
        is.na(sp_name) ||
        !nzchar(trimws(sp_name))) {
      stop(
        "Catalog entry ",
        i,
        " has an invalid `sp_name`.",
        call. = FALSE
      )
    }
    if (length(entry$tkey) != 1L || is.na(entry$tkey)) {
      stop(
        "Catalog entry for ",
        sp_name,
        " has an invalid `tkey`.",
        call. = FALSE
      )
    }
    if (!is.list(entry$paths_list_bin)) {
      stop(
        "Catalog entry for ",
        sp_name,
        " has no list-valued `paths_list_bin`.",
        call. = FALSE
      )
    }

    for (j in seq_len(nrow(per_scn))) {
      period <- per_scn$period[[j]]
      scenario <- per_scn$scenario[[j]]
      period_paths <- entry$paths_list_bin[[period]]
      if (!is.list(period_paths)) {
        stop(
          "Catalog entry for ",
          sp_name,
          " has no list-valued binary paths for period ",
          period,
          ".",
          call. = FALSE
        )
      }
      path <- period_paths[[scenario]]
      if (!is.character(path) ||
          length(path) != 1L ||
          is.na(path) ||
          !nzchar(trimws(path))) {
        stop(
          "Catalog entry for ",
          sp_name,
          " has no valid binary raster path for ",
          period,
          "/",
          scenario,
          ".",
          call. = FALSE
        )
      }
    }

    normalized_names[[i]] <- tolower(trimws(sp_name))
  }

  if (anyDuplicated(normalized_names)) {
    stop(
      "The best-model cache contains duplicated species names.",
      call. = FALSE
    )
  }

  invisible(TRUE)
}

potential_maps_model_type_disclaimer <- function(entry) {
  path_field <- if (!is.null(entry$paths_bin)) {
    "paths_bin"
  } else if (!is.null(entry$paths_list_bin)) {
    "paths_list_bin"
  } else {
    stop(
      "Catalog entry has neither `paths_bin` nor `paths_list_bin`.",
      call. = FALSE
    )
  }
  paths <- unlist(entry[[path_field]], use.names = FALSE)
  if (length(paths) == 0L ||
      !is.character(paths) ||
      anyNA(paths) ||
      any(!nzchar(trimws(paths)))) {
    stop(
      "Catalog entry has invalid paths in `",
      path_field,
      "`.",
      call. = FALSE
    )
  }

  normalized_paths <- tolower(gsub("\\\\", "/", paths))
  is_climate <- grepl("/climate/", normalized_paths, fixed = TRUE)
  if (all(is_climate)) {
    "Climate-only model"
  } else {
    "Climate-habitat ensemble"
  }
}

build_potential_map_jobs <- function(
    catalog,
    output_dir = potential_maps_config$output_dir,
    occurrence_dir = potential_maps_config$occurrence_dir) {
  potential_maps_validate_catalog(catalog)
  potential_maps_assert_path(output_dir, "output_dir")
  potential_maps_assert_path(occurrence_dir, "occurrence_dir")

  rows <- vector("list", length(catalog) * nrow(per_scn))
  row_index <- 0L

  for (species_index in seq_along(catalog)) {
    entry <- catalog[[species_index]]
    display_name <- potential_maps_species_display_name(entry$sp_name)
    species_slug <- potential_maps_species_slug(entry$sp_name)
    mod_type_disclaimer <- potential_maps_model_type_disclaimer(entry)
    point_path <- potential_maps_occurrence_path(
      entry$sp_name,
      occurrence_dir = occurrence_dir
    )

    for (scenario_index in seq_len(nrow(per_scn))) {
      scenario_row <- per_scn[scenario_index, , drop = FALSE]
      period <- scenario_row$period[[1L]]
      scenario <- scenario_row$scenario[[1L]]
      raster_path <- entry$paths_list_bin[[period]][[scenario]]
      output_filename <- potential_maps_output_filename(
        entry$sp_name,
        scenario_row$proj_code[[1L]]
      )

      row_index <- row_index + 1L
      rows[[row_index]] <- data.frame(
        job_id = paste(entry$tkey, period, scenario, sep = "__"),
        species_index = species_index,
        sp_name = entry$sp_name,
        display_name = display_name,
        species_slug = species_slug,
        tkey = entry$tkey,
        period = period,
        scenario = scenario,
        ps_name = scenario_row$ps_name[[1L]],
        proj_code = scenario_row$proj_code[[1L]],
        proj_name = scenario_row$proj_name[[1L]],
        mod_type_disclaimer = mod_type_disclaimer,
        raster_path = raster_path,
        point_path = point_path,
        output_filename = output_filename,
        output_path = file.path(output_dir, output_filename),
        stringsAsFactors = FALSE,
        check.names = FALSE
      )
    }
  }

  jobs <- do.call(rbind, rows)
  rownames(jobs) <- NULL

  if (anyDuplicated(jobs$job_id)) {
    stop("Potential-map job IDs are not unique.", call. = FALSE)
  }

  normalized_outputs <- tolower(normalizePath(
    jobs$output_path,
    winslash = "/",
    mustWork = FALSE
  ))
  if (anyDuplicated(normalized_outputs)) {
    duplicated_outputs <- unique(
      jobs$output_filename[duplicated(normalized_outputs)]
    )
    stop(
      "Legacy binomial naming produces duplicate map output(s): ",
      paste(duplicated_outputs, collapse = ", "),
      ".",
      call. = FALSE
    )
  }

  jobs
}

# ==============================================================================
# Spatial preparation
# ==============================================================================

potential_maps_validate_binary_metadata <- function(raster, context) {
  if (!inherits(raster, "SpatRaster")) {
    stop(context, " is not a terra SpatRaster.", call. = FALSE)
  }
  if (terra::nlyr(raster) != 1L) {
    stop(context, " must contain exactly one raster layer.", call. = FALSE)
  }

  value_range <- terra::minmax(raster)
  finite_range <- value_range[is.finite(value_range)]
  if (length(finite_range) > 0L &&
      (min(finite_range) < 0 || max(finite_range) > 1)) {
    stop(
      context,
      " has values outside the binary 0/1 range.",
      call. = FALSE
    )
  }

  invisible(raster)
}

potential_maps_build_mask <- function(bounds_vec, template) {
  if (!inherits(bounds_vec, "SpatVector")) {
    stop("`bounds_vec` must be a terra SpatVector.", call. = FALSE)
  }
  if (!inherits(template, "SpatRaster")) {
    stop("`template` must be a terra SpatRaster.", call. = FALSE)
  }

  if (!terra::same.crs(bounds_vec, template)) {
    bounds_vec <- terra::project(bounds_vec, terra::crs(template))
  }

  rst_mask <- terra::rasterize(bounds_vec, template, field = 1)
  if (!isTRUE(terra::compareGeom(
    rst_mask,
    template,
    lyrs = FALSE,
    stopOnError = FALSE,
    messages = FALSE
  ))) {
    stop(
      "The European raster mask does not match the model template geometry.",
      call. = FALSE
    )
  }

  rst_mask
}

potential_maps_load_plot_boundaries <- function(
    bounds_vec,
    source_path,
    cache_path,
    tolerance_m,
    target_crs,
    layer_label = "GADM plotting boundaries") {
  if (!inherits(bounds_vec, "SpatVector")) {
    stop("`bounds_vec` must be a terra SpatVector.", call. = FALSE)
  }
  potential_maps_assert_path(source_path, "source_path", must_exist = TRUE)
  potential_maps_assert_path(cache_path, "cache_path")
  potential_maps_assert_positive_number(tolerance_m, "tolerance_m")
  if (!is.character(layer_label) ||
      length(layer_label) != 1L ||
      is.na(layer_label) ||
      !nzchar(layer_label)) {
    stop("`layer_label` must be one non-empty string.", call. = FALSE)
  }

  source_normalized <- normalizePath(
    source_path,
    winslash = "/",
    mustWork = TRUE
  )
  cache_normalized <- normalizePath(
    cache_path,
    winslash = "/",
    mustWork = FALSE
  )
  if (tolower(source_normalized) == tolower(cache_normalized)) {
    stop(
      "The simplified plotting cache must not overwrite the source GADM file.",
      call. = FALSE
    )
  }

  cache_is_current <- file.exists(cache_path) &&
    file.info(cache_path)$mtime >= file.info(source_path)$mtime

  if (!cache_is_current) {
    cli::cli_alert_info(
      paste0(
        "Creating cached ",
        layer_label,
        " with a ",
        format(tolerance_m, trim = TRUE),
        " m simplification tolerance."
      )
    )
    simplified_vec <- terra::simplifyGeom(
      bounds_vec,
      tolerance = tolerance_m
    )
    cache_dir <- dirname(cache_path)
    dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
    if (!dir.exists(cache_dir)) {
      stop(
        "Could not create plotting-boundary cache directory: ",
        cache_dir,
        call. = FALSE
      )
    }
    terra::writeVector(
      simplified_vec,
      cache_path,
      overwrite = TRUE,
      filetype = "GPKG"
    )
    cli::cli_alert_success(
      "Cached simplified {layer_label} at {.path {cache_path}}."
    )
  } else {
    cli::cli_alert_info(
      "Using cached simplified {layer_label}: {.path {cache_path}}."
    )
  }

  plot_bounds_vec <- terra::vect(cache_path)
  if (terra::nrow(plot_bounds_vec) != terra::nrow(bounds_vec)) {
    stop(
      "The simplified ",
      layer_label,
      " cache has ",
      terra::nrow(plot_bounds_vec),
      " features but the source has ",
      terra::nrow(bounds_vec),
      ". Delete the cache and rerun to rebuild it.",
      call. = FALSE
    )
  }

  plot_boundaries <- sf::st_as_sf(plot_bounds_vec)
  valid <- sf::st_is_valid(plot_boundaries)
  if (anyNA(valid) || !all(valid)) {
    stop(
      "The simplified ",
      layer_label,
      " cache contains invalid geometry. ",
      "Delete the cache and rerun to rebuild it.",
      call. = FALSE
    )
  }
  if (is.na(sf::st_crs(plot_boundaries))) {
    stop("The simplified ", layer_label, " cache has no CRS.", call. = FALSE)
  }
  if (sf::st_crs(plot_boundaries) != target_crs) {
    plot_boundaries <- sf::st_transform(plot_boundaries, target_crs)
  }

  plot_boundaries
}

potential_maps_validate_basemap_countries <- function(
    study_bounds_vec,
    basemap_vec) {
  if (!inherits(study_bounds_vec, "SpatVector") ||
      !inherits(basemap_vec, "SpatVector")) {
    stop("Study boundaries and basemap must be terra SpatVectors.", call. = FALSE)
  }
  boundary_sets <- list(
    "study-area boundaries" = study_bounds_vec,
    "GADM 4.1 basemap" = basemap_vec
  )
  for (label in names(boundary_sets)) {
    item <- boundary_sets[[label]]
    if (!"NAME_0" %in% names(item)) {
      stop(
        label,
        " must contain a `NAME_0` country field.",
        call. = FALSE
      )
    }
  }

  study_names <- unique(trimws(as.character(study_bounds_vec$NAME_0)))
  basemap_names <- unique(trimws(as.character(basemap_vec$NAME_0)))
  if (anyNA(study_names) || any(!nzchar(study_names)) ||
      anyNA(basemap_names) || any(!nzchar(basemap_names))) {
    stop("Country names in `NAME_0` must be complete.", call. = FALSE)
  }
  missing_names <- setdiff(study_names, basemap_names)
  if (length(missing_names) > 0L) {
    stop(
      "The GADM 4.1 basemap is missing study-area countries: ",
      paste(missing_names, collapse = ", "),
      ".",
      call. = FALSE
    )
  }

  invisible(TRUE)
}

potential_maps_grid_source_mtime <- function(grid_path) {
  grid_dir <- dirname(grid_path)
  grid_stem <- tools::file_path_sans_ext(basename(grid_path))
  candidates <- list.files(grid_dir, full.names = TRUE)
  components <- candidates[
    tools::file_path_sans_ext(basename(candidates)) == grid_stem
  ]
  if (length(components) == 0L) {
    stop("No source-grid files found for: ", grid_path, call. = FALSE)
  }
  max(file.info(components)$mtime, na.rm = TRUE)
}

potential_maps_load_study_grid <- function(
    grid_ref,
    bounds_vec,
    grid_path,
    bounds_path,
    cache_path,
    grid_id_field = potential_maps_config$grid_id_field) {
  if (!inherits(grid_ref, "sf")) {
    stop("`grid_ref` must be an sf object.", call. = FALSE)
  }
  if (!inherits(bounds_vec, "SpatVector")) {
    stop("`bounds_vec` must be a terra SpatVector.", call. = FALSE)
  }
  potential_maps_assert_path(grid_path, "grid_path", must_exist = TRUE)
  potential_maps_assert_path(bounds_path, "bounds_path", must_exist = TRUE)
  potential_maps_assert_path(cache_path, "study_grid_gpkg")
  if (!is.character(grid_id_field) ||
      length(grid_id_field) != 1L ||
      is.na(grid_id_field) ||
      !nzchar(grid_id_field)) {
    stop("`grid_id_field` must be one non-empty field name.", call. = FALSE)
  }
  if (!grid_id_field %in% names(grid_ref)) {
    stop(
      "The 20 km grid is missing ID field `",
      grid_id_field,
      "`.",
      call. = FALSE
    )
  }
  grid_ids <- as.character(grid_ref[[grid_id_field]])
  if (anyNA(grid_ids) || anyDuplicated(grid_ids)) {
    stop(
      "The 20 km grid ID field must be complete and unique.",
      call. = FALSE
    )
  }

  source_mtime <- max(
    potential_maps_grid_source_mtime(grid_path),
    file.info(bounds_path)$mtime
  )
  cache_is_current <- file.exists(cache_path) &&
    file.info(cache_path)$mtime >= source_mtime

  if (!cache_is_current) {
    cli::cli_alert_info(
      "Creating the cached 20 km grid intersecting the exact GADM study area."
    )
    bounds_sf <- sf::st_as_sf(bounds_vec)
    if (is.na(sf::st_crs(bounds_sf))) {
      stop("The GADM study boundaries have no CRS.", call. = FALSE)
    }
    if (sf::st_crs(bounds_sf) != sf::st_crs(grid_ref)) {
      bounds_sf <- sf::st_transform(bounds_sf, sf::st_crs(grid_ref))
    }
    study_grid <- sf::st_filter(
      grid_ref,
      bounds_sf,
      .predicate = sf::st_intersects
    )
    cache_dir <- dirname(cache_path)
    dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
    if (!dir.exists(cache_dir)) {
      stop(
        "Could not create study-grid cache directory: ",
        cache_dir,
        call. = FALSE
      )
    }
    suppressMessages(sf::st_write(
      study_grid,
      cache_path,
      delete_dsn = file.exists(cache_path),
      quiet = TRUE
    ))
    cli::cli_alert_success(
      "Cached {nrow(study_grid)} study-area grid cells at {.path {cache_path}}."
    )
  } else {
    cli::cli_alert_info(
      "Using cached study-area 20 km grid: {.path {cache_path}}."
    )
    study_grid <- sf::read_sf(cache_path)
  }

  if (!grid_id_field %in% names(study_grid)) {
    stop(
      "The cached study grid is missing ID field `",
      grid_id_field,
      "`. Delete it and rerun to rebuild the cache.",
      call. = FALSE
    )
  }
  study_ids <- as.character(study_grid[[grid_id_field]])
  if (anyNA(study_ids) ||
      anyDuplicated(study_ids) ||
      !all(study_ids %in% grid_ids)) {
    stop(
      "The cached study-grid IDs are invalid. Delete the cache and rerun.",
      call. = FALSE
    )
  }
  if (sf::st_crs(study_grid) != sf::st_crs(grid_ref)) {
    study_grid <- sf::st_transform(study_grid, sf::st_crs(grid_ref))
  }

  study_grid
}

potential_maps_partition_species_grid <- function(
    species_grid,
    study_grid,
    grid_id_field = potential_maps_config$grid_id_field) {
  if (!inherits(species_grid, "sf") || !inherits(study_grid, "sf")) {
    stop("Species and study grids must both be sf objects.", call. = FALSE)
  }
  if (!grid_id_field %in% names(species_grid) ||
      !grid_id_field %in% names(study_grid)) {
    stop(
      "Species and study grids must both contain ID field `",
      grid_id_field,
      "`.",
      call. = FALSE
    )
  }

  species_ids <- as.character(species_grid[[grid_id_field]])
  study_ids <- as.character(study_grid[[grid_id_field]])
  inside_ids <- species_ids[species_ids %in% study_ids]

  list(
    inside = study_grid[
      as.character(study_grid[[grid_id_field]]) %in% inside_ids,
      ,
      drop = FALSE
    ],
    outside = species_grid[
      !species_ids %in% study_ids,
      ,
      drop = FALSE
    ]
  )
}

potential_maps_prepare_raster <- function(
    raster,
    scenario,
    rst_mask,
    context = "Raster") {
  if (!inherits(raster, "SpatRaster")) {
    stop(context, " is not a terra SpatRaster.", call. = FALSE)
  }
  if (!scenario %in% per_scn$scenario) {
    stop(context, " has an unknown scenario: ", scenario, call. = FALSE)
  }
  if (!isTRUE(terra::compareGeom(
    raster,
    rst_mask,
    lyrs = FALSE,
    stopOnError = FALSE,
    messages = FALSE
  ))) {
    stop(context, " does not match the European mask geometry.", call. = FALSE)
  }

  if (identical(scenario, "hist")) {
    return(raster)
  }

  terra::mask(raster, rst_mask)
}

potential_maps_prepare_display_raster <- function(
    raster,
    factor = potential_maps_config$display_aggregate_factor,
    fun = potential_maps_config$display_aggregate_fun) {
  if (!inherits(raster, "SpatRaster")) {
    stop("`raster` must be a terra SpatRaster.", call. = FALSE)
  }
  potential_maps_assert_positive_number(
    factor,
    "display_aggregate_factor",
    integer = TRUE
  )
  if (!is.character(fun) ||
      length(fun) != 1L ||
      is.na(fun) ||
      !identical(fun, "modal")) {
    stop(
      "`display_aggregate_fun` must be \"modal\" for binary maps.",
      call. = FALSE
    )
  }

  factor <- as.integer(factor)
  if (factor == 1L) {
    return(raster)
  }

  display_raster <- terra::aggregate(
    raster,
    fact = factor,
    fun = fun,
    na.rm = TRUE
  )
  names(display_raster) <- names(raster)
  display_raster
}

potential_maps_intersect_grid <- function(grid_ref, points, context = "Species") {
  if (!inherits(grid_ref, "sf")) {
    stop("`grid_ref` must be an sf object.", call. = FALSE)
  }
  if (!inherits(points, c("sf", "sfc"))) {
    stop("`points` must be an sf or sfc object.", call. = FALSE)
  }
  if (is.na(sf::st_crs(grid_ref))) {
    stop("The 20 km grid has no CRS.", call. = FALSE)
  }
  if (is.na(sf::st_crs(points))) {
    stop(context, " occurrence points have no CRS.", call. = FALSE)
  }

  non_empty <- !sf::st_is_empty(points)
  if (inherits(points, "sf")) {
    points <- points[non_empty, , drop = FALSE]
  } else {
    points <- points[non_empty]
  }
  if (length(points) == 0L ||
      (inherits(points, "sf") && nrow(points) == 0L)) {
    stop(context, " occurrence dataset contains no valid points.", call. = FALSE)
  }

  if (sf::st_crs(points) != sf::st_crs(grid_ref)) {
    points <- sf::st_transform(points, sf::st_crs(grid_ref))
  }

  sf::st_filter(
    grid_ref,
    sf::st_geometry(points),
    .predicate = sf::st_intersects
  )
}

potential_maps_load_species_grid <- function(point_path, grid_ref, sp_name) {
  potential_maps_assert_path(point_path, "point_path", must_exist = TRUE)
  points <- sf::st_read(point_path, quiet = TRUE)
  geometry_types <- unique(as.character(sf::st_geometry_type(points)))
  if (!all(geometry_types %in% c("POINT", "MULTIPOINT"))) {
    stop(
      "Filtered occurrence file for ",
      sp_name,
      " contains non-point geometry: ",
      paste(geometry_types, collapse = ", "),
      ".",
      call. = FALSE
    )
  }

  potential_maps_intersect_grid(
    grid_ref,
    points,
    context = paste0("Species ", sp_name)
  )
}

potential_maps_validate_occurrence_metadata <- function(
    point_path,
    sp_name,
    target_crs) {
  potential_maps_assert_path(point_path, "point_path", must_exist = TRUE)
  layers <- sf::st_layers(point_path)
  if (length(layers$name) != 1L) {
    stop(
      "Filtered occurrence GPKG for ",
      sp_name,
      " must contain exactly one layer; found ",
      length(layers$name),
      ".",
      call. = FALSE
    )
  }

  geometry_type <- as.character(layers$geomtype[[1L]])
  if (!tolower(geometry_type) %in% c("point", "multipoint")) {
    stop(
      "Filtered occurrence GPKG for ",
      sp_name,
      " has non-point layer geometry: ",
      geometry_type,
      ".",
      call. = FALSE
    )
  }
  if (is.na(layers$features[[1L]]) || layers$features[[1L]] < 1) {
    stop(
      "Filtered occurrence GPKG for ",
      sp_name,
      " contains no features.",
      call. = FALSE
    )
  }

  escaped_layer <- gsub(
    "\"",
    "\"\"",
    layers$name[[1L]],
    fixed = TRUE
  )
  sample_query <- paste0(
    "SELECT * FROM \"",
    escaped_layer,
    "\" LIMIT 1"
  )
  point_sample <- sf::st_read(
    point_path,
    query = sample_query,
    quiet = TRUE
  )
  if (nrow(point_sample) != 1L || is.na(sf::st_crs(point_sample))) {
    stop(
      "Filtered occurrence GPKG for ",
      sp_name,
      " has no readable feature with a valid CRS.",
      call. = FALSE
    )
  }

  # Transforming one feature verifies that the stored CRS is compatible with
  # the common raster/grid CRS without loading the full occurrence dataset.
  tryCatch(
    sf::st_transform(point_sample, target_crs),
    error = function(e) {
      stop(
        "Filtered occurrence GPKG for ",
        sp_name,
        " cannot be transformed to the model CRS: ",
        conditionMessage(e),
        call. = FALSE
      )
    }
  )

  invisible(TRUE)
}

potential_maps_preflight <- function(
    jobs,
    bounds_gpkg = potential_maps_config$bounds_gpkg,
    plot_bounds_gpkg = potential_maps_config$plot_bounds_gpkg,
    plot_bounds_tolerance_m =
      potential_maps_config$plot_bounds_tolerance_m,
    basemap_gpkg = potential_maps_config$basemap_gpkg,
    basemap_cache_gpkg = potential_maps_config$basemap_cache_gpkg,
    basemap_tolerance_m = potential_maps_config$basemap_tolerance_m,
    grid_path = potential_maps_config$grid_path,
    study_grid_gpkg = potential_maps_config$study_grid_gpkg,
    grid_id_field = potential_maps_config$grid_id_field,
    expected_species = 120L) {
  potential_maps_require_packages()
  potential_maps_assert_path(bounds_gpkg, "bounds_gpkg", must_exist = TRUE)
  potential_maps_assert_path(plot_bounds_gpkg, "plot_bounds_gpkg")
  potential_maps_assert_positive_number(
    plot_bounds_tolerance_m,
    "plot_bounds_tolerance_m"
  )
  potential_maps_assert_path(basemap_gpkg, "basemap_gpkg", must_exist = TRUE)
  potential_maps_assert_path(basemap_cache_gpkg, "basemap_cache_gpkg")
  potential_maps_assert_positive_number(
    basemap_tolerance_m,
    "basemap_tolerance_m"
  )
  potential_maps_assert_path(grid_path, "grid_path", must_exist = TRUE)
  potential_maps_assert_path(study_grid_gpkg, "study_grid_gpkg")

  if (!is.data.frame(jobs) || nrow(jobs) == 0L) {
    stop("`jobs` must be a non-empty potential-map job table.", call. = FALSE)
  }
  if (!is.null(expected_species)) {
    if (!is.numeric(expected_species) ||
        length(expected_species) != 1L ||
        is.na(expected_species) ||
        expected_species < 1) {
      stop("`expected_species` must be NULL or one positive number.", call. = FALSE)
    }
    observed_species <- length(unique(jobs$species_index))
    if (observed_species != expected_species) {
      stop(
        "Expected ",
        expected_species,
        " species but the cache produced ",
        observed_species,
        ".",
        call. = FALSE
      )
    }
  }
  if (nrow(jobs) != length(unique(jobs$species_index)) * nrow(per_scn)) {
    stop("Each species must produce exactly seven map jobs.", call. = FALSE)
  }

  missing_rasters <- unique(jobs$raster_path[!file.exists(jobs$raster_path)])
  if (length(missing_rasters) > 0L) {
    stop(
      "Missing model raster(s):\n- ",
      paste(utils::head(missing_rasters, 20L), collapse = "\n- "),
      call. = FALSE
    )
  }
  point_paths <- unique(jobs$point_path)
  missing_points <- point_paths[!file.exists(point_paths)]
  if (length(missing_points) > 0L) {
    stop(
      "Missing filtered occurrence GPKG(s):\n- ",
      paste(utils::head(missing_points, 20L), collapse = "\n- "),
      call. = FALSE
    )
  }

  template <- terra::rast(jobs$raster_path[[1L]])
  potential_maps_validate_binary_metadata(template, "Model template")

  raster_progress <- potential_maps_progress_bar(
    name = "Raster metadata preflight",
    total = nrow(jobs),
    .envir = environment()
  )
  raster_progress_done <- FALSE
  on.exit({
    if (!raster_progress_done) {
      try(cli::cli_progress_done(id = raster_progress), silent = TRUE)
    }
  }, add = TRUE)

  for (i in seq_len(nrow(jobs))) {
    candidate <- terra::rast(jobs$raster_path[[i]])
    context <- paste0(
      jobs$sp_name[[i]],
      " [",
      jobs$period[[i]],
      "/",
      jobs$scenario[[i]],
      "]"
    )
    potential_maps_progress_update(
      id = raster_progress,
      set = i - 1L,
      status = context
    )
    potential_maps_validate_binary_metadata(candidate, context)
    if (!isTRUE(terra::compareGeom(
      candidate,
      template,
      lyrs = FALSE,
      stopOnError = FALSE,
      messages = FALSE
    ))) {
      stop(context, " does not match the shared model geometry.", call. = FALSE)
    }
    potential_maps_progress_update(
      id = raster_progress,
      set = i,
      status = paste("Validated", context)
    )
  }
  cli::cli_progress_done(id = raster_progress)
  raster_progress_done <- TRUE

  cli::cli_alert_info("Building the shared GADM Europe raster mask.")
  bounds_vec <- terra::vect(bounds_gpkg)
  if (terra::nrow(bounds_vec) == 0L) {
    stop("The GADM boundary dataset is empty.", call. = FALSE)
  }
  rst_mask <- potential_maps_build_mask(bounds_vec, template)

  target_crs <- sf::st_crs(terra::crs(template))
  if (is.na(target_crs)) {
    stop("The model template has no valid CRS.", call. = FALSE)
  }

  cli::cli_alert_info("Loading simplified GADM study boundaries and the 20 km grid.")
  boundaries <- potential_maps_load_plot_boundaries(
    bounds_vec = bounds_vec,
    source_path = bounds_gpkg,
    cache_path = plot_bounds_gpkg,
    tolerance_m = plot_bounds_tolerance_m,
    target_crs = target_crs,
    layer_label = "GADM study-area boundaries"
  )

  cli::cli_alert_info("Loading the complete GADM 4.1 bounding-box basemap.")
  basemap_vec <- terra::vect(basemap_gpkg)
  if (terra::nrow(basemap_vec) == 0L) {
    stop("The GADM 4.1 basemap is empty.", call. = FALSE)
  }
  potential_maps_validate_basemap_countries(bounds_vec, basemap_vec)
  basemap_boundaries <- potential_maps_load_plot_boundaries(
    bounds_vec = basemap_vec,
    source_path = basemap_gpkg,
    cache_path = basemap_cache_gpkg,
    tolerance_m = basemap_tolerance_m,
    target_crs = target_crs,
    layer_label = "GADM 4.1 basemap boundaries"
  )

  grid_ref <- sf::read_sf(grid_path)
  if (nrow(grid_ref) == 0L) {
    stop("The 20 km reference grid is empty.", call. = FALSE)
  }
  if (is.na(sf::st_crs(grid_ref))) {
    stop("The 20 km reference grid has no CRS.", call. = FALSE)
  }
  if (sf::st_crs(grid_ref) != target_crs) {
    grid_ref <- sf::st_transform(grid_ref, target_crs)
  }

  study_grid <- potential_maps_load_study_grid(
    grid_ref = grid_ref,
    bounds_vec = bounds_vec,
    grid_path = grid_path,
    bounds_path = bounds_gpkg,
    cache_path = study_grid_gpkg,
    grid_id_field = grid_id_field
  )

  species_rows <- jobs[!duplicated(jobs$species_index), , drop = FALSE]
  occurrence_progress <- potential_maps_progress_bar(
    name = "Filtered occurrence preflight",
    total = nrow(species_rows),
    .envir = environment()
  )
  occurrence_progress_done <- FALSE
  on.exit({
    if (!occurrence_progress_done) {
      try(cli::cli_progress_done(id = occurrence_progress), silent = TRUE)
    }
  }, add = TRUE)

  for (i in seq_len(nrow(species_rows))) {
    potential_maps_progress_update(
      id = occurrence_progress,
      set = i - 1L,
      status = species_rows$sp_name[[i]]
    )
    potential_maps_validate_occurrence_metadata(
      point_path = species_rows$point_path[[i]],
      sp_name = species_rows$sp_name[[i]],
      target_crs = target_crs
    )
    potential_maps_progress_update(
      id = occurrence_progress,
      set = i,
      status = paste("Validated", species_rows$sp_name[[i]])
    )
  }
  cli::cli_progress_done(id = occurrence_progress)
  occurrence_progress_done <- TRUE

  list(
    template = template,
    rst_mask = rst_mask,
    basemap_boundaries = basemap_boundaries,
    boundaries = boundaries,
    grid_ref = grid_ref,
    study_grid = study_grid
  )
}

# ==============================================================================
# Map rendering
# ==============================================================================

potential_maps_assert_binary_vector <- function(values, context) {
  chunk_size <- 1000000L
  starts <- seq.int(1L, length(values), by = chunk_size)
  for (start in starts) {
    end <- min(length(values), start + chunk_size - 1L)
    chunk <- values[start:end]
    invalid <- !is.na(chunk) & chunk != 0 & chunk != 1
    if (any(invalid)) {
      stop(context, " contains non-binary cell values.", call. = FALSE)
    }
  }
  invisible(values)
}

potential_maps_raster_to_df <- function(raster, context = "Raster") {
  raster_df <- terra::as.data.frame(
    raster,
    xy = TRUE,
    na.rm = FALSE
  )
  if (ncol(raster_df) != 3L) {
    stop(context, " did not convert to one x/y/value table.", call. = FALSE)
  }
  names(raster_df) <- c("x", "y", "suitability")
  potential_maps_assert_binary_vector(raster_df$suitability, context)
  raw_suitability <- raster_df$suitability
  suitability_class <- rep(
    potential_maps_condition_levels[[3L]],
    length(raw_suitability)
  )
  suitability_class[!is.na(raw_suitability) & raw_suitability == 0] <-
    potential_maps_condition_levels[[2L]]
  suitability_class[!is.na(raw_suitability) & raw_suitability == 1] <-
    potential_maps_condition_levels[[1L]]
  raster_df$suitability <- factor(
    suitability_class,
    levels = potential_maps_condition_levels
  )
  raster_df
}

potential_maps_full_disclaimer <- function(
    mod_type_disclaimer,
    config = potential_maps_config) {
  fields <- c(
    "disclaimer_text",
    "basemap_disclaimer",
    "version",
    "ref_date"
  )
  missing_fields <- setdiff(fields, names(config))
  if (length(missing_fields) > 0L) {
    stop(
      "Map configuration is missing disclaimer field(s): ",
      paste(missing_fields, collapse = ", "),
      ".",
      call. = FALSE
    )
  }
  invalid <- !vapply(
    config[fields],
    function(value) {
      is.character(value) &&
        length(value) == 1L &&
        !is.na(value) &&
        nzchar(value)
    },
    logical(1)
  )
  if (any(invalid)) {
    stop(
      paste0(
        "Disclaimer text, basemap information, version, and reference date ",
        "must be non-empty strings."
      ),
      call. = FALSE
    )
  }
  allowed_model_types <- c(
    "Climate-only model",
    "Climate-habitat ensemble"
  )
  if (!is.character(mod_type_disclaimer) ||
      length(mod_type_disclaimer) != 1L ||
      is.na(mod_type_disclaimer) ||
      !mod_type_disclaimer %in% allowed_model_types) {
    stop(
      "`mod_type_disclaimer` must be one of: ",
      paste(allowed_model_types, collapse = ", "),
      ".",
      call. = FALSE
    )
  }

  paste0(
    config$disclaimer_text,
    " | Model represented: ",
    mod_type_disclaimer,
    " | ",
    config$basemap_disclaimer,
    " | Version: ",
    config$version,
    " | Reference date: ",
    config$ref_date
  )
}

potential_maps_prepare_occurrence_grid <- function(
    grid_inside,
    grid_outside) {
  if (!inherits(grid_inside, "sf") || !inherits(grid_outside, "sf")) {
    stop("Inside and outside occurrence grids must be sf objects.", call. = FALSE)
  }
  if (is.na(sf::st_crs(grid_inside)) || is.na(sf::st_crs(grid_outside))) {
    stop("Occurrence-grid layers must have valid CRSs.", call. = FALSE)
  }
  if (sf::st_crs(grid_outside) != sf::st_crs(grid_inside)) {
    grid_outside <- sf::st_transform(grid_outside, sf::st_crs(grid_inside))
  }

  make_layer <- function(grid, occurrence_class, linewidth, alpha) {
    sf::st_sf(
      occurrence_class = factor(
        rep(occurrence_class, nrow(grid)),
        levels = potential_maps_occurrence_levels
      ),
      occurrence_linewidth = rep(linewidth, nrow(grid)),
      occurrence_alpha = rep(alpha, nrow(grid)),
      geometry = sf::st_geometry(grid)
    )
  }

  # Draw outside-study cells first so the emphasized blue cells remain on top.
  rbind(
    make_layer(
      grid_outside,
      potential_maps_occurrence_levels[[2L]],
      linewidth = 0.12,
      alpha = 0.65
    ),
    make_layer(
      grid_inside,
      potential_maps_occurrence_levels[[1L]],
      linewidth = 0.15,
      alpha = 0.8
    )
  )
}

create_sdm_map_v2 <- function(
    raster_df,
    basemap_boundaries,
    boundaries,
    grid_inside,
    grid_outside,
    title,
    subtitle,
    mod_type_disclaimer,
    config = potential_maps_config,
    extent_buffer = 1.025) {
  x_range <- range(raster_df$x, na.rm = TRUE)
  y_range <- range(raster_df$y, na.rm = TRUE)
  x_center <- mean(x_range)
  y_center <- mean(y_range)
  x_span <- diff(x_range) * extent_buffer
  y_span <- diff(y_range) * extent_buffer
  xlim <- c(x_center - x_span / 2, x_center + x_span / 2)
  ylim <- c(y_center - y_span / 2, y_center + y_span / 2)
  if (!is.numeric(config$legend_position_inside) ||
      length(config$legend_position_inside) != 2L ||
      anyNA(config$legend_position_inside) ||
      any(config$legend_position_inside < 0) ||
      any(config$legend_position_inside > 1)) {
    stop(
      "`config$legend_position_inside` must contain two values from 0 to 1.",
      call. = FALSE
    )
  }
  occurrence_grid <- potential_maps_prepare_occurrence_grid(
    grid_inside,
    grid_outside
  )
  full_disclaimer <- potential_maps_full_disclaimer(
    mod_type_disclaimer = mod_type_disclaimer,
    config = config
  )

  ggplot2::ggplot() +
    ggplot2::geom_raster(
      data = raster_df,
      ggplot2::aes(x = x, y = y, fill = suitability)
    ) +

    ggplot2::geom_sf(
      data = basemap_boundaries,
      fill = NA,
      color = config$color_basemap_border,
      linewidth = 0.08,
      alpha = 0.65
    ) +
    ggplot2::geom_sf(
      data = boundaries,
      fill = NA,
      color = config$color_border,
      linewidth = 0.1,
      alpha = 0.6
    ) +
    
    ggplot2::geom_sf(
      data = occurrence_grid,
      ggplot2::aes(
        colour = occurrence_class,
        linewidth = occurrence_linewidth,
        alpha = occurrence_alpha
      ),
      fill = NA,
      show.legend = TRUE
    ) +
    ggplot2::scale_fill_manual(
      values = stats::setNames(
        c(
          config$color_suitable,
          config$color_unsuitable,
          config$color_outside_area
        ),
        potential_maps_condition_levels
      ),
      limits = potential_maps_condition_levels,
      breaks = potential_maps_condition_levels,
      name = "Modelled conditions",
      na.value = config$color_outside_area,
      drop = FALSE
    ) +
    ggplot2::scale_colour_manual(
      values = stats::setNames(
        c(config$color_grid, config$color_grid_outside),
        potential_maps_occurrence_levels
      ),
      limits = potential_maps_occurrence_levels,
      breaks = potential_maps_occurrence_levels,
      name = "Occurrence data (20×20 km cells)",
      na.value = config$color_grid_outside,
      drop = FALSE
    ) +
    ggplot2::scale_linewidth_identity(guide = "none") +
    ggplot2::scale_alpha_identity(guide = "none") +
    ggplot2::guides(
      fill = ggplot2::guide_legend(
        order = 1,
        title.position = "top",
        override.aes = list(
          colour = config$color_legend_border,
          linewidth = 0.25
        )
      ),
      colour = ggplot2::guide_legend(
        order = 2,
        title.position = "top",
        override.aes = list(
          fill = NA,
          alpha = 1,
          linewidth = 0.65
        )
      )
    ) +
    ggspatial::annotation_scale(
      location = "br",
      width_hint = 0.25,
      pad_x = grid::unit(0.2, "in"),
      pad_y = grid::unit(0.2, "in"),
      style = "ticks",
      line_col = "grey40",
      text_col = "grey40"
    ) +
    ggplot2::coord_sf(xlim = xlim, ylim = ylim, expand = FALSE) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(
        hjust = 0.5,
        face = "italic",
        size = 16
      ),
      plot.subtitle = ggplot2::element_text(
        hjust = 0.5,
        size = 10,
        color = "grey40"
      ),
      legend.position = "inside",
      legend.position.inside = config$legend_position_inside,
      legend.justification = c(0, 0.5),
      legend.box = "vertical",
      legend.box.just = "left",
      legend.background = ggplot2::element_blank(),
      legend.box.background = ggplot2::element_rect(
        fill = grDevices::adjustcolor("white", alpha.f = 0.92),
        colour = config$color_legend_border,
        linewidth = 0.3
      ),
      legend.margin = ggplot2::margin(4, 5, 4, 5),
      legend.spacing.y = grid::unit(2, "pt"),
      legend.title = ggplot2::element_text(
        face = "bold",
        size = 8.2,
        lineheight = 0.95,
        hjust = 0
      ),
      legend.text = ggplot2::element_text(size = 7.3),
      legend.key = ggplot2::element_rect(
        fill = "white",
        colour = config$color_legend_border,
        linewidth = 0.25
      ),
      legend.key.height = grid::unit(0.38, "cm"),
      legend.key.width = grid::unit(0.48, "cm"),
      panel.grid.major = ggplot2::element_line(
        color = "grey90",
        linewidth = 0.2
      ),
      panel.grid.minor = ggplot2::element_blank(),
      panel.background = ggplot2::element_rect(fill = "white", color = NA),
      plot.background = ggplot2::element_rect(fill = "white", color = NA),
      axis.title = ggplot2::element_text(size = 10),
      plot.caption.position = "plot",
      plot.caption = ggtext::element_textbox_simple(
        fill = "black",
        colour = "white",
        size = 7.3,
        halign = 0.5,
        padding = ggplot2::margin(2.5, 5, 2.5, 5),
        margin = ggplot2::margin(4, 0, 0, 0),
        width = grid::unit(1, "npc")
      ),
      plot.margin = ggplot2::margin(10, 10, 3, 10)
    ) +
    ggplot2::labs(
      title = title,
      subtitle = subtitle,
      x = "Longitude",
      y = "Latitude",
      caption = full_disclaimer
    )
}

save_sdm_map_v2 <- function(
    plot,
    output_path,
    config = potential_maps_config) {
  potential_maps_assert_path(output_path, "output_path")
  ggplot2::ggsave(
    filename = output_path,
    plot = plot,
    width = round(config$map_width * config$map_dpi),
    height = round(config$map_height * config$map_dpi),
    units = "px",
    dpi = config$map_dpi,
    bg = "white"
  )
  invisible(output_path)
}

# ==============================================================================
# Summary-statistics persistence
# ==============================================================================

potential_maps_empty_summary <- function() {
  data.frame(
    sp_name = character(0),
    proj_code = character(0),
    total_cells = numeric(0),
    na_cells = numeric(0),
    valid_cells = numeric(0),
    suitable_cells = numeric(0),
    unsuitable_cells = numeric(0),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
}

potential_maps_load_summary <- function(path) {
  if (!file.exists(path)) {
    return(potential_maps_empty_summary())
  }

  summary <- utils::read.csv(
    path,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  missing_columns <- setdiff(potential_maps_summary_columns, names(summary))
  if (length(missing_columns) > 0L) {
    stop(
      "Existing summary CSV is missing column(s): ",
      paste(missing_columns, collapse = ", "),
      ".",
      call. = FALSE
    )
  }
  summary <- summary[, potential_maps_summary_columns, drop = FALSE]

  key <- paste(summary$sp_name, summary$proj_code, sep = "\r")
  summary[!duplicated(key, fromLast = TRUE), , drop = FALSE]
}

potential_maps_has_summary <- function(summary, sp_name, proj_name) {
  any(summary$sp_name == sp_name & summary$proj_code == proj_name)
}

potential_maps_upsert_summary <- function(summary, row) {
  if (!is.data.frame(row) || nrow(row) != 1L) {
    stop("`row` must contain exactly one summary record.", call. = FALSE)
  }
  missing_columns <- setdiff(potential_maps_summary_columns, names(row))
  if (length(missing_columns) > 0L) {
    stop(
      "Summary record is missing column(s): ",
      paste(missing_columns, collapse = ", "),
      ".",
      call. = FALSE
    )
  }

  keep <- !(
    summary$sp_name == row$sp_name[[1L]] &
      summary$proj_code == row$proj_code[[1L]]
  )
  result <- rbind(
    summary[keep, potential_maps_summary_columns, drop = FALSE],
    row[, potential_maps_summary_columns, drop = FALSE]
  )
  rownames(result) <- NULL
  result
}

potential_maps_write_summary <- function(summary, path) {
  potential_maps_assert_path(path, "summary_path")
  utils::write.csv(summary, path, row.names = FALSE)
  invisible(path)
}

potential_maps_summary_stats <- function(
    raster,
    sp_name,
    proj_name,
    report = TRUE) {
  potential_maps_assert_scalar_logical(report, "report")
  counts <- terra::global(
    raster,
    fun = c("isNA", "notNA", "sum"),
    na.rm = TRUE
  )
  total_cells <- terra::ncell(raster)
  na_cells <- as.numeric(counts$isNA[[1L]])
  valid_cells <- as.numeric(counts$notNA[[1L]])
  suitable_cells <- as.numeric(counts$sum[[1L]])
  unsuitable_cells <- valid_cells - suitable_cells

  valid_percent <- if (valid_cells > 0) {
    100 * valid_cells / total_cells
  } else {
    0
  }
  suitable_percent <- if (valid_cells > 0) {
    100 * suitable_cells / valid_cells
  } else {
    0
  }
  unsuitable_percent <- if (valid_cells > 0) {
    100 * unsuitable_cells / valid_cells
  } else {
    0
  }

  cell_area <- prod(terra::res(raster))
  if (report) {
    cli::cli_h3("Summary statistics")
    cli::cli_bullets(c(
      "*" = sprintf("Total cells: %d", total_cells),
      "*" = sprintf(
        "Valid cells: %d (%.1f%%)",
        valid_cells,
        valid_percent
      ),
      "*" = sprintf(
        "Suitable habitat: %d cells (%.1f%% of valid)",
        suitable_cells,
        suitable_percent
      ),
      "*" = sprintf(
        "Unsuitable habitat: %d cells (%.1f%% of valid)",
        unsuitable_cells,
        unsuitable_percent
      ),
      "*" = sprintf(
        "Suitable area: %.2f km² (assuming meters)",
        suitable_cells * cell_area / 1e6
      )
    ))
  }

  # Keep v1's column schema. Historically `proj_code` stored the printable
  # projection label supplied to print_summary_stats(), not the filename code.
  data.frame(
    sp_name = sp_name,
    proj_code = proj_name,
    total_cells = total_cells,
    na_cells = na_cells,
    valid_cells = valid_cells,
    suitable_cells = suitable_cells,
    unsuitable_cells = unsuitable_cells,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
}

# ==============================================================================
# Batch runner
# ==============================================================================

run_potential_distribution_maps <- function(
    cache_path = potential_maps_config$cache_path,
    bounds_gpkg = potential_maps_config$bounds_gpkg,
    plot_bounds_gpkg = potential_maps_config$plot_bounds_gpkg,
    plot_bounds_tolerance_m =
      potential_maps_config$plot_bounds_tolerance_m,
    basemap_gpkg = potential_maps_config$basemap_gpkg,
    basemap_cache_gpkg = potential_maps_config$basemap_cache_gpkg,
    basemap_tolerance_m = potential_maps_config$basemap_tolerance_m,
    grid_path = potential_maps_config$grid_path,
    study_grid_gpkg = potential_maps_config$study_grid_gpkg,
    grid_id_field = potential_maps_config$grid_id_field,
    occurrence_dir = potential_maps_config$occurrence_dir,
    output_dir = potential_maps_config$output_dir,
    overwrite = FALSE,
    expected_species = 120L,
    config = potential_maps_config) {
  potential_maps_require_packages()
  potential_maps_assert_scalar_logical(overwrite, "overwrite")
  potential_maps_assert_positive_number(
    config$display_aggregate_factor,
    "config$display_aggregate_factor",
    integer = TRUE
  )
  if (!is.character(config$display_aggregate_fun) ||
      length(config$display_aggregate_fun) != 1L ||
      is.na(config$display_aggregate_fun) ||
      !identical(config$display_aggregate_fun, "modal")) {
    stop(
      "`config$display_aggregate_fun` must be \"modal\".",
      call. = FALSE
    )
  }
  potential_maps_assert_path(cache_path, "cache_path", must_exist = TRUE)
  potential_maps_assert_path(occurrence_dir, "occurrence_dir", must_exist = TRUE)
  potential_maps_assert_path(output_dir, "output_dir")

  cli::cli_h1("wiSDM 2.0 species distribution mapping")
  cli::cli_alert_info("Loading the cached best-model raster catalog.")
  catalog <- readRDS(cache_path)
  jobs <- build_potential_map_jobs(
    catalog = catalog,
    output_dir = output_dir,
    occurrence_dir = occurrence_dir
  )
  cli::cli_alert_success(
    "Prepared {nrow(jobs)} map jobs for {length(unique(jobs$species_index))} species."
  )

  cli::cli_h2("Input preflight")
  spatial <- potential_maps_preflight(
    jobs,
    bounds_gpkg = bounds_gpkg,
    plot_bounds_gpkg = plot_bounds_gpkg,
    plot_bounds_tolerance_m = plot_bounds_tolerance_m,
    basemap_gpkg = basemap_gpkg,
    basemap_cache_gpkg = basemap_cache_gpkg,
    basemap_tolerance_m = basemap_tolerance_m,
    grid_path = grid_path,
    study_grid_gpkg = study_grid_gpkg,
    grid_id_field = grid_id_field,
    expected_species = expected_species
  )
  cli::cli_alert_success("Raster, boundary, grid, and occurrence preflight passed.")

  cli::cli_h2("Output preparation")
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  if (!dir.exists(output_dir)) {
    stop("Could not create output directory: ", output_dir, call. = FALSE)
  }
  summary_path <- file.path(output_dir, "summary_stats_by_scn.csv")
  all_sum_stats <- potential_maps_load_summary(summary_path)
  cli::cli_alert_info("Maps will be written to {.path {output_dir}}.")
  cli::cli_alert_info("Summary rows will be written to {.path {summary_path}}.")

  cli::cli_h2("Rendering maps")
  render_progress <- potential_maps_progress_bar(
    name = "Potential-distribution maps",
    total = nrow(jobs),
    .envir = environment()
  )
  render_progress_done <- FALSE
  on.exit({
    if (!render_progress_done) {
      try(cli::cli_progress_done(id = render_progress), silent = TRUE)
    }
  }, add = TRUE)

  species_indices <- unique(jobs$species_index)
  completed <- 0L
  for (species_index in species_indices) {
    species_jobs <- jobs[
      jobs$species_index == species_index,
      ,
      drop = FALSE
    ]

    needs_render <- overwrite | !file.exists(species_jobs$output_path)
    species_grid_layers <- NULL
    if (any(needs_render)) {
      potential_maps_progress_update(
        id = render_progress,
        set = completed,
        status = paste(
          "Selecting blue 20 km grid cells for",
          species_jobs$sp_name[[1L]]
        )
      )
      species_grid <- potential_maps_load_species_grid(
        point_path = species_jobs$point_path[[1L]],
        grid_ref = spatial$grid_ref,
        sp_name = species_jobs$sp_name[[1L]]
      )
      species_grid_layers <- potential_maps_partition_species_grid(
        species_grid = species_grid,
        study_grid = spatial$study_grid,
        grid_id_field = grid_id_field
      )
    }

    for (job_index in seq_len(nrow(species_jobs))) {
      job <- species_jobs[job_index, , drop = FALSE]
      completed <- completed + 1L
      output_exists <- file.exists(job$output_path[[1L]])
      summary_exists <- potential_maps_has_summary(
        all_sum_stats,
        job$display_name[[1L]],
        job$proj_name[[1L]]
      )

      job_label <- paste0(
        job$sp_name[[1L]],
        " | ",
        job$period[[1L]],
        "/",
        job$scenario[[1L]]
      )
      potential_maps_progress_update(
        id = render_progress,
        set = completed - 1L,
        status = paste("Checking", job_label)
      )

      if (!overwrite && output_exists && summary_exists) {
        potential_maps_progress_update(
          id = render_progress,
          set = completed,
          status = paste("Skipped existing", job_label)
        )
        next
      }

      context <- paste0(
        job$sp_name[[1L]],
        " [",
        job$period[[1L]],
        "/",
        job$scenario[[1L]],
        "]"
      )
      potential_maps_progress_update(
        id = render_progress,
        set = completed - 1L,
        status = paste(
          if (identical(job$scenario[[1L]], "hist")) {
            "Loading historical raster for"
          } else {
            "Loading and masking future raster for"
          },
          job_label
        )
      )
      sdm_raster <- tryCatch(
        {
          loaded <- terra::rast(job$raster_path[[1L]])
          potential_maps_validate_binary_metadata(loaded, context)
          potential_maps_prepare_raster(
            loaded,
            scenario = job$scenario[[1L]],
            rst_mask = spatial$rst_mask,
            context = context
          )
        },
        error = function(e) {
          stop(context, ": ", conditionMessage(e), call. = FALSE)
        }
      )

      if (overwrite || !summary_exists) {
        potential_maps_progress_update(
          id = render_progress,
          set = completed - 1L,
          status = paste("Calculating summary statistics for", job_label)
        )
        sum_stats <- potential_maps_summary_stats(
          sdm_raster,
          sp_name = job$display_name[[1L]],
          proj_name = job$proj_name[[1L]],
          report = FALSE
        )
        all_sum_stats <- potential_maps_upsert_summary(
          all_sum_stats,
          sum_stats
        )
        potential_maps_write_summary(all_sum_stats, summary_path)
      }

      if (overwrite || !output_exists) {
        potential_maps_progress_update(
          id = render_progress,
          set = completed - 1L,
          status = paste0(
            "Aggregating display raster ",
            config$display_aggregate_factor,
            "x with ",
            config$display_aggregate_fun,
            " for ",
            job_label
          )
        )
        display_raster <- potential_maps_prepare_display_raster(
          sdm_raster,
          factor = config$display_aggregate_factor,
          fun = config$display_aggregate_fun
        )
        potential_maps_progress_update(
          id = render_progress,
          set = completed - 1L,
          status = paste("Preparing display pixels for", job_label)
        )
        raster_df <- potential_maps_raster_to_df(display_raster, context)
        potential_maps_progress_update(
          id = render_progress,
          set = completed - 1L,
          status = paste("Building map layers for", job_label)
        )
        sdm_map <- create_sdm_map_v2(
          raster_df = raster_df,
          basemap_boundaries = spatial$basemap_boundaries,
          boundaries = spatial$boundaries,
          grid_inside = species_grid_layers$inside,
          grid_outside = species_grid_layers$outside,
          title = job$display_name[[1L]],
          subtitle = job$proj_name[[1L]],
          mod_type_disclaimer = job$mod_type_disclaimer[[1L]],
          config = config
        )
        potential_maps_progress_update(
          id = render_progress,
          set = completed - 1L,
          status = paste("Saving PNG for", job_label)
        )
        save_sdm_map_v2(
          plot = sdm_map,
          output_path = job$output_path[[1L]],
          config = config
        )
        rm(display_raster, raster_df, sdm_map)
      }

      rm(sdm_raster)
      invisible(gc(verbose = FALSE))
      potential_maps_progress_update(
        id = render_progress,
        set = completed,
        status = paste(
          if (output_exists && !overwrite) {
            "Restored summary for"
          } else {
            "Finished"
          },
          job_label
        )
      )
    }
  }

  cli::cli_progress_done(id = render_progress)
  render_progress_done <- TRUE
  cli::cli_alert_success(
    "Completed {nrow(jobs)} map jobs. Summary: {.path {summary_path}}"
  )

  invisible(list(
    jobs = jobs,
    summary = all_sum_stats,
    summary_path = summary_path
  ))
}

potential_maps_autorun_on_source <- getOption(
  "wisdm.run_potential_maps_on_source",
  interactive()
)

if (sys.nframe() == 0L || isTRUE(potential_maps_autorun_on_source)) {
  run_potential_distribution_maps()
}
