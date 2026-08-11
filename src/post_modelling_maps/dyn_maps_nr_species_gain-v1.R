# Species gain-overlap maps from historical-to-future dynamics rasters.
#
# Direct execution runs the complete batch:
# Rscript src/post_modelling_maps/dyn_maps_nr_species_gain-v1.R
#
# Sourcing loads the functions without starting the batch:
# source("src/post_modelling_maps/dyn_maps_nr_species_gain-v1.R")
#
# To run while sourcing, opt in explicitly:
# options(wisdm.run_species_gain_aggregations_on_source = TRUE)
# source("src/post_modelling_maps/dyn_maps_nr_species_gain-v1.R")

gain_overlap_expected_combinations <- data.frame(
  future_period = rep(
    c("mid_2041_2070", "late_2071_2100"),
    each = 3L
  ),
  future_scenario = rep(
    c("ssp126", "ssp370", "ssp585"),
    times = 2L
  ),
  stringsAsFactors = FALSE
)

gain_overlap_expected_levels <- data.frame(
  value = 1:4,
  dynamics = c(
    "Stable unsuitable (0->0)",
    "Gain (0->1)",
    "Loss (1->0)",
    "Stable suitable (1->1)"
  ),
  stringsAsFactors = FALSE
)

gain_overlap_timestamp <- function() {
  format(Sys.time(), tz = "UTC", usetz = TRUE)
}

gain_overlap_message <- function(show_progress, ...) {
  if (isTRUE(show_progress)) {
    message(...)
  }
  invisible(NULL)
}

gain_overlap_assert_scalar_logical <- function(value, argument) {
  if (!is.logical(value) || length(value) != 1L || is.na(value)) {
    stop("`", argument, "` must be TRUE or FALSE.", call. = FALSE)
  }
  invisible(value)
}

gain_overlap_assert_path <- function(path, argument) {
  if (!is.character(path) ||
      length(path) != 1L ||
      is.na(path) ||
      !nzchar(trimws(path))) {
    stop("`", argument, "` must be one non-empty path.", call. = FALSE)
  }
  invisible(path)
}

gain_overlap_require_packages <- function() {
  required <- c("terra", "readxl")
  missing <- required[
    !vapply(required, requireNamespace, logical(1), quietly = TRUE)
  ]
  if (length(missing) > 0L) {
    stop(
      "Required R package(s) are not installed: ",
      paste(missing, collapse = ", "),
      ".",
      call. = FALSE
    )
  }
  invisible(TRUE)
}

gain_overlap_clean_slug <- function(x) {
  ascii <- iconv(
    enc2utf8(as.character(x)),
    from = "UTF-8",
    to = "ASCII//TRANSLIT",
    sub = ""
  )
  ascii[is.na(ascii)] <- as.character(x)[is.na(ascii)]
  slug <- tolower(ascii)
  slug <- gsub("[^a-z0-9]+", "_", slug)
  slug <- gsub("_+", "_", slug)
  slug <- gsub("^_|_$", "", slug)
  slug
}

gain_overlap_resolve_source_path <- function(path, dynamics_dir) {
  path <- path.expand(as.character(path))
  candidates <- unique(c(
    path,
    file.path(dynamics_dir, basename(path))
  ))
  existing <- candidates[file.exists(candidates)]
  if (length(existing) == 0L) {
    return(path)
  }
  normalizePath(existing[[1L]], winslash = "/", mustWork = TRUE)
}

gain_overlap_read_manifest <- function(dynamics_dir) {
  gain_overlap_assert_path(dynamics_dir, "dynamics_dir")
  manifest_path <- file.path(
    dynamics_dir,
    "sp_hist_fut_dynamics_manifest.rds"
  )
  if (!file.exists(manifest_path)) {
    stop(
      "Dynamics manifest does not exist: ",
      manifest_path,
      call. = FALSE
    )
  }

  manifest <- readRDS(manifest_path)
  required <- c(
    "sp_name",
    "tkey",
    "species_stem",
    "future_period",
    "future_scenario",
    "output_path"
  )
  missing_columns <- setdiff(required, names(manifest))
  if (length(missing_columns) > 0L) {
    stop(
      "Dynamics manifest is missing required column(s): ",
      paste(missing_columns, collapse = ", "),
      ".",
      call. = FALSE
    )
  }
  if (nrow(manifest) == 0L) {
    stop("Dynamics manifest is empty.", call. = FALSE)
  }

  manifest$sp_name <- trimws(as.character(manifest$sp_name))
  manifest$tkey <- as.character(manifest$tkey)
  manifest$species_stem <- as.character(manifest$species_stem)
  manifest$future_period <- as.character(manifest$future_period)
  manifest$future_scenario <- as.character(manifest$future_scenario)
  manifest$source_path <- vapply(
    manifest$output_path,
    gain_overlap_resolve_source_path,
    character(1),
    dynamics_dir = dynamics_dir
  )
  manifest$combo_id <- paste(
    manifest$future_period,
    manifest$future_scenario,
    sep = "__"
  )

  missing_files <- unique(
    manifest$source_path[!file.exists(manifest$source_path)]
  )
  if (length(missing_files) > 0L) {
    stop(
      "Dynamics manifest references missing raster(s):\n- ",
      paste(missing_files, collapse = "\n- "),
      call. = FALSE
    )
  }
  if (any(!nzchar(manifest$sp_name)) ||
      any(!nzchar(manifest$tkey)) ||
      any(!nzchar(manifest$species_stem))) {
    stop(
      "Dynamics manifest contains empty species identifiers.",
      call. = FALSE
    )
  }

  expected_ids <- paste(
    gain_overlap_expected_combinations$future_period,
    gain_overlap_expected_combinations$future_scenario,
    sep = "__"
  )
  actual_ids <- unique(manifest$combo_id)
  missing_combinations <- setdiff(expected_ids, actual_ids)
  unexpected_combinations <- setdiff(actual_ids, expected_ids)
  if (length(missing_combinations) > 0L ||
      length(unexpected_combinations) > 0L) {
    stop(
      "Dynamics manifest period/scenario combinations differ from the ",
      "expected six. Missing: ",
      if (length(missing_combinations)) {
        paste(missing_combinations, collapse = ", ")
      } else {
        "none"
      },
      "; unexpected: ",
      if (length(unexpected_combinations)) {
        paste(unexpected_combinations, collapse = ", ")
      } else {
        "none"
      },
      ".",
      call. = FALSE
    )
  }

  species <- sort(unique(manifest$sp_name))
  if (length(species) != 120L) {
    stop(
      "Expected 120 unique species in the dynamics manifest; found ",
      length(species),
      ".",
      call. = FALSE
    )
  }
  for (combo_id in expected_ids) {
    rows <- manifest[manifest$combo_id == combo_id, , drop = FALSE]
    if (nrow(rows) != length(species) ||
        anyDuplicated(rows$sp_name) ||
        !setequal(rows$sp_name, species)) {
      stop(
        "Combination ",
        combo_id,
        " does not contain exactly one raster for each of the 120 species.",
        call. = FALSE
      )
    }
  }
  if (anyDuplicated(manifest$source_path)) {
    stop(
      "Dynamics manifest reuses one or more raster paths across jobs.",
      call. = FALSE
    )
  }

  manifest <- manifest[
    order(
      match(manifest$combo_id, expected_ids),
      manifest$sp_name
    ),
    ,
    drop = FALSE
  ]
  rownames(manifest) <- NULL
  attr(manifest, "manifest_path") <- normalizePath(
    manifest_path,
    winslash = "/",
    mustWork = TRUE
  )
  manifest
}

gain_overlap_read_groups <- function(groups_path) {
  gain_overlap_assert_path(groups_path, "groups_path")
  if (!file.exists(groups_path)) {
    stop(
      "Species-group workbook does not exist: ",
      groups_path,
      call. = FALSE
    )
  }

  groups <- as.data.frame(
    readxl::read_excel(
      groups_path,
      sheet = "sp_names"
    ),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  required <- c("Sp_name", "Organism_Type")
  missing_columns <- setdiff(required, names(groups))
  if (length(missing_columns) > 0L) {
    stop(
      "Species-group workbook is missing required column(s): ",
      paste(missing_columns, collapse = ", "),
      ".",
      call. = FALSE
    )
  }

  groups$Sp_name <- trimws(as.character(groups$Sp_name))
  groups$Organism_Type <- trimws(as.character(groups$Organism_Type))
  if (nrow(groups) == 0L ||
      any(is.na(groups$Sp_name)) ||
      any(!nzchar(groups$Sp_name)) ||
      any(is.na(groups$Organism_Type)) ||
      any(!nzchar(groups$Organism_Type))) {
    stop(
      "Species-group workbook contains empty species or Organism_Type values.",
      call. = FALSE
    )
  }
  if (anyDuplicated(groups$Sp_name)) {
    duplicated_names <- unique(
      groups$Sp_name[
        duplicated(groups$Sp_name) |
          duplicated(groups$Sp_name, fromLast = TRUE)
      ]
    )
    stop(
      "Species-group workbook contains duplicated Sp_name value(s): ",
      paste(duplicated_names, collapse = ", "),
      ".",
      call. = FALSE
    )
  }
  groups
}

gain_overlap_build_membership <- function(manifest, groups) {
  manifest_species <- unique(manifest$sp_name)
  workbook_species <- unique(groups$Sp_name)
  workbook_only <- setdiff(workbook_species, manifest_species)
  manifest_only <- setdiff(manifest_species, workbook_species)
  if (length(workbook_only) > 0L ||
      length(manifest_only) > 0L) {
    stop(
      "Sp_name does not match the dynamics manifest exactly. Workbook-only: ",
      if (length(workbook_only)) {
        paste(workbook_only, collapse = ", ")
      } else {
        "none"
      },
      "; manifest-only: ",
      if (length(manifest_only)) {
        paste(manifest_only, collapse = ", ")
      } else {
        "none"
      },
      ".",
      call. = FALSE
    )
  }

  species_rows <- manifest[
    !duplicated(manifest$sp_name),
    c("sp_name", "tkey", "species_stem"),
    drop = FALSE
  ]
  species_rows <- species_rows[
    order(species_rows$sp_name),
    ,
    drop = FALSE
  ]
  group_index <- match(species_rows$sp_name, groups$Sp_name)
  membership <- data.frame(
    sp_name = species_rows$sp_name,
    tkey = species_rows$tkey,
    species_stem = species_rows$species_stem,
    organism_type = groups$Organism_Type[group_index],
    stringsAsFactors = FALSE
  )
  membership$group_slug <- gain_overlap_clean_slug(
    membership$organism_type
  )

  if (any(!nzchar(membership$group_slug))) {
    stop(
      "At least one Organism_Type cannot be converted to a path-safe slug.",
      call. = FALSE
    )
  }
  group_map <- unique(
    membership[c("organism_type", "group_slug")]
  )
  if (anyDuplicated(group_map$group_slug)) {
    stop(
      "Path-safe Organism_Type slugs are not unique.",
      call. = FALSE
    )
  }
  if (nrow(group_map) != 14L) {
    stop(
      "Expected 14 Organism_Type groups; found ",
      nrow(group_map),
      ".",
      call. = FALSE
    )
  }
  rownames(membership) <- NULL
  membership
}

gain_overlap_group_summary <- function(membership) {
  counts <- table(membership$organism_type)
  summary <- data.frame(
    organism_type = names(counts),
    species_count = as.integer(counts),
    stringsAsFactors = FALSE
  )
  slug_lookup <- unique(
    membership[c("organism_type", "group_slug")]
  )
  summary$group_slug <- slug_lookup$group_slug[
    match(summary$organism_type, slug_lookup$organism_type)
  ]
  summary <- summary[
    order(summary$group_slug),
    c("organism_type", "group_slug", "species_count"),
    drop = FALSE
  ]
  rownames(summary) <- NULL
  summary
}

gain_overlap_build_jobs <- function(manifest, membership, out_dir) {
  gain_overlap_assert_path(out_dir, "out_dir")
  group_summary <- gain_overlap_group_summary(membership)
  rows <- vector(
    "list",
    nrow(gain_overlap_expected_combinations) *
      (nrow(group_summary) + 1L)
  )
  row_index <- 0L

  for (combo_index in seq_len(nrow(gain_overlap_expected_combinations))) {
    period <- gain_overlap_expected_combinations$future_period[[combo_index]]
    scenario <- gain_overlap_expected_combinations$future_scenario[[combo_index]]
    combo_id <- paste(period, scenario, sep = "__")

    row_index <- row_index + 1L
    all_filename <- paste0(
      "gain_count_all_species_hist_to_",
      period,
      "_",
      scenario,
      "_n",
      nrow(membership),
      ".tif"
    )
    rows[[row_index]] <- data.frame(
      job_id = paste("all_species", combo_id, sep = "__"),
      combo_id = combo_id,
      scope = "all_species",
      organism_type = "All species",
      group_slug = "all_species",
      species_count = nrow(membership),
      future_period = period,
      future_scenario = scenario,
      input_count = nrow(membership),
      output_filename = all_filename,
      output_path = file.path(
        out_dir,
        "all_species",
        all_filename
      ),
      status = "pending",
      message = "",
      started_at = "",
      finished_at = "",
      output_size_bytes = NA_real_,
      valid_cell_count = NA_real_,
      nonzero_cell_count = NA_real_,
      minimum = NA_real_,
      maximum = NA_real_,
      stringsAsFactors = FALSE,
      check.names = FALSE
    )

    for (group_index in seq_len(nrow(group_summary))) {
      group <- group_summary[group_index, , drop = FALSE]
      row_index <- row_index + 1L
      filename <- paste0(
        "gain_count_",
        group$group_slug,
        "_hist_to_",
        period,
        "_",
        scenario,
        "_n",
        group$species_count,
        ".tif"
      )
      rows[[row_index]] <- data.frame(
        job_id = paste(group$group_slug, combo_id, sep = "__"),
        combo_id = combo_id,
        scope = "organism_group",
        organism_type = group$organism_type,
        group_slug = group$group_slug,
        species_count = group$species_count,
        future_period = period,
        future_scenario = scenario,
        input_count = group$species_count,
        output_filename = filename,
        output_path = file.path(
          out_dir,
          "by_organism_type",
          group$group_slug,
          filename
        ),
        status = "pending",
        message = "",
        started_at = "",
        finished_at = "",
        output_size_bytes = NA_real_,
        valid_cell_count = NA_real_,
        nonzero_cell_count = NA_real_,
        minimum = NA_real_,
        maximum = NA_real_,
        stringsAsFactors = FALSE,
        check.names = FALSE
      )
    }
  }

  jobs <- do.call(rbind, rows)
  rownames(jobs) <- NULL
  if (nrow(jobs) != 90L) {
    stop(
      "Expected 90 aggregation jobs; built ",
      nrow(jobs),
      ".",
      call. = FALSE
    )
  }
  if (anyDuplicated(jobs$job_id) ||
      anyDuplicated(tolower(jobs$output_path))) {
    stop(
      "Aggregation job IDs or output paths are not unique.",
      call. = FALSE
    )
  }
  jobs
}

gain_overlap_extract_levels <- function(raster, path) {
  levels_list <- terra::levels(raster)
  if (length(levels_list) != 1L ||
      is.null(levels_list[[1L]]) ||
      nrow(levels_list[[1L]]) == 0L) {
    stop(
      "Dynamics raster does not contain a categorical levels table: ",
      path,
      call. = FALSE
    )
  }
  levels_table <- as.data.frame(levels_list[[1L]])
  value_candidates <- names(levels_table)[
    tolower(names(levels_table)) %in% c("value", "id")
  ]
  label_candidates <- names(levels_table)[
    tolower(names(levels_table)) == "dynamics"
  ]
  if (length(value_candidates) == 0L ||
      length(label_candidates) == 0L) {
    stop(
      "Dynamics raster levels table lacks value/dynamics columns: ",
      path,
      call. = FALSE
    )
  }
  data.frame(
    value = as.integer(levels_table[[value_candidates[[1L]]]]),
    dynamics = as.character(levels_table[[label_candidates[[1L]]]]),
    stringsAsFactors = FALSE
  )
}

gain_overlap_validate_input_rasters <- function(manifest,
                                                show_progress = TRUE) {
  gain_overlap_assert_scalar_logical(show_progress, "show_progress")
  paths <- manifest$source_path
  template <- terra::rast(paths[[1L]])
  if (terra::nlyr(template) != 1L) {
    stop(
      "Dynamics raster must contain exactly one layer: ",
      paths[[1L]],
      call. = FALSE
    )
  }

  for (i in seq_along(paths)) {
    raster <- tryCatch(
      terra::rast(paths[[i]]),
      error = function(e) {
        stop(
          "Cannot open dynamics raster ",
          paths[[i]],
          ": ",
          conditionMessage(e),
          call. = FALSE
        )
      }
    )
    if (terra::nlyr(raster) != 1L) {
      stop(
        "Dynamics raster must contain exactly one layer: ",
        paths[[i]],
        call. = FALSE
      )
    }
    if (!isTRUE(terra::compareGeom(
      raster,
      template,
      lyrs = FALSE,
      stopOnError = FALSE
    ))) {
      stop(
        "Dynamics raster geometry differs from the common template: ",
        paths[[i]],
        call. = FALSE
      )
    }

    levels_table <- gain_overlap_extract_levels(raster, paths[[i]])
    if (!identical(
      levels_table,
      gain_overlap_expected_levels
    )) {
      stop(
        "Dynamics raster categorical levels differ from the expected ",
        "1-4 transition codebook: ",
        paths[[i]],
        call. = FALSE
      )
    }
    stored_range <- terra::minmax(raster)
    stored_min <- as.numeric(stored_range[1L, 1L])
    stored_max <- as.numeric(stored_range[2L, 1L])
    if (!is.finite(stored_min) ||
        !is.finite(stored_max) ||
        stored_min < 1 ||
        stored_max > 4) {
      stop(
        "Dynamics raster stored min/max falls outside class codes 1-4: ",
        paths[[i]],
        call. = FALSE
      )
    }

    if (isTRUE(show_progress) &&
        (i %% 100L == 0L || i == length(paths))) {
      message(
        "Validated ",
        i,
        "/",
        length(paths),
        " source rasters."
      )
    }
  }
  template
}

gain_overlap_aggregate_stack <- function(raster_stack,
                                         group_index,
                                         work_dir = NULL) {
  if (!inherits(raster_stack, "SpatRaster")) {
    stop("`raster_stack` must be a terra SpatRaster.", call. = FALSE)
  }
  if (terra::nlyr(raster_stack) != length(group_index)) {
    stop(
      "The raster layer count does not match `group_index`.",
      call. = FALSE
    )
  }
  group_index <- as.character(group_index)
  if (any(is.na(group_index)) || any(!nzchar(group_index))) {
    stop("`group_index` contains empty values.", call. = FALSE)
  }

  group_levels <- sort(unique(group_index))
  group_factor <- factor(
    group_index,
    levels = group_levels
  )
  cleanup_paths <- character()

  if (!is.null(work_dir)) {
    gain_overlap_assert_path(work_dir, "work_dir")
    dir.create(work_dir, recursive = TRUE, showWarnings = FALSE)
    gain_path <- tempfile(
      pattern = ".group_gain_counts_",
      tmpdir = work_dir,
      fileext = ".tif"
    )
    valid_path <- tempfile(
      pattern = ".group_valid_counts_",
      tmpdir = work_dir,
      fileext = ".tif"
    )
    cleanup_paths <- c(
      gain_path,
      paste0(gain_path, ".aux.xml"),
      valid_path,
      paste0(valid_path, ".aux.xml")
    )
  } else {
    gain_path <- ""
    valid_path <- ""
  }

  gain_indicator <- raster_stack == 2L
  valid_indicator <- !is.na(raster_stack)
  group_gain <- terra::tapp(
    gain_indicator,
    index = group_factor,
    fun = "sum",
    na.rm = TRUE,
    filename = gain_path,
    overwrite = nzchar(gain_path),
    wopt = list(
      datatype = "INT1U",
      NAflag = 255,
      gdal = c(
        "COMPRESS=LZW",
        "TILED=YES",
        "BIGTIFF=IF_SAFER"
      )
    )
  )
  group_valid <- terra::tapp(
    valid_indicator,
    index = group_factor,
    fun = "sum",
    na.rm = TRUE,
    filename = valid_path,
    overwrite = nzchar(valid_path),
    wopt = list(
      datatype = "INT1U",
      NAflag = 255,
      gdal = c(
        "COMPRESS=LZW",
        "TILED=YES",
        "BIGTIFF=IF_SAFER"
      )
    )
  )
  names(group_gain) <- group_levels
  names(group_valid) <- group_levels

  group_counts <- terra::ifel(
    group_valid > 0,
    group_gain,
    NA
  )
  names(group_counts) <- group_levels
  all_gain <- terra::app(
    group_gain,
    fun = "sum",
    na.rm = TRUE
  )
  all_valid <- terra::app(
    group_valid,
    fun = "sum",
    na.rm = TRUE
  )
  all_count <- terra::ifel(
    all_valid > 0,
    all_gain,
    NA
  )
  names(all_count) <- "species_gain_count"

  list(
    group_counts = group_counts,
    group_valid_counts = group_valid,
    all_count = all_count,
    all_valid_count = all_valid,
    cleanup_paths = cleanup_paths
  )
}

gain_overlap_cleanup_products <- function(products) {
  paths <- products$cleanup_paths
  rm(products)
  invisible(gc(verbose = FALSE))
  paths <- paths[file.exists(paths)]
  if (length(paths) > 0L) {
    unlink(paths, force = TRUE)
  }
  invisible(TRUE)
}

gain_overlap_atomic_replace <- function(target, writer) {
  dir.create(dirname(target), recursive = TRUE, showWarnings = FALSE)
  extension <- tools::file_ext(target)
  extension <- if (nzchar(extension)) paste0(".", extension) else ""
  stage <- tempfile(
    pattern = ".gain_overlap_metadata_",
    tmpdir = dirname(target),
    fileext = extension
  )
  backup <- NULL
  promoted <- FALSE
  on.exit({
    if (file.exists(stage)) {
      unlink(stage, force = TRUE)
    }
    if (!is.null(backup) && file.exists(backup)) {
      if (!promoted && !file.exists(target)) {
        file.rename(backup, target)
      } else if (promoted) {
        unlink(backup, force = TRUE)
      }
    }
  }, add = TRUE)

  writer(stage)
  if (!file.exists(stage) || is.na(file.info(stage)$size) ||
      file.info(stage)$size <= 0) {
    stop(
      "Metadata staging write failed for ",
      target,
      ".",
      call. = FALSE
    )
  }

  if (file.exists(target)) {
    backup <- tempfile(
      pattern = ".gain_overlap_backup_",
      tmpdir = dirname(target),
      fileext = extension
    )
    if (!isTRUE(file.rename(target, backup))) {
      stop(
        "Could not back up existing metadata file: ",
        target,
        call. = FALSE
      )
    }
  }
  if (!isTRUE(file.rename(stage, target))) {
    stop(
      "Could not promote metadata file to ",
      target,
      ".",
      call. = FALSE
    )
  }
  promoted <- TRUE
  invisible(target)
}

gain_overlap_write_manifests <- function(jobs, metadata_dir) {
  rds_path <- file.path(
    metadata_dir,
    "sp_hist_fut_dyn_spr_gains_manifest.rds"
  )
  csv_path <- file.path(
    metadata_dir,
    "sp_hist_fut_dyn_spr_gains_manifest.csv"
  )
  gain_overlap_atomic_replace(
    rds_path,
    function(path) saveRDS(jobs, path)
  )
  gain_overlap_atomic_replace(
    csv_path,
    function(path) {
      utils::write.csv(
        jobs,
        path,
        row.names = FALSE,
        na = "",
        fileEncoding = "UTF-8"
      )
    }
  )
  invisible(list(rds = rds_path, csv = csv_path))
}

gain_overlap_write_membership <- function(membership, metadata_dir) {
  path <- file.path(metadata_dir, "species_group_membership.csv")
  gain_overlap_atomic_replace(
    path,
    function(staged) {
      utils::write.csv(
        membership,
        staged,
        row.names = FALSE,
        na = "",
        fileEncoding = "UTF-8"
      )
    }
  )
  invisible(path)
}

gain_overlap_readme_lines <- function(manifest,
                                      membership,
                                      template,
                                      dynamics_dir,
                                      groups_path) {
  group_summary <- gain_overlap_group_summary(membership)
  group_lines <- paste0(
    "- ",
    group_summary$organism_type,
    " (`",
    group_summary$group_slug,
    "`): ",
    group_summary$species_count,
    " species"
  )
  resolution <- terra::res(template)
  extent <- as.vector(terra::ext(template))
  c(
    "# Species gain-overlap rasters",
    "",
    "## Contents",
    "",
    paste0(
      "This folder contains ",
      nrow(manifest),
      " single-band GeoTIFF maps: six for all 120 species and ",
      "six for each of 14 organism groups."
    ),
    "",
    "## Method",
    "",
    paste0(
      "For each future period/scenario, each input dynamics raster was ",
      "converted to a gain indicator using `dynamics == 2`, where class 2 ",
      "is `Gain (0->1)`. Indicators were summed by pixel."
    ),
    "",
    paste0(
      "Partial NoData is ignored when counting gains. A result cell is ",
      "NoData only when every species raster contributing to that ",
      "aggregation is NoData."
    ),
    "",
    paste0(
      "All-species results are derived from the 14 mutually exclusive ",
      "group totals and are checked for exact cellwise reconciliation."
    ),
    "",
    "No intermediate Boolean rasters are retained.",
    "",
    "## Inputs",
    "",
    paste0("- Dynamics directory: `", dynamics_dir, "`"),
    paste0("- Species-group workbook: `", groups_path, "`"),
    paste0(
      "- Species join: exact match between dynamics-manifest `sp_name` ",
      "and workbook `Sp_name`"
    ),
    "",
    "## Raster properties",
    "",
    "- Band name: `species_gain_count`",
    "- Data type: unsigned 8-bit integer (`INT1U`), NoData 255",
    "- Compression: LZW, tiled GeoTIFF",
    paste0(
      "- Dimensions: ",
      terra::nrow(template),
      " rows x ",
      terra::ncol(template),
      " columns"
    ),
    paste0(
      "- Resolution: ",
      format(resolution[[1L]], scientific = FALSE),
      " x ",
      format(resolution[[2L]], scientific = FALSE),
      " map units"
    ),
    paste0(
      "- Extent (xmin, xmax, ymin, ymax): ",
      paste(format(extent, scientific = FALSE), collapse = ", ")
    ),
    paste0("- CRS: `", terra::crs(template, proj = TRUE), "`"),
    "",
    "## Naming",
    "",
    paste0(
      "- All species: ",
      "`all_species/gain_count_all_species_hist_to_",
      "{period}_{scenario}_n120.tif`"
    ),
    paste0(
      "- Organism groups: ",
      "`by_organism_type/{group_slug}/gain_count_{group_slug}_",
      "hist_to_{period}_{scenario}_n{N}.tif`"
    ),
    "",
    "The `n{N}` suffix records the maximum species denominator.",
    "",
    "## Organism groups",
    "",
    group_lines,
    "",
    "## Metadata",
    "",
    paste0(
      "`metadata/sp_hist_fut_dyn_spr_gains_manifest.csv` records one row ",
      "per output, including processing status, species denominator, ",
      "file size, valid-cell count, nonzero-cell count, and observed range."
    ),
    "",
    paste0(
      "`metadata/species_group_membership.csv` records the exact species ",
      "membership used for each group."
    ),
    "",
    paste0("Generated: ", gain_overlap_timestamp())
  )
}

gain_overlap_write_readme <- function(manifest,
                                      membership,
                                      template,
                                      dynamics_dir,
                                      groups_path,
                                      out_dir) {
  path <- file.path(out_dir, "metadata", "README.md")
  lines <- gain_overlap_readme_lines(
    manifest,
    membership,
    template,
    dynamics_dir,
    groups_path
  )
  gain_overlap_atomic_replace(
    path,
    function(staged) {
      writeLines(
        enc2utf8(lines),
        con = staged,
        useBytes = TRUE
      )
    }
  )
  invisible(path)
}

gain_overlap_verify_written_raster <- function(path,
                                               template,
                                               max_count) {
  if (!file.exists(path)) {
    stop("Output raster does not exist: ", path, call. = FALSE)
  }
  if (is.na(file.info(path)$size) || file.info(path)$size <= 0) {
    stop("Output raster is empty: ", path, call. = FALSE)
  }
  raster <- terra::rast(path)
  if (terra::nlyr(raster) != 1L ||
      !isTRUE(terra::compareGeom(
        raster,
        template,
        lyrs = FALSE,
        stopOnError = FALSE
      ))) {
    stop(
      "Output raster geometry or layer count is invalid: ",
      path,
      call. = FALSE
    )
  }
  if (!identical(terra::datatype(raster), "INT1U")) {
    stop(
      "Output raster datatype is not INT1U: ",
      path,
      call. = FALSE
    )
  }
  if (!identical(names(raster), "species_gain_count")) {
    stop(
      "Output raster band name is not species_gain_count: ",
      path,
      call. = FALSE
    )
  }
  description <- terra::describe(path)
  if (!any(grepl("COMPRESSION=LZW", description, fixed = TRUE))) {
    stop(
      "Output raster does not report LZW compression: ",
      path,
      call. = FALSE
    )
  }

  indicators <- c(!is.na(raster), raster > 0)
  indicator_sums <- terra::global(
    indicators,
    fun = "sum",
    na.rm = TRUE
  )
  range_stats <- terra::global(
    raster,
    fun = c("min", "max"),
    na.rm = TRUE
  )
  valid_cell_count <- as.numeric(indicator_sums[1L, 1L])
  nonzero_cell_count <- as.numeric(indicator_sums[2L, 1L])
  minimum <- as.numeric(range_stats[1L, "min"])
  maximum <- as.numeric(range_stats[1L, "max"])
  if (!is.finite(valid_cell_count) ||
      valid_cell_count < 1 ||
      !is.finite(nonzero_cell_count) ||
      nonzero_cell_count < 0 ||
      !is.finite(minimum) ||
      !is.finite(maximum) ||
      minimum < 0 ||
      maximum > max_count ||
      minimum != floor(minimum) ||
      maximum != floor(maximum)) {
    stop(
      "Output raster statistics are invalid for a 0-",
      max_count,
      " species count: ",
      path,
      ". Observed valid_cell_count=",
      valid_cell_count,
      ", nonzero_cell_count=",
      nonzero_cell_count,
      ", minimum=",
      minimum,
      ", maximum=",
      maximum,
      call. = FALSE
    )
  }

  list(
    valid_cell_count = valid_cell_count,
    nonzero_cell_count = nonzero_cell_count,
    minimum = minimum,
    maximum = maximum,
    output_size_bytes = as.numeric(file.info(path)$size)
  )
}

gain_overlap_write_raster_atomic <- function(raster,
                                             path,
                                             template,
                                             max_count,
                                             overwrite = FALSE) {
  gain_overlap_assert_scalar_logical(overwrite, "overwrite")
  if (!inherits(raster, "SpatRaster") ||
      terra::nlyr(raster) != 1L) {
    stop(
      "`raster` must be a single-layer terra SpatRaster.",
      call. = FALSE
    )
  }
  if (!isTRUE(terra::compareGeom(
    raster,
    template,
    lyrs = FALSE,
    stopOnError = FALSE
  ))) {
    stop(
      "Raster to write does not match the output template.",
      call. = FALSE
    )
  }

  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  if (file.exists(path) && !isTRUE(overwrite)) {
    stop(
      "Output already exists and overwrite=FALSE: ",
      path,
      call. = FALSE
    )
  }
  stage <- tempfile(
    pattern = ".gain_overlap_raster_",
    tmpdir = dirname(path),
    fileext = ".tif"
  )
  stage_aux <- paste0(stage, ".aux.xml")
  target_aux <- paste0(path, ".aux.xml")
  backup <- NULL
  backup_aux <- NULL
  promoted <- FALSE
  on.exit({
    unlink(c(stage, stage_aux), force = TRUE)
    if (!promoted && !is.null(backup) && file.exists(backup)) {
      if (file.exists(path)) {
        unlink(path, force = TRUE)
      }
      file.rename(backup, path)
      if (!is.null(backup_aux) && file.exists(backup_aux)) {
        if (file.exists(target_aux)) {
          unlink(target_aux, force = TRUE)
        }
        file.rename(backup_aux, target_aux)
      }
    }
    if (promoted) {
      unlink(c(backup, backup_aux), force = TRUE)
    }
  }, add = TRUE)

  names(raster) <- "species_gain_count"
  terra::writeRaster(
    raster,
    filename = stage,
    overwrite = TRUE,
    wopt = list(
      datatype = "INT1U",
      NAflag = 255,
      gdal = c(
        "COMPRESS=LZW",
        "TILED=YES",
        "BIGTIFF=IF_SAFER"
      )
    )
  )
  gain_overlap_verify_written_raster(
    stage,
    template,
    max_count
  )

  if (file.exists(path)) {
    backup <- tempfile(
      pattern = ".gain_overlap_raster_backup_",
      tmpdir = dirname(path),
      fileext = ".tif"
    )
    if (!isTRUE(file.rename(path, backup))) {
      stop(
        "Could not back up existing output raster: ",
        path,
        call. = FALSE
      )
    }
    if (file.exists(target_aux)) {
      backup_aux <- paste0(backup, ".aux.xml")
      if (!isTRUE(file.rename(target_aux, backup_aux))) {
        stop(
          "Could not back up existing raster sidecar: ",
          target_aux,
          call. = FALSE
        )
      }
    }
  }

  if (!isTRUE(file.rename(stage, path))) {
    stop(
      "Could not promote staged raster to ",
      path,
      ".",
      call. = FALSE
    )
  }
  if (file.exists(stage_aux) &&
      !isTRUE(file.rename(stage_aux, target_aux))) {
    stop(
      "Could not promote staged raster sidecar to ",
      target_aux,
      ".",
      call. = FALSE
    )
  }
  promoted <- TRUE
  stats <- gain_overlap_verify_written_raster(
    path,
    template,
    max_count
  )
  invisible(stats)
}

gain_overlap_existing_output_state <- function(path,
                                               template,
                                               max_count,
                                               latest_input_mtime) {
  if (!file.exists(path)) {
    return(list(state = "missing", message = "", stats = NULL))
  }
  validation <- tryCatch(
    gain_overlap_verify_written_raster(
      path,
      template,
      max_count
    ),
    error = function(e) e
  )
  if (inherits(validation, "error")) {
    return(list(
      state = "invalid",
      message = conditionMessage(validation),
      stats = NULL
    ))
  }
  stale <- file.info(path)$mtime < latest_input_mtime
  list(
    state = if (isTRUE(stale)) "stale" else "current",
    message = if (isTRUE(stale)) {
      "Output is older than at least one input."
    } else {
      ""
    },
    stats = validation
  )
}

gain_overlap_apply_stats <- function(jobs, row_index, stats) {
  jobs$output_size_bytes[[row_index]] <- stats$output_size_bytes
  jobs$valid_cell_count[[row_index]] <- stats$valid_cell_count
  jobs$nonzero_cell_count[[row_index]] <- stats$nonzero_cell_count
  jobs$minimum[[row_index]] <- stats$minimum
  jobs$maximum[[row_index]] <- stats$maximum
  jobs
}

gain_overlap_validate_reconciliation <- function(jobs,
                                                 show_progress = TRUE) {
  for (combo_index in seq_len(nrow(gain_overlap_expected_combinations))) {
    combo_id <- paste(
      gain_overlap_expected_combinations$future_period[[combo_index]],
      gain_overlap_expected_combinations$future_scenario[[combo_index]],
      sep = "__"
    )
    combo_jobs <- jobs[jobs$combo_id == combo_id, , drop = FALSE]
    all_path <- combo_jobs$output_path[
      combo_jobs$scope == "all_species"
    ]
    group_paths <- combo_jobs$output_path[
      combo_jobs$scope == "organism_group"
    ]
    all_raster <- terra::rast(all_path)
    group_stack <- terra::rast(group_paths)
    if (!isTRUE(terra::compareGeom(
      group_stack,
      all_raster,
      lyrs = FALSE,
      stopOnError = FALSE
    ))) {
      stop(
        "Group/all output geometry differs for ",
        combo_id,
        ".",
        call. = FALSE
      )
    }
    group_sum <- terra::app(
      group_stack,
      fun = "sum",
      na.rm = TRUE
    )
    valid_groups <- terra::app(
      !is.na(group_stack),
      fun = "sum",
      na.rm = TRUE
    )
    reconstructed <- terra::ifel(
      valid_groups > 0,
      group_sum,
      NA
    )
    mask_differences <- terra::global(
      is.na(reconstructed) != is.na(all_raster),
      fun = "sum",
      na.rm = TRUE
    )
    maximum_difference <- terra::global(
      abs(reconstructed - all_raster),
      fun = "max",
      na.rm = TRUE
    )
    if (as.numeric(mask_differences[1L, 1L]) != 0 ||
        as.numeric(maximum_difference[1L, 1L]) != 0) {
      stop(
        "Organism-group outputs do not reconcile exactly with the ",
        "all-species output for ",
        combo_id,
        ".",
        call. = FALSE
      )
    }
    gain_overlap_message(
      show_progress,
      "Reconciled group and all-species outputs for ",
      combo_id,
      "."
    )
  }
  invisible(TRUE)
}

run_species_gain_aggregations <- function(
    dynamics_dir = file.path(
      "data",
      "post_model_maps",
      "sp_hist_fut_dynamics"
    ),
    groups_path = file.path(
      "data",
      "sel_species-standardised_organism_types_v2.xlsx"
    ),
    out_dir = file.path(
      "data",
      "post_model_maps",
      "sp_hist_fut_dyn_spr_gains"
    ),
    overwrite = FALSE,
    show_progress = TRUE) {
  gain_overlap_assert_scalar_logical(overwrite, "overwrite")
  gain_overlap_assert_scalar_logical(show_progress, "show_progress")
  gain_overlap_assert_path(dynamics_dir, "dynamics_dir")
  gain_overlap_assert_path(groups_path, "groups_path")
  gain_overlap_assert_path(out_dir, "out_dir")
  gain_overlap_require_packages()

  old_terra_options <- terra::terraOptions()
  on.exit({
    terra::terraOptions(progress = old_terra_options$progress)
  }, add = TRUE)
  terra::terraOptions(progress = if (isTRUE(show_progress)) 1 else 0)

  gain_overlap_message(show_progress, "Loading input catalogues.")
  manifest <- gain_overlap_read_manifest(dynamics_dir)
  groups <- gain_overlap_read_groups(groups_path)
  membership <- gain_overlap_build_membership(manifest, groups)
  jobs <- gain_overlap_build_jobs(manifest, membership, out_dir)

  gain_overlap_message(
    show_progress,
    "Preflighting ",
    nrow(manifest),
    " dynamics rasters."
  )
  template <- gain_overlap_validate_input_rasters(
    manifest,
    show_progress = show_progress
  )

  metadata_dir <- file.path(out_dir, "metadata")
  work_dir <- file.path(out_dir, ".working")
  dir.create(metadata_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(work_dir, recursive = TRUE, showWarnings = FALSE)
  on.exit({
    remaining <- list.files(
      work_dir,
      all.files = TRUE,
      full.names = TRUE,
      no.. = TRUE
    )
    if (length(remaining) == 0L && dir.exists(work_dir)) {
      unlink(work_dir, recursive = TRUE, force = TRUE)
    }
  }, add = TRUE)

  membership_path <- gain_overlap_write_membership(
    membership,
    metadata_dir
  )
  manifest_path <- attr(manifest, "manifest_path")
  source_metadata_mtime <- max(
    file.info(c(manifest_path, groups_path))$mtime,
    na.rm = TRUE
  )

  gain_overlap_message(show_progress, "Checking for resumable outputs.")
  for (row_index in seq_len(nrow(jobs))) {
    combo_rows <- manifest[
      manifest$combo_id == jobs$combo_id[[row_index]],
      ,
      drop = FALSE
    ]
    latest_input_mtime <- max(
      c(
        source_metadata_mtime,
        file.info(combo_rows$source_path)$mtime
      ),
      na.rm = TRUE
    )
    state <- gain_overlap_existing_output_state(
      jobs$output_path[[row_index]],
      template,
      jobs$species_count[[row_index]],
      latest_input_mtime
    )
    if (identical(state$state, "current") && !isTRUE(overwrite)) {
      jobs$status[[row_index]] <- "skipped_existing"
      jobs$message[[row_index]] <- "Existing output is valid and current."
      jobs$finished_at[[row_index]] <- gain_overlap_timestamp()
      jobs <- gain_overlap_apply_stats(jobs, row_index, state$stats)
    } else if (state$state %in% c("invalid", "stale") &&
               !isTRUE(overwrite)) {
      stop(
        "Cannot resume because output is ",
        state$state,
        ": ",
        jobs$output_path[[row_index]],
        ". ",
        state$message,
        " Re-run with overwrite=TRUE to regenerate it.",
        call. = FALSE
      )
    }
  }

  gain_overlap_write_readme(
    jobs,
    membership,
    template,
    dynamics_dir,
    groups_path,
    out_dir
  )
  gain_overlap_write_manifests(jobs, metadata_dir)

  for (combo_index in seq_len(nrow(gain_overlap_expected_combinations))) {
    period <- gain_overlap_expected_combinations$future_period[[combo_index]]
    scenario <- gain_overlap_expected_combinations$future_scenario[[combo_index]]
    combo_id <- paste(period, scenario, sep = "__")
    combo_job_indices <- which(jobs$combo_id == combo_id)
    pending_indices <- combo_job_indices[
      jobs$status[combo_job_indices] == "pending" |
        isTRUE(overwrite)
    ]
    if (length(pending_indices) == 0L) {
      gain_overlap_message(
        show_progress,
        "All outputs are current for ",
        combo_id,
        "; skipping aggregation."
      )
      next
    }

    gain_overlap_message(
      show_progress,
      "Aggregating ",
      combo_id,
      " (",
      combo_index,
      "/",
      nrow(gain_overlap_expected_combinations),
      ")."
    )
    combo_manifest <- manifest[
      manifest$combo_id == combo_id,
      ,
      drop = FALSE
    ]
    combo_manifest <- combo_manifest[
      order(combo_manifest$sp_name),
      ,
      drop = FALSE
    ]
    membership_index <- match(
      combo_manifest$sp_name,
      membership$sp_name
    )
    raster_stack <- terra::rast(combo_manifest$source_path)
    products <- NULL

    combo_error <- tryCatch(
      {
        products <- gain_overlap_aggregate_stack(
          raster_stack,
          membership$group_slug[membership_index],
          work_dir = work_dir
        )
        for (job_index in pending_indices) {
          jobs$status[[job_index]] <- "processing"
          jobs$started_at[[job_index]] <- gain_overlap_timestamp()
          gain_overlap_write_manifests(jobs, metadata_dir)

          if (jobs$scope[[job_index]] == "all_species") {
            output <- products$all_count
          } else {
            output <- products$group_counts[[
              jobs$group_slug[[job_index]]
            ]]
          }
          names(output) <- "species_gain_count"
          stats <- gain_overlap_write_raster_atomic(
            output,
            jobs$output_path[[job_index]],
            template,
            jobs$species_count[[job_index]],
            overwrite = isTRUE(overwrite)
          )
          jobs$status[[job_index]] <- "written"
          jobs$message[[job_index]] <- "Output written and verified."
          jobs$finished_at[[job_index]] <- gain_overlap_timestamp()
          jobs <- gain_overlap_apply_stats(jobs, job_index, stats)
          gain_overlap_write_manifests(jobs, metadata_dir)
          gain_overlap_message(
            show_progress,
            "  Wrote ",
            jobs$output_filename[[job_index]],
            "."
          )
        }
        NULL
      },
      error = function(e) e
    )

    if (!is.null(products)) {
      cleanup_paths <- products$cleanup_paths
      rm(products)
      invisible(gc(verbose = FALSE))
      unlink(cleanup_paths[file.exists(cleanup_paths)], force = TRUE)
    }
    if (inherits(combo_error, "error")) {
      failed_indices <- combo_job_indices[
        jobs$status[combo_job_indices] %in% c("pending", "processing")
      ]
      jobs$status[failed_indices] <- "failed"
      jobs$message[failed_indices] <- conditionMessage(combo_error)
      jobs$finished_at[failed_indices] <- gain_overlap_timestamp()
      gain_overlap_write_manifests(jobs, metadata_dir)
      stop(conditionMessage(combo_error), call. = FALSE)
    }
  }

  unfinished <- which(
    !jobs$status %in% c("written", "skipped_existing")
  )
  if (length(unfinished) > 0L) {
    stop(
      length(unfinished),
      " aggregation job(s) did not complete. See the output manifest.",
      call. = FALSE
    )
  }

  gain_overlap_message(
    show_progress,
    "Validating exact group-to-all reconciliation."
  )
  gain_overlap_validate_reconciliation(
    jobs,
    show_progress = show_progress
  )
  gain_overlap_write_manifests(jobs, metadata_dir)
  gain_overlap_write_readme(
    jobs,
    membership,
    template,
    dynamics_dir,
    groups_path,
    out_dir
  )

  staging_files <- list.files(
    out_dir,
    pattern = "^[.](gain_overlap|group_gain|group_valid)",
    all.files = TRUE,
    recursive = TRUE,
    full.names = TRUE
  )
  if (length(staging_files) > 0L) {
    stop(
      "One or more staging files remain after completion:\n- ",
      paste(staging_files, collapse = "\n- "),
      call. = FALSE
    )
  }
  if (file.exists(membership_path) &&
      nrow(jobs) == 90L) {
    gain_overlap_message(
      show_progress,
      "Completed ",
      nrow(jobs),
      " species gain-overlap rasters."
    )
  }
  invisible(jobs)
}

gain_overlap_run_on_source <- isTRUE(getOption(
  "wisdm.run_species_gain_aggregations_on_source",
  FALSE
))

if (sys.nframe() == 0L || gain_overlap_run_on_source) {
  run_species_gain_aggregations()
}
