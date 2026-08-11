# Species-specific historical-to-future binary suitability dynamics.
#
# From an interactive R or RStudio console, source the script to run the
# complete 120-species batch:
# source("./src/post_modelling_maps/dyn_maps_by_sp_scn-v2.R")
#
# To source only the function definitions without starting the batch:
# options(wisdm.run_species_dynamics_on_source = FALSE)
# source("./src/post_modelling_maps/dyn_maps_by_sp_scn-v2.R")
#
# After loading functions only, run a subset with:
# run_species_dynamics(
#   species = "Acacia mearnsii De Wild.",
#   overwrite = FALSE
# )
#
# The default fast preflight checks raster metadata and stored binary ranges.
# Use input_validation = "full" only when an exhaustive cell-value audit is
# required.
#
# Direct execution with Rscript remains supported:
# Rscript src/post_modelling_maps/dyn_maps_by_sp_scn-v2.R

dynamics_catalog_path <- file.path(
  "data",
  "processed",
  "best_mod_set_paths.rds"
)

if (file.exists(dynamics_catalog_path)) {
  best_mod_set_paths <- readRDS(dynamics_catalog_path)
} else {
  source(file.path(
    "src",
    "post_modelling_analyses",
    "check_best_model-v1.R"
  ))
}

source(file.path(
  "src",
  "post_modelling_maps",
  "rast_files_crawler-v1.R"
))
source(file.path("src", "helper_functions.R"))

if (!requireNamespace("terra", quietly = TRUE)) {
  stop("Package 'terra' is required to calculate dynamics maps.", call. = FALSE)
}
if (!requireNamespace("cli", quietly = TRUE)) {
  stop("Package 'cli' is required to report batch progress.", call. = FALSE)
}

dyn_levels <- data.frame(
  id = 1:4,
  dynamics = c(
    "Stable unsuitable (0->0)",
    "Gain (0->1)",
    "Loss (1->0)",
    "Stable suitable (1->1)"
  ),
  stringsAsFactors = FALSE
)

species_dynamics_future_scenarios <- data.frame(
  period = rep(c("mid_2041_2070", "late_2071_2100"), each = 3L),
  scenario = rep(c("ssp126", "ssp370", "ssp585"), times = 2L),
  stringsAsFactors = FALSE
)

species_dynamics_timestamp <- function() {
  format(Sys.time(), tz = "UTC", usetz = TRUE)
}

species_dynamics_progress_bar <- function(name, total,
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

species_dynamics_assert_scalar_logical <- function(value, argument) {
  if (!is.logical(value) || length(value) != 1L || is.na(value)) {
    stop("`", argument, "` must be TRUE or FALSE.", call. = FALSE)
  }
  invisible(value)
}

species_dynamics_assert_output_dir <- function(out_dir) {
  if (!is.character(out_dir) ||
      length(out_dir) != 1L ||
      is.na(out_dir) ||
      !nzchar(trimws(out_dir))) {
    stop("`out_dir` must be one non-empty path.", call. = FALSE)
  }
  invisible(out_dir)
}

species_dynamics_make_stem <- function(sp_name) {
  suppressMessages(species_output_stem(sp_name))
}

build_species_dynamics_jobs <- function(
    species = "all",
    out_dir = file.path(
      "data",
      "post_model_maps",
      "sp_hist_fut_dynamics"
    ),
    fuzzy = FALSE,
    max_fuzzy_distance = 0.10,
    catalog = best_mod_set_paths) {
  species_dynamics_assert_output_dir(out_dir)

  historical <- crawl_raster_paths(
    species = species,
    map_type = "binary",
    period = "hist",
    scenario = "hist",
    fuzzy = fuzzy,
    max_fuzzy_distance = max_fuzzy_distance,
    catalog = catalog
  )

  if (anyDuplicated(as.character(historical$tkey))) {
    stop(
      "Historical crawler results contain duplicated taxon keys.",
      call. = FALSE
    )
  }

  future_by_scenario <- lapply(
    seq_len(nrow(species_dynamics_future_scenarios)),
    function(i) {
      period <- species_dynamics_future_scenarios$period[[i]]
      scenario <- species_dynamics_future_scenarios$scenario[[i]]
      result <- crawl_raster_paths(
        species = species,
        map_type = "binary",
        period = period,
        scenario = scenario,
        fuzzy = fuzzy,
        max_fuzzy_distance = max_fuzzy_distance,
        catalog = catalog
      )
      result$future_period <- period
      result$future_scenario <- scenario
      result
    }
  )

  species_stems <- vapply(
    historical$sp_name,
    species_dynamics_make_stem,
    character(1)
  )
  taxon_keys <- as.character(historical$tkey)

  rows <- vector(
    "list",
    nrow(historical) * nrow(species_dynamics_future_scenarios)
  )
  row_index <- 0L

  for (species_index in seq_len(nrow(historical))) {
    for (scenario_index in seq_along(future_by_scenario)) {
      future <- future_by_scenario[[scenario_index]]
      future_index <- match(
        taxon_keys[[species_index]],
        as.character(future$tkey)
      )
      if (is.na(future_index)) {
        stop(
          "Future crawler results are missing taxon key ",
          taxon_keys[[species_index]],
          " for ",
          future$future_period[[1L]],
          "/",
          future$future_scenario[[1L]],
          ".",
          call. = FALSE
        )
      }

      period <- future$future_period[[future_index]]
      scenario <- future$future_scenario[[future_index]]
      output_filename <- paste0(
        "dyn_",
        species_stems[[species_index]],
        "_",
        taxon_keys[[species_index]],
        "_hist_to_",
        period,
        "_",
        scenario,
        "_binary10pct.tif"
      )
      output_path <- fortify_output_path(file.path(
        out_dir,
        output_filename
      ))

      row_index <- row_index + 1L
      rows[[row_index]] <- data.frame(
        job_id = paste(
          taxon_keys[[species_index]],
          period,
          scenario,
          sep = "__"
        ),
        sp_name = historical$sp_name[[species_index]],
        tkey = historical$tkey[[species_index]],
        species_stem = species_stems[[species_index]],
        map_type = "binary",
        historical_period = "hist",
        historical_scenario = "hist",
        future_period = period,
        future_scenario = scenario,
        historical_path = historical$path[[species_index]],
        future_path = future$path[[future_index]],
        output_filename = basename(output_path),
        output_path = output_path,
        status = "pending",
        message = "",
        started_at = "",
        finished_at = "",
        output_size_bytes = NA_real_,
        stringsAsFactors = FALSE,
        check.names = FALSE
      )
    }
  }

  jobs <- do.call(rbind, rows)
  rownames(jobs) <- NULL

  if (anyDuplicated(jobs$job_id)) {
    stop("Dynamics job IDs are not unique.", call. = FALSE)
  }
  normalized_outputs <- tolower(normalizePath(
    jobs$output_path,
    winslash = "/",
    mustWork = FALSE
  ))
  if (anyDuplicated(normalized_outputs)) {
    stop("Dynamics output paths are not unique.", call. = FALSE)
  }

  expected_jobs <- nrow(historical) *
    nrow(species_dynamics_future_scenarios)
  if (nrow(jobs) != expected_jobs) {
    stop(
      "Expected ",
      expected_jobs,
      " dynamics jobs but built ",
      nrow(jobs),
      ".",
      call. = FALSE
    )
  }

  jobs
}

species_dynamics_binary_summary <- function(raster, label = "Raster") {
  if (!inherits(raster, "SpatRaster")) {
    stop(label, " must be a terra SpatRaster.", call. = FALSE)
  }
  if (terra::nlyr(raster) != 1L) {
    stop(label, " must contain exactly one raster layer.", call. = FALSE)
  }

  valid_cells <- !is.na(raster)
  invalid_cells <- (raster != 0) & (raster != 1)
  validation_layers <- c(valid_cells, invalid_cells)
  names(validation_layers) <- c("valid", "invalid")
  counts <- terra::global(
    validation_layers,
    fun = "sum",
    na.rm = TRUE
  )
  valid_count <- as.numeric(counts["valid", 1L])
  invalid_count <- as.numeric(counts["invalid", 1L])

  if (!is.finite(valid_count) || valid_count < 1) {
    stop(label, " contains no non-missing cells.", call. = FALSE)
  }
  if (!is.finite(invalid_count) || invalid_count > 0) {
    stop(
      label,
      " contains ",
      if (is.finite(invalid_count)) invalid_count else "unknown",
      " non-binary cell(s); only 0, 1, and NA are allowed.",
      call. = FALSE
    )
  }

  invisible(list(
    valid_cells = valid_count,
    invalid_cells = invalid_count
  ))
}

species_dynamics_validation_mode <- function(validation_mode) {
  match.arg(validation_mode, c("metadata", "full"))
}

species_dynamics_validate_input_file <- function(
    path,
    validation_mode = c("metadata", "full")) {
  validation_mode <- species_dynamics_validation_mode(validation_mode)
  if (!is.character(path) ||
      length(path) != 1L ||
      is.na(path) ||
      !nzchar(path)) {
    stop("Input raster path must be one non-empty value.", call. = FALSE)
  }
  if (!file.exists(path)) {
    stop("Input raster does not exist: ", path, call. = FALSE)
  }

  raster <- tryCatch(
    terra::rast(path),
    error = function(e) {
      stop(
        "Input raster cannot be opened: ",
        path,
        ". ",
        conditionMessage(e),
        call. = FALSE
      )
    }
  )
  label <- paste0("Input raster '", path, "'")
  if (terra::nlyr(raster) != 1L) {
    stop(label, " must contain exactly one raster layer.", call. = FALSE)
  }
  if (terra::nrow(raster) < 1L ||
      terra::ncol(raster) < 1L ||
      terra::ncell(raster) < 1) {
    stop(label, " has invalid or empty dimensions.", call. = FALSE)
  }
  raster_crs <- terra::crs(raster, proj = TRUE)
  if (is.na(raster_crs) || !nzchar(raster_crs)) {
    stop(label, " does not define a coordinate reference system.", call. = FALSE)
  }
  raster_resolution <- terra::res(raster)
  if (any(!is.finite(raster_resolution)) ||
      any(raster_resolution <= 0)) {
    stop(label, " has an invalid raster resolution.", call. = FALSE)
  }

  if (isTRUE(terra::hasMinMax(raster))) {
    stored_range <- as.numeric(terra::minmax(raster))
    if (length(stored_range) != 2L ||
        any(!is.finite(stored_range)) ||
        any(!stored_range %in% c(0, 1))) {
      stop(
        label,
        " has stored min/max values incompatible with a binary raster: ",
        paste(stored_range, collapse = ", "),
        ".",
        call. = FALSE
      )
    }
  }

  if (identical(validation_mode, "full")) {
    species_dynamics_binary_summary(raster, label)
  }
  invisible(list(
    validation_mode = validation_mode,
    stored_minmax_checked = isTRUE(terra::hasMinMax(raster))
  ))
}

species_dynamics_normalize_path <- function(path) {
  normalized <- normalizePath(
    raster_paths_mark_utf8(path),
    winslash = "/",
    mustWork = FALSE
  )
  if (.Platform$OS.type == "windows") {
    normalized <- tolower(normalized)
  }
  normalized
}

species_dynamics_input_metadata <- function(paths) {
  paths <- unique(raster_paths_mark_utf8(paths))
  info <- file.info(paths)
  data.frame(
    path = paths,
    normalized_path = species_dynamics_normalize_path(paths),
    size_bytes = as.numeric(info$size),
    mtime_numeric = as.numeric(info$mtime),
    mtime_utc = ifelse(
      is.na(info$mtime),
      "",
      format(info$mtime, tz = "UTC", usetz = TRUE)
    ),
    validated_at = "",
    validation_method = "",
    validation_mode = "",
    validated = FALSE,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
}

species_dynamics_input_cache_paths <- function(cache_dir) {
  list(
    rds = file.path(
      cache_dir,
      "sp_hist_fut_dynamics_input_validation.rds"
    ),
    csv = file.path(
      cache_dir,
      "sp_hist_fut_dynamics_input_validation.csv"
    )
  )
}

load_species_dynamics_input_cache <- function(cache_dir) {
  if (is.null(cache_dir)) {
    return(NULL)
  }
  paths <- species_dynamics_input_cache_paths(cache_dir)
  if (!file.exists(paths$rds)) {
    return(NULL)
  }

  cache <- tryCatch(
    readRDS(paths$rds),
    error = function(e) NULL
  )
  required_columns <- c(
    "path",
    "normalized_path",
    "size_bytes",
    "mtime_numeric",
    "mtime_utc",
    "validated_at",
    "validation_method"
  )
  if (!is.data.frame(cache) ||
      length(setdiff(required_columns, names(cache))) > 0L) {
    cli::cli_alert_warning(
      "Ignoring an unreadable or incompatible binary-input validation cache."
    )
    return(NULL)
  }

  if (!"validation_mode" %in% names(cache)) {
    cache$validation_mode <- ifelse(
      cache$validation_method == "binary_scan",
      "full",
      "metadata"
    )
  }
  cache[, c(required_columns, "validation_mode"), drop = FALSE]
}

write_species_dynamics_input_cache <- function(records, cache_dir) {
  if (is.null(cache_dir)) {
    return(invisible(NULL))
  }
  if (!dir.exists(cache_dir)) {
    stop(
      "Input-validation cache directory does not exist: ",
      cache_dir,
      call. = FALSE
    )
  }

  cache <- records[
    records$validated,
    c(
      "path",
      "normalized_path",
      "size_bytes",
      "mtime_numeric",
      "mtime_utc",
      "validated_at",
      "validation_method",
      "validation_mode"
    ),
    drop = FALSE
  ]
  rownames(cache) <- NULL
  paths <- species_dynamics_input_cache_paths(cache_dir)
  species_dynamics_write_atomic(
    paths$rds,
    function(path) saveRDS(cache, path)
  )
  species_dynamics_write_atomic(
    paths$csv,
    function(path) utils::write.csv(
      cache,
      file = path,
      row.names = FALSE,
      na = "",
      fileEncoding = "UTF-8"
    )
  )
  invisible(paths)
}

species_dynamics_same_path <- function(first, second) {
  identical(
    species_dynamics_normalize_path(first),
    species_dynamics_normalize_path(second)
  )
}

confirm_species_dynamics_outputs <- function(jobs,
                                             out_dir,
                                             show_progress = TRUE) {
  species_dynamics_assert_scalar_logical(show_progress, "show_progress")
  manifest_path <- file.path(
    out_dir,
    "sp_hist_fut_dynamics_manifest.rds"
  )
  if (!file.exists(manifest_path)) {
    return(list(job_ids = character(0), input_paths = character(0)))
  }

  previous <- tryCatch(
    readRDS(manifest_path),
    error = function(e) NULL
  )
  required_columns <- c(
    "job_id",
    "historical_path",
    "future_path",
    "output_path",
    "status"
  )
  if (!is.data.frame(previous) ||
      length(setdiff(required_columns, names(previous))) > 0L ||
      anyDuplicated(previous$job_id)) {
    return(list(job_ids = character(0), input_paths = character(0)))
  }

  previous_index <- match(jobs$job_id, previous$job_id)
  candidates <- which(!is.na(previous_index))
  candidates <- candidates[vapply(candidates, function(i) {
    old <- previous[previous_index[[i]], , drop = FALSE]
    old$status[[1L]] %in% c("written", "skipped_existing") &&
      species_dynamics_same_path(
        jobs$historical_path[[i]],
        old$historical_path[[1L]]
      ) &&
      species_dynamics_same_path(
        jobs$future_path[[i]],
        old$future_path[[1L]]
      ) &&
      species_dynamics_same_path(
        jobs$output_path[[i]],
        old$output_path[[1L]]
      ) &&
      file.exists(jobs$output_path[[i]])
  }, logical(1))]

  if (length(candidates) == 0L) {
    return(list(job_ids = character(0), input_paths = character(0)))
  }

  progress_id <- NULL
  progress_done <- FALSE
  if (show_progress) {
    progress_id <- species_dynamics_progress_bar(
      "Confirming existing outputs",
      length(candidates)
    )
    on.exit({
      if (!progress_done) {
        cli::cli_progress_done(id = progress_id)
      }
    }, add = TRUE)
  }

  confirmed_jobs <- character(0)
  confirmed_inputs <- character(0)
  for (position in seq_along(candidates)) {
    i <- candidates[[position]]
    if (show_progress) {
      cli::cli_progress_update(
        id = progress_id,
        set = position - 1L,
        inc = 0,
        status = jobs$job_id[[i]],
        force = TRUE
      )
    }
    valid <- tryCatch(
      {
        validate_species_dynamics_output(
          jobs$output_path[[i]],
          jobs$historical_path[[i]],
          jobs$future_path[[i]]
        )
        TRUE
      },
      error = function(e) FALSE
    )
    if (valid) {
      confirmed_jobs <- c(confirmed_jobs, jobs$job_id[[i]])
      confirmed_inputs <- c(
        confirmed_inputs,
        jobs$historical_path[[i]],
        jobs$future_path[[i]]
      )
    }
    if (show_progress) {
      cli::cli_progress_update(
        id = progress_id,
        set = position,
        inc = 0,
        status = if (valid) "Confirmed" else "Not reusable",
        force = TRUE
      )
    }
  }

  if (show_progress) {
    cli::cli_progress_done(id = progress_id)
    progress_done <- TRUE
  }

  list(
    job_ids = unique(confirmed_jobs),
    input_paths = unique(confirmed_inputs)
  )
}

validate_species_dynamics_inputs <- function(
    jobs,
    show_progress = TRUE,
    cache_dir = NULL,
    confirm_existing_outputs = FALSE,
    validation_mode = c("metadata", "full"),
    cache_checkpoint_every = 25L) {
  validation_mode <- species_dynamics_validation_mode(validation_mode)
  species_dynamics_assert_scalar_logical(show_progress, "show_progress")
  species_dynamics_assert_scalar_logical(
    confirm_existing_outputs,
    "confirm_existing_outputs"
  )
  if (confirm_existing_outputs && is.null(cache_dir)) {
    stop(
      "`cache_dir` is required when `confirm_existing_outputs = TRUE`.",
      call. = FALSE
    )
  }
  if (length(cache_checkpoint_every) != 1L ||
      is.na(cache_checkpoint_every) ||
      cache_checkpoint_every != as.integer(cache_checkpoint_every) ||
      cache_checkpoint_every < 1L) {
    stop(
      "`cache_checkpoint_every` must be one positive integer.",
      call. = FALSE
    )
  }
  cache_checkpoint_every <- as.integer(cache_checkpoint_every)
  required_columns <- c("historical_path", "future_path")
  missing_columns <- setdiff(required_columns, names(jobs))
  if (length(missing_columns) > 0L) {
    stop(
      "Jobs are missing input path column(s): ",
      paste(missing_columns, collapse = ", "),
      ".",
      call. = FALSE
    )
  }

  input_paths <- unique(c(jobs$historical_path, jobs$future_path))
  records <- species_dynamics_input_metadata(input_paths)
  cache <- load_species_dynamics_input_cache(cache_dir)
  cache_hits <- 0L

  if (!is.null(cache) && nrow(cache) > 0L) {
    cache_index <- match(records$normalized_path, cache$normalized_path)
    reusable <- !is.na(cache_index)
    reusable[reusable] <- vapply(which(reusable), function(i) {
      old <- cache[cache_index[[i]], , drop = FALSE]
      identical(records$size_bytes[[i]], old$size_bytes[[1L]]) &&
        identical(
          records$mtime_numeric[[i]],
          old$mtime_numeric[[1L]]
        ) &&
        if (identical(validation_mode, "full")) {
          identical(old$validation_mode[[1L]], "full")
        } else {
          old$validation_mode[[1L]] %in% c("metadata", "full")
        }
    }, logical(1))
    records$validated[reusable] <- TRUE
    records$validated_at[reusable] <- cache$validated_at[
      cache_index[reusable]
    ]
    records$validation_method[reusable] <- "cache"
    records$validation_mode[reusable] <- cache$validation_mode[
      cache_index[reusable]
    ]
    cache_hits <- sum(reusable)
  }

  confirmed <- list(job_ids = character(0), input_paths = character(0))
  if (confirm_existing_outputs &&
      identical(validation_mode, "metadata") &&
      any(!records$validated)) {
    confirmed <- confirm_species_dynamics_outputs(
      jobs,
      out_dir = cache_dir,
      show_progress = show_progress
    )
    confirmed_normalized <- species_dynamics_normalize_path(
      confirmed$input_paths
    )
    confirmed_rows <- records$normalized_path %in% confirmed_normalized
    confirmed_rows <- confirmed_rows & !records$validated
    if (any(confirmed_rows)) {
      records$validated[confirmed_rows] <- TRUE
      records$validated_at[confirmed_rows] <-
        species_dynamics_timestamp()
      records$validation_method[confirmed_rows] <-
        "existing_output_manifest"
      records$validation_mode[confirmed_rows] <- "metadata"
      write_species_dynamics_input_cache(records, cache_dir)
    }
  }

  unresolved <- which(!records$validated)
  failures <- character(0)
  progress_id <- NULL
  progress_done <- FALSE

  if (length(unresolved) == 0L) {
    if (show_progress) {
      cli::cli_alert_success(
        "Reused prior validation for all {nrow(records)} unchanged binary inputs."
      )
    }
    result <- records$path
    attr(result, "validation_records") <- records
    attr(result, "cache_hits") <- cache_hits
    attr(result, "validated_count") <- 0L
    attr(result, "prevalidated_job_ids") <- confirmed$job_ids
    return(invisible(result))
  }

  if (show_progress) {
    progress_id <- species_dynamics_progress_bar(
      if (identical(validation_mode, "metadata")) {
        "Checking binary input metadata"
      } else {
        "Auditing binary input values"
      },
      length(unresolved)
    )
    on.exit({
      if (!progress_done) {
        cli::cli_progress_done(id = progress_id)
      }
    }, add = TRUE)
  }

  validated_count <- 0L
  checkpointed_count <- 0L
  for (position in seq_along(unresolved)) {
    i <- unresolved[[position]]
    path <- records$path[[i]]
    if (show_progress) {
      cli::cli_progress_update(
        id = progress_id,
        set = position - 1L,
        inc = 0,
        status = basename(path),
        force = TRUE
      )
    }

    error <- tryCatch(
      {
        species_dynamics_validate_input_file(
          path,
          validation_mode = validation_mode
        )
        NULL
      },
      error = function(e) conditionMessage(e)
    )
    if (!is.null(error)) {
      failures <- c(failures, error)
    } else {
      records$validated[[i]] <- TRUE
      records$validated_at[[i]] <- species_dynamics_timestamp()
      records$validation_method[[i]] <- if (
        identical(validation_mode, "metadata")
      ) {
        "metadata_check"
      } else {
        "binary_scan"
      }
      records$validation_mode[[i]] <- validation_mode
      validated_count <- validated_count + 1L
      if (validated_count %% cache_checkpoint_every == 0L ||
          position == length(unresolved)) {
        write_species_dynamics_input_cache(records, cache_dir)
        checkpointed_count <- validated_count
      }
    }

    if (show_progress) {
      cli::cli_progress_update(
        id = progress_id,
        set = position,
        inc = 0,
        status = if (is.null(error)) "Validated" else "Validation failed",
        force = TRUE
      )
    }
  }

  if (show_progress) {
    cli::cli_progress_done(id = progress_id)
    progress_done <- TRUE
  }

  if (validated_count > checkpointed_count) {
    write_species_dynamics_input_cache(records, cache_dir)
  }

  if (length(failures) > 0L) {
    stop(
      paste0(
        length(failures),
        " input raster validation failure(s):\n- ",
        paste(failures, collapse = "\n- ")
      ),
      call. = FALSE
    )
  }

  result <- records$path
  attr(result, "validation_records") <- records
  attr(result, "cache_hits") <- cache_hits
  attr(result, "validated_count") <- validated_count
  attr(result, "prevalidated_job_ids") <- confirmed$job_ids
  invisible(result)
}

classify_species_dynamics <- function(historical,
                                      future,
                                      validate_values = TRUE) {
  species_dynamics_assert_scalar_logical(
    validate_values,
    "validate_values"
  )
  if (!inherits(historical, "SpatRaster") ||
      !inherits(future, "SpatRaster")) {
    stop(
      "`historical` and `future` must both be terra SpatRaster objects.",
      call. = FALSE
    )
  }
  if (terra::nlyr(historical) != 1L ||
      terra::nlyr(future) != 1L) {
    stop(
      "Historical and future inputs must each contain exactly one layer.",
      call. = FALSE
    )
  }
  if (!isTRUE(terra::compareGeom(
    historical,
    future,
    lyrs = FALSE,
    stopOnError = FALSE
  ))) {
    stop(
      "Historical and future raster geometries differ.",
      call. = FALSE
    )
  }

  if (validate_values) {
    species_dynamics_binary_summary(historical, "Historical raster")
    species_dynamics_binary_summary(future, "Future raster")
  }

  dynamics <- terra::ifel(
    (historical == 0) & (future == 0),
    1L,
    terra::ifel(
      (historical == 0) & (future == 1),
      2L,
      terra::ifel(
        (historical == 1) & (future == 0),
        3L,
        terra::ifel(
          (historical == 1) & (future == 1),
          4L,
          NA
        )
      )
    )
  )
  dynamics <- terra::as.factor(dynamics)
  levels(dynamics) <- dyn_levels
  names(dynamics) <- "dynamics"
  dynamics
}

species_dynamics_output_levels <- function(raster) {
  level_tables <- terra::levels(raster)
  if (length(level_tables) != 1L) {
    stop("Output raster does not contain the required RAT.", call. = FALSE)
  }
  level_table <- level_tables[[1L]]
  if (!is.data.frame(level_table) ||
      ncol(level_table) < 2L ||
      nrow(level_table) < 1L) {
    stop("Output raster does not contain the required RAT.", call. = FALSE)
  }
  data.frame(
    id = suppressWarnings(as.integer(level_table[[1L]])),
    dynamics = as.character(level_table[[2L]]),
    stringsAsFactors = FALSE
  )
}

validate_species_dynamics_output <- function(output_path,
                                             historical_path,
                                             future_path) {
  if (!file.exists(output_path)) {
    stop("Dynamics output does not exist: ", output_path, call. = FALSE)
  }
  output_info <- file.info(output_path)
  if (is.na(output_info$size) || output_info$size <= 0) {
    stop("Dynamics output is empty: ", output_path, call. = FALSE)
  }

  input_info <- file.info(c(historical_path, future_path))
  if (any(is.na(input_info$mtime))) {
    stop(
      "Cannot determine source raster modification times for output validation.",
      call. = FALSE
    )
  }
  if (output_info$mtime < max(input_info$mtime)) {
    stop(
      "Existing dynamics output is older than one or both source rasters: ",
      output_path,
      ". Rerun with `overwrite = TRUE`.",
      call. = FALSE
    )
  }

  output <- tryCatch(
    terra::rast(output_path),
    error = function(e) {
      stop(
        "Dynamics output cannot be opened: ",
        output_path,
        ". ",
        conditionMessage(e),
        call. = FALSE
      )
    }
  )
  historical <- terra::rast(historical_path)
  future <- terra::rast(future_path)

  if (terra::nlyr(output) != 1L) {
    stop("Dynamics output must contain exactly one layer.", call. = FALSE)
  }
  if (!isTRUE(terra::compareGeom(
    historical,
    future,
    lyrs = FALSE,
    stopOnError = FALSE
  )) ||
      !isTRUE(terra::compareGeom(
        output,
        historical,
        lyrs = FALSE,
        stopOnError = FALSE
      ))) {
    stop(
      "Dynamics output geometry does not match both source rasters.",
      call. = FALSE
    )
  }
  if (!identical(terra::datatype(output), "INT1U")) {
    stop(
      "Dynamics output datatype must be INT1U, not ",
      terra::datatype(output),
      ".",
      call. = FALSE
    )
  }

  actual_levels <- species_dynamics_output_levels(output)
  if (!identical(actual_levels, dyn_levels)) {
    stop(
      "Dynamics output RAT does not match the required four-class codebook.",
      call. = FALSE
    )
  }

  invisible(list(
    output_path = output_path,
    output_size_bytes = as.numeric(output_info$size)
  ))
}

species_dynamics_move_file <- function(source, target) {
  moved <- suppressWarnings(file.rename(source, target))
  if (!isTRUE(moved)) {
    moved <- file.copy(
      source,
      target,
      overwrite = FALSE,
      copy.mode = TRUE
    )
    if (isTRUE(moved)) {
      unlink(source, force = TRUE)
    }
  }
  isTRUE(moved) && file.exists(target)
}

write_species_dynamics_raster_safely <- function(dynamics,
                                                  output_path,
                                                  historical_path,
                                                  future_path,
                                                  overwrite = FALSE) {
  species_dynamics_assert_scalar_logical(overwrite, "overwrite")
  target <- fortify_output_path(output_path)
  target_sidecar <- paste0(target, ".aux.xml")
  target_directory <- dirname(target)

  if (!dir.exists(target_directory)) {
    created <- dir.create(
      target_directory,
      recursive = TRUE,
      showWarnings = FALSE
    )
    if (!isTRUE(created) && !dir.exists(target_directory)) {
      stop(
        "Could not create dynamics output directory: ",
        target_directory,
        call. = FALSE
      )
    }
  }
  if (file.exists(target) && !overwrite) {
    stop(
      "Dynamics output already exists and `overwrite = FALSE`: ",
      target,
      call. = FALSE
    )
  }

  staging <- tempfile(
    pattern = ".dynamics_write_",
    tmpdir = target_directory,
    fileext = ".tif"
  )
  staging_sidecar <- paste0(staging, ".aux.xml")
  had_target <- file.exists(target)
  had_sidecar <- file.exists(target_sidecar)
  raster_backup <- NULL
  sidecar_backup <- NULL
  promotion_started <- FALSE
  promoted <- FALSE

  on.exit({
    for (path in c(staging, staging_sidecar)) {
      if (file.exists(path)) {
        unlink(path, force = TRUE)
      }
    }
    if (!promoted) {
      if (!is.null(raster_backup) && file.exists(raster_backup)) {
        if (file.exists(target)) {
          unlink(target, force = TRUE)
        }
        species_dynamics_move_file(raster_backup, target)
      } else if (promotion_started && !had_target && file.exists(target)) {
        unlink(target, force = TRUE)
      }
      if (!is.null(sidecar_backup) && file.exists(sidecar_backup)) {
        if (file.exists(target_sidecar)) {
          unlink(target_sidecar, force = TRUE)
        }
        species_dynamics_move_file(sidecar_backup, target_sidecar)
      } else if (promotion_started &&
                 !had_sidecar &&
                 file.exists(target_sidecar)) {
        unlink(target_sidecar, force = TRUE)
      }
    } else {
      for (path in c(raster_backup, sidecar_backup)) {
        if (!is.null(path) && file.exists(path)) {
          unlink(path, force = TRUE)
        }
      }
    }
  }, add = TRUE)

  terra::writeRaster(
    dynamics,
    filename = staging,
    overwrite = TRUE,
    wopt = list(
      datatype = "INT1U",
      gdal = c("COMPRESS=LZW", "TILED=YES")
    )
  )
  if (!file.exists(staging) ||
      is.na(file.info(staging)$size) ||
      file.info(staging)$size <= 0) {
    stop("Dynamics staging raster was not written correctly.", call. = FALSE)
  }
  if (!file.exists(staging_sidecar)) {
    stop(
      "Dynamics staging raster is missing its categorical RAT sidecar.",
      call. = FALSE
    )
  }
  staging_raster <- terra::rast(staging)
  if (!isTRUE(terra::compareGeom(
    staging_raster,
    dynamics,
    lyrs = FALSE,
    stopOnError = FALSE
  )) ||
      !identical(
        species_dynamics_output_levels(staging_raster),
        dyn_levels
      )) {
    stop(
      "Dynamics staging raster failed geometry or RAT verification.",
      call. = FALSE
    )
  }

  if (file.exists(target)) {
    raster_backup <- tempfile(
      pattern = ".dynamics_backup_",
      tmpdir = target_directory,
      fileext = ".tif"
    )
    if (!species_dynamics_move_file(target, raster_backup)) {
      stop("Could not back up existing raster: ", target, call. = FALSE)
    }
  }
  if (file.exists(target_sidecar)) {
    sidecar_backup <- tempfile(
      pattern = ".dynamics_backup_",
      tmpdir = target_directory,
      fileext = ".tif.aux.xml"
    )
    if (!species_dynamics_move_file(target_sidecar, sidecar_backup)) {
      stop(
        "Could not back up existing RAT sidecar: ",
        target_sidecar,
        call. = FALSE
      )
    }
  }

  promotion_started <- TRUE
  if (!species_dynamics_move_file(staging, target)) {
    stop("Could not promote dynamics raster to: ", target, call. = FALSE)
  }
  if (!species_dynamics_move_file(staging_sidecar, target_sidecar)) {
    stop(
      "Could not promote dynamics RAT sidecar to: ",
      target_sidecar,
      call. = FALSE
    )
  }

  validation <- validate_species_dynamics_output(
    target,
    historical_path,
    future_path
  )
  promoted <- TRUE
  for (path in c(raster_backup, sidecar_backup)) {
    if (!is.null(path) && file.exists(path)) {
      unlink(path, force = TRUE)
    }
  }
  raster_backup <- NULL
  sidecar_backup <- NULL

  invisible(validation$output_path)
}

process_species_dynamics_job <- function(job,
                                         overwrite = FALSE,
                                         validate_values = TRUE) {
  species_dynamics_assert_scalar_logical(overwrite, "overwrite")
  species_dynamics_assert_scalar_logical(
    validate_values,
    "validate_values"
  )
  required_columns <- c(
    "historical_path",
    "future_path",
    "output_path"
  )
  missing_columns <- setdiff(required_columns, names(job))
  if (nrow(job) != 1L || length(missing_columns) > 0L) {
    stop(
      "`job` must contain one row and historical, future, and output paths.",
      call. = FALSE
    )
  }

  historical_path <- job$historical_path[[1L]]
  future_path <- job$future_path[[1L]]
  output_path <- job$output_path[[1L]]

  if (file.exists(output_path) && !overwrite) {
    validation <- validate_species_dynamics_output(
      output_path,
      historical_path,
      future_path
    )
    return(list(
      status = "skipped_existing",
      message = "Existing output passed geometry, RAT, and freshness checks.",
      output_size_bytes = validation$output_size_bytes
    ))
  }

  historical <- terra::rast(historical_path)
  future <- terra::rast(future_path)
  dynamics <- classify_species_dynamics(
    historical,
    future,
    validate_values = validate_values
  )

  resolved_output <- write_species_dynamics_raster_safely(
    dynamics,
    output_path = output_path,
    historical_path = historical_path,
    future_path = future_path,
    overwrite = overwrite
  )
  validation <- validate_species_dynamics_output(
    resolved_output,
    historical_path,
    future_path
  )

  list(
    status = "written",
    message = "Dynamics raster written and verified.",
    output_size_bytes = validation$output_size_bytes
  )
}

species_dynamics_replace_file <- function(staging_path, target_path) {
  backup_path <- NULL
  promoted <- FALSE
  on.exit({
    if (file.exists(staging_path)) {
      unlink(staging_path, force = TRUE)
    }
    if (!is.null(backup_path) && file.exists(backup_path)) {
      if (!promoted && !file.exists(target_path)) {
        file.rename(backup_path, target_path)
      } else if (promoted) {
        unlink(backup_path, force = TRUE)
      }
    }
  }, add = TRUE)

  if (file.exists(target_path)) {
    backup_path <- tempfile(
      pattern = ".manifest_backup_",
      tmpdir = dirname(target_path),
      fileext = paste0(".", tools::file_ext(target_path))
    )
    backed_up <- suppressWarnings(file.rename(target_path, backup_path))
    if (!isTRUE(backed_up)) {
      backed_up <- file.copy(
        target_path,
        backup_path,
        overwrite = FALSE,
        copy.mode = TRUE
      )
      if (isTRUE(backed_up)) {
        unlink(target_path, force = TRUE)
      }
    }
    if (!isTRUE(backed_up) || file.exists(target_path)) {
      stop(
        "Could not back up existing manifest: ",
        target_path,
        call. = FALSE
      )
    }
  }

  renamed <- suppressWarnings(file.rename(staging_path, target_path))
  if (!isTRUE(renamed)) {
    copied <- file.copy(
      staging_path,
      target_path,
      overwrite = FALSE,
      copy.mode = TRUE
    )
    if (isTRUE(copied)) {
      unlink(staging_path, force = TRUE)
    }
  }
  if (!file.exists(target_path)) {
    if (!is.null(backup_path) && file.exists(backup_path)) {
      file.rename(backup_path, target_path)
      backup_path <- NULL
    }
    stop("Could not promote manifest to: ", target_path, call. = FALSE)
  }

  promoted <- TRUE
  if (!is.null(backup_path) && file.exists(backup_path)) {
    unlink(backup_path, force = TRUE)
  }
  backup_path <- NULL
  invisible(target_path)
}

species_dynamics_write_atomic <- function(target_path, writer) {
  staging_path <- tempfile(
    pattern = ".manifest_write_",
    tmpdir = dirname(target_path),
    fileext = paste0(".", tools::file_ext(target_path))
  )
  on.exit({
    if (file.exists(staging_path)) {
      unlink(staging_path, force = TRUE)
    }
  }, add = TRUE)
  writer(staging_path)
  if (!file.exists(staging_path) ||
      is.na(file.info(staging_path)$size) ||
      file.info(staging_path)$size <= 0) {
    stop(
      "Manifest staging write produced no usable file for: ",
      target_path,
      call. = FALSE
    )
  }
  species_dynamics_replace_file(staging_path, target_path)
}

write_species_dynamics_manifests <- function(manifest, out_dir) {
  species_dynamics_assert_output_dir(out_dir)
  if (!dir.exists(out_dir)) {
    stop("Manifest output directory does not exist: ", out_dir, call. = FALSE)
  }

  rds_path <- file.path(
    out_dir,
    "sp_hist_fut_dynamics_manifest.rds"
  )
  csv_path <- file.path(
    out_dir,
    "sp_hist_fut_dynamics_manifest.csv"
  )

  species_dynamics_write_atomic(
    rds_path,
    function(path) saveRDS(manifest, path)
  )
  species_dynamics_write_atomic(
    csv_path,
    function(path) utils::write.csv(
      manifest,
      file = path,
      row.names = FALSE,
      na = "",
      fileEncoding = "UTF-8"
    )
  )

  invisible(list(rds = rds_path, csv = csv_path))
}

write_species_dynamics_levels <- function(out_dir) {
  levels_path <- file.path(
    out_dir,
    "sp_hist_fut_dynamics_levels.csv"
  )
  species_dynamics_write_atomic(
    levels_path,
    function(path) utils::write.csv(
      dyn_levels,
      file = path,
      row.names = FALSE,
      fileEncoding = "UTF-8"
    )
  )
  invisible(levels_path)
}

run_species_dynamics <- function(
    species = "all",
    out_dir = file.path(
      "data",
      "post_model_maps",
      "sp_hist_fut_dynamics"
    ),
    overwrite = FALSE,
    fuzzy = FALSE,
    max_fuzzy_distance = 0.10,
    input_validation = c("metadata", "full"),
    catalog = best_mod_set_paths) {
  species_dynamics_assert_output_dir(out_dir)
  species_dynamics_assert_scalar_logical(overwrite, "overwrite")
  input_validation <- species_dynamics_validation_mode(input_validation)

  if (!dir.exists(out_dir)) {
    created <- dir.create(
      out_dir,
      recursive = TRUE,
      showWarnings = FALSE
    )
    if (!isTRUE(created) && !dir.exists(out_dir)) {
      stop("Could not create output directory: ", out_dir, call. = FALSE)
    }
  }

  terra::terraOptions(progress = 0)
  jobs <- build_species_dynamics_jobs(
    species = species,
    out_dir = out_dir,
    fuzzy = fuzzy,
    max_fuzzy_distance = max_fuzzy_distance,
    catalog = catalog
  )

  cli::cli_alert_info(
    "Prepared {nrow(jobs)} dynamics map job{?s} for {length(unique(jobs$tkey))} species."
  )
  cli::cli_alert_info(
    if (identical(input_validation, "metadata")) {
      paste(
        "Using fast metadata validation;",
        "set `input_validation = \"full\"` for an exhaustive cell audit."
      )
    } else {
      "Using exhaustive cell-value validation for every changed input."
    }
  )
  validation_result <- validate_species_dynamics_inputs(
    jobs,
    show_progress = TRUE,
    cache_dir = out_dir,
    confirm_existing_outputs = !overwrite,
    validation_mode = input_validation
  )
  prevalidated_job_ids <- attr(
    validation_result,
    "prevalidated_job_ids",
    exact = TRUE
  )
  if (is.null(prevalidated_job_ids)) {
    prevalidated_job_ids <- character(0)
  }
  write_species_dynamics_levels(out_dir)
  write_species_dynamics_manifests(jobs, out_dir)

  progress_id <- species_dynamics_progress_bar(
    "Calculating maps",
    nrow(jobs)
  )
  progress_done <- FALSE
  on.exit({
    if (!progress_done) {
      try(cli::cli_progress_done(id = progress_id), silent = TRUE)
    }
    try(write_species_dynamics_manifests(jobs, out_dir), silent = TRUE)
  }, add = TRUE)

  failures <- character(0)
  for (i in seq_len(nrow(jobs))) {
    status_label <- paste(
      jobs$sp_name[[i]],
      jobs$future_period[[i]],
      jobs$future_scenario[[i]]
    )
    cli::cli_progress_update(
      id = progress_id,
      set = i - 1L,
      inc = 0,
      status = paste("Processing", status_label),
      force = TRUE
    )

    jobs$started_at[[i]] <- species_dynamics_timestamp()
    result <- tryCatch(
      if (!overwrite && jobs$job_id[[i]] %in% prevalidated_job_ids) {
        list(
          status = "skipped_existing",
          message = paste(
            "Existing output was verified while rebuilding the",
            "binary-input validation cache."
          ),
          output_size_bytes = as.numeric(
            file.info(jobs$output_path[[i]])$size
          )
        )
      } else {
        process_species_dynamics_job(
          jobs[i, , drop = FALSE],
          overwrite = overwrite,
          validate_values = FALSE
        )
      },
      error = function(e) {
        list(
          status = "failed",
          message = conditionMessage(e),
          output_size_bytes = if (file.exists(jobs$output_path[[i]])) {
            as.numeric(file.info(jobs$output_path[[i]])$size)
          } else {
            NA_real_
          }
        )
      }
    )

    jobs$status[[i]] <- result$status
    jobs$message[[i]] <- result$message
    jobs$finished_at[[i]] <- species_dynamics_timestamp()
    jobs$output_size_bytes[[i]] <- result$output_size_bytes
    if (identical(result$status, "failed")) {
      failures <- c(
        failures,
        paste0(jobs$job_id[[i]], ": ", result$message)
      )
    }

    write_species_dynamics_manifests(jobs, out_dir)
    cli::cli_progress_update(
      id = progress_id,
      set = i,
      inc = 0,
      status = paste(
        if (identical(result$status, "failed")) "Failed" else "Finished",
        status_label
      ),
      force = TRUE
    )
  }

  cli::cli_progress_done(id = progress_id)
  progress_done <- TRUE
  write_species_dynamics_manifests(jobs, out_dir)

  if (length(failures) > 0L) {
    shown <- utils::head(failures, 20L)
    stop(
      paste0(
        length(failures),
        " dynamics job(s) failed. See the manifest for full details:\n- ",
        paste(shown, collapse = "\n- "),
        if (length(failures) > length(shown)) {
          paste0(
            "\n- ... and ",
            length(failures) - length(shown),
            " more."
          )
        } else {
          ""
        }
      ),
      call. = FALSE
    )
  }

  cli::cli_alert_success(
    paste0(
      "Completed {nrow(jobs)} dynamics map job{?s}; ",
      "{sum(jobs$status == 'written')} written and ",
      "{sum(jobs$status == 'skipped_existing')} validated/skipped."
    )
  )
  invisible(jobs)
}

species_dynamics_autorun_on_source <- getOption(
  "wisdm.run_species_dynamics_on_source",
  interactive()
)

if (sys.nframe() == 0L || isTRUE(species_dynamics_autorun_on_source)) {
  run_species_dynamics()
}
