#-------------------------------------------------------------------------------
# Manual RStudio background chunk runner for wiSDM
#-------------------------------------------------------------------------------
# Usage:
#   1. Edit the user settings below, especially `nr_active_block`.
#   2. In RStudio, run this file as a background job.
#   3. Start other chunk jobs by changing `nr_active_block` and running again.
#
# Parameter guide:
#   `n_blocks`
#      Total number of species chunks. Use the same value for every background
#      job that belongs to the same batch.
#
#   `nr_active_block`
#      The chunk number this R process should run. Manually change only this
#      value when starting each RStudio background job. Must be between 1 and
#      `n_blocks`.
#
#   `project_prefix`
#      Base name shared by this batch of chunk jobs. Each species gets its own
#      project name built as:
#        <project_prefix>_<chunk_number>_<sanitized_species_name>
#      Example: onestop_01_arthurdendyus_triangulatus
#
#   `species_list_path`
#      Path to the Excel file containing the candidate species list.
#
#   `species_column`
#      Column in `species_list_path` that contains the species names.
#
#   `filter_column` and `filter_value`
#      Optional row filter used to select species from the list before chunking.
#      Set `filter_column <- NULL` or `filter_column <- ""` to use all rows.
#      `filter_value` may be a single value or a vector of accepted values.
#
#   `retry_failed`
#      If TRUE, a species previously logged as failed can be claimed again on a
#      later run. If FALSE, failed species are skipped on reruns.
#
#   `retry_skipped`
#      If TRUE, species previously logged as skipped can be claimed again. This
#      is usually FALSE because skipped species are commonly data-limited
#      species without enough information to build a climate model.
#
#   `verify_completed_outputs`
#      If TRUE, a registry row marked completed is skipped only when the expected
#      project folder and terminal model file(s) still exist. This protects reruns
#      after project folders were moved, deleted, or partially restored.
#
#   `reuse_successful_stages`
#      If TRUE, retry attempts can reuse earlier stage outputs instead of running
#      the stage again. Reuse is opt-in, only applies to attempt 2+, and requires
#      both a prior successful/reused stage event and the downstream files that
#      prove the stage output is still available.
#
#   `run_setup_01`
#      If TRUE, stage 01 is run once per `project_prefix` under a shared setup
#      lock. Other chunks skip it after the setup marker exists.
#
#   `force_setup_01`
#      If TRUE, stage 01 is rerun even if the setup marker already exists.
#
#   `dry_run`
#      If TRUE, the runner reports chunk membership, project names, and skip or
#      claim decisions without sourcing stages 02-05.
#
#   `config_guard`
#      Controls how preflight problems are handled. Use "stop" for the safest
#      default, "warn" to continue after warnings, or "ignore" to suppress
#      guard failures.
#
# Advanced parameter guide:
#   `max_project_name_chars`
#      Maximum generated project-name length. Long species names are shortened
#      with a stable digest suffix so project folders remain filesystem-safe.
#
#   `registry_lock_timeout` and `event_lock_timeout`
#      How long to wait for shared registry/event log locks. `Inf` is safest for
#      normal background jobs because writes to these files are brief.
#
#   `species_lock_timeout`
#      How long to wait for a per-species lock. The default of 5 seconds avoids
#      rare false skips during simultaneous chunk-worker claim attempts.
#
# Control files and verbose logs are written outside species projects under:
#   data/projects/<project_prefix>_chunk_control/
#
# This runner intentionally leaves 00_configurations.R and scripts 01-05 unchanged.
# It injects `project` and `species_to_model` before each stage is sourced.


#-------------------------------------------------------------------------------
# User settings
#-------------------------------------------------------------------------------
n_blocks <- 6
nr_active_block <- 1
project_prefix <- "onestop"

species_list_path <- file.path("data", "external", "Species_list_v5.xlsx")
species_column <- "Species"
filter_column <- "Modelling"
filter_value <- "Yes"

retry_failed <- TRUE
retry_skipped <- FALSE
verify_completed_outputs <- TRUE
reuse_successful_stages <- FALSE

run_setup_01 <- TRUE
force_setup_01 <- FALSE

dry_run <- FALSE
config_guard <- "stop" # one of "stop", "warn", or "ignore"


#-------------------------------------------------------------------------------
# Advanced settings
#-------------------------------------------------------------------------------
max_project_name_chars <- 90L
registry_lock_timeout <- Inf
event_lock_timeout <- Inf
species_lock_timeout <- 5


#-------------------------------------------------------------------------------
# Configuration preflight
#-------------------------------------------------------------------------------
check_config_guard_value <- function(config_guard) {
  valid_values <- c("stop", "warn", "ignore")
  if (!identical(config_guard, tolower(config_guard)) || !config_guard %in% valid_values) {
    stop(
      "config_guard must be one of: ",
      paste(valid_values, collapse = ", "),
      call. = FALSE
    )
  }
}

find_protected_config_assignments <- function(config_path,
                                              protected_names = c("project", "species_to_model")) {
  if (!file.exists(config_path)) {
    stop("Configuration file not found: ", config_path, call. = FALSE)
  }

  parsed_config <- parse(config_path, keep.source = TRUE)
  found <- character()

  inspect_expr <- function(expr) {
    if (!is.call(expr)) {
      return(invisible(NULL))
    }

    operator <- as.character(expr[[1]])

    if (operator %in% c("<-", "<<-", "=")) {
      lhs <- expr[[2]]
      if (is.symbol(lhs)) {
        lhs_name <- as.character(lhs)
        if (lhs_name %in% protected_names) {
          found <<- unique(c(found, lhs_name))
        }
      }
    }

    if (operator %in% c("->", "->>")) {
      rhs <- expr[[3]]
      if (is.symbol(rhs)) {
        rhs_name <- as.character(rhs)
        if (rhs_name %in% protected_names) {
          found <<- unique(c(found, rhs_name))
        }
      }
    }

    for (i in seq_along(expr)[-1]) {
      inspect_expr(expr[[i]])
    }

    invisible(NULL)
  }

  for (expr in parsed_config) {
    inspect_expr(expr)
  }

  found
}

format_preflight_value <- function(value) {
  if (is.null(value) || length(value) == 0) {
    return("")
  }
  value <- paste(as.character(value), collapse = "; ")
  value <- gsub("[\r\n\t]+", " ", value)
  value <- gsub("\\s+", " ", value)
  trimws(value)
}

preflight_pass_message <- function() {
  message("\033[32m✅ All wiSDM chunk-runner preflight checks passed.\033[39m")
}

handle_preflight_issue <- function(message_text, config_guard) {
  if (identical(config_guard, "stop")) {
    stop(message_text, call. = FALSE)
  }

  warning(message_text, call. = FALSE)
  invisible(FALSE)
}

run_configuration_preflight <- function(config_path, config_guard) {
  check_config_guard_value(config_guard)
  if (identical(config_guard, "ignore")) {
    return(invisible(FALSE))
  }

  active_assignments <- find_protected_config_assignments(config_path)
  if (length(active_assignments) == 0) {
    return(invisible(TRUE))
  }

  message_text <- paste0(
    "The chunk runner controls these variable(s), but active assignment(s) ",
    "were found in ", config_path, ": ",
    paste(active_assignments, collapse = ", "), ".\n",
    "Please comment out these assignment lines in 00_configurations.R before ",
    "running src/_run_wisdm_by_chunks.R. Keep the other shared configuration ",
    "parameters active."
  )

  if (identical(config_guard, "stop")) {
    stop(message_text, call. = FALSE)
  }

  warning(message_text, call. = FALSE)
  invisible(FALSE)
}

run_input_path_preflight <- function(config_path,
                                     config_guard,
                                     species_list_path,
                                     path_parameters = c(
                                       "user_specific_climate_data",
                                       "user_specific_landcover_data",
                                       "custom_eu_boundary_path",
                                       "custom_country_boundary_path"
                                     )) {
  if (identical(config_guard, "ignore")) {
    return(invisible(FALSE))
  }

  config_env <- new.env(parent = baseenv())
  sys.source(config_path, envir = config_env)

  missing_paths <- data.frame(
    parameter = character(),
    path = character(),
    stringsAsFactors = FALSE
  )

  if (is.null(species_list_path) ||
      length(species_list_path) != 1 ||
      is.na(species_list_path) ||
      !nzchar(species_list_path) ||
      !file.exists(species_list_path)) {
    missing_paths <- rbind(
      missing_paths,
      data.frame(
        parameter = "species_list_path",
        path = format_preflight_value(species_list_path),
        stringsAsFactors = FALSE
      )
    )
  }

  for (parameter in path_parameters) {
    if (!exists(parameter, envir = config_env, inherits = FALSE)) {
      next
    }

    value <- get(parameter, envir = config_env, inherits = FALSE)
    if (is.null(value)) {
      next
    }

    value <- as.character(value)
    value <- value[!is.na(value) & nzchar(value)]
    if (length(value) == 0) {
      next
    }

    missing_value <- value[!file.exists(value)]
    if (length(missing_value) > 0) {
      missing_paths <- rbind(
        missing_paths,
        data.frame(
          parameter = rep(parameter, length(missing_value)),
          path = missing_value,
          stringsAsFactors = FALSE
        )
      )
    }
  }

  if (nrow(missing_paths) == 0) {
    return(invisible(TRUE))
  }

  missing_lines <- paste0("  - ", missing_paths$parameter, ": ", missing_paths$path)
  message_text <- paste0(
    "The following path parameter(s) in ", config_path,
    " point to files that do not exist:\n",
    paste(missing_lines, collapse = "\n"), "\n",
    "Set the parameter to NULL if it should not be used, or fix the path before ",
    "running src/_run_wisdm_by_chunks.R."
  )

  handle_preflight_issue(message_text, config_guard)
}

preflight_config_ok <- run_configuration_preflight(
  config_path = file.path("src", "00_configurations.R"),
  config_guard = config_guard
)

preflight_paths_ok <- run_input_path_preflight(
  config_path = file.path("src", "00_configurations.R"),
  config_guard = config_guard,
  species_list_path = species_list_path
)


#-------------------------------------------------------------------------------
# Package setup
#-------------------------------------------------------------------------------
required_packages <- c("readxl", "dplyr", "stringr", "readr", "filelock", "digest")
installed_packages <- rownames(installed.packages())

for (package in required_packages) {
  if (!package %in% installed_packages) {
    install.packages(package)
  }
  suppressPackageStartupMessages(library(package, character.only = TRUE))
}


#-------------------------------------------------------------------------------
# Constants
#-------------------------------------------------------------------------------
stage_scripts <- list(
  "01" = file.path("src", "01_prepare_files_and_folders.R"),
  "02" = file.path("src", "02_global_occurrence_download.R"),
  "03" = file.path("src", "03_fit_climate_model.R"),
  "04" = file.path("src", "04_fit_habitat_model.R"),
  "05" = file.path("src", "05_cross_validation.R")
)

stage_labels <- c(
  "01" = "prepare_files_and_folders",
  "02" = "global_occurrence_download",
  "03" = "fit_climate_model",
  "04" = "fit_habitat_model",
  "05" = "cross_validation"
)

registry_columns <- c(
  "species_requested", "species_key", "chunk_id", "project", "attempt",
  "status", "outcome", "stage", "started_at", "updated_at", "completed_at",
  "pid", "run_id", "accepted_species", "accepted_taxonkey", "error_message"
)

event_columns <- c(
  "timestamp", "run_id", "pid", "chunk_id", "species_requested",
  "species_key", "project", "attempt", "stage", "state", "message", "log_file"
)

completed_columns <- c(
  "species_requested", "species_key", "project", "outcome", "completed_at", "run_id"
)


#-------------------------------------------------------------------------------
# General helpers
#-------------------------------------------------------------------------------
timestamp_now <- function() {
  format(Sys.time(), "%Y-%m-%d %H:%M:%S")
}

timestamp_file <- function() {
  format(Sys.time(), "%Y%m%d_%H%M%S")
}

clean_one_line <- function(x, max_chars = 2000L) {
  if (length(x) == 0 || is.null(x) || is.na(x)) {
    return("")
  }
  x <- paste(as.character(x), collapse = " ")
  x <- gsub("[\r\n\t]+", " ", x)
  x <- gsub("\\s+", " ", x)
  x <- trimws(x)
  if (nchar(x) > max_chars) {
    x <- paste0(substr(x, 1, max_chars - 3L), "...")
  }
  x
}

sanitize_token <- function(x, fallback = "item") {
  original <- x
  x <- iconv(as.character(x), from = "", to = "ASCII//TRANSLIT", sub = "")
  x <- tolower(x)
  x <- gsub("[^a-z0-9]+", "_", x)
  x <- gsub("_+", "_", x)
  x <- gsub("^_|_$", "", x)
  if (!nzchar(x)) {
    x <- paste0(fallback, "_", substr(digest::digest(original, algo = "xxhash32"), 1, 8))
  }
  x
}

trim_with_hash <- function(prefix, slug, original, max_chars) {
  project <- paste0(prefix, "_", slug)
  if (nchar(project) <= max_chars) {
    return(project)
  }

  hash <- substr(digest::digest(original, algo = "xxhash32"), 1, 8)
  available <- max_chars - nchar(prefix) - nchar(hash) - 2L
  if (available < 10L) {
    stop("max_project_name_chars is too short for the configured project prefix.", call. = FALSE)
  }

  slug <- substr(slug, 1, available)
  slug <- gsub("_+$", "", slug)
  paste0(prefix, "_", slug, "_", hash)
}

make_unique_slugs <- function(species) {
  slugs <- vapply(species, sanitize_token, character(1), fallback = "species")
  duplicate_slugs <- slugs %in% slugs[duplicated(slugs)]

  if (any(duplicate_slugs)) {
    hashes <- substr(vapply(species, digest::digest, character(1), algo = "xxhash32"), 1, 8)
    slugs[duplicate_slugs] <- paste0(slugs[duplicate_slugs], "_", hashes[duplicate_slugs])
  }

  slugs
}

species_key_from_slug <- function(slug, species_requested) {
  paste0(substr(slug, 1, 70), "_", substr(digest::digest(species_requested, algo = "xxhash32"), 1, 8))
}

safe_dir_create <- function(path) {
  if (!dir.exists(path)) {
    dir.create(path, recursive = TRUE, showWarnings = FALSE)
  }
}

lock_file <- function(path, timeout = Inf) {
  safe_dir_create(dirname(path))
  tryCatch(filelock::lock(path, timeout = timeout), error = function(e) NULL)
}

unlock_file <- function(lock) {
  if (!is.null(lock)) {
    try(filelock::unlock(lock), silent = TRUE)
  }
}

with_locked_file <- function(path, expr, timeout = Inf) {
  lock <- lock_file(path, timeout = timeout)
  if (is.null(lock)) {
    stop("Could not acquire lock: ", path, call. = FALSE)
  }
  on.exit(unlock_file(lock), add = TRUE)
  force(expr)
}

csv_has_content <- function(path) {
  file.exists(path) && isTRUE(file.info(path)$size > 0)
}

runner_message <- function(state,
                           species_requested,
                           project,
                           stage,
                           index,
                           total,
                           detail = NULL) {
  species_label <- if (length(species_requested) == 0 || is.na(species_requested)) {
    "setup"
  } else {
    species_requested
  }
  index_label <- if (is.na(index)) "setup" else paste0(index, "/", total)
  detail <- clean_one_line(detail, max_chars = 500L)
  detail <- if (nzchar(detail)) paste0(" - ", detail) else ""

  message(
    "[", timestamp_now(), "] ",
    "chunk ", chunk_id, " | ",
    index_label, " | ",
    species_label, " | ",
    project, " | ",
    stage, " | ",
    state,
    detail
  )
}


#-------------------------------------------------------------------------------
# Control file helpers
#-------------------------------------------------------------------------------
empty_registry <- function() {
  out <- as.data.frame(setNames(rep(list(character()), length(registry_columns)), registry_columns))
  out$attempt <- integer()
  out
}

empty_completed <- function() {
  as.data.frame(setNames(rep(list(character()), length(completed_columns)), completed_columns))
}

empty_events <- function() {
  as.data.frame(setNames(rep(list(character()), length(event_columns)), event_columns))
}

read_registry_file <- function(path) {
  if (!csv_has_content(path)) {
    return(empty_registry())
  }

  registry <- suppressWarnings(readr::read_csv(
    path,
    col_types = readr::cols(.default = readr::col_character(), attempt = readr::col_integer()),
    show_col_types = FALSE
  ))
  registry <- as.data.frame(registry, stringsAsFactors = FALSE)

  for (column in setdiff(registry_columns, names(registry))) {
    registry[[column]] <- NA_character_
  }

  registry <- registry[, registry_columns, drop = FALSE]
  registry$attempt <- as.integer(registry$attempt)
  registry$attempt[is.na(registry$attempt)] <- 0L
  registry
}

write_registry_file <- function(path, registry) {
  safe_dir_create(dirname(path))
  registry <- registry[, registry_columns, drop = FALSE]
  readr::write_csv(registry, path, na = "")
}

read_completed_file <- function(path) {
  if (!csv_has_content(path)) {
    return(empty_completed())
  }

  completed <- suppressWarnings(readr::read_csv(
    path,
    col_types = readr::cols(.default = readr::col_character()),
    show_col_types = FALSE
  ))
  completed <- as.data.frame(completed, stringsAsFactors = FALSE)

  for (column in setdiff(completed_columns, names(completed))) {
    completed[[column]] <- NA_character_
  }

  completed[, completed_columns, drop = FALSE]
}

read_event_file <- function(path) {
  if (!csv_has_content(path)) {
    return(empty_events())
  }

  events <- suppressWarnings(readr::read_csv(
    path,
    col_types = readr::cols(.default = readr::col_character()),
    show_col_types = FALSE
  ))
  events <- as.data.frame(events, stringsAsFactors = FALSE)

  for (column in setdiff(event_columns, names(events))) {
    events[[column]] <- NA_character_
  }

  events[, event_columns, drop = FALSE]
}

read_all_event_files <- function(control_dir) {
  files <- list.files(
    control_dir,
    pattern = "^chunk_[0-9]+_events\\.csv$",
    full.names = TRUE
  )
  files <- files[file.exists(files)]
  if (length(files) == 0) {
    return(empty_events())
  }

  dplyr::bind_rows(lapply(files, read_event_file))
}

modify_registry_row <- function(species_key, values) {
  with_locked_file(registry_lock_path, {
    registry <- read_registry_file(registry_path)
    row_index <- which(registry$species_key == species_key)
    if (length(row_index) == 0) {
      stop("Cannot update registry; species is not registered: ", species_key, call. = FALSE)
    }

    row_index <- row_index[length(row_index)]
    for (name in names(values)) {
      registry[[name]][row_index] <- values[[name]]
    }

    registry$updated_at[row_index] <- timestamp_now()
    write_registry_file(registry_path, registry)
    registry[row_index, , drop = FALSE]
  }, timeout = registry_lock_timeout)
}

append_event <- function(species_requested,
                         species_key,
                         project,
                         attempt,
                         stage,
                         state,
                         message_text = "",
                         log_file = "") {
  row <- data.frame(
    timestamp = timestamp_now(),
    run_id = run_id,
    pid = as.character(Sys.getpid()),
    chunk_id = chunk_id,
    species_requested = clean_one_line(species_requested),
    species_key = clean_one_line(species_key),
    project = clean_one_line(project),
    attempt = as.character(attempt),
    stage = clean_one_line(stage),
    state = clean_one_line(state),
    message = clean_one_line(message_text),
    log_file = clean_one_line(log_file, max_chars = 1000L),
    stringsAsFactors = FALSE
  )
  row <- row[, event_columns, drop = FALSE]

  with_locked_file(chunk_events_lock_path, {
    safe_dir_create(dirname(chunk_events_path))
    readr::write_csv(row, chunk_events_path, append = csv_has_content(chunk_events_path), na = "")
    invisible(row)
  }, timeout = event_lock_timeout)
}

append_completed_species <- function(species_requested, species_key, project, outcome) {
  row <- data.frame(
    species_requested = clean_one_line(species_requested),
    species_key = species_key,
    project = project,
    outcome = outcome,
    completed_at = timestamp_now(),
    run_id = run_id,
    stringsAsFactors = FALSE
  )
  row <- row[, completed_columns, drop = FALSE]

  with_locked_file(chunk_completed_lock_path, {
    completed <- read_completed_file(chunk_completed_path)
    if (!species_key %in% completed$species_key) {
      completed <- dplyr::bind_rows(completed, row)
      safe_dir_create(dirname(chunk_completed_path))
      readr::write_csv(completed, chunk_completed_path, na = "")
    }
    invisible(row)
  }, timeout = event_lock_timeout)
}


#-------------------------------------------------------------------------------
# Species list and expected output helpers
#-------------------------------------------------------------------------------
load_species_list <- function(path, species_column, filter_column, filter_value) {
  if (!file.exists(path)) {
    stop("Species list file not found: ", path, call. = FALSE)
  }

  species_table <- readxl::read_excel(path)

  if (!species_column %in% names(species_table)) {
    stop("Species column not found in species list: ", species_column, call. = FALSE)
  }

  if (!is.null(filter_column) && nzchar(filter_column)) {
    if (!filter_column %in% names(species_table)) {
      stop("Filter column not found in species list: ", filter_column, call. = FALSE)
    }
    species_table <- species_table[as.character(species_table[[filter_column]]) %in% filter_value, , drop = FALSE]
  }

  species <- trimws(as.character(species_table[[species_column]]))
  species <- species[!is.na(species) & nzchar(species)]
  species <- unique(species)

  if (length(species) == 0) {
    stop("No species were selected from the species list.", call. = FALSE)
  }

  species
}

make_chunk_assignments <- function(n_species, n_blocks) {
  if (n_species == 0) {
    return(integer())
  }
  if (n_species <= n_blocks) {
    return(seq_len(n_species))
  }
  as.integer(cut(seq_len(n_species), breaks = n_blocks, labels = FALSE))
}

species_model_name <- function(species) {
  sub("^(\\w+)\\s+(\\w+).*", "\\1_\\2", species)
}

read_taxa_metadata <- function(project) {
  taxa_info_path <- file.path("data", "projects", project, paste0(project, "_taxa_info.csv"))
  if (!file.exists(taxa_info_path)) {
    stop("Expected taxa info file was not created: ", taxa_info_path, call. = FALSE)
  }

  taxa_info <- read.csv2(taxa_info_path, stringsAsFactors = FALSE)
  required <- c("acceptedScientificName", "acceptedTaxonKey")
  missing <- setdiff(required, names(taxa_info))
  if (length(missing) > 0) {
    stop("Taxa info file is missing column(s): ", paste(missing, collapse = ", "), call. = FALSE)
  }

  taxa_info <- unique(taxa_info[, required, drop = FALSE])
  taxa_info <- taxa_info[!is.na(taxa_info$acceptedScientificName) &
                           !is.na(taxa_info$acceptedTaxonKey), , drop = FALSE]

  if (nrow(taxa_info) == 0) {
    stop("Taxa info file contains no usable accepted species/taxon key rows.", call. = FALSE)
  }

  taxa_info$speciesName <- vapply(taxa_info$acceptedScientificName, species_model_name, character(1))
  taxa_info$base_dir <- file.path(
    "data", "projects", project,
    paste0(taxa_info$speciesName, "_", taxa_info$acceptedTaxonKey)
  )
  taxa_info$climate_qs <- file.path(
    taxa_info$base_dir, "Climate",
    paste0("Climate_model_", taxa_info$speciesName, "_", taxa_info$acceptedTaxonKey, ".qs")
  )
  taxa_info$habitat_qs <- file.path(
    taxa_info$base_dir, "Habitat",
    paste0("Habitat_model_", taxa_info$speciesName, "_", taxa_info$acceptedTaxonKey, ".qs")
  )

  taxa_info
}

taxa_summary_values <- function(taxa_info) {
  list(
    accepted_species = paste(unique(taxa_info$acceptedScientificName), collapse = "; "),
    accepted_taxonkey = paste(unique(taxa_info$acceptedTaxonKey), collapse = "; ")
  )
}

completed_outputs_ready <- function(registry_row) {
  if (!isTRUE(verify_completed_outputs)) {
    return(TRUE)
  }

  project_value <- registry_row$project[[1]]
  if (length(project_value) == 0 || is.na(project_value) || !nzchar(project_value)) {
    return(FALSE)
  }

  project_dir <- file.path("data", "projects", project_value)
  if (!dir.exists(project_dir)) {
    return(FALSE)
  }

  taxa_info <- tryCatch(read_taxa_metadata(project_value), error = function(e) NULL)
  if (is.null(taxa_info)) {
    return(FALSE)
  }

  climate_exists <- any(file.exists(taxa_info$climate_qs))
  habitat_exists <- any(file.exists(taxa_info$habitat_qs))
  outcome_value <- registry_row$outcome[[1]]
  outcome_value <- ifelse(is.na(outcome_value), "", outcome_value)

  if (identical(outcome_value, "full_model")) {
    return(climate_exists && habitat_exists)
  }
  if (identical(outcome_value, "climate_only")) {
    return(climate_exists)
  }

  climate_exists || habitat_exists
}

file_ready <- function(path) {
  if (length(path) == 0 || is.null(path) || is.na(path) || !nzchar(path)) {
    return(FALSE)
  }
  file.exists(path) && isTRUE(file.info(path)$size > 0)
}

all_files_ready <- function(paths) {
  paths <- unique(paths[!is.na(paths) & nzchar(paths)])
  length(paths) > 0 && all(vapply(paths, file_ready, logical(1)))
}

project_taxa_info_path <- function(project) {
  file.path("data", "projects", project, paste0(project, "_taxa_info.csv"))
}

project_occurrences_path <- function(project) {
  file.path("data", "projects", project, paste0(project, "_processed_occurrences.qs"))
}

climate_core_paths <- function(taxa_info) {
  paths <- character()
  for (row_index in seq_len(nrow(taxa_info))) {
    species_name <- taxa_info$speciesName[[row_index]]
    taxon_key <- taxa_info$acceptedTaxonKey[[row_index]]
    base_dir <- taxa_info$base_dir[[row_index]]
    base_file <- paste0(species_name, "_Climate_")

    paths <- c(
      paths,
      file.path(base_dir, "Climate", paste0("Climate_model_", species_name, "_", taxon_key, ".qs")),
      file.path(base_dir, "Climate", "Current", "Predictions", "Rasters", paste0(base_file, "current_ensemble.tif")),
      file.path(base_dir, "Climate", "Current", "Interim", paste0(base_file, "current_ensemble_mean.tif")),
      file.path(base_dir, "Climate", "Current", "Diagnostics", "Confidence_maps", "Rasters", paste0(base_file, "current_ensemble_SD.tif"))
    )

    for (period in c("2041-2070", "2071-2100")) {
      for (scenario in c("ssp126", "ssp370", "ssp585")) {
        paths <- c(
          paths,
          file.path(base_dir, "Climate", period, scenario, "Predictions", "Rasters", paste0(base_file, period, "_", scenario, "_ensemble.tif")),
          file.path(base_dir, "Climate", "Current", "Interim", paste0(base_file, period, "_", scenario, "_ensemble_mean.tif")),
          file.path(base_dir, "Climate", period, scenario, "Diagnostics", "Confidence_maps", "Rasters", paste0(base_file, period, "_", scenario, "_ensemble_SD.tif"))
        )
      }
    }
  }

  unique(paths)
}

stage_02_outputs_ready <- function(project) {
  all_files_ready(c(project_occurrences_path(project), project_taxa_info_path(project))) &&
    !is.null(tryCatch(read_taxa_metadata(project), error = function(e) NULL))
}

stage_03_outputs_ready <- function(project) {
  if (!stage_02_outputs_ready(project)) {
    return(FALSE)
  }

  taxa_info <- tryCatch(read_taxa_metadata(project), error = function(e) NULL)
  !is.null(taxa_info) && all_files_ready(climate_core_paths(taxa_info))
}

stage_04_outputs_ready <- function(project) {
  taxa_info <- tryCatch(read_taxa_metadata(project), error = function(e) NULL)
  !is.null(taxa_info) && all_files_ready(taxa_info$habitat_qs)
}

stage_04_climate_only_ready <- function(project) {
  taxa_info <- tryCatch(read_taxa_metadata(project), error = function(e) NULL)
  !is.null(taxa_info) && !any(file.exists(taxa_info$habitat_qs))
}

stage_05_outputs_ready <- function(project) {
  file_ready(file.path("data", "projects", project, "Model_validation", "Validation_summary.csv"))
}

stage_outputs_ready <- function(stage_id, species_row) {
  project <- species_row$project[[1]]
  if (identical(stage_id, "02")) {
    return(stage_02_outputs_ready(project))
  }
  if (identical(stage_id, "03")) {
    return(stage_03_outputs_ready(project))
  }
  if (identical(stage_id, "04")) {
    return(stage_04_outputs_ready(project))
  }
  if (identical(stage_id, "05")) {
    return(stage_05_outputs_ready(project))
  }

  FALSE
}

latest_stage_event <- function(species_key, stage_id) {
  events <- read_all_event_files(control_dir)
  if (nrow(events) == 0) {
    return(NULL)
  }

  events <- events[events$species_key == species_key & events$stage == stage_id, , drop = FALSE]
  if (nrow(events) == 0) {
    return(NULL)
  }

  events$row_order <- seq_len(nrow(events))
  events$timestamp_sort <- suppressWarnings(as.POSIXct(events$timestamp, tz = Sys.timezone()))
  events$timestamp_sort[is.na(events$timestamp_sort)] <- as.POSIXct("1970-01-01", tz = "UTC")
  events <- events[order(events$timestamp_sort, events$row_order), , drop = FALSE]
  events[nrow(events), event_columns, drop = FALSE]
}

stage_has_success_event <- function(stage_event) {
  !is.null(stage_event) && stage_event$state[[1]] %in% c("OK", "REUSE")
}

stage_has_climate_only_skip_event <- function(stage_event) {
  !is.null(stage_event) &&
    identical(stage_event$stage[[1]], "04") &&
    identical(stage_event$state[[1]], "SKIP") &&
    isTRUE(grepl(
      "no habitat model file; stage 05 will run climate validation only",
      stage_event$message[[1]],
      fixed = TRUE
    ))
}

can_reuse_stage <- function(stage_id, species_row, attempt) {
  if (!isTRUE(reuse_successful_stages)) {
    return(list(reuse = FALSE, reason = "reuse_successful_stages is FALSE.", mode = "none"))
  }
  if (is.na(attempt) || attempt <= 1L) {
    return(list(reuse = FALSE, reason = "Stage reuse only applies to retry attempts.", mode = "none"))
  }

  stage_event <- latest_stage_event(species_row$species_key[[1]], stage_id)
  if (identical(stage_id, "04") &&
      stage_has_climate_only_skip_event(stage_event) &&
      stage_04_climate_only_ready(species_row$project[[1]])) {
    return(list(
      reuse = TRUE,
      reason = "Reusing previous stage 04 climate-only decision; no habitat model file exists.",
      mode = "climate_only"
    ))
  }

  if (!stage_has_success_event(stage_event)) {
    return(list(reuse = FALSE, reason = "No prior OK/REUSE event for this stage.", mode = "none"))
  }
  if (!stage_outputs_ready(stage_id, species_row)) {
    return(list(reuse = FALSE, reason = "Required output file(s) for this stage are missing.", mode = "none"))
  }

  list(
    reuse = TRUE,
    reason = paste0("Reusing previous stage ", stage_id, " outputs from prior attempt."),
    mode = "outputs"
  )
}


#-------------------------------------------------------------------------------
# Stage execution
#-------------------------------------------------------------------------------
run_stage <- function(stage_id,
                      script_path,
                      project,
                      species_to_model,
                      species_slug,
                      attempt,
                      log_root) {
  if (!file.exists(script_path)) {
    stop("Stage script not found: ", script_path, call. = FALSE)
  }

  stage_label <- paste0(stage_id, "_", stage_labels[[stage_id]])
  log_dir <- file.path(log_root, species_slug)
  safe_dir_create(log_dir)
  log_file <- file.path(
    log_dir,
    paste0(stage_label, "_attempt_", attempt, "_", timestamp_file(), ".log")
  )

  assign("project", project, envir = .GlobalEnv)
  assign("species_to_model", species_to_model, envir = .GlobalEnv)

  stage_env <- new.env(parent = .GlobalEnv)
  assign("project", project, envir = stage_env)
  assign("species_to_model", species_to_model, envir = stage_env)

  output_sink_count <- sink.number()
  message_sink_count <- sink.number(type = "message")
  log_connection <- file(log_file, open = "at", encoding = "UTF-8")

  result <- tryCatch(
    {
      sink(log_connection, split = FALSE)
      sink(log_connection, type = "message")

      cat("wiSDM chunk runner stage log\n")
      cat("Timestamp: ", timestamp_now(), "\n", sep = "")
      cat("Run id: ", run_id, "\n", sep = "")
      cat("Process id: ", Sys.getpid(), "\n", sep = "")
      cat("Chunk: ", chunk_id, "\n", sep = "")
      cat("Project: ", project, "\n", sep = "")
      cat("Species: ", paste(species_to_model, collapse = "; "), "\n", sep = "")
      cat("Stage: ", stage_id, " - ", basename(script_path), "\n\n", sep = "")

      source(script_path, local = stage_env, echo = FALSE, chdir = FALSE)

      cat("\nStage completed at ", timestamp_now(), "\n", sep = "")
      list(ok = TRUE, error_message = "", log_file = log_file)
    },
    error = function(e) {
      error_message <- conditionMessage(e)
      cat("\nStage failed at ", timestamp_now(), "\n", sep = "")
      cat("Error: ", error_message, "\n", sep = "")
      call <- conditionCall(e)
      if (!is.null(call)) {
        cat("Call: ", paste(deparse(call), collapse = " "), "\n", sep = "")
      }
      list(ok = FALSE, error_message = error_message, log_file = log_file)
    },
    finally = {
      while (sink.number(type = "message") > message_sink_count) {
        sink(type = "message")
      }
      while (sink.number() > output_sink_count) {
        sink()
      }
      close(log_connection)
    }
  )

  result$error_message <- clean_one_line(result$error_message)
  result
}


#-------------------------------------------------------------------------------
# Registry state transitions
#-------------------------------------------------------------------------------
claim_species <- function(species_row) {
  species_lock_path <- file.path(lock_dir, paste0("species_", species_row$species_key, ".lock"))
  species_lock <- lock_file(species_lock_path, timeout = species_lock_timeout)

  if (is.null(species_lock)) {
    return(list(
      claimed = FALSE,
      reason = "Species lock is held by another process.",
      species_lock = NULL,
      attempt = NA_integer_
    ))
  }

  decision <- with_locked_file(registry_lock_path, {
    registry <- read_registry_file(registry_path)
    current <- registry[registry$species_key == species_row$species_key, , drop = FALSE]
    current <- if (nrow(current) > 0) current[nrow(current), , drop = FALSE] else NULL

    if (!is.null(current) &&
        identical(current$status[[1]], "completed") &&
        completed_outputs_ready(current)) {
      list(claimed = FALSE, reason = "Species is already completed.", attempt = current$attempt[[1]])
    } else if (!is.null(current) && identical(current$status[[1]], "skipped") && !isTRUE(retry_skipped)) {
      list(claimed = FALSE, reason = "Species was previously skipped and retry_skipped is FALSE.", attempt = current$attempt[[1]])
    } else if (!is.null(current) && identical(current$status[[1]], "failed") && !isTRUE(retry_failed)) {
      list(claimed = FALSE, reason = "Species previously failed and retry_failed is FALSE.", attempt = current$attempt[[1]])
    } else {
      previous_attempt <- if (is.null(current)) 0L else as.integer(current$attempt[[1]])
      previous_attempt <- if (is.na(previous_attempt)) 0L else previous_attempt
      attempt <- previous_attempt + 1L

      row <- data.frame(
        species_requested = species_row$species_requested,
        species_key = species_row$species_key,
        chunk_id = chunk_id,
        project = species_row$project,
        attempt = attempt,
        status = "running",
        outcome = "",
        stage = "claim",
        started_at = timestamp_now(),
        updated_at = timestamp_now(),
        completed_at = "",
        pid = as.character(Sys.getpid()),
        run_id = run_id,
        accepted_species = "",
        accepted_taxonkey = "",
        error_message = "",
        stringsAsFactors = FALSE
      )
      row <- row[, registry_columns, drop = FALSE]

      registry <- registry[registry$species_key != species_row$species_key, , drop = FALSE]
      registry <- dplyr::bind_rows(registry, row)
      write_registry_file(registry_path, registry)

      reason <- ""
      if (!is.null(current) && identical(current$status[[1]], "running")) {
        reason <- "Previous registry state was running, but no active species lock was held; starting a new attempt."
      } else if (!is.null(current) && identical(current$status[[1]], "completed")) {
        reason <- "Previous registry state was completed, but expected output file(s) are missing; starting a new attempt."
      }

      list(claimed = TRUE, reason = reason, attempt = attempt)
    }
  }, timeout = registry_lock_timeout)

  if (!isTRUE(decision$claimed)) {
    unlock_file(species_lock)
    species_lock <- NULL
  }

  list(
    claimed = isTRUE(decision$claimed),
    reason = decision$reason,
    species_lock = species_lock,
    attempt = decision$attempt
  )
}

mark_stage_ok <- function(species_row, attempt, stage_id, message_text = "") {
  modify_registry_row(species_row$species_key, list(stage = stage_id, error_message = ""))
  append_event(
    species_requested = species_row$species_requested,
    species_key = species_row$species_key,
    project = species_row$project,
    attempt = attempt,
    stage = stage_id,
    state = "OK",
    message_text = message_text
  )
}

mark_stage_reused <- function(species_row, attempt, stage_id, message_text = "") {
  modify_registry_row(species_row$species_key, list(stage = stage_id, error_message = ""))
  append_event(
    species_requested = species_row$species_requested,
    species_key = species_row$species_key,
    project = species_row$project,
    attempt = attempt,
    stage = stage_id,
    state = "REUSE",
    message_text = message_text
  )
}

mark_species_failed <- function(species_row, attempt, stage_id, error_message, log_file = "") {
  modify_registry_row(
    species_row$species_key,
    list(
      status = "failed",
      outcome = "stage_error",
      stage = stage_id,
      completed_at = timestamp_now(),
      error_message = clean_one_line(error_message)
    )
  )
  append_event(
    species_requested = species_row$species_requested,
    species_key = species_row$species_key,
    project = species_row$project,
    attempt = attempt,
    stage = stage_id,
    state = "ERROR",
    message_text = error_message,
    log_file = log_file
  )
}

mark_species_skipped <- function(species_row, attempt, stage_id, outcome, message_text) {
  modify_registry_row(
    species_row$species_key,
    list(
      status = "skipped",
      outcome = outcome,
      stage = stage_id,
      completed_at = timestamp_now(),
      error_message = clean_one_line(message_text)
    )
  )
  append_event(
    species_requested = species_row$species_requested,
    species_key = species_row$species_key,
    project = species_row$project,
    attempt = attempt,
    stage = stage_id,
    state = "SKIP",
    message_text = message_text
  )
}

mark_species_completed <- function(species_row, attempt, outcome, taxa_info) {
  taxa_values <- taxa_summary_values(taxa_info)
  modify_registry_row(
    species_row$species_key,
    list(
      status = "completed",
      outcome = outcome,
      stage = "05",
      completed_at = timestamp_now(),
      accepted_species = taxa_values$accepted_species,
      accepted_taxonkey = taxa_values$accepted_taxonkey,
      error_message = ""
    )
  )
  append_event(
    species_requested = species_row$species_requested,
    species_key = species_row$species_key,
    project = species_row$project,
    attempt = attempt,
    stage = "complete",
    state = "DONE",
    message_text = outcome
  )
  append_completed_species(
    species_requested = species_row$species_requested,
    species_key = species_row$species_key,
    project = species_row$project,
    outcome = outcome
  )
}

handle_successful_stage_outputs <- function(species_row,
                                            attempt,
                                            stage_id,
                                            index,
                                            total,
                                            reused = FALSE,
                                            reuse_mode = "outputs") {
  if (identical(stage_id, "02")) {
    taxa_info <- read_taxa_metadata(species_row$project)
    taxa_values <- taxa_summary_values(taxa_info)
    modify_registry_row(
      species_row$species_key,
      list(
        accepted_species = taxa_values$accepted_species,
        accepted_taxonkey = taxa_values$accepted_taxonkey
      )
    )
  }

  if (identical(stage_id, "03")) {
    taxa_info <- read_taxa_metadata(species_row$project)
    if (!any(file.exists(taxa_info$climate_qs))) {
      reason <- "Stage 03 finished without creating a climate model file."
      mark_species_skipped(species_row, attempt, "03", "no_climate_model", reason)
      runner_message("SKIP", species_row$species_requested, species_row$project, "03", index, total, reason)
      return(FALSE)
    }
  }

  if (identical(stage_id, "04")) {
    taxa_info <- read_taxa_metadata(species_row$project)
    if (!any(file.exists(taxa_info$habitat_qs))) {
      if (isTRUE(reused) && identical(reuse_mode, "climate_only")) {
        return(TRUE)
      }

      reason <- "Stage 04 finished without creating a habitat model file; stage 05 will run climate validation only."
      append_event(
        species_requested = species_row$species_requested,
        species_key = species_row$species_key,
        project = species_row$project,
        attempt = attempt,
        stage = "04",
        state = "SKIP",
        message_text = reason
      )
      runner_message("SKIP", species_row$species_requested, species_row$project, "04", index, total, reason)
    }
  }

  TRUE
}


#-------------------------------------------------------------------------------
# Setup and species processing
#-------------------------------------------------------------------------------
config_value <- function(config_env, name, default = NULL) {
  if (exists(name, envir = config_env, inherits = FALSE)) {
    get(name, envir = config_env, inherits = FALSE)
  } else {
    default
  }
}

expected_setup_output_paths <- function(config_path = file.path("src", "00_configurations.R")) {
  config_env <- new.env(parent = baseenv())
  sys.source(config_path, envir = config_env)

  use_user_specific_climate <- !is.null(config_value(config_env, "user_specific_climate_data", NULL))
  use_user_specific_landcover <- !is.null(config_value(config_env, "user_specific_landcover_data", NULL))
  country_of_interest_value <- config_value(config_env, "country_of_interest", "Europe")
  custom_country_boundary_path_value <- config_value(config_env, "custom_country_boundary_path", NULL)

  paths <- character()

  if (!isTRUE(use_user_specific_climate)) {
    processed_folder <- file.path("data", "external", "climate", "chelsa_current", "processed")
    paths <- c(
      paths,
      file.path(processed_folder, "globalclimpreds.tif"),
      file.path(processed_folder, "globalclim_5k.tif"),
      file.path(processed_folder, "euclimpreds.tif")
    )

    if (tolower(as.character(country_of_interest_value)) != "europe" ||
        !is.null(custom_country_boundary_path_value)) {
      paths <- c(paths, file.path(processed_folder, "country_climpreds.tif"))
    }

    for (period in c("2041-2070", "2071-2100")) {
      for (scenario in c("ssp126", "ssp370", "ssp585")) {
        paths <- c(
          paths,
          file.path(
            "data", "external", "climate", "chelsa_future", "country",
            period, scenario, paste0(period, "_", scenario, "_masked.tif")
          )
        )
      }
    }
  }

  if (!isTRUE(use_user_specific_landcover)) {
    paths <- c(paths, file.path("data", "external", "habitat", "processed", "habitat_stack.tif"))
  }

  unique(paths)
}

missing_setup_output_paths <- function(paths) {
  if (length(paths) == 0) {
    return(character())
  }

  info <- file.info(paths)
  paths[!file.exists(paths) | is.na(info$size) | info$size <= 0]
}

format_setup_output_summary <- function(paths, max_paths = 8L) {
  if (length(paths) == 0) {
    return("none")
  }

  shown <- head(paths, max_paths)
  suffix <- if (length(paths) > max_paths) {
    paste0(" ... and ", length(paths) - max_paths, " more")
  } else {
    ""
  }

  paste0(paste(shown, collapse = "; "), suffix)
}

run_setup_once <- function() {
  if (!isTRUE(run_setup_01)) {
    runner_message("SKIP", NA_character_, project_prefix_safe, "01", NA_integer_, NA_integer_, "run_setup_01 is FALSE.")
    return(invisible(FALSE))
  }

  setup_marker <- file.path(control_dir, "setup_01_completed.txt")
  setup_outputs <- expected_setup_output_paths()
  missing_setup_outputs <- missing_setup_output_paths(setup_outputs)

  if (file.exists(setup_marker) && !isTRUE(force_setup_01) && length(missing_setup_outputs) == 0) {
    runner_message("SKIP", NA_character_, project_prefix_safe, "01", NA_integer_, NA_integer_, "Setup marker already exists.")
    return(invisible(FALSE))
  } else if (file.exists(setup_marker) && !isTRUE(force_setup_01)) {
    runner_message(
      "START",
      NA_character_,
      project_prefix_safe,
      "01",
      NA_integer_,
      NA_integer_,
      paste0(
        "Setup marker exists, but required setup output(s) are missing. Rerunning setup. Missing: ",
        format_setup_output_summary(missing_setup_outputs)
      )
    )
  }

  with_locked_file(setup_lock_path, {
    missing_setup_outputs <- missing_setup_output_paths(setup_outputs)

    if (file.exists(setup_marker) && !isTRUE(force_setup_01) && length(missing_setup_outputs) == 0) {
      runner_message(
        "SKIP",
        NA_character_,
        project_prefix_safe,
        "01",
        NA_integer_,
        NA_integer_,
        "Another job completed setup while this job was waiting."
      )
    } else {
      runner_message("START", NA_character_, project_prefix_safe, "01", NA_integer_, NA_integer_, "Running shared setup.")

      result <- run_stage(
        stage_id = "01",
        script_path = stage_scripts[["01"]],
        project = project_prefix_safe,
        species_to_model = character(0),
        species_slug = "setup",
        attempt = 1L,
        log_root = file.path(log_root, "setup")
      )

      if (!isTRUE(result$ok)) {
        runner_message("ERROR", NA_character_, project_prefix_safe, "01", NA_integer_, NA_integer_, result$error_message)
        stop("Stage 01 setup failed. See log: ", result$log_file, call. = FALSE)
      }

      missing_after_setup <- missing_setup_output_paths(setup_outputs)
      if (length(missing_after_setup) > 0) {
        message_text <- paste0(
          "Stage 01 setup finished, but required setup output(s) are still missing: ",
          format_setup_output_summary(missing_after_setup)
        )
        runner_message("ERROR", NA_character_, project_prefix_safe, "01", NA_integer_, NA_integer_, message_text)
        stop(message_text, call. = FALSE)
      }

      writeLines(
        c(
          paste0("completed_at=", timestamp_now()),
          paste0("run_id=", run_id),
          paste0("pid=", Sys.getpid()),
          paste0("log_file=", result$log_file)
        ),
        setup_marker
      )

      runner_message("OK", NA_character_, project_prefix_safe, "01", NA_integer_, NA_integer_, "Shared setup completed.")
    }

    invisible(TRUE)
  }, timeout = Inf)
}

process_species <- function(species_row, index, total) {
  claim <- claim_species(species_row)

  if (!isTRUE(claim$claimed)) {
    runner_message("SKIP", species_row$species_requested, species_row$project, "claim", index, total, claim$reason)
    append_event(
      species_requested = species_row$species_requested,
      species_key = species_row$species_key,
      project = species_row$project,
      attempt = ifelse(is.na(claim$attempt), "", claim$attempt),
      stage = "claim",
      state = "SKIP",
      message_text = claim$reason
    )
    return(invisible(FALSE))
  }

  on.exit(unlock_file(claim$species_lock), add = TRUE)
  attempt <- claim$attempt
  current_stage <- "claim"

  append_event(
    species_requested = species_row$species_requested,
    species_key = species_row$species_key,
    project = species_row$project,
    attempt = attempt,
    stage = "claim",
    state = "START",
    message_text = claim$reason
  )
  runner_message("START", species_row$species_requested, species_row$project, "claim", index, total, paste0("Attempt ", attempt, "."))

  tryCatch(
    {
      upstream_reuse_chain <- TRUE
      for (stage_id in c("02", "03", "04", "05")) {
        current_stage <- stage_id

        reuse_decision <- if (isTRUE(upstream_reuse_chain)) {
          can_reuse_stage(stage_id, species_row, attempt)
        } else {
          list(
            reuse = FALSE,
            reason = "An upstream stage was rerun in this attempt.",
            mode = "none"
          )
        }
        if (isTRUE(reuse_decision$reuse)) {
          mark_stage_reused(species_row, attempt, stage_id, reuse_decision$reason)
          runner_message(
            "REUSE",
            species_row$species_requested,
            species_row$project,
            stage_id,
            index,
            total,
            reuse_decision$reason
          )

          if (!handle_successful_stage_outputs(
            species_row = species_row,
            attempt = attempt,
            stage_id = stage_id,
            index = index,
            total = total,
            reused = TRUE,
            reuse_mode = reuse_decision$mode
          )) {
            return(invisible(FALSE))
          }

          next
        }

        upstream_reuse_chain <- FALSE

        runner_message(
          "START",
          species_row$species_requested,
          species_row$project,
          stage_id,
          index,
          total,
          basename(stage_scripts[[stage_id]])
        )
        append_event(
          species_requested = species_row$species_requested,
          species_key = species_row$species_key,
          project = species_row$project,
          attempt = attempt,
          stage = stage_id,
          state = "START",
          message_text = basename(stage_scripts[[stage_id]])
        )

        result <- run_stage(
          stage_id = stage_id,
          script_path = stage_scripts[[stage_id]],
          project = species_row$project,
          species_to_model = species_row$species_requested,
          species_slug = species_row$species_slug,
          attempt = attempt,
          log_root = file.path(log_root, paste0("chunk_", chunk_id))
        )

        if (!isTRUE(result$ok)) {
          mark_species_failed(species_row, attempt, stage_id, result$error_message, result$log_file)
          runner_message(
            "ERROR",
            species_row$species_requested,
            species_row$project,
            stage_id,
            index,
            total,
            paste0(result$error_message, " Log: ", result$log_file)
          )
          return(invisible(FALSE))
        }

        mark_stage_ok(species_row, attempt, stage_id, paste0("Log: ", result$log_file))
        runner_message(
          "OK",
          species_row$species_requested,
          species_row$project,
          stage_id,
          index,
          total,
          paste0("Log: ", result$log_file)
        )

        if (!handle_successful_stage_outputs(
          species_row = species_row,
          attempt = attempt,
          stage_id = stage_id,
          index = index,
          total = total
        )) {
          return(invisible(FALSE))
        }
      }

      taxa_info <- read_taxa_metadata(species_row$project)
      outcome <- if (any(file.exists(taxa_info$habitat_qs))) "full_model" else "climate_only"
      mark_species_completed(species_row, attempt, outcome, taxa_info)
      runner_message("DONE", species_row$species_requested, species_row$project, "complete", index, total, outcome)
      invisible(TRUE)
    },
    error = function(e) {
      mark_species_failed(species_row, attempt, current_stage, conditionMessage(e))
      runner_message("ERROR", species_row$species_requested, species_row$project, current_stage, index, total, conditionMessage(e))
      invisible(FALSE)
    }
  )
}

dry_run_report <- function(active_species_table) {
  registry <- read_registry_file(registry_path)
  runner_message(
    "START",
    NA_character_,
    chunk_prefix,
    "dry_run",
    NA_integer_,
    NA_integer_,
    paste0("Selected ", nrow(active_species_table), " species for chunk ", chunk_id, ".")
  )

  if (isTRUE(run_setup_01)) {
    setup_marker <- file.path(control_dir, "setup_01_completed.txt")
    setup_detail <- if (file.exists(setup_marker) && !force_setup_01) {
      "Stage 01 would be skipped because setup marker exists."
    } else {
      "Stage 01 would run under setup lock."
    }
    runner_message("SKIP", NA_character_, project_prefix_safe, "01", NA_integer_, NA_integer_, setup_detail)
  }

  for (i in seq_len(nrow(active_species_table))) {
    row <- active_species_table[i, , drop = FALSE]
    previous <- registry[registry$species_key == row$species_key, , drop = FALSE]
    previous <- if (nrow(previous) > 0) previous[nrow(previous), , drop = FALSE] else NULL

    decision <- "would run"
    if (!is.null(previous) &&
        identical(previous$status[[1]], "completed") &&
        completed_outputs_ready(previous)) {
      decision <- "would skip: completed"
    } else if (!is.null(previous) && identical(previous$status[[1]], "completed")) {
      decision <- "would rerun: completed registry row has missing output files"
    } else if (!is.null(previous) && identical(previous$status[[1]], "skipped") && !retry_skipped) {
      decision <- "would skip: previously skipped"
    } else if (!is.null(previous) && identical(previous$status[[1]], "failed") && retry_failed) {
      decision <- "would retry: previously failed"
    } else if (!is.null(previous) && identical(previous$status[[1]], "failed") && !retry_failed) {
      decision <- "would skip: previously failed"
    }

    previous_attempt <- if (is.null(previous)) 0L else suppressWarnings(as.integer(previous$attempt[[1]]))
    previous_attempt <- ifelse(is.na(previous_attempt), 0L, previous_attempt)
    if (isTRUE(reuse_successful_stages) && previous_attempt > 0L && grepl("^would (retry|rerun|run)", decision)) {
      candidate_attempt <- previous_attempt + 1L
      upstream_reuse_chain <- TRUE
      reusable_stages <- character()
      for (stage_id in c("02", "03", "04", "05")) {
        if (isTRUE(upstream_reuse_chain) && isTRUE(can_reuse_stage(stage_id, row, candidate_attempt)$reuse)) {
          reusable_stages <- c(reusable_stages, stage_id)
        } else {
          upstream_reuse_chain <- FALSE
        }
      }
      if (length(reusable_stages) > 0) {
        decision <- paste0(decision, "; reusable stages: ", paste(reusable_stages, collapse = ", "))
      }
    }

    runner_message("SKIP", row$species_requested, row$project, "dry_run", i, nrow(active_species_table), decision)
  }

  invisible(active_species_table)
}

run_runner_settings_preflight <- function() {
  if (length(n_blocks) != 1 || is.na(n_blocks) || n_blocks != as.integer(n_blocks) || n_blocks < 1) {
    stop("n_blocks must be a single positive integer.", call. = FALSE)
  }
  if (length(nr_active_block) != 1 ||
      is.na(nr_active_block) ||
      nr_active_block != as.integer(nr_active_block) ||
      nr_active_block < 1 ||
      nr_active_block > n_blocks) {
    stop("nr_active_block must be a single integer between 1 and n_blocks.", call. = FALSE)
  }
  if (length(project_prefix) != 1 || is.na(project_prefix) || !nzchar(project_prefix)) {
    stop("project_prefix must be a single non-empty string.", call. = FALSE)
  }
  if (length(max_project_name_chars) != 1 ||
      is.na(max_project_name_chars) ||
      max_project_name_chars != as.integer(max_project_name_chars) ||
      max_project_name_chars < 30) {
    stop("max_project_name_chars must be a single integer >= 30.", call. = FALSE)
  }
  if (!is.logical(retry_failed) || length(retry_failed) != 1 || is.na(retry_failed)) {
    stop("retry_failed must be TRUE or FALSE.", call. = FALSE)
  }
  if (!is.logical(retry_skipped) || length(retry_skipped) != 1 || is.na(retry_skipped)) {
    stop("retry_skipped must be TRUE or FALSE.", call. = FALSE)
  }
  if (!is.logical(verify_completed_outputs) ||
      length(verify_completed_outputs) != 1 ||
      is.na(verify_completed_outputs)) {
    stop("verify_completed_outputs must be TRUE or FALSE.", call. = FALSE)
  }
  if (!is.logical(reuse_successful_stages) ||
      length(reuse_successful_stages) != 1 ||
      is.na(reuse_successful_stages)) {
    stop("reuse_successful_stages must be TRUE or FALSE.", call. = FALSE)
  }
  if (!is.logical(run_setup_01) || length(run_setup_01) != 1 || is.na(run_setup_01)) {
    stop("run_setup_01 must be TRUE or FALSE.", call. = FALSE)
  }
  if (!is.logical(force_setup_01) || length(force_setup_01) != 1 || is.na(force_setup_01)) {
    stop("force_setup_01 must be TRUE or FALSE.", call. = FALSE)
  }
  if (!is.logical(dry_run) || length(dry_run) != 1 || is.na(dry_run)) {
    stop("dry_run must be TRUE or FALSE.", call. = FALSE)
  }

  invisible(TRUE)
}

run_stage_scripts_preflight <- function(stage_scripts) {
  missing_stages <- names(stage_scripts)[!file.exists(unlist(stage_scripts, use.names = FALSE))]
  if (length(missing_stages) > 0) {
    missing_lines <- paste0("  - ", missing_stages, ": ", unlist(stage_scripts[missing_stages], use.names = FALSE))
    stop(
      "Required wiSDM stage script(s) are missing:\n",
      paste(missing_lines, collapse = "\n"),
      call. = FALSE
    )
  }

  invisible(TRUE)
}

run_package_preflight <- function(required_packages) {
  still_missing <- setdiff(required_packages, rownames(installed.packages()))
  if (length(still_missing) > 0) {
    stop(
      "The following required package(s) are still missing after installation attempt: ",
      paste(still_missing, collapse = ", "),
      call. = FALSE
    )
  }

  invisible(TRUE)
}

run_config_value_preflight <- function(config_path, config_guard, dry_run) {
  if (identical(config_guard, "ignore")) {
    return(invisible(FALSE))
  }

  config_env <- new.env(parent = baseenv())
  sys.source(config_path, envir = config_env)
  problems <- character()

  add_problem <- function(text) {
    problems <<- c(problems, text)
  }

  check_choice <- function(name, allowed) {
    if (!exists(name, envir = config_env, inherits = FALSE)) {
      return(invisible(NULL))
    }
    value <- get(name, envir = config_env, inherits = FALSE)
    if (length(value) != 1 || is.na(value) || !tolower(as.character(value)) %in% allowed) {
      add_problem(paste0(name, " must be one of: ", paste(allowed, collapse = ", "), "."))
    }
  }

  check_choice("occurrence_thinning_method", c("random", "kmeans_clustering"))
  check_choice("pseudoabsence_thinning_method", c("random", "kmeans_clustering"))
  check_choice("update_files", c("ask", "yes", "no"))
  check_choice("workflow", c("single_step", "two_step"))

  if (exists("mtp_probabilities", envir = config_env, inherits = FALSE)) {
    mtp_probabilities_value <- get("mtp_probabilities", envir = config_env, inherits = FALSE)
    if (!is.numeric(mtp_probabilities_value) ||
        any(is.na(mtp_probabilities_value)) ||
        any(mtp_probabilities_value <= 0 | mtp_probabilities_value >= 1)) {
      add_problem("mtp_probabilities must be numeric values greater than 0 and lower than 1.")
    }
  }

  if (exists("boyce_background_size", envir = config_env, inherits = FALSE)) {
    boyce_background_size_value <- get("boyce_background_size", envir = config_env, inherits = FALSE)
    if (length(boyce_background_size_value) != 1 ||
        is.na(boyce_background_size_value) ||
        boyce_background_size_value != as.integer(boyce_background_size_value) ||
        boyce_background_size_value < 1) {
      add_problem("boyce_background_size must be a single positive integer.")
    }
  }

  if (exists("habitat_filter_near_zero_variance_predictors", envir = config_env, inherits = FALSE)) {
    habitat_nzv_value <- get("habitat_filter_near_zero_variance_predictors", envir = config_env, inherits = FALSE)
    if (!is.logical(habitat_nzv_value) || length(habitat_nzv_value) != 1 || is.na(habitat_nzv_value)) {
      add_problem("habitat_filter_near_zero_variance_predictors must be TRUE or FALSE.")
    }
  }

  if (!isTRUE(dry_run) &&
      exists("update_files", envir = config_env, inherits = FALSE) &&
      identical(tolower(as.character(get("update_files", envir = config_env, inherits = FALSE))), "ask")) {
    add_problem("update_files is set to 'ask', which can block unattended RStudio background jobs. Use 'yes' or 'no'.")
  }

  if (!isTRUE(dry_run)) {
    required_vars <- c("GBIF_USER", "GBIF_PWD", "GBIF_EMAIL")
    missing_vars <- required_vars[!nzchar(Sys.getenv(required_vars))]
    if (length(missing_vars) > 0) {
      add_problem(paste0(
        "Missing GBIF environment variable(s): ",
        paste(missing_vars, collapse = ", "),
        ". Add them to ~/.Renviron before running real chunk jobs."
      ))
    }
  }

  if (length(problems) == 0) {
    return(invisible(TRUE))
  }

  message_text <- paste0(
    "The following configuration value issue(s) were found in ",
    config_path,
    ":\n  - ",
    paste(problems, collapse = "\n  - ")
  )

  handle_preflight_issue(message_text, config_guard)
}

run_species_catalog_preflight <- function(species_catalog, active_species_table) {
  if (anyDuplicated(species_catalog$species_key)) {
    duplicated_keys <- unique(species_catalog$species_key[duplicated(species_catalog$species_key)])
    stop(
      "Duplicate species keys were generated after sanitizing species names: ",
      paste(duplicated_keys, collapse = ", "),
      call. = FALSE
    )
  }

  if (anyDuplicated(species_catalog$project)) {
    duplicated_projects <- unique(species_catalog$project[duplicated(species_catalog$project)])
    stop(
      "Duplicate project names were generated after sanitizing species names: ",
      paste(duplicated_projects, collapse = ", "),
      call. = FALSE
    )
  }

  too_long_projects <- species_catalog$project[nchar(species_catalog$project) > max_project_name_chars]
  if (length(too_long_projects) > 0) {
    stop(
      "Generated project name(s) exceed max_project_name_chars: ",
      paste(too_long_projects, collapse = ", "),
      call. = FALSE
    )
  }

  if (nrow(active_species_table) == 0) {
    warning(
      "No species are assigned to chunk ",
      chunk_id,
      ". The job will finish without running models.",
      call. = FALSE
    )
    return(invisible(FALSE))
  }

  invisible(TRUE)
}

run_control_dir_preflight <- function(paths) {
  for (path in paths) {
    safe_dir_create(path)
    if (!dir.exists(path)) {
      stop("Could not create required control directory: ", path, call. = FALSE)
    }
    if (file.access(path, mode = 2) != 0) {
      stop("Required control directory is not writable: ", path, call. = FALSE)
    }
  }

  invisible(TRUE)
}


#-------------------------------------------------------------------------------
# Main
#-------------------------------------------------------------------------------
preflight_settings_ok <- run_runner_settings_preflight()
preflight_packages_ok <- run_package_preflight(required_packages)
preflight_stage_scripts_ok <- run_stage_scripts_preflight(stage_scripts)
preflight_config_values_ok <- run_config_value_preflight(
  config_path = file.path("src", "00_configurations.R"),
  config_guard = config_guard,
  dry_run = dry_run
)

project_prefix_safe <- sanitize_token(project_prefix, fallback = "project")
chunk_id <- stringr::str_pad(nr_active_block, width = 2, pad = "0")
chunk_prefix <- paste0(project_prefix_safe, "_", chunk_id)
run_id <- paste0("chunk_", chunk_id, "_pid_", Sys.getpid(), "_", timestamp_file())

control_dir <- file.path("data", "projects", paste0(project_prefix_safe, "_chunk_control"))
lock_dir <- file.path(control_dir, "locks")
log_root <- file.path(control_dir, "logs")
preflight_control_dir_ok <- run_control_dir_preflight(c(control_dir, lock_dir, log_root))

registry_path <- file.path(control_dir, "species_registry.csv")
chunk_events_path <- file.path(control_dir, paste0("chunk_", chunk_id, "_events.csv"))
chunk_completed_path <- file.path(control_dir, paste0("chunk_", chunk_id, "_completed_species.csv"))

setup_lock_path <- file.path(lock_dir, "setup_01.lock")
registry_lock_path <- file.path(lock_dir, "species_registry.lock")
chunk_events_lock_path <- file.path(lock_dir, paste0("chunk_", chunk_id, "_events.lock"))
chunk_completed_lock_path <- file.path(lock_dir, paste0("chunk_", chunk_id, "_completed.lock"))

if (!identical(project_prefix_safe, project_prefix)) {
  message("Project prefix was sanitized from '", project_prefix, "' to '", project_prefix_safe, "'.")
}

all_species <- load_species_list(
  path = species_list_path,
  species_column = species_column,
  filter_column = filter_column,
  filter_value = filter_value
)

chunk_assignments <- make_chunk_assignments(length(all_species), n_blocks)
species_slugs <- make_unique_slugs(all_species)
species_keys <- mapply(species_key_from_slug, species_slugs, all_species, USE.NAMES = FALSE)
projects <- mapply(
  trim_with_hash,
  prefix = chunk_prefix,
  slug = species_slugs,
  original = all_species,
  MoreArgs = list(max_chars = max_project_name_chars),
  USE.NAMES = FALSE
)

species_catalog <- data.frame(
  species_requested = all_species,
  species_slug = species_slugs,
  species_key = species_keys,
  assigned_chunk = chunk_assignments,
  project = projects,
  stringsAsFactors = FALSE
)

active_species_table <- species_catalog[species_catalog$assigned_chunk == nr_active_block, , drop = FALSE]
preflight_species_catalog_ok <- run_species_catalog_preflight(species_catalog, active_species_table)

preflight_all_ok <- all(vapply(
  list(
    preflight_config_ok,
    preflight_paths_ok,
    preflight_settings_ok,
    preflight_packages_ok,
    preflight_stage_scripts_ok,
    preflight_config_values_ok,
    preflight_control_dir_ok,
    preflight_species_catalog_ok
  ),
  isTRUE,
  logical(1)
))

if (preflight_all_ok) {
  preflight_pass_message()
}

message(
  "\nwiSDM chunk runner\n",
  "Run id: ", run_id, "\n",
  "Chunk: ", chunk_id, " of ", n_blocks, "\n",
  "Species selected for this chunk: ", nrow(active_species_table), "\n",
  "Control directory: ", normalizePath(control_dir, winslash = "/", mustWork = FALSE), "\n"
)

if (nrow(active_species_table) == 0) {
  runner_message("DONE", NA_character_, chunk_prefix, "chunk", NA_integer_, NA_integer_, "No species assigned to this chunk.")
} else if (isTRUE(dry_run)) {
  dry_run_report(active_species_table)
} else {
  run_setup_once()

  for (i in seq_len(nrow(active_species_table))) {
    process_species(active_species_table[i, , drop = FALSE], index = i, total = nrow(active_species_table))
  }

  runner_message("DONE", NA_character_, chunk_prefix, "chunk", NA_integer_, NA_integer_, "Chunk job finished.")
}
