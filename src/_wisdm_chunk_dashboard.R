#-------------------------------------------------------------------------------
# wiSDM chunk-runner dashboard
#-------------------------------------------------------------------------------
# Prefer launching with:
#   source("src/_launch_wisdm_chunk_dashboard.R")
#
# The dashboard is designed to be open while chunk-worker jobs are running. Use
# the detached launcher above for live monitoring; avoid launching this file
# directly with RStudio's "Run App" button during live runs because RStudio's
# console Stop action can affect unrelated background jobs on some setups.
#
# The app reads the control folder created by src/_run_wisdm_by_chunks.R:
#   data/projects/<project_prefix>_chunk_control/
#
# It does not run models or acquire runner lock files. The progress tabs only
# read registry, event, completed-species, and log files, plus local PID state.
# The Imports / Merge tab can explicitly build a separate master control folder
# and adopt imported species project folders after source jobs have stopped.

required_packages <- c(
  "shiny", "DT", "dplyr", "readr", "tidyr", "stringr", "ggplot2",
  "ps", "filelock"
)

installed_packages <- rownames(installed.packages())
for (package in required_packages) {
  if (!package %in% installed_packages) {
    install.packages(package)
  }
  suppressPackageStartupMessages(library(package, character.only = TRUE))
}

`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}


#-------------------------------------------------------------------------------
# Paths and file readers
#-------------------------------------------------------------------------------
find_repo_root <- function(start = getwd()) {
  current <- normalizePath(start, winslash = "/", mustWork = FALSE)

  repeat {
    runner_path <- file.path(current, "src", "_run_wisdm_by_chunks.R")
    if (file.exists(runner_path)) {
      return(current)
    }

    parent <- dirname(current)
    if (identical(parent, current)) {
      break
    }
    current <- parent
  }

  normalizePath(getwd(), winslash = "/", mustWork = FALSE)
}

repo_root <- find_repo_root()

registry_columns <- c(
  "species_requested", "species_key", "chunk_id", "project", "attempt",
  "status", "outcome", "stage", "started_at", "updated_at", "completed_at",
  "pid", "run_id", "accepted_species", "accepted_taxonkey", "error_message",
  "source_snapshot", "source_file"
)

event_columns <- c(
  "timestamp", "run_id", "pid", "chunk_id", "species_requested",
  "species_key", "project", "attempt", "stage", "state", "message", "log_file",
  "source_snapshot", "source_file"
)

completed_columns <- c(
  "species_requested", "species_key", "project", "outcome", "completed_at",
  "run_id", "source_snapshot", "source_file"
)

empty_table <- function(columns) {
  as.data.frame(setNames(rep(list(character()), length(columns)), columns))
}

empty_registry <- function() {
  out <- empty_table(registry_columns)
  out$attempt <- integer()
  out
}

empty_events <- function() {
  out <- empty_table(event_columns)
  out$log_file_resolved <- character()
  out
}

empty_completed <- function() {
  empty_table(completed_columns)
}

has_content <- function(path) {
  file.exists(path) && isTRUE(file.info(path)$size > 0)
}

read_csv_safe <- function(path, columns, integer_columns = character()) {
  if (!has_content(path)) {
    out <- empty_table(columns)
    for (column in integer_columns) {
      out[[column]] <- integer()
    }
    return(out)
  }

  out <- suppressWarnings(readr::read_csv(
    path,
    col_types = readr::cols(.default = readr::col_character()),
    show_col_types = FALSE
  ))
  out <- as.data.frame(out, stringsAsFactors = FALSE)

  for (column in setdiff(columns, names(out))) {
    out[[column]] <- NA_character_
  }
  out <- out[, columns, drop = FALSE]

  for (column in integer_columns) {
    out[[column]] <- suppressWarnings(as.integer(out[[column]]))
  }

  out
}

list_control_dirs <- function(root = repo_root) {
  projects_dir <- file.path(root, "data", "projects")
  if (!dir.exists(projects_dir)) {
    return(character())
  }

  dirs <- list.dirs(projects_dir, full.names = TRUE, recursive = FALSE)
  dirs <- dirs[grepl("_chunk_control$", basename(dirs))]
  normalizePath(dirs, winslash = "/", mustWork = FALSE)
}

path_label <- function(path) {
  if (!nzchar(path)) {
    return("")
  }
  basename(path)
}

named_choices <- function(values, labels = values) {
  if (length(values) == 0) {
    return(character())
  }

  stats::setNames(values, labels)
}

resolve_existing_path <- function(path, control_dir = NULL, root = repo_root) {
  if (length(path) == 0 || is.na(path) || !nzchar(path)) {
    return(NA_character_)
  }

  candidates <- unique(c(
    path,
    file.path(root, path),
    if (!is.null(control_dir)) file.path(control_dir, path) else character()
  ))

  existing <- candidates[file.exists(candidates)]
  if (length(existing) == 0) {
    return(normalizePath(candidates[[1]], winslash = "/", mustWork = FALSE))
  }

  normalizePath(existing[[1]], winslash = "/", mustWork = TRUE)
}

read_registry <- function(control_dir) {
  path <- file.path(control_dir, "species_registry.csv")
  read_csv_safe(path, registry_columns, integer_columns = "attempt")
}

read_events <- function(control_dir) {
  files <- list.files(
    control_dir,
    pattern = "^chunk_[0-9]+_events\\.csv$",
    full.names = TRUE
  )

  if (length(files) == 0) {
    return(empty_events())
  }

  events <- lapply(files, read_csv_safe, columns = event_columns)
  events <- dplyr::bind_rows(events)

  events$log_file_resolved <- vapply(
    events$log_file,
    resolve_existing_path,
    character(1),
    control_dir = control_dir
  )

  events
}

read_completed <- function(control_dir) {
  files <- list.files(
    control_dir,
    pattern = "^chunk_[0-9]+_completed_species\\.csv$",
    full.names = TRUE
  )

  if (length(files) == 0) {
    return(empty_completed())
  }

  completed <- lapply(files, read_csv_safe, columns = completed_columns)
  dplyr::bind_rows(completed)
}

empty_log_info <- function() {
  data.frame(
    path = character(),
    chunk = character(),
    species_slug = character(),
    folder = character(),
    stage = character(),
    attempt = integer(),
    modified = as.POSIXct(character()),
    size_mb = numeric(),
    name = character(),
    stringsAsFactors = FALSE
  )
}

build_log_info <- function(log_files) {
  if (length(log_files) == 0) {
    return(empty_log_info())
  }

  normalized_paths <- normalizePath(log_files, winslash = "/", mustWork = FALSE)
  log_stat <- file.info(log_files)
  log_names <- basename(log_files)
  log_stems <- sub("\\.log$", "", log_names)
  stages <- ifelse(
    grepl("_attempt_[0-9]+_", log_stems),
    sub("_attempt_[0-9]+_.*$", "", log_stems),
    log_stems
  )
  attempts <- ifelse(
    grepl("_attempt_[0-9]+_", log_stems),
    sub("^.*_attempt_([0-9]+)_.*$", "\\1", log_stems),
    NA_character_
  )

  data.frame(
    path = normalized_paths,
    chunk = basename(dirname(dirname(normalized_paths))),
    species_slug = basename(dirname(normalized_paths)),
    folder = basename(dirname(normalized_paths)),
    stage = stages,
    attempt = suppressWarnings(as.integer(attempts)),
    modified = log_stat$mtime,
    size_mb = round(log_stat$size / 1024^2, 3),
    name = log_names,
    stringsAsFactors = FALSE
  )
}

parse_dashboard_time <- function(x) {
  if (inherits(x, "POSIXct")) {
    return(x)
  }
  if (length(x) == 0 || is.null(x) || is.na(x) || !nzchar(as.character(x))) {
    return(as.POSIXct(NA))
  }

  out <- suppressWarnings(as.POSIXct(as.character(x), tz = Sys.timezone()))
  if (is.na(out)) {
    out <- suppressWarnings(readr::parse_datetime(as.character(x)))
  }
  out
}

pid_state <- function(pid, started_at) {
  if (length(pid) == 0 || is.null(pid) || is.na(pid) || !nzchar(as.character(pid))) {
    return("missing")
  }

  pid_number <- suppressWarnings(as.integer(pid))
  if (is.na(pid_number) || pid_number <= 0) {
    return("missing")
  }

  handle <- tryCatch(ps::ps_handle(pid_number), error = function(e) NULL)
  if (is.null(handle)) {
    return("inactive")
  }

  is_running <- tryCatch(ps::ps_is_running(handle), error = function(e) FALSE)
  if (!isTRUE(is_running)) {
    return("inactive")
  }

  process_started_at <- tryCatch(ps::ps_create_time(handle), error = function(e) as.POSIXct(NA))
  registry_started_at <- parse_dashboard_time(started_at)

  if (!is.na(process_started_at) &&
      !is.na(registry_started_at) &&
      as.numeric(process_started_at) > as.numeric(registry_started_at)) {
    return("reused")
  }

  "active"
}

species_lock_path <- function(control_dir, species_key) {
  file.path(control_dir, "locks", paste0("species_", species_key, ".lock"))
}

probe_species_lock_state <- function(control_dir, species_key) {
  "not_checked"
}

stale_reason <- function(pid_state_value, lock_state_value) {
  reasons <- character()

  if (pid_state_value %in% c("inactive", "missing", "reused")) {
    reasons <- c(reasons, paste0("pid_", pid_state_value))
  }
  if (lock_state_value %in% c("free", "missing")) {
    reasons <- c(reasons, paste0("lock_", lock_state_value))
  }

  paste(reasons, collapse = "; ")
}

add_effective_status <- function(registry, control_dir) {
  row_count <- nrow(registry)
  if (!"status" %in% names(registry)) {
    registry$status <- rep(NA_character_, row_count)
  }

  registry$registry_status <- as.character(registry$status)
  registry$effective_status <- as.character(registry$status)
  registry$pid_state <- rep("", row_count)
  registry$lock_state <- rep("", row_count)
  registry$stale_reason <- rep("", row_count)

  if (row_count == 0) {
    return(registry)
  }

  running_rows <- which(!is.na(registry$registry_status) & registry$registry_status == "running")
  for (row_index in running_rows) {
    current_pid_state <- pid_state(registry$pid[[row_index]], registry$started_at[[row_index]])
    current_lock_state <- probe_species_lock_state(control_dir, registry$species_key[[row_index]])
    current_stale_reason <- stale_reason(current_pid_state, current_lock_state)

    registry$pid_state[[row_index]] <- current_pid_state
    registry$lock_state[[row_index]] <- current_lock_state
    registry$stale_reason[[row_index]] <- current_stale_reason

    if (nzchar(current_stale_reason)) {
      registry$effective_status[[row_index]] <- "stale_running"
    } else if (identical(current_pid_state, "active")) {
      registry$effective_status[[row_index]] <- "active_running"
    } else {
      registry$effective_status[[row_index]] <- "running_unknown"
    }
  }

  registry
}

add_model_status <- function(registry) {
  row_count <- nrow(registry)
  if (!"outcome" %in% names(registry)) {
    registry$outcome <- rep(NA_character_, row_count)
  }
  if (!"effective_status" %in% names(registry)) {
    registry$effective_status <- rep(NA_character_, row_count)
  }
  if (!"registry_status" %in% names(registry)) {
    registry$registry_status <- as.character(registry$status)
  }

  registry$model_status <- as.character(registry$effective_status)
  registry$model_status_note <- rep("", row_count)

  if (row_count == 0) {
    return(registry)
  }

  registry_status <- as.character(registry$registry_status)
  outcome <- as.character(registry$outcome)

  full_rows <- registry_status %in% "completed" & outcome %in% "full_model"
  partial_rows <- registry_status %in% "completed" & outcome %in% "climate_only"
  unknown_completed_rows <- registry_status %in% "completed" & !(outcome %in% c("full_model", "climate_only"))
  no_climate_rows <- registry_status %in% "skipped" & outcome %in% "no_climate_model"

  registry$model_status[full_rows] <- "full_model_completed"
  registry$model_status_note[full_rows] <- "All stages completed and both climate and habitat model files were found."

  registry$model_status[partial_rows] <- "partial_climate_only"
  registry$model_status_note[partial_rows] <- "Stage 05 completed, but no habitat model file was found after stage 04."

  registry$model_status[unknown_completed_rows] <- "completed_unknown"
  registry$model_status_note[unknown_completed_rows] <- "Registry status is completed, but outcome does not identify full_model or climate_only."

  registry$model_status[no_climate_rows] <- "no_climate_model"
  registry$model_status_note[no_climate_rows] <- "Stage 03 finished without creating a climate model file."

  registry
}

read_control_data <- function(control_dir) {
  if (is.null(control_dir) || !nzchar(control_dir) || !dir.exists(control_dir)) {
    return(list(
      control_dir = control_dir,
      registry = add_model_status(add_effective_status(empty_registry(), control_dir)),
      events = empty_events(),
      completed = empty_completed(),
      log_files = empty_log_info(),
      loaded_at = Sys.time(),
      valid = FALSE
    ))
  }

  log_files <- list.files(
    file.path(control_dir, "logs"),
    pattern = "\\.log$",
    full.names = TRUE,
    recursive = TRUE
  )

  log_info <- build_log_info(log_files)

  list(
    control_dir = control_dir,
    registry = add_model_status(add_effective_status(read_registry(control_dir), control_dir)),
    events = read_events(control_dir),
    completed = read_completed(control_dir),
    log_files = log_info,
    loaded_at = Sys.time(),
    valid = TRUE
  )
}


#-------------------------------------------------------------------------------
# Summary helpers
#-------------------------------------------------------------------------------
latest_event_states <- function(events) {
  if (nrow(events) == 0) {
    return(data.frame())
  }

  events %>%
    dplyr::mutate(
      timestamp_sort = suppressWarnings(as.POSIXct(timestamp)),
      timestamp_sort = dplyr::if_else(
        is.na(timestamp_sort),
        as.POSIXct("1970-01-01", tz = "UTC"),
        timestamp_sort
      )
    ) %>%
    dplyr::arrange(species_key, stage, timestamp_sort) %>%
    dplyr::group_by(species_key, stage) %>%
    dplyr::slice_tail(n = 1) %>%
    dplyr::ungroup() %>%
    dplyr::select(species_key, stage, state)
}

stage_matrix <- function(registry, events) {
  if (nrow(registry) == 0) {
    return(data.frame())
  }

  states <- latest_event_states(events)
  if (nrow(states) > 0) {
    states <- states %>%
      tidyr::pivot_wider(
        names_from = stage,
        values_from = state,
        values_fill = ""
      )
  }

  out <- registry %>%
    dplyr::select(
      species_requested, species_key, chunk_id, project, attempt,
      registry_status, effective_status, model_status, outcome, current_stage = stage
    )

  if (nrow(states) > 0) {
    out <- dplyr::left_join(out, states, by = "species_key")
  }

  stage_columns <- c("claim", "02", "03", "04", "05", "complete")
  for (column in setdiff(stage_columns, names(out))) {
    out[[column]] <- ""
  }

  out[, c(
    "chunk_id", "species_requested", "project", "attempt", "registry_status",
    "effective_status", "model_status", "outcome", "current_stage", stage_columns
  ), drop = FALSE]
}

filter_by_chunk <- function(data, chunk_id) {
  if (is.null(chunk_id) || !nzchar(chunk_id) || identical(chunk_id, "all")) {
    return(data)
  }

  data[data$chunk_id == chunk_id, , drop = FALSE]
}

summary_counts <- function(registry, events, completed) {
  claimed <- nrow(registry)
  full_model_n <- sum(registry$model_status == "full_model_completed", na.rm = TRUE)
  partial_n <- sum(registry$model_status == "partial_climate_only", na.rm = TRUE)
  completed_unknown_n <- sum(registry$model_status == "completed_unknown", na.rm = TRUE)
  no_climate_n <- sum(registry$model_status == "no_climate_model", na.rm = TRUE)
  failed_n <- sum(registry$model_status == "failed", na.rm = TRUE)
  active_running_n <- sum(registry$model_status == "active_running", na.rm = TRUE)
  stale_running_n <- sum(registry$model_status == "stale_running", na.rm = TRUE)
  running_unknown_n <- sum(registry$model_status == "running_unknown", na.rm = TRUE)
  skipped_n <- sum(registry$model_status == "skipped", na.rm = TRUE)
  event_n <- nrow(events)
  log_n <- if ("log_file_resolved" %in% names(events)) {
    sum(!is.na(events$log_file_resolved) & file.exists(events$log_file_resolved))
  } else {
    0L
  }

  out <- data.frame(
    label = c(
      "Claimed species", "Full models", "Partial: climate only", "Completed unknown",
      "No climate model", "Failed", "Active running", "Stale running",
      "Running unknown", "Other skipped", "Events", "Linked logs"
    ),
    value = c(
      claimed, full_model_n, partial_n, completed_unknown_n, no_climate_n, failed_n,
      active_running_n, stale_running_n, running_unknown_n, skipped_n, event_n, log_n
    ),
    tone = c(
      "neutral", "completed", "partial", "unknown", "blocked", "failed",
      "running", "stale", "unknown", "skipped", "neutral", "neutral"
    ),
    show_when_zero = c(
      TRUE, TRUE, TRUE, FALSE, TRUE, TRUE, TRUE, TRUE,
      FALSE, FALSE, TRUE, TRUE
    ),
    stringsAsFactors = FALSE
  )

  out <- out[out$show_when_zero | out$value > 0, , drop = FALSE]
  out$show_when_zero <- NULL
  out
}

status_palette <- c(
  completed = "#198754",
  full_model_completed = "#198754",
  partial_climate_only = "#b7791f",
  completed_unknown = "#20c997",
  no_climate_model = "#f59f00",
  running = "#0d6efd",
  active_running = "#0d6efd",
  stale_running = "#7c3aed",
  running_unknown = "#6c757d",
  failed = "#dc3545",
  skipped = "#f59f00",
  unknown = "#6c757d"
)

state_palette <- c(
  START = "#0d6efd",
  OK = "#198754",
  REUSE = "#0f766e",
  SKIP = "#f59f00",
  ERROR = "#dc3545",
  DONE = "#198754",
  unknown = "#6c757d"
)

state_badge <- function(x) {
  ifelse(
    is.na(x) | !nzchar(x),
    "",
    paste0("<span class='state-badge state-", tolower(x), "'>", x, "</span>")
  )
}

status_label <- function(x) {
  labels <- c(
    active_running = "active running",
    stale_running = "stale running",
    running_unknown = "running unknown",
    full_model_completed = "full model completed",
    partial_climate_only = "partial: climate only",
    completed_unknown = "completed unknown",
    no_climate_model = "no climate model"
  )
  out <- ifelse(x %in% names(labels), labels[x], gsub("_", " ", x))
  unname(out)
}

status_badge <- function(x) {
  class_name <- gsub("[^a-z0-9]+", "_", tolower(x))
  ifelse(
    is.na(x) | !nzchar(x),
    "",
    paste0("<span class='status-badge status-", class_name, "'>", status_label(x), "</span>")
  )
}

clean_log_text <- function(path, n_lines = 1000L) {
  if (length(path) == 0 || is.na(path) || !file.exists(path)) {
    return("Select an existing log file to preview it here.")
  }

  lines <- readLines(path, warn = FALSE, encoding = "UTF-8")
  if (length(lines) > n_lines) {
    lines <- tail(lines, n_lines)
    lines <- c(
      paste0("Showing last ", n_lines, " lines from: ", path),
      strrep("-", 80),
      lines
    )
  } else {
    lines <- c(paste0("Showing full log: ", path), strrep("-", 80), lines)
  }

  paste(lines, collapse = "\n")
}

open_local_path <- function(path) {
  if (length(path) == 0 || is.na(path) || !file.exists(path)) {
    return(FALSE)
  }

  normalized <- normalizePath(path, winslash = "\\", mustWork = TRUE)
  if (identical(.Platform$OS.type, "windows")) {
    shell.exec(normalized)
  } else {
    utils::browseURL(paste0("file://", normalizePath(path, mustWork = TRUE)))
  }

  TRUE
}


#-------------------------------------------------------------------------------
# Import snapshots and master merge helpers
#-------------------------------------------------------------------------------
projects_root <- file.path(repo_root, "data", "projects")

merge_manifest_columns <- c(
  "merge_id", "timestamp", "snapshot_id", "chunk_id", "species_key", "project",
  "source_control_dir", "source_project_dir", "destination_project_dir",
  "artifact_state", "action"
)

import_summary_columns <- c(
  "snapshot_id", "status", "chunks", "registry_rows", "event_rows",
  "completed_rows", "log_files", "projects_expected", "projects_in_import",
  "projects_missing", "projects_already_adopted", "project_collisions",
  "warnings", "snapshot_dir", "control_dir", "projects_dir"
)

import_warning_columns <- c("snapshot_id", "severity", "message")

merge_result_columns <- c("item", "value")

safe_dir_create <- function(path) {
  if (!dir.exists(path)) {
    dir.create(path, recursive = TRUE, showWarnings = FALSE)
  }
  invisible(path)
}

safe_write_csv <- function(data, path) {
  safe_dir_create(dirname(path))
  readr::write_csv(data, path, na = "")
  invisible(path)
}

timestamp_file <- function() {
  format(Sys.time(), "%Y%m%d_%H%M%S")
}

timestamp_now <- function() {
  format(Sys.time(), "%Y-%m-%d %H:%M:%S")
}

nonempty_unique <- function(x) {
  x <- as.character(x)
  sort(unique(x[!is.na(x) & nzchar(x)]))
}

collapse_values <- function(x, empty = "") {
  x <- nonempty_unique(x)
  if (length(x) == 0) {
    return(empty)
  }
  paste(x, collapse = ", ")
}

empty_import_scan <- function() {
  list(
    summary = as.data.frame(setNames(rep(list(character()), length(import_summary_columns)), import_summary_columns)),
    warnings = as.data.frame(setNames(rep(list(character()), length(import_warning_columns)), import_warning_columns))
  )
}

empty_merge_result <- function() {
  as.data.frame(setNames(rep(list(character()), length(merge_result_columns)), merge_result_columns))
}

control_prefix <- function(control_dir) {
  if (is.null(control_dir) || !nzchar(control_dir)) {
    return("")
  }
  base <- basename(normalizePath(control_dir, winslash = "/", mustWork = FALSE))
  base <- sub("_master_chunk_control$", "", base)
  sub("_chunk_control$", "", base)
}

runner_control_dir_for <- function(control_dir) {
  prefix <- control_prefix(control_dir)
  if (!nzchar(prefix)) {
    return("")
  }
  file.path(dirname(normalizePath(control_dir, winslash = "/", mustWork = FALSE)), paste0(prefix, "_chunk_control"))
}

master_control_dir_for <- function(control_dir) {
  prefix <- control_prefix(control_dir)
  if (!nzchar(prefix)) {
    return("")
  }
  file.path(dirname(normalizePath(control_dir, winslash = "/", mustWork = FALSE)), paste0(prefix, "_master_chunk_control"))
}

imports_root_for <- function(control_dir) {
  base_control <- runner_control_dir_for(control_dir)
  if (!nzchar(base_control)) {
    return("")
  }
  file.path(base_control, "imports")
}

snapshot_control_dir <- function(snapshot_dir) {
  candidates <- c(
    file.path(snapshot_dir, "control"),
    snapshot_dir
  )

  child_dirs <- list.dirs(snapshot_dir, full.names = TRUE, recursive = FALSE)
  candidates <- unique(c(
    candidates,
    child_dirs[grepl("_chunk_control$", basename(child_dirs))]
  ))

  matches <- candidates[file.exists(file.path(candidates, "species_registry.csv"))]
  if (length(matches) == 0) {
    return("")
  }
  matches[[1]]
}

snapshot_projects_dir <- function(snapshot_dir) {
  file.path(snapshot_dir, "projects")
}

snapshot_adoption_path <- function(snapshot_dir) {
  file.path(snapshot_dir, "merge_adopted_projects.csv")
}

read_generic_csv <- function(path, columns) {
  read_csv_safe(path, columns = columns)
}

read_merge_manifest <- function(master_dir) {
  read_generic_csv(file.path(master_dir, "merge_manifest.csv"), merge_manifest_columns)
}

read_snapshot_adoptions <- function(snapshot_dir) {
  read_generic_csv(snapshot_adoption_path(snapshot_dir), merge_manifest_columns)
}

project_already_adopted <- function(project, snapshot_id, master_dir, snapshot_dir) {
  if (!nzchar(project)) {
    return(FALSE)
  }

  manifest <- read_merge_manifest(master_dir)
  snapshot_manifest <- read_snapshot_adoptions(snapshot_dir)
  records <- dplyr::bind_rows(manifest, snapshot_manifest)
  if (nrow(records) == 0) {
    return(FALSE)
  }

  any(
    records$project == project &
      records$snapshot_id == snapshot_id &
      records$artifact_state %in% c("adopted", "already_adopted"),
    na.rm = TRUE
  )
}

read_import_events <- function(control_dir, snapshot_id) {
  files <- list.files(
    control_dir,
    pattern = "^chunk_[0-9]+_events\\.csv$",
    full.names = TRUE
  )
  if (length(files) == 0) {
    return(empty_events()[, event_columns, drop = FALSE])
  }

  events <- lapply(files, function(path) {
    out <- read_csv_safe(path, event_columns)
    out$source_snapshot <- snapshot_id
    out$source_file <- normalizePath(path, winslash = "/", mustWork = FALSE)
    out
  })

  dplyr::bind_rows(events)
}

read_import_completed <- function(control_dir, snapshot_id) {
  files <- list.files(
    control_dir,
    pattern = "^chunk_[0-9]+_completed_species\\.csv$",
    full.names = TRUE
  )
  if (length(files) == 0) {
    return(empty_completed())
  }

  completed <- lapply(files, function(path) {
    out <- read_csv_safe(path, completed_columns)
    out$source_snapshot <- snapshot_id
    out$source_file <- normalizePath(path, winslash = "/", mustWork = FALSE)
    out
  })

  dplyr::bind_rows(completed)
}

read_snapshot_data <- function(snapshot_dir) {
  snapshot_dir <- normalizePath(snapshot_dir, winslash = "/", mustWork = FALSE)
  snapshot_id <- basename(snapshot_dir)
  control_dir <- snapshot_control_dir(snapshot_dir)
  projects_dir <- snapshot_projects_dir(snapshot_dir)

  if (!nzchar(control_dir)) {
    return(list(
      snapshot_id = snapshot_id,
      snapshot_dir = snapshot_dir,
      control_dir = "",
      projects_dir = projects_dir,
      registry = empty_registry(),
      events = empty_events()[, event_columns, drop = FALSE],
      completed = empty_completed(),
      log_files = character(),
      valid = FALSE,
      error = "No control folder with species_registry.csv was found."
    ))
  }

  registry_path <- file.path(control_dir, "species_registry.csv")
  registry <- read_csv_safe(registry_path, registry_columns, integer_columns = "attempt")
  registry$source_snapshot <- snapshot_id
  registry$source_file <- normalizePath(registry_path, winslash = "/", mustWork = FALSE)

  log_files <- list.files(
    file.path(control_dir, "logs"),
    pattern = "\\.log$",
    full.names = TRUE,
    recursive = TRUE
  )

  list(
    snapshot_id = snapshot_id,
    snapshot_dir = snapshot_dir,
    control_dir = normalizePath(control_dir, winslash = "/", mustWork = FALSE),
    projects_dir = normalizePath(projects_dir, winslash = "/", mustWork = FALSE),
    registry = registry,
    events = read_import_events(control_dir, snapshot_id),
    completed = read_import_completed(control_dir, snapshot_id),
    log_files = normalizePath(log_files, winslash = "/", mustWork = FALSE),
    valid = TRUE,
    error = ""
  )
}

scan_import_snapshots <- function(control_dir) {
  imports_root <- imports_root_for(control_dir)
  master_dir <- master_control_dir_for(control_dir)
  if (!nzchar(imports_root) || !dir.exists(imports_root)) {
    return(empty_import_scan())
  }

  snapshot_dirs <- list.dirs(imports_root, full.names = TRUE, recursive = FALSE)
  if (length(snapshot_dirs) == 0) {
    return(empty_import_scan())
  }

  summaries <- list()
  warnings <- list()
  duplicate_tracker <- data.frame(
    snapshot_id = character(),
    chunk_id = character(),
    species_key = character(),
    project = character(),
    stringsAsFactors = FALSE
  )

  for (snapshot_dir in snapshot_dirs) {
    snapshot <- read_snapshot_data(snapshot_dir)
    snapshot_warnings <- character()

    if (!isTRUE(snapshot$valid)) {
      snapshot_warnings <- c(snapshot_warnings, snapshot$error)
      warnings[[length(warnings) + 1L]] <- data.frame(
        snapshot_id = snapshot$snapshot_id,
        severity = "error",
        message = snapshot$error,
        stringsAsFactors = FALSE
      )
    }

    registry <- snapshot$registry
    expected_projects <- nonempty_unique(registry$project)
    import_projects <- if (dir.exists(snapshot$projects_dir)) {
      basename(list.dirs(snapshot$projects_dir, full.names = TRUE, recursive = FALSE))
    } else {
      character()
    }

    destination_exists <- file.exists(file.path(projects_root, expected_projects))
    already_adopted <- vapply(
      expected_projects,
      project_already_adopted,
      logical(1),
      snapshot_id = snapshot$snapshot_id,
      master_dir = master_dir,
      snapshot_dir = snapshot$snapshot_dir
    )
    in_import <- expected_projects %in% import_projects
    missing_projects <- expected_projects[!in_import & !destination_exists]
    project_collisions <- expected_projects[destination_exists & !already_adopted]

    if (length(missing_projects) > 0) {
      msg <- paste("Missing project folder(s):", collapse_values(missing_projects))
      snapshot_warnings <- c(snapshot_warnings, msg)
      warnings[[length(warnings) + 1L]] <- data.frame(
        snapshot_id = snapshot$snapshot_id,
        severity = "warning",
        message = msg,
        stringsAsFactors = FALSE
      )
    }

    if (length(project_collisions) > 0) {
      msg <- paste("Destination project folder already exists and is not recorded as adopted:", collapse_values(project_collisions))
      snapshot_warnings <- c(snapshot_warnings, msg)
      warnings[[length(warnings) + 1L]] <- data.frame(
        snapshot_id = snapshot$snapshot_id,
        severity = "error",
        message = msg,
        stringsAsFactors = FALSE
      )
    }

    if (nrow(registry) > 0) {
      duplicate_tracker <- dplyr::bind_rows(
        duplicate_tracker,
        registry[, c("source_snapshot", "chunk_id", "species_key", "project"), drop = FALSE] %>%
          dplyr::rename(snapshot_id = source_snapshot)
      )
    }

    status <- if (!isTRUE(snapshot$valid) || length(project_collisions) > 0) {
      "Blocked"
    } else if (length(missing_projects) > 0) {
      "Warning"
    } else {
      "Ready"
    }

    summaries[[length(summaries) + 1L]] <- data.frame(
      snapshot_id = snapshot$snapshot_id,
      status = status,
      chunks = collapse_values(registry$chunk_id),
      registry_rows = nrow(registry),
      event_rows = nrow(snapshot$events),
      completed_rows = nrow(snapshot$completed),
      log_files = length(snapshot$log_files),
      projects_expected = length(expected_projects),
      projects_in_import = sum(expected_projects %in% import_projects),
      projects_missing = length(missing_projects),
      projects_already_adopted = sum(already_adopted),
      project_collisions = length(project_collisions),
      warnings = paste(snapshot_warnings, collapse = " | "),
      snapshot_dir = snapshot$snapshot_dir,
      control_dir = snapshot$control_dir,
      projects_dir = snapshot$projects_dir,
      stringsAsFactors = FALSE
    )
  }

  if (nrow(duplicate_tracker) > 0) {
    for (field in c("chunk_id", "species_key", "project")) {
      field_data <- duplicate_tracker[!is.na(duplicate_tracker[[field]]) & nzchar(duplicate_tracker[[field]]), , drop = FALSE]
      if (identical(field, "chunk_id")) {
        duplicated_values <- field_data %>%
          dplyr::distinct(snapshot_id, chunk_id) %>%
          dplyr::count(chunk_id, name = "sources") %>%
          dplyr::filter(sources > 1) %>%
          dplyr::pull(chunk_id)
      } else {
        values <- field_data[[field]]
        duplicated_values <- unique(values[duplicated(values)])
      }
      if (length(duplicated_values) > 0) {
        warnings[[length(warnings) + 1L]] <- data.frame(
          snapshot_id = "multiple",
          severity = "error",
          message = paste("Duplicate", field, "across import snapshots:", collapse_values(duplicated_values)),
          stringsAsFactors = FALSE
        )
      }
    }
  }

  summary <- dplyr::bind_rows(summaries)
  warnings_df <- if (length(warnings) == 0) {
    as.data.frame(setNames(rep(list(character()), length(import_warning_columns)), import_warning_columns))
  } else {
    dplyr::bind_rows(warnings)
  }

  list(summary = summary[, import_summary_columns, drop = FALSE], warnings = warnings_df)
}

copy_import_logs <- function(snapshot, master_dir) {
  source_logs <- file.path(snapshot$control_dir, "logs")
  if (!dir.exists(source_logs)) {
    return(setNames(character(), character()))
  }

  log_files <- list.files(source_logs, pattern = "\\.log$", full.names = TRUE, recursive = TRUE)
  if (length(log_files) == 0) {
    return(setNames(character(), character()))
  }

  source_logs_norm <- normalizePath(source_logs, winslash = "/", mustWork = TRUE)
  copied <- character()
  names(copied) <- character()

  for (source_file in log_files) {
    source_norm <- normalizePath(source_file, winslash = "/", mustWork = TRUE)
    rel <- if (startsWith(source_norm, paste0(source_logs_norm, "/"))) {
      substring(source_norm, nchar(source_logs_norm) + 2L)
    } else {
      basename(source_norm)
    }
    dest_rel <- file.path("logs", "imports", snapshot$snapshot_id, rel)
    dest_file <- file.path(master_dir, dest_rel)
    safe_dir_create(dirname(dest_file))
    ok <- file.copy(source_norm, dest_file, overwrite = TRUE, copy.date = TRUE)
    if (isTRUE(ok)) {
      copied[[source_norm]] <- normalizePath(dest_rel, winslash = "/", mustWork = FALSE)
    }
  }

  copied
}

rewrite_event_log_paths <- function(events, snapshot, copied_logs) {
  if (nrow(events) == 0 || length(copied_logs) == 0) {
    return(events)
  }

  for (row_index in seq_len(nrow(events))) {
    raw_path <- events$log_file[[row_index]]
    if (is.na(raw_path) || !nzchar(raw_path)) {
      next
    }
    resolved <- resolve_existing_path(raw_path, control_dir = snapshot$control_dir)
    if (!file.exists(resolved)) {
      next
    }
    resolved <- normalizePath(resolved, winslash = "/", mustWork = TRUE)
    if (resolved %in% names(copied_logs)) {
      events$log_file[[row_index]] <- copied_logs[[resolved]]
    }
  }

  events
}

existing_rows_from_same_snapshot <- function(existing, snapshot_id) {
  if (nrow(existing) == 0 || !"source_snapshot" %in% names(existing)) {
    return(rep(FALSE, nrow(existing)))
  }
  !is.na(existing$source_snapshot) &
    nzchar(existing$source_snapshot) &
    existing$source_snapshot == snapshot_id
}

validate_merge_conflicts <- function(imports, master_dir) {
  imported_registry <- dplyr::bind_rows(lapply(imports, `[[`, "registry"))
  if (nrow(imported_registry) == 0) {
    stop("Selected imports do not contain registry rows.", call. = FALSE)
  }

  for (field in c("chunk_id", "species_key", "project")) {
    if (identical(field, "chunk_id")) {
      duplicates <- imported_registry %>%
        dplyr::filter(!is.na(chunk_id) & nzchar(chunk_id)) %>%
        dplyr::distinct(source_snapshot, chunk_id) %>%
        dplyr::count(chunk_id, name = "sources") %>%
        dplyr::filter(sources > 1) %>%
        dplyr::pull(chunk_id)
    } else {
      values <- imported_registry[[field]]
      values <- values[!is.na(values) & nzchar(values)]
      duplicates <- unique(values[duplicated(values)])
    }
    if (length(duplicates) > 0) {
      stop("Duplicate ", field, " in selected imports: ", collapse_values(duplicates), call. = FALSE)
    }
  }

  master_registry <- read_registry(master_dir)
  if (nrow(master_registry) == 0) {
    return(invisible(TRUE))
  }

  for (field in c("chunk_id", "species_key", "project")) {
    imported_pairs <- unique(imported_registry[, c("source_snapshot", field), drop = FALSE])
    imported_pairs <- imported_pairs[!is.na(imported_pairs[[field]]) & nzchar(imported_pairs[[field]]), , drop = FALSE]
    for (row_index in seq_len(nrow(imported_pairs))) {
      value <- imported_pairs[[field]][[row_index]]
      snapshot_id <- imported_pairs$source_snapshot[[row_index]]
      existing <- master_registry[master_registry[[field]] == value, , drop = FALSE]
      if (nrow(existing) == 0) {
        next
      }
      if (!all(existing_rows_from_same_snapshot(existing, snapshot_id))) {
        stop(
          "Master registry already contains ", field, " '", value,
          "' from a different or unknown source. Merge stopped.",
          call. = FALSE
        )
      }
    }
  }

  invisible(TRUE)
}

adopt_project_folder <- function(snapshot, project, master_dir) {
  source_dir <- file.path(snapshot$projects_dir, project)
  dest_dir <- file.path(projects_root, project)
  already_adopted <- project_already_adopted(project, snapshot$snapshot_id, master_dir, snapshot$snapshot_dir)

  if (dir.exists(dest_dir)) {
    if (already_adopted && !dir.exists(source_dir)) {
      return(list(state = "already_adopted", action = "none"))
    }
    if (already_adopted && dir.exists(source_dir)) {
      stop(
        "Project '", project, "' is already adopted, but the import staging folder still contains a duplicate copy. ",
        "Remove the staged duplicate before rerunning the merge.",
        call. = FALSE
      )
    }
    stop(
      "Destination project folder already exists and is not recorded as adopted: ",
      normalizePath(dest_dir, winslash = "/", mustWork = FALSE),
      call. = FALSE
    )
  }

  if (!dir.exists(source_dir)) {
    return(list(state = "missing", action = "none"))
  }

  safe_dir_create(dirname(dest_dir))
  moved <- file.rename(source_dir, dest_dir)
  if (!isTRUE(moved)) {
    stop(
      "Could not move project folder with file.rename(). This often means the import and destination are on different filesystems, ",
      "or the folder is locked. No copy/delete fallback was attempted for large model artifacts. Source: ",
      normalizePath(source_dir, winslash = "/", mustWork = FALSE),
      " Destination: ",
      normalizePath(dest_dir, winslash = "/", mustWork = FALSE),
      call. = FALSE
    )
  }

  list(state = "adopted", action = "moved")
}

write_chunked_events <- function(events, master_dir) {
  if (nrow(events) == 0) {
    return(invisible(NULL))
  }

  events <- events[, event_columns, drop = FALSE]
  for (chunk_id in nonempty_unique(events$chunk_id)) {
    chunk_events <- events[events$chunk_id == chunk_id, event_columns, drop = FALSE]
    safe_write_csv(chunk_events, file.path(master_dir, paste0("chunk_", chunk_id, "_events.csv")))
  }

  invisible(NULL)
}

write_chunked_completed <- function(completed, registry, master_dir) {
  if (nrow(completed) == 0) {
    return(invisible(NULL))
  }

  registry_chunks <- unique(registry[, c("species_key", "chunk_id"), drop = FALSE])
  completed_chunks <- dplyr::left_join(completed, registry_chunks, by = "species_key")
  for (chunk_id in nonempty_unique(completed_chunks$chunk_id)) {
    chunk_completed <- completed_chunks[completed_chunks$chunk_id == chunk_id, completed_columns, drop = FALSE]
    safe_write_csv(chunk_completed, file.path(master_dir, paste0("chunk_", chunk_id, "_completed_species.csv")))
  }

  invisible(NULL)
}

merge_import_snapshots <- function(snapshot_dirs, control_dir) {
  if (length(snapshot_dirs) == 0) {
    stop("Select at least one import snapshot to merge.", call. = FALSE)
  }

  master_dir <- master_control_dir_for(control_dir)
  if (!nzchar(master_dir)) {
    stop("Could not resolve the master control folder for the selected control folder.", call. = FALSE)
  }

  safe_dir_create(master_dir)
  safe_dir_create(file.path(master_dir, "locks"))
  merge_lock <- filelock::lock(file.path(master_dir, "locks", "merge.lock"), timeout = 10)
  if (is.null(merge_lock)) {
    stop("Another merge is already running for this master control folder.", call. = FALSE)
  }
  on.exit(filelock::unlock(merge_lock), add = TRUE)

  imports <- lapply(snapshot_dirs, read_snapshot_data)
  invalid <- vapply(imports, function(x) !isTRUE(x$valid), logical(1))
  if (any(invalid)) {
    stop(
      "Selected import(s) are invalid: ",
      paste(vapply(imports[invalid], `[[`, character(1), "snapshot_id"), collapse = ", "),
      call. = FALSE
    )
  }

  validate_merge_conflicts(imports, master_dir)

  merge_id <- paste0("merge_", timestamp_file())
  registry_master <- read_registry(master_dir)
  events_master <- read_events(master_dir)
  if ("log_file_resolved" %in% names(events_master)) {
    events_master$log_file_resolved <- NULL
  }
  completed_master <- read_completed(master_dir)
  manifest_master <- read_merge_manifest(master_dir)

  registry_to_add <- list()
  events_to_add <- list()
  completed_to_add <- list()
  manifest_to_add <- list()
  warnings <- list()
  current_artifact_states <- character()

  for (snapshot in imports) {
    copied_logs <- copy_import_logs(snapshot, master_dir)
    events <- rewrite_event_log_paths(snapshot$events, snapshot, copied_logs)

    existing_same_snapshot <- if (nrow(registry_master) > 0) {
      !is.na(registry_master$source_snapshot) &
        registry_master$source_snapshot == snapshot$snapshot_id
    } else {
      logical()
    }

    snapshot_registry <- snapshot$registry
    if (nrow(registry_master) > 0 && any(existing_same_snapshot, na.rm = TRUE)) {
      existing_species <- registry_master$species_key[existing_same_snapshot]
      snapshot_registry <- snapshot_registry[!snapshot_registry$species_key %in% existing_species, , drop = FALSE]
    }

    registry_to_add[[length(registry_to_add) + 1L]] <- snapshot_registry
    events_to_add[[length(events_to_add) + 1L]] <- events
    completed_to_add[[length(completed_to_add) + 1L]] <- snapshot$completed

    for (row_index in seq_len(nrow(snapshot$registry))) {
      project <- snapshot$registry$project[[row_index]]
      if (is.na(project) || !nzchar(project)) {
        next
      }

      adoption <- adopt_project_folder(snapshot, project, master_dir)
      current_artifact_states <- c(current_artifact_states, adoption$state)
      source_project_dir <- file.path(snapshot$projects_dir, project)
      destination_project_dir <- file.path(projects_root, project)

      if (identical(adoption$state, "missing")) {
        warnings[[length(warnings) + 1L]] <- data.frame(
          snapshot_id = snapshot$snapshot_id,
          severity = "warning",
          message = paste("Missing project folder for registry project:", project),
          stringsAsFactors = FALSE
        )
      }

      if (!identical(adoption$state, "already_adopted")) {
        manifest_to_add[[length(manifest_to_add) + 1L]] <- data.frame(
          merge_id = merge_id,
          timestamp = timestamp_now(),
          snapshot_id = snapshot$snapshot_id,
          chunk_id = snapshot$registry$chunk_id[[row_index]],
          species_key = snapshot$registry$species_key[[row_index]],
          project = project,
          source_control_dir = snapshot$control_dir,
          source_project_dir = normalizePath(source_project_dir, winslash = "/", mustWork = FALSE),
          destination_project_dir = normalizePath(destination_project_dir, winslash = "/", mustWork = FALSE),
          artifact_state = adoption$state,
          action = adoption$action,
          stringsAsFactors = FALSE
        )
      }
    }
  }

  registry_all <- dplyr::bind_rows(registry_master, dplyr::bind_rows(registry_to_add))
  events_all <- dplyr::bind_rows(events_master, dplyr::bind_rows(events_to_add))
  completed_all <- dplyr::bind_rows(completed_master, dplyr::bind_rows(completed_to_add))
  manifest_all <- dplyr::bind_rows(manifest_master, dplyr::bind_rows(manifest_to_add))

  registry_all <- dplyr::distinct(registry_all[, registry_columns, drop = FALSE])
  events_all <- dplyr::distinct(events_all[, event_columns, drop = FALSE])
  completed_all <- dplyr::distinct(completed_all[, completed_columns, drop = FALSE])
  manifest_all <- dplyr::distinct(manifest_all[, merge_manifest_columns, drop = FALSE])

  safe_write_csv(registry_all, file.path(master_dir, "species_registry.csv"))
  write_chunked_events(events_all, master_dir)
  write_chunked_completed(completed_all, registry_all, master_dir)
  safe_write_csv(manifest_all, file.path(master_dir, "merge_manifest.csv"))

  warning_df <- if (length(warnings) == 0) {
    as.data.frame(setNames(rep(list(character()), length(import_warning_columns)), import_warning_columns))
  } else {
    dplyr::bind_rows(warnings)
  }
  safe_write_csv(warning_df, file.path(master_dir, "merge_warnings.csv"))

  snapshot_manifest <- dplyr::bind_rows(manifest_to_add)
  if (nrow(snapshot_manifest) > 0) {
    for (snapshot_id in nonempty_unique(snapshot_manifest$snapshot_id)) {
      snapshot_dir <- imports[[which(vapply(imports, `[[`, character(1), "snapshot_id") == snapshot_id)[[1]]]]$snapshot_dir
      existing_snapshot_manifest <- read_snapshot_adoptions(snapshot_dir)
      out <- dplyr::distinct(dplyr::bind_rows(existing_snapshot_manifest, snapshot_manifest[snapshot_manifest$snapshot_id == snapshot_id, , drop = FALSE]))
      safe_write_csv(out[, merge_manifest_columns, drop = FALSE], snapshot_adoption_path(snapshot_dir))
    }
  }

  data.frame(
    item = c(
      "Master control folder", "Snapshots merged", "Registry rows added",
      "Total master registry rows", "Event rows in master", "Completed rows in master",
      "Artifacts adopted or already adopted", "Missing artifact folders"
    ),
    value = c(
      normalizePath(master_dir, winslash = "/", mustWork = FALSE),
      length(imports),
      sum(vapply(registry_to_add, nrow, integer(1))),
      nrow(registry_all),
      nrow(events_all),
      nrow(completed_all),
      sum(current_artifact_states %in% c("adopted", "already_adopted"), na.rm = TRUE),
      sum(current_artifact_states == "missing", na.rm = TRUE)
    ),
    stringsAsFactors = FALSE
  )
}


#-------------------------------------------------------------------------------
# UI
#-------------------------------------------------------------------------------
app_css <- "
body {
  background: #f5f7fa;
  color: #17212b;
}
.container-fluid {
  max-width: 1680px;
}
.app-title {
  margin-top: 12px;
  margin-bottom: 16px;
}
.app-subtitle {
  color: #5a6673;
  font-size: 0.96rem;
}
.dashboard-shell {
  display: flex;
  align-items: flex-start;
  gap: 16px;
  width: 100%;
}
.dashboard-controls {
  flex: 0 0 calc(25% - 8px);
  max-width: 420px;
  min-width: 280px;
  background: #ffffff;
  border: 1px solid #dce3ea;
  border-radius: 8px;
  padding: 16px;
  box-shadow: 0 1px 2px rgba(16, 24, 40, 0.04);
}
.dashboard-main {
  flex: 1 1 calc(75% - 8px);
  min-width: 0;
}
.dashboard-shell.controls-collapsed .dashboard-controls {
  display: none;
}
.dashboard-shell.controls-collapsed .dashboard-main {
  flex-basis: 100%;
  max-width: 100%;
}
.controls-toolbar {
  display: flex;
  justify-content: flex-end;
  margin-bottom: 10px;
}
.controls-toggle {
  min-width: 116px;
}
.section-card {
  background: #ffffff;
  border: 1px solid #dce3ea;
  border-radius: 8px;
  padding: 16px;
  margin-bottom: 16px;
  box-shadow: 0 1px 2px rgba(16, 24, 40, 0.04);
}
.summary-grid {
  display: grid;
  grid-template-columns: repeat(auto-fit, minmax(150px, 1fr));
  gap: 12px;
  margin-bottom: 16px;
}
.summary-card {
  background: #ffffff;
  border: 1px solid #dce3ea;
  border-left: 5px solid #6c757d;
  border-radius: 8px;
  padding: 14px 16px;
  min-height: 96px;
}
.summary-card.completed { border-left-color: #198754; }
.summary-card.partial { border-left-color: #b7791f; }
.summary-card.blocked { border-left-color: #f59f00; }
.summary-card.running { border-left-color: #0d6efd; }
.summary-card.stale { border-left-color: #7c3aed; }
.summary-card.failed { border-left-color: #dc3545; }
.summary-card.skipped { border-left-color: #f59f00; }
.summary-card.unknown { border-left-color: #6c757d; }
.summary-label {
  color: #5a6673;
  font-size: 0.83rem;
  text-transform: uppercase;
  letter-spacing: 0;
  margin-bottom: 6px;
}
.summary-value {
  font-size: 2rem;
  font-weight: 700;
  line-height: 1;
}
.status-badge,
.state-badge {
  display: inline-block;
  min-width: 72px;
  padding: 4px 8px;
  border-radius: 999px;
  color: #ffffff;
  font-weight: 700;
  text-align: center;
  font-size: 0.78rem;
}
.status-completed,
.status-full_model_completed,
.state-ok,
.state-done { background: #198754; }
.state-reuse { background: #0f766e; }
.status-partial_climate_only { background: #b7791f; }
.status-completed_unknown { background: #20c997; color: #102a43; }
.status-running,
.status-active_running,
.state-start { background: #0d6efd; }
.status-stale_running { background: #7c3aed; }
.status-running_unknown { background: #6c757d; }
.status-failed,
.state-error { background: #dc3545; }
.status-no_climate_model,
.status-skipped,
.state-skip { background: #f59f00; color: #1f2933; }
.stale-note {
  background: #f8f5ff;
  border: 1px solid #ded0ff;
  border-left: 5px solid #7c3aed;
  border-radius: 8px;
  color: #3b2f56;
  padding: 12px 14px;
  margin-bottom: 16px;
}
.completion-note {
  background: #f7fbf8;
  border: 1px solid #cfe8d6;
  border-left: 5px solid #198754;
  border-radius: 8px;
  color: #1f3d2b;
  padding: 12px 14px;
  margin-bottom: 16px;
}
.log-preview {
  white-space: pre-wrap;
  background: #101820;
  color: #e8edf2;
  border-radius: 8px;
  padding: 16px;
  min-height: 560px;
  max-height: calc(100vh - 260px);
  overflow-y: auto;
  font-family: Consolas, 'Liberation Mono', monospace;
  font-size: 0.86rem;
}
.log-card-header {
  display: flex;
  align-items: flex-start;
  justify-content: space-between;
  gap: 16px;
  flex-wrap: wrap;
  margin-bottom: 12px;
}
.log-card-header h4 {
  margin-top: 0;
  margin-bottom: 4px;
}
.log-actions {
  display: flex;
  align-items: center;
  gap: 8px;
  flex-wrap: wrap;
}
.logs-table-wrap .dataTables_wrapper {
  width: 100%;
}
.small-muted {
  color: #6c757d;
  font-size: 0.86rem;
}
.sidebar-note {
  color: #5a6673;
  font-size: 0.86rem;
  margin-top: 8px;
}
.merge-note {
  background: #f7fafc;
  border: 1px solid #dce3ea;
  border-left: 5px solid #2a6f97;
  border-radius: 8px;
  color: #26323f;
  padding: 12px 14px;
  margin-bottom: 16px;
}
.merge-actions {
  display: flex;
  align-items: center;
  gap: 8px;
  flex-wrap: wrap;
  margin-bottom: 12px;
}
.merge-result {
  margin-top: 12px;
}
.import-help-code {
  background: #101820;
  color: #e8edf2;
  border-radius: 8px;
  padding: 12px;
  white-space: pre;
  overflow-x: auto;
  font-family: Consolas, 'Liberation Mono', monospace;
  font-size: 0.86rem;
}
@media (max-width: 900px) {
  .dashboard-shell {
    display: block;
  }
  .dashboard-controls {
    max-width: none;
    min-width: 0;
    width: 100%;
    margin-bottom: 16px;
  }
  .dashboard-main {
    width: 100%;
  }
  .controls-toolbar {
    justify-content: stretch;
  }
  .controls-toggle {
    width: 100%;
  }
  .log-preview {
    min-height: 420px;
    max-height: 70vh;
  }
}
"

ui <- fluidPage(
  tags$head(
    tags$title("wiSDM Chunk Progress Dashboard"),
    tags$style(HTML(app_css)),
    tags$script(HTML("
      $(document).on('click', '#controls_toggle', function() {
        var shell = $('#dashboard_shell');
        shell.toggleClass('controls-expanded controls-collapsed');
        var isCollapsed = shell.hasClass('controls-collapsed');
        $('#controls_toggle').text(isCollapsed ? 'Show controls' : 'Hide controls');
        setTimeout(function() {
          $(window).trigger('resize');
          if ($.fn.dataTable) {
            $.fn.dataTable.tables({ visible: true, api: true }).columns.adjust();
          }
        }, 250);
      });
    "))
  ),
  div(
    class = "app-title",
    h2("wiSDM Chunk Progress Dashboard"),
    div(
      class = "app-subtitle",
      "Read-only view of chunk registries, stage events, completed species, and per-stage log files."
    )
  ),
  div(
    id = "dashboard_shell",
    class = "dashboard-shell controls-expanded",
    div(
      class = "dashboard-controls",
      actionButton("refresh_dirs", "Refresh control folders", class = "btn-primary"),
      selectInput("control_dir_choice", "Control folder", choices = character()),
      textInput(
        "manual_control_dir",
        "Manual control folder",
        value = "",
        placeholder = "Optional full path to *_chunk_control"
      ),
      actionButton("refresh", "Reload logs", class = "btn-success"),
      actionButton("open_control_dir", "Open control folder"),
      tags$hr(),
      selectInput("chunk_filter", "Chunk", choices = c("All chunks" = "all")),
      checkboxInput("auto_refresh", "Auto-refresh every 60 seconds", value = TRUE),
      numericInput("log_lines", "Log preview lines", value = 1000, min = 100, max = 10000, step = 100),
      tags$hr(),
      uiOutput("control_folder_status"),
      div(
        class = "sidebar-note",
        "Tip: select a species in the Species tab, then open the Logs tab to inspect its stage log files."
      )
    ),
    div(
      class = "dashboard-main",
      div(
        class = "controls-toolbar",
        tags$button(
          id = "controls_toggle",
          type = "button",
          class = "btn btn-default controls-toggle",
          "Hide controls"
        )
      ),
      tabsetPanel(
        id = "main_tabs",
        tabPanel(
          "Overview",
          br(),
          uiOutput("summary_cards"),
          div(
            class = "stale-note",
            "Stale running rows usually mean an RStudio background job was killed, the R session crashed, or the process ended before the registry could be marked failed. The dashboard does not touch runner locks; it classifies running rows from the recorded PID only."
          ),
          div(
            class = "completion-note",
            "Full models have both climate and habitat model files. Partial climate-only runs reached cross-validation but did not produce a habitat model file after stage 04."
          ),
          fluidRow(
            column(6, div(class = "section-card", h4("Species Status"), plotOutput("status_plot", height = 280))),
            column(6, div(class = "section-card", h4("Outcomes"), plotOutput("outcome_plot", height = 280)))
          ),
          fluidRow(
            column(12, div(class = "section-card", h4("Latest Activity"), DTOutput("latest_events_table")))
          )
        ),
        tabPanel(
          "Chunks",
          br(),
          div(class = "section-card", h4("Chunk Summary"), DTOutput("chunk_summary_table")),
          div(class = "section-card", h4("Chunk Stage Errors"), DTOutput("chunk_error_table"))
        ),
        tabPanel(
          "Species",
          br(),
          div(class = "section-card", h4("Species Registry"), DTOutput("species_table")),
          div(class = "section-card", h4("Stage Matrix"), DTOutput("stage_matrix_table"))
        ),
        tabPanel(
          "Events",
          br(),
          checkboxInput("events_selected_only", "Show events for selected species only", value = FALSE),
          div(class = "section-card", h4("Event Log"), DTOutput("events_table"))
        ),
        tabPanel(
          "Imports / Merge",
          br(),
          div(
            class = "merge-note",
            tags$strong("No-duplication merge: "),
            "control CSVs and logs are merged into a separate master control folder, while large species project folders are adopted with file.rename() from the import staging area into data/projects/. No overwrite and no copy/delete fallback are used for model artifacts."
          ),
          div(
            class = "section-card",
            h4("Import Snapshots"),
            div(class = "small-muted", textOutput("imports_root_label")),
            div(class = "small-muted", textOutput("master_control_label")),
            div(
              class = "merge-actions",
              actionButton("show_import_help", "Import instructions", class = "btn-default"),
              actionButton("scan_imports", "Scan imports", class = "btn-primary"),
              actionButton("merge_imports", "Merge selected imports", class = "btn-danger")
            ),
            DTOutput("import_scan_table")
          ),
          div(
            class = "section-card",
            h4("Import Warnings"),
            DTOutput("import_warnings_table")
          ),
          div(
            class = "section-card merge-result",
            h4("Latest Merge Result"),
            DTOutput("merge_result_table")
          )
        ),
        tabPanel(
          "Logs",
          br(),
          div(
            class = "section-card",
            div(
              class = "log-card-header",
              div(
                h4("Log Files"),
                div(class = "small-muted", textOutput("selected_species_label"))
              ),
              div(
                class = "log-actions",
                actionButton("open_log_file", "Open selected log"),
                actionButton("open_log_folder", "Open log folder"),
                downloadButton("download_log", "Download selected log")
              )
            ),
            div(class = "logs-table-wrap", DTOutput("logs_table"))
          ),
          div(
            class = "section-card",
            h4("Log Preview"),
            div(class = "log-preview", textOutput("log_preview"))
          )
        )
      )
    )
  )
)


#-------------------------------------------------------------------------------
# Server
#-------------------------------------------------------------------------------
server <- function(input, output, session) {
  auto_timer <- reactiveTimer(60000, session = session)

  observe({
    dirs <- list_control_dirs()
    choices <- named_choices(dirs, basename(dirs))
    selected <- isolate(input$control_dir_choice)

    if (length(choices) == 0) {
      choices <- c("No *_chunk_control folders found" = "")
    }

    updateSelectInput(
      session,
      "control_dir_choice",
      choices = choices,
      selected = if (length(selected) == 1 && selected %in% choices) selected else choices[[1]]
    )
  })

  observeEvent(input$refresh_dirs, {
    dirs <- list_control_dirs()
    choices <- named_choices(dirs, basename(dirs))

    if (length(choices) == 0) {
      choices <- c("No *_chunk_control folders found" = "")
    }

    updateSelectInput(session, "control_dir_choice", choices = choices, selected = choices[[1]])
  })

  selected_control_dir <- reactive({
    manual <- trimws(input$manual_control_dir %||% "")
    if (nzchar(manual)) {
      return(normalizePath(manual, winslash = "/", mustWork = FALSE))
    }

    choice <- input$control_dir_choice %||% ""
    if (!nzchar(choice)) {
      return("")
    }

    normalizePath(choice, winslash = "/", mustWork = FALSE)
  })

  import_scan_state <- reactiveVal(empty_import_scan())
  merge_result_state <- reactiveVal(empty_merge_result())

  observeEvent(selected_control_dir(), {
    import_scan_state(empty_import_scan())
    merge_result_state(empty_merge_result())
  })

  output$imports_root_label <- renderText({
    paste("Imports folder:", normalizePath(imports_root_for(selected_control_dir()), winslash = "/", mustWork = FALSE))
  })

  output$master_control_label <- renderText({
    paste("Master control folder:", normalizePath(master_control_dir_for(selected_control_dir()), winslash = "/", mustWork = FALSE))
  })

  observeEvent(input$show_import_help, {
    prefix <- control_prefix(selected_control_dir())
    prefix <- if (nzchar(prefix)) prefix else "<project_prefix>"

    showModal(modalDialog(
      title = "Import Folder Structure",
      size = "l",
      easyClose = TRUE,
      footer = tagList(
        modalButton("Close")
      ),
      tags$p(
        "Use one import snapshot folder per machine, server, or transfer batch. ",
        "A snapshot can contain one chunk or multiple chunks."
      ),
      tags$pre(
        class = "import-help-code",
        paste0(
          "data/projects/", prefix, "_chunk_control/imports/serverA_20260516/\n",
          "  control/\n",
          "    species_registry.csv\n",
          "    chunk_01_events.csv\n",
          "    chunk_02_events.csv\n",
          "    chunk_01_completed_species.csv\n",
          "    chunk_02_completed_species.csv\n",
          "    logs/\n",
          "  projects/\n",
          "    ", prefix, "_01_species_a/\n",
          "    ", prefix, "_02_species_b/\n"
        )
      ),
      tags$h4("What To Include"),
      tags$ul(
        tags$li(tags$strong("Required: "), "control/species_registry.csv."),
        tags$li(tags$strong("Recommended: "), "control/chunk_*_events.csv and control/chunk_*_completed_species.csv."),
        tags$li(tags$strong("Recommended for log viewing: "), "control/logs/."),
        tags$li(tags$strong("Required for model artifacts: "), "projects/<species_project_folder>/ for each imported project.")
      ),
      tags$h4("What To Exclude"),
      tags$ul(
        tags$li("Do not include the source control folder's locks/ directory. Runtime locks are machine-local and are ignored after jobs stop."),
        tags$li("Do not paste imported files directly over the active runner control folder.")
      ),
      tags$h4("Merge Behavior"),
      tags$ul(
        tags$li("Control CSVs and logs are merged into a separate master control folder."),
        tags$li("Large species project folders are moved from imports/<snapshot>/projects/ into data/projects/ with file.rename()."),
        tags$li("No project folder is overwritten. If a destination already exists unexpectedly, the merge stops."),
        tags$li("If file.rename() fails, the merge stops. It does not fall back to copy/delete for large model folders."),
        tags$li("The same snapshot folder name can be reused for later imports from the same machine. Add only new project folders under projects/.")
      ),
      tags$p(
        tags$strong("After merge: "),
        "select ", tags$code(paste0(prefix, "_master_chunk_control")),
        " in the control folder dropdown to inspect the combined results."
      )
    ))
  })

  observeEvent(input$scan_imports, {
    scan <- scan_import_snapshots(selected_control_dir())
    import_scan_state(scan)
    showNotification(
      paste("Scanned", nrow(scan$summary), "import snapshot(s)."),
      type = if (nrow(scan$warnings) > 0 && any(scan$warnings$severity == "error")) "warning" else "message"
    )
  })

  selected_import_dirs <- reactive({
    scan <- import_scan_state()$summary
    selected <- input$import_scan_table_rows_selected
    if (length(selected) == 0 || nrow(scan) == 0) {
      return(character())
    }
    scan$snapshot_dir[selected]
  })

  observeEvent(input$merge_imports, {
    selected_dirs <- selected_import_dirs()
    if (length(selected_dirs) == 0) {
      showNotification("Select at least one import snapshot first.", type = "error")
      return(invisible(NULL))
    }

    result <- tryCatch(
      merge_import_snapshots(selected_dirs, selected_control_dir()),
      error = function(e) e
    )

    if (inherits(result, "error")) {
      merge_result_state(data.frame(
        item = "Merge failed",
        value = conditionMessage(result),
        stringsAsFactors = FALSE
      ))
      showNotification(conditionMessage(result), type = "error", duration = 12)
      return(invisible(NULL))
    }

    merge_result_state(result)
    import_scan_state(scan_import_snapshots(selected_control_dir()))

    dirs <- list_control_dirs()
    choices <- named_choices(dirs, basename(dirs))
    master_dir <- master_control_dir_for(selected_control_dir())
    if (length(choices) > 0 && master_dir %in% choices) {
      updateSelectInput(session, "control_dir_choice", choices = choices, selected = master_dir)
    }

    showNotification("Merge completed. The master control folder is ready for inspection.", type = "message")
  })

  control_data <- reactive({
    input$refresh
    if (isTRUE(input$auto_refresh)) {
      auto_timer()
    }

    read_control_data(selected_control_dir())
  })

  filtered_registry <- reactive({
    filter_by_chunk(control_data()$registry, input$chunk_filter)
  })

  filtered_events <- reactive({
    filter_by_chunk(control_data()$events, input$chunk_filter)
  })

  filtered_completed <- reactive({
    completed <- control_data()$completed
    if (nrow(completed) == 0 || identical(input$chunk_filter, "all")) {
      return(completed)
    }

    registry <- control_data()$registry[, c("species_key", "chunk_id"), drop = FALSE]
    completed <- dplyr::left_join(completed, registry, by = "species_key")
    completed[completed$chunk_id == input$chunk_filter, , drop = FALSE]
  })

  observe({
    data <- control_data()
    chunks <- sort(unique(data$registry$chunk_id[!is.na(data$registry$chunk_id)]))
    event_chunks <- sort(unique(data$events$chunk_id[!is.na(data$events$chunk_id)]))
    chunks <- sort(unique(c(chunks, event_chunks)))

    choices <- c("All chunks" = "all", named_choices(chunks, paste("Chunk", chunks)))
    selected <- isolate(input$chunk_filter)
    if (length(selected) != 1 || !selected %in% choices) {
      selected <- "all"
    }
    updateSelectInput(session, "chunk_filter", choices = choices, selected = selected)
  })

  output$control_folder_status <- renderUI({
    data <- control_data()
    status <- if (isTRUE(data$valid)) "Readable" else "Not found"
    registry_n <- nrow(data$registry)
    events_n <- nrow(data$events)
    logs_n <- nrow(data$log_files)

    tags$div(
      tags$strong("Current folder"),
      tags$br(),
      tags$code(data$control_dir %||% ""),
      tags$br(),
      tags$span(class = "small-muted", paste("Status:", status)),
      tags$br(),
      tags$span(class = "small-muted", paste("Registry rows:", registry_n)),
      tags$br(),
      tags$span(class = "small-muted", paste("Events:", events_n)),
      tags$br(),
      tags$span(class = "small-muted", paste("Log files:", logs_n)),
      tags$br(),
      tags$span(class = "small-muted", paste("Loaded:", format(data$loaded_at, "%Y-%m-%d %H:%M:%S")))
    )
  })

  output$import_scan_table <- renderDT({
    scan <- import_scan_state()$summary
    if (nrow(scan) == 0) {
      return(DT::datatable(data.frame(), rownames = FALSE))
    }

    display <- scan
    DT::datatable(
      display,
      rownames = FALSE,
      selection = "multiple",
      filter = "top",
      options = list(
        pageLength = 10,
        scrollX = TRUE,
        columnDefs = list(list(visible = FALSE, targets = c(13, 14, 15)))
      )
    )
  })

  output$import_warnings_table <- renderDT({
    warnings <- import_scan_state()$warnings
    DT::datatable(
      warnings,
      rownames = FALSE,
      filter = if (nrow(warnings) > 0) "top" else "none",
      options = list(pageLength = 10, scrollX = TRUE)
    )
  })

  output$merge_result_table <- renderDT({
    DT::datatable(
      merge_result_state(),
      rownames = FALSE,
      options = list(pageLength = 10, scrollX = TRUE, dom = "tip")
    )
  })

  output$summary_cards <- renderUI({
    counts <- summary_counts(filtered_registry(), filtered_events(), filtered_completed())

    div(
      class = "summary-grid",
      lapply(seq_len(nrow(counts)), function(i) {
        div(
          class = paste("summary-card", counts$tone[[i]]),
          div(class = "summary-label", counts$label[[i]]),
          div(class = "summary-value", counts$value[[i]])
        )
      })
    )
  })

  output$status_plot <- renderPlot({
    registry <- filtered_registry()
    if (nrow(registry) == 0) {
      plot.new()
      text(0.5, 0.5, "No registry rows found")
      return(invisible(NULL))
    }

    plot_data <- registry %>%
      dplyr::mutate(
        model_status = dplyr::if_else(
          is.na(model_status) | !nzchar(model_status),
          "unknown",
          model_status
        )
      ) %>%
      dplyr::count(model_status, name = "species") %>%
      dplyr::mutate(status_label = status_label(model_status))

    ggplot(plot_data, aes(x = reorder(status_label, species), y = species, fill = model_status)) +
      geom_col(width = 0.72) +
      coord_flip() +
      scale_fill_manual(values = status_palette, na.value = "#6c757d") +
      labs(x = NULL, y = "Species", fill = NULL) +
      theme_minimal(base_size = 12) +
      theme(legend.position = "none", panel.grid.major.y = element_blank())
  })

  output$outcome_plot <- renderPlot({
    registry <- filtered_registry()
    if (nrow(registry) == 0 || all(is.na(registry$outcome) | !nzchar(registry$outcome))) {
      plot.new()
      text(0.5, 0.5, "No outcomes recorded yet")
      return(invisible(NULL))
    }

    plot_data <- registry %>%
      dplyr::mutate(outcome = dplyr::if_else(is.na(outcome) | !nzchar(outcome), "not_final", outcome)) %>%
      dplyr::count(outcome, name = "species")

    ggplot(plot_data, aes(x = reorder(outcome, species), y = species, fill = outcome)) +
      geom_col(width = 0.72, fill = "#2a6f97") +
      coord_flip() +
      labs(x = NULL, y = "Species") +
      theme_minimal(base_size = 12) +
      theme(panel.grid.major.y = element_blank())
  })

  output$latest_events_table <- renderDT({
    events <- filtered_events()
    if (nrow(events) > 0) {
      events <- events %>%
        dplyr::arrange(dplyr::desc(timestamp)) %>%
        dplyr::select(timestamp, chunk_id, species_requested, project, stage, state, message) %>%
        head(25)
      events$state <- state_badge(events$state)
    }

    DT::datatable(
      events,
      escape = FALSE,
      rownames = FALSE,
      options = list(pageLength = 10, scrollX = TRUE, dom = "tip")
    )
  })

  output$chunk_summary_table <- renderDT({
    registry <- control_data()$registry
    if (nrow(registry) == 0) {
      return(DT::datatable(data.frame(), rownames = FALSE))
    }

    summary <- registry %>%
      dplyr::mutate(
        model_status = dplyr::if_else(
          is.na(model_status) | !nzchar(model_status),
          "unknown",
          model_status
        )
      ) %>%
      dplyr::count(chunk_id, model_status, name = "species") %>%
      tidyr::pivot_wider(names_from = model_status, values_from = species, values_fill = 0) %>%
      dplyr::arrange(chunk_id)

    DT::datatable(summary, rownames = FALSE, options = list(pageLength = 20, scrollX = TRUE))
  })

  output$chunk_error_table <- renderDT({
    errors <- control_data()$events %>%
      dplyr::filter(state == "ERROR") %>%
      dplyr::count(chunk_id, stage, name = "errors") %>%
      dplyr::arrange(chunk_id, stage)

    DT::datatable(errors, rownames = FALSE, options = list(pageLength = 20, scrollX = TRUE))
  })

  output$species_table <- renderDT({
    registry <- filtered_registry()
    if (nrow(registry) > 0) {
      registry <- registry %>%
        dplyr::select(
          chunk_id, species_requested, project, attempt, registry_status,
          effective_status, model_status, model_status_note, pid_state, lock_state, stale_reason, outcome, stage,
          started_at, updated_at, completed_at, pid, run_id, source_snapshot, source_file, error_message
        )
      registry$registry_status <- status_badge(registry$registry_status)
      registry$effective_status <- status_badge(registry$effective_status)
      registry$model_status <- status_badge(registry$model_status)
    }

    DT::datatable(
      registry,
      escape = FALSE,
      rownames = FALSE,
      selection = "single",
      filter = "top",
      options = list(pageLength = 25, scrollX = TRUE)
    )
  })

  selected_species <- reactive({
    registry <- filtered_registry()
    selected <- input$species_table_rows_selected
    if (length(selected) == 0 || nrow(registry) == 0) {
      return(NULL)
    }

    registry[selected[[1]], , drop = FALSE]
  })

  output$stage_matrix_table <- renderDT({
    matrix <- stage_matrix(filtered_registry(), filtered_events())
    if (nrow(matrix) > 0) {
      state_columns <- intersect(c("claim", "02", "03", "04", "05", "complete"), names(matrix))
      for (column in state_columns) {
        matrix[[column]] <- state_badge(matrix[[column]])
      }
      matrix$registry_status <- status_badge(matrix$registry_status)
      matrix$effective_status <- status_badge(matrix$effective_status)
      matrix$model_status <- status_badge(matrix$model_status)
    }

    DT::datatable(
      matrix,
      escape = FALSE,
      rownames = FALSE,
      filter = "top",
      options = list(pageLength = 25, scrollX = TRUE)
    )
  })

  output$events_table <- renderDT({
    events <- filtered_events()
    selected <- selected_species()

    if (isTRUE(input$events_selected_only) && !is.null(selected) && nrow(events) > 0) {
      events <- events[events$species_key == selected$species_key[[1]], , drop = FALSE]
    }

    if (nrow(events) > 0) {
      events <- events %>%
        dplyr::arrange(dplyr::desc(timestamp)) %>%
        dplyr::select(
          timestamp, chunk_id, species_requested, project, attempt, stage, state,
          message, log_file, source_snapshot, source_file
        )
      events$state <- state_badge(events$state)
    }

    DT::datatable(
      events,
      escape = FALSE,
      rownames = FALSE,
      filter = "top",
      options = list(pageLength = 30, scrollX = TRUE)
    )
  })

  output$selected_species_label <- renderText({
    selected <- selected_species()
    if (is.null(selected)) {
      return("No species selected. Showing all log files for the current chunk filter.")
    }

    paste("Selected species:", selected$species_requested[[1]])
  })

  logs_for_selection <- reactive({
    data <- control_data()
    logs <- data$log_files
    events <- filtered_events()
    selected <- selected_species()

    if (nrow(logs) == 0) {
      return(logs)
    }

    if (!is.null(selected)) {
      species_events <- events[events$species_key == selected$species_key[[1]], , drop = FALSE]
      event_paths <- unique(species_events$log_file_resolved)
      event_paths <- event_paths[!is.na(event_paths) & nzchar(event_paths)]

      if (length(event_paths) > 0) {
        logs <- logs[logs$path %in% event_paths, , drop = FALSE]
      } else {
        project_slug <- sub("^.*_[0-9]{2}_", "", selected$project[[1]])
        logs <- logs[logs$folder == project_slug, , drop = FALSE]
      }
    } else if (!identical(input$chunk_filter, "all")) {
      chunk_folder <- paste0("chunk_", input$chunk_filter)
      logs <- logs[grepl(paste0("/", chunk_folder, "/"), logs$path, fixed = TRUE), , drop = FALSE]
    }

    logs[order(logs$modified, decreasing = TRUE), , drop = FALSE]
  })

  output$logs_table <- renderDT({
    logs <- logs_for_selection()
    display <- logs %>%
      dplyr::mutate(
        attempt = ifelse(is.na(attempt), "", as.character(attempt)),
        modified = ifelse(is.na(modified), "", format(modified, "%Y-%m-%d %H:%M:%S"))
      ) %>%
      dplyr::select(chunk, species_slug, stage, attempt, modified, size_mb, name, path)

    DT::datatable(
      display,
      rownames = FALSE,
      selection = "single",
      filter = "top",
      options = list(
        pageLength = 10,
        scrollX = TRUE,
        order = list(list(4, "desc")),
        columnDefs = list(list(visible = FALSE, targets = 7))
      )
    )
  })

  selected_log_path <- reactive({
    logs <- logs_for_selection()
    selected <- input$logs_table_rows_selected
    if (length(selected) == 0 || nrow(logs) == 0) {
      return(NA_character_)
    }

    logs$path[[selected[[1]]]]
  })

  output$log_preview <- renderText({
    clean_log_text(selected_log_path(), n_lines = input$log_lines %||% 1000L)
  })

  observeEvent(input$open_control_dir, {
    opened <- open_local_path(selected_control_dir())
    if (!opened) {
      showNotification("Control folder does not exist.", type = "error")
    }
  })

  observeEvent(input$open_log_file, {
    opened <- open_local_path(selected_log_path())
    if (!opened) {
      showNotification("Select an existing log file first.", type = "error")
    }
  })

  observeEvent(input$open_log_folder, {
    path <- selected_log_path()
    if (length(path) == 0 || is.na(path) || !file.exists(path)) {
      showNotification("Select an existing log file first.", type = "error")
      return(invisible(NULL))
    }

    opened <- open_local_path(dirname(path))
    if (!opened) {
      showNotification("Could not open the log folder.", type = "error")
    }
  })

  output$download_log <- downloadHandler(
    filename = function() {
      path <- selected_log_path()
      if (length(path) == 0 || is.na(path) || !file.exists(path)) {
        return("wisdm_log.txt")
      }
      basename(path)
    },
    content = function(file) {
      path <- selected_log_path()
      if (length(path) == 0 || is.na(path) || !file.exists(path)) {
        writeLines("No log file selected.", file)
      } else {
        file.copy(path, file, overwrite = TRUE)
      }
    }
  )
}

shinyApp(ui, server)
