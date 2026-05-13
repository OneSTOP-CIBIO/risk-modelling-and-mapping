#-------------------------------------------------------------------------------
# Detached launcher for the wiSDM chunk-runner dashboard
#-------------------------------------------------------------------------------
# Usage from the RStudio main console:
#   source("src/_launch_wisdm_chunk_dashboard.R")
#
# This file intentionally does not run shiny::runApp() in the foreground console.
# It starts the dashboard in a separate Rscript process, opens the browser, and
# immediately returns control to RStudio. This avoids RStudio console Stop events
# being coupled to running chunk-worker background jobs.
#
# Optional environment variables:
#   WISDM_DASHBOARD_HOST: default "127.0.0.1"
#   WISDM_DASHBOARD_PORT: default "3840"
#
# Optional command-line actions:
#   Rscript src/_launch_wisdm_chunk_dashboard.R --status
#   Rscript src/_launch_wisdm_chunk_dashboard.R --stop
#   Rscript src/_launch_wisdm_chunk_dashboard.R --foreground

`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}

trailing_args <- commandArgs(trailingOnly = TRUE)
host <- Sys.getenv("WISDM_DASHBOARD_HOST", "127.0.0.1")
port <- as.integer(Sys.getenv("WISDM_DASHBOARD_PORT", "3840"))

if (is.na(port) || port <= 0 || port > 65535) {
  stop("WISDM_DASHBOARD_PORT must be an integer between 1 and 65535.", call. = FALSE)
}

find_script_path <- function() {
  command_args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", command_args, value = TRUE)
  if (length(file_arg) > 0) {
    path <- sub("^--file=", "", file_arg[[1]])
    return(normalizePath(path, winslash = "/", mustWork = TRUE))
  }

  frame_files <- vapply(
    sys.frames(),
    function(frame) {
      file <- frame$ofile %||% ""
      if (length(file) == 0 || is.null(file) || is.na(file)) "" else file
    },
    character(1)
  )
  frame_files <- frame_files[nzchar(frame_files)]
  if (length(frame_files) > 0) {
    return(normalizePath(frame_files[[length(frame_files)]], winslash = "/", mustWork = TRUE))
  }

  normalizePath(file.path("src", "_launch_wisdm_chunk_dashboard.R"), winslash = "/", mustWork = TRUE)
}

find_repo_root <- function(start) {
  current <- normalizePath(start, winslash = "/", mustWork = FALSE)

  repeat {
    runner_path <- file.path(current, "src", "_run_wisdm_by_chunks.R")
    dashboard_path <- file.path(current, "src", "_wisdm_chunk_dashboard.R")
    if (file.exists(runner_path) && file.exists(dashboard_path)) {
      return(current)
    }

    parent <- dirname(current)
    if (identical(parent, current)) {
      break
    }
    current <- parent
  }

  stop("Could not find the wiSDM repository root from: ", start, call. = FALSE)
}

script_path <- find_script_path()
repo_root <- find_repo_root(dirname(script_path))
app_path <- file.path(repo_root, "src", "_wisdm_chunk_dashboard.R")
dashboard_url <- paste0("http://", host, ":", port)
dashboard_state_dir <- file.path(repo_root, "data", "projects")
dashboard_pid_file <- file.path(dashboard_state_dir, "wisdm_chunk_dashboard.pid")
dashboard_out_log <- file.path(dashboard_state_dir, "wisdm_chunk_dashboard.out.log")
dashboard_err_log <- file.path(dashboard_state_dir, "wisdm_chunk_dashboard.err.log")

dir.create(dashboard_state_dir, recursive = TRUE, showWarnings = FALSE)

timestamp_now <- function() {
  format(Sys.time(), "%Y-%m-%d %H:%M:%S")
}

rscript_path <- function() {
  exe <- if (.Platform$OS.type == "windows") "Rscript.exe" else "Rscript"
  path <- file.path(R.home("bin"), exe)
  if (file.exists(path)) {
    return(path)
  }
  exe
}

port_open <- function(host, port, timeout = 0.25) {
  con <- tryCatch(
    socketConnection(host = host, port = port, open = "r+", blocking = TRUE, timeout = timeout),
    error = function(e) NULL
  )
  if (is.null(con)) {
    return(FALSE)
  }
  close(con)
  TRUE
}

wait_for_dashboard <- function(seconds = 15) {
  deadline <- Sys.time() + seconds
  while (Sys.time() < deadline) {
    if (port_open(host, port)) {
      return(TRUE)
    }
    Sys.sleep(0.5)
  }
  FALSE
}

write_pid_file <- function() {
  writeLines(
    c(
      paste0("pid=", Sys.getpid()),
      paste0("started_at=", timestamp_now()),
      paste0("host=", host),
      paste0("port=", port),
      paste0("url=", dashboard_url)
    ),
    dashboard_pid_file
  )
}

read_pid_file <- function() {
  if (!file.exists(dashboard_pid_file)) {
    return(list())
  }

  lines <- readLines(dashboard_pid_file, warn = FALSE)
  parts <- strsplit(lines, "=", fixed = TRUE)
  parts <- parts[lengths(parts) >= 2]
  if (length(parts) == 0) {
    return(list())
  }

  out <- stats::setNames(
    lapply(parts, function(x) paste(x[-1], collapse = "=")),
    vapply(parts, `[[`, character(1), 1)
  )
  as.list(out)
}

pid_state <- function(info) {
  pid <- suppressWarnings(as.integer(info$pid %||% NA_integer_))
  if (is.na(pid) || pid <= 0) {
    return("missing")
  }
  if (!requireNamespace("ps", quietly = TRUE)) {
    return("unknown")
  }

  handle <- tryCatch(ps::ps_handle(pid), error = function(e) NULL)
  if (is.null(handle) || !isTRUE(tryCatch(ps::ps_is_running(handle), error = function(e) FALSE))) {
    return("inactive")
  }

  started_at <- suppressWarnings(as.POSIXct(info$started_at %||% NA_character_, tz = Sys.timezone()))
  process_started_at <- tryCatch(ps::ps_create_time(handle), error = function(e) as.POSIXct(NA))
  if (!is.na(started_at) &&
      !is.na(process_started_at) &&
      as.numeric(process_started_at) > as.numeric(started_at) + 2) {
    return("reused")
  }

  "active"
}

message_status <- function() {
  info <- read_pid_file()
  state <- pid_state(info)
  message("wiSDM dashboard status: ", state)
  if (length(info) > 0) {
    message("PID file: ", dashboard_pid_file)
    message("URL: ", info$url %||% dashboard_url)
  }
  invisible(state)
}

run_dashboard_foreground <- function() {
  setwd(repo_root)
  shiny::runApp(
    shiny::shinyAppFile(app_path),
    host = host,
    port = port,
    launch.browser = TRUE
  )
}

run_dashboard_child <- function() {
  setwd(repo_root)

  out_connection <- file(dashboard_out_log, open = "at")
  err_connection <- file(dashboard_err_log, open = "at")
  sink(out_connection, split = FALSE)
  sink(err_connection, type = "message")
  on.exit({
    while (sink.number(type = "message") > 0) sink(type = "message")
    while (sink.number() > 0) sink()
    close(out_connection)
    close(err_connection)
  }, add = TRUE)

  write_pid_file()
  on.exit(unlink(dashboard_pid_file), add = TRUE)

  message("Starting wiSDM dashboard at ", dashboard_url, " at ", timestamp_now())
  shiny::runApp(
    shiny::shinyAppFile(app_path),
    host = host,
    port = port,
    launch.browser = FALSE
  )
}

start_dashboard_detached <- function() {
  if (port_open(host, port)) {
    message("A dashboard server is already responding at ", dashboard_url)
    utils::browseURL(dashboard_url)
    return(invisible(TRUE))
  }

  info <- read_pid_file()
  state <- pid_state(info)
  if (identical(state, "active")) {
    message("The dashboard process is already running at ", info$url %||% dashboard_url)
    utils::browseURL(info$url %||% dashboard_url)
    return(invisible(TRUE))
  }

  if (state %in% c("inactive", "reused")) {
    unlink(dashboard_pid_file)
  }

  args <- c(shQuote(script_path), "--child")
  system2(rscript_path(), args = args, wait = FALSE, invisible = TRUE)

  if (wait_for_dashboard()) {
    utils::browseURL(dashboard_url)
    message("wiSDM dashboard started in a detached R process: ", dashboard_url)
    message("RStudio console Stop is no longer used to stop this dashboard.")
    message("To stop only the dashboard, run: Rscript src/_launch_wisdm_chunk_dashboard.R --stop")
    return(invisible(TRUE))
  }

  warning(
    "The dashboard process was started, but the port did not respond within the wait window. ",
    "Check logs:\n  ",
    dashboard_out_log,
    "\n  ",
    dashboard_err_log,
    call. = FALSE
  )
  invisible(FALSE)
}

stop_dashboard <- function() {
  info <- read_pid_file()
  state <- pid_state(info)

  if (!identical(state, "active")) {
    unlink(dashboard_pid_file)
    message("No active dashboard process was found.")
    return(invisible(FALSE))
  }

  if (!requireNamespace("ps", quietly = TRUE)) {
    stop("Package 'ps' is required to stop the dashboard safely by PID.", call. = FALSE)
  }

  pid <- as.integer(info$pid)
  handle <- ps::ps_handle(pid)
  ps::ps_kill(handle)
  unlink(dashboard_pid_file)
  message("Stopped dashboard process PID ", pid, ".")
  invisible(TRUE)
}

if ("--help" %in% trailing_args) {
  cat(
    "Usage:\n",
    "  source('src/_launch_wisdm_chunk_dashboard.R')\n",
    "  Rscript src/_launch_wisdm_chunk_dashboard.R --status\n",
    "  Rscript src/_launch_wisdm_chunk_dashboard.R --stop\n",
    "  Rscript src/_launch_wisdm_chunk_dashboard.R --foreground\n",
    sep = ""
  )
} else if ("--child" %in% trailing_args) {
  run_dashboard_child()
} else if ("--foreground" %in% trailing_args) {
  run_dashboard_foreground()
} else if ("--status" %in% trailing_args) {
  message_status()
} else if ("--stop" %in% trailing_args) {
  stop_dashboard()
} else {
  start_dashboard_detached()
}
