# ==========================================================================
# Zonal suitability sweep  ---  parallel, instrumented        (v3, fixed)
#
# Fixes over v2:
#   * tryCatch evaluates its expr in the CALLING frame, so `<<-` skipped the
#     local `tmid` and wrote to globalenv -- read/count times were always NA.
#     Now uses an environment (reference semantics, scoping-proof).
#   * terra::rasterize() takes write options via wopt=, not via `...`.
#   * terraOptions(tempdir=) requires the directory to already exist.
#   * The performance report pooled the sequential baseline groups into the
#     parallel statistics. Groups are now tagged by phase.
#   * The baseline is skipped unless enough groups remain to parallelise.
#   * Every conditional wrapped in isTRUE() so NaN cannot abort the report.
# ==========================================================================

suppressPackageStartupMessages({
  library(terra); library(dplyr); library(tibble); library(purrr)
  library(future); library(furrr); library(progressr); library(cli)
})

source("./src/post_modelling_analyses/check_best_model-v1.R")

# ---- 0. configuration ----------------------------------------------------

CFG <- list(
  bounds_gpkg = "./data/external/gadm/europe_selected_countries_wgs84cea_v3-1.gpkg",
  zone_dir    = "./data/interim/zone_rasters",
  out_file    = "./data/processed/suitability_by_country.rds",
  err_file    = "./data/processed/suitability_errors.rds",
  perf_file   = "./data/processed/suitability_perf.rds",
  map_used    = "binary10pct",
  
  workers        = max(1L, parallel::detectCores(logical = FALSE) - 1L),
  mem_cell_limit = 2e7,
  
  # Sequential calibration groups. Automatically skipped when too few groups
  # remain to make the parallel measurement meaningful (see min_par_groups).
  benchmark_groups = 3L,
  
  check_workers = TRUE,   # verify multisession really forks distinct PIDs
  max_tasks     = Inf     # e.g. 140 for a smoke test; needs >= ~2x workers groups
)

dir.create(CFG$zone_dir,          recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(CFG$out_file), recursive = TRUE, showWarnings = FALSE)
terra::tmpFiles(remove = TRUE, orphan = TRUE)

# ---- helpers -------------------------------------------------------------

PHASE <- new.env(parent = emptyenv())
phase <- function(label, expr) {
  t0 <- Sys.time(); val <- force(expr)
  PHASE[[label]] <- as.numeric(difftime(Sys.time(), t0, units = "secs")); val
}
t_script_start <- Sys.time()

fmt_s <- function(s) {
  if (!is.finite(s)) return("n/a")
  if (s < 60) sprintf("%.1fs", s) else
    sprintf("%dm %02ds", as.integer(s) %/% 60L, as.integer(s) %% 60L)
}
pct_of <- function(a, b) if (isTRUE(is.finite(a) && is.finite(b) && b > 0))
  sprintf("%.0f%%", 100 * a / b) else "n/a"
ratio <- function(a, b, digits = 2) if (isTRUE(is.finite(a) && is.finite(b) && b > 0))
  sprintf(paste0("%.", digits, "fx"), a / b) else "n/a"

# ---- progressr / cli handler --------------------------------------------
# cli_progress_bar() cannot be driven from a worker: the bar lives in the
# main session. progressr is the transport; handler_cli renders it through
# cli, so {cli::pb_*} fields still apply.

setup_handlers <- function() {
  ok <- tryCatch({
    progressr::handlers(progressr::handler_cli(
      format      = paste0("{cli::pb_spin} {cli::pb_name} ",
                           "{cli::pb_current}/{cli::pb_total} {cli::pb_bar} ",
                           "{cli::pb_percent} | {cli::pb_rate} | ETA {cli::pb_eta}"),
      format_done = paste0("{cli::col_green(cli::symbol$tick)} {cli::pb_name}: ",
                           "{cli::pb_total} rasters in {cli::pb_elapsed}"),
      clear = FALSE
    )); TRUE
  }, error = function(e) FALSE)
  if (!ok) {
    tryCatch(progressr::handlers("cli"),
             error = function(e) progressr::handlers("txtprogressbar"))
    cli::cli_alert_warning("Custom cli format unsupported by progressr {utils::packageVersion('progressr')}; using default.")
  }
}
setup_handlers()

# ---- 1. zones ------------------------------------------------------------

cli::cli_h1("Setup")
eu_bounds      <- phase("load_bounds", vect(CFG$bounds_gpkg))
eu_bounds$ID   <- seq_len(nrow(eu_bounds))
countries_list <- as.data.frame(eu_bounds)
n_zone         <- nrow(eu_bounds)
cli::cli_alert_info("{n_zone} zone{?s} from {.file {basename(CFG$bounds_gpkg)}}")

# ---- 2. task table -------------------------------------------------------

per_scn <- tibble(
  period   = c(rep(c("mid_2041_2070", "late_2071_2100"), each = 3), "hist"),
  scenario = c(rep(c("ssp126", "ssp370", "ssp585"), times = 2), "hist")
) |> mutate(ps_name = paste(period, scenario, sep = "_"))

get_path <- function(m, p, s) {
  x <- tryCatch(m$paths_list_bin[[p]][[s]], error = function(e) NULL)
  if (length(x) == 0L) NA_character_ else as.character(x)[[1L]]
}

tasks <- phase("build_tasks", {
  map(best_mod_set_paths, function(m) {
    per_scn |> mutate(sp_name = m$sp_name, tkey = m$tkey,
                      fpath = map2_chr(period, scenario, \(p, s) get_path(m, p, s)))
  }) |> list_rbind()
})

bad <- which(is.na(tasks$fpath) | !file.exists(tasks$fpath))
if (length(bad)) {
  cli::cli_alert_warning("{length(bad)}/{nrow(tasks)} raster{?s} missing -- skipping.")
  print(utils::head(tasks[bad, c("sp_name", "ps_name", "fpath")], 10L))
  tasks <- tasks[-bad, , drop = FALSE]
}
# Truncate on whole species so groups stay intact.
if (is.finite(CFG$max_tasks)) {
  keep_sp <- unique(tasks$sp_name)[seq_len(max(1L, ceiling(CFG$max_tasks / nrow(per_scn))))]
  tasks <- tasks[tasks$sp_name %in% keep_sp, , drop = FALSE]
}
stopifnot(nrow(tasks) > 0L)

# ---- 3. zone raster per distinct grid ------------------------------------

geom_key <- function(f) {
  r <- rast(f)
  paste(c(dim(r)[1:2], as.vector(ext(r)), res(r), crs(r, proj = TRUE)), collapse = "|")
}

tasks$geom <- phase("geom_scan", {
  pb <- cli::cli_progress_bar(
    name  = "Scanning geometries", total = nrow(tasks),
    format = "{cli::pb_spin} {cli::pb_name} {cli::pb_current}/{cli::pb_total} {cli::pb_bar} {cli::pb_percent} | ETA {cli::pb_eta}",
    clear = TRUE)
  out <- character(nrow(tasks))
  for (i in seq_len(nrow(tasks))) {
    out[i] <- geom_key(tasks$fpath[[i]]); cli::cli_progress_update(id = pb, set = i)
  }
  cli::cli_progress_done(id = pb); out
})

zone_lookup <- tasks |> distinct(geom, .keep_all = TRUE) |>
  transmute(geom, template = fpath,
            zone_file = file.path(CFG$zone_dir,
                                  sprintf("zones_%s.tif",
                                          substr(vapply(geom, rlang::hash, character(1L)), 1L, 10L))))

cli::cli_alert_info("{nrow(zone_lookup)} distinct geometr{?y/ies} across {nrow(tasks)} file{?s}.")

phase("rasterize_zones", {
  todo <- which(!file.exists(zone_lookup$zone_file))
  if (length(todo)) {
    for (i in todo) {
      cli::cli_alert_info("Rasterizing {basename(zone_lookup$zone_file[[i]])}")
      # NOTE: write options go in wopt=, NOT in `...`  (terra warns otherwise)
      rasterize(eu_bounds, rast(zone_lookup$template[[i]]), field = "ID",
                filename = zone_lookup$zone_file[[i]], overwrite = TRUE,
                wopt = list(datatype = "INT2U",
                            gdal = c("COMPRESS=DEFLATE", "TILED=YES", "PREDICTOR=2")))
    }
  } else cli::cli_alert_success("Zone grids cached.")
})

tasks <- left_join(tasks, select(zone_lookup, geom, zone_file), by = "geom")

# ---- 4. counting kernel (runs in workers) --------------------------------

count_group <- function(df, n_zone, mem_cell_limit, n_workers, group_id,
                        exec = "parallel", pgr = NULL) {
  
  # A worker does not own the machine. nbrOfWorkers() returns 1 inside a
  # worker, so the count is passed in. tempdir must exist before it is set.
  td <- file.path(tempdir(), paste0("terra_", Sys.getpid()))
  dir.create(td, recursive = TRUE, showWarnings = FALSE)
  terra::terraOptions(memfrac = max(0.1, 0.6 / max(1L, n_workers)),
                      tempdir = td, progress = 0)
  Sys.setenv(GDAL_NUM_THREADS = "1")
  
  pid <- Sys.getpid(); g_t0 <- as.numeric(Sys.time())
  zr <- terra::rast(df$zone_file[[1L]])
  use_mem <- terra::ncell(zr) <= mem_cell_limit
  
  z_t0 <- as.numeric(Sys.time())
  if (use_mem) {
    z <- as.integer(terra::values(zr, mat = FALSE))
    keep_z <- !is.na(z); zz <- z[keep_z]; rm(z); gc(verbose = FALSE)
  }
  secs_zone <- as.numeric(Sys.time()) - z_t0
  
  rows <- lapply(seq_len(nrow(df)), function(i) {
    row <- df[i, , drop = FALSE]
    tk  <- new.env(parent = emptyenv())   # scoping-proof timing slot
    tk$mid <- NA_real_
    t0  <- as.numeric(Sys.time())
    
    counts <- tryCatch({
      xr <- terra::rast(row$fpath[[1L]])
      if (terra::nlyr(xr) != 1L) stop("expected 1 layer, got ", terra::nlyr(xr))
      
      if (use_mem) {
        x <- terra::values(xr, mat = FALSE)[keep_z]
        tk$mid <- as.numeric(Sys.time())        # read finished
        ok <- !is.na(x); zk <- zz[ok]; xk <- as.integer(x[ok]); rm(x, ok)
        list(total = tabulate(zk,           nbins = n_zone),
             n1    = tabulate(zk[xk == 1L], nbins = n_zone),
             n0    = tabulate(zk[xk == 0L], nbins = n_zone))
      } else {
        s <- terra::zonal(xr, zr, fun = "sum", na.rm = TRUE)
        tk$mid <- as.numeric(Sys.time())
        nn <- terra::zonal(xr, zr, fun = "notNA")
        total <- n1 <- integer(n_zone); idx <- as.integer(s[[1L]])
        n1[idx] <- as.integer(s[[2L]]); total[idx] <- as.integer(nn[[2L]])
        list(total = total, n1 = n1, n0 = total - n1)
      }
    }, error = function(e) conditionMessage(e))
    
    t1 <- as.numeric(Sys.time()); good <- is.list(counts)
    if (!is.null(pgr)) pgr(message = sprintf("%s / %s", row$sp_name[[1L]], row$ps_name[[1L]]))
    
    list(
      cnt = data.frame(sp_name = row$sp_name, tkey = row$tkey, period = row$period,
                       scenario = row$scenario, ps_name = row$ps_name, id = seq_len(n_zone),
                       n0 = if (good) counts$n0 else NA_integer_,
                       n1 = if (good) counts$n1 else NA_integer_,
                       total = if (good) counts$total else NA_integer_,
                       error = if (good) NA_character_ else counts, stringsAsFactors = FALSE),
      tim = data.frame(group_id = group_id, pid = pid, exec = exec,
                       sp_name = row$sp_name, ps_name = row$ps_name,
                       mode = if (use_mem) "memory" else "out-of-core",
                       t_start = t0, t_end = t1,
                       secs_read  = tk$mid - t0,
                       secs_count = t1 - tk$mid,
                       secs_total = t1 - t0,
                       ncell = terra::ncell(zr),
                       bytes = suppressWarnings(file.size(row$fpath[[1L]])),
                       ok = good, stringsAsFactors = FALSE))
  })
  
  g_t1 <- as.numeric(Sys.time())
  list(counts  = do.call(rbind, lapply(rows, `[[`, "cnt")),
       timings = do.call(rbind, lapply(rows, `[[`, "tim")),
       group   = data.frame(group_id = group_id, pid = pid, exec = exec,
                            n_rasters = nrow(df), t_start = g_t0, t_end = g_t1,
                            secs_zone = secs_zone, secs_group = g_t1 - g_t0,
                            stringsAsFactors = FALSE))
}

# ---- 5. run --------------------------------------------------------------

groups <- tasks |> group_by(sp_name, geom) |> group_split()

# The baseline steals groups from the parallel phase. Only worth it when
# enough remain to saturate the workers -- otherwise you measure noise.
min_par_groups <- 2L * CFG$workers
n_bench <- if (length(groups) >= min_par_groups + CFG$benchmark_groups)
  CFG$benchmark_groups else 0L
bench_i <- seq_len(n_bench); par_i <- setdiff(seq_along(groups), bench_i)

cli::cli_h1("Sweep")
cli::cli_alert_info("{nrow(tasks)} raster{?s} / {length(groups)} group{?s} / {CFG$workers} worker{?s}")
if (n_bench == 0L && CFG$benchmark_groups > 0L)
  cli::cli_alert_warning(
    "Skipping the sequential baseline: needs >= {min_par_groups + CFG$benchmark_groups} groups, have {length(groups)}. Speedup will be reported as effective concurrency only.")

## 5a. sequential calibration ----------------------------------------------
res_bench <- list(); bench_rasters <- 0L; PHASE[["benchmark"]] <- 0
if (n_bench > 0L) {
  bench_rasters <- sum(vapply(groups[bench_i], nrow, integer(1)))
  pb <- cli::cli_progress_bar(
    name = "Sequential baseline", total = bench_rasters,
    format = "{cli::pb_spin} {cli::pb_name} {cli::pb_current}/{cli::pb_total} {cli::pb_bar} {cli::pb_percent} | ETA {cli::pb_eta}",
    clear = TRUE)
  done <- 0L
  PHASE[["benchmark"]] <- as.numeric(system.time({
    res_bench <- lapply(bench_i, function(i) {
      r <- count_group(groups[[i]], n_zone, CFG$mem_cell_limit, n_workers = 1L,
                       group_id = i, exec = "sequential", pgr = NULL)
      done <<- done + nrow(groups[[i]]); cli::cli_progress_update(id = pb, set = done); r
    })
  })[["elapsed"]])
  cli::cli_progress_done(id = pb)
  cli::cli_alert_success("Baseline: {bench_rasters} rasters single-threaded in {fmt_s(PHASE[['benchmark']])}")
}

## 5b. parallel sweep ------------------------------------------------------
plan(multisession, workers = CFG$workers)

if (isTRUE(CFG$check_workers)) {
  probe <- unlist(future_map(seq_len(CFG$workers), function(i) {
    Sys.sleep(0.05); Sys.getpid()
  }, .options = furrr_options(seed = NULL, chunk_size = 1L)))
  wp <- sort(unique(probe))
  cli::cli_alert_info("Main PID {Sys.getpid()}; worker PIDs {paste(wp, collapse = ', ')}")
  if (isTRUE(any(probe == Sys.getpid())))
    cli::cli_alert_danger("Futures are executing in the MAIN session -- plan(multisession) did not take effect.")
  else if (isTRUE(length(wp) < CFG$workers))
    cli::cli_alert_warning("Only {length(wp)} distinct worker{?s} seen; probe may not have reached all of them.")
}

res_par <- list(); n_par_rasters <- 0L; PHASE[["parallel_sweep"]] <- 0
if (length(par_i)) {
  n_par_rasters <- sum(vapply(groups[par_i], nrow, integer(1)))
  PHASE[["parallel_sweep"]] <- as.numeric(system.time({
    res_par <- with_progress({
      pgr <- progressor(steps = n_par_rasters)
      future_map(par_i,
                 function(i) count_group(groups[[i]], n_zone, CFG$mem_cell_limit,
                                         n_workers = CFG$workers, group_id = i,
                                         exec = "parallel", pgr = pgr),
                 # one group per future: best load balance and finest progress relay
                 .options = furrr_options(seed = NULL, packages = "terra", chunk_size = 1L))
    })
  })[["elapsed"]])
}
plan(sequential)

res <- c(res_bench, res_par)

# ---- 6. assemble ---------------------------------------------------------

suitability <- phase("assemble", {
  map(res, "counts") |> list_rbind() |>
    mutate(map_used = CFG$map_used,
           pct_suitb = if_else(total > 0L, round(100 * n1 / total, 3), NA_real_)) |>
    left_join(countries_list, by = c("id" = "ID")) |>
    relocate(sp_name, tkey, map_used, period, scenario, ps_name, id)
})
timings <- map(res, "timings") |> list_rbind()
gtimes  <- map(res, "group")   |> list_rbind()

odd <- suitability |> filter(!is.na(total), n0 + n1 != total)
if (nrow(odd)) cli::cli_alert_warning("{nrow(odd)} row{?s} where n0+n1 != total: rasters hold values other than 0/1.")

failed <- suitability |> filter(!is.na(error)) |> distinct(sp_name, ps_name, error)
if (nrow(failed)) {
  cli::cli_alert_danger("{nrow(failed)} raster{?s} failed -> {.file {CFG$err_file}}")
  saveRDS(failed, CFG$err_file)
}
saveRDS(suitability, CFG$out_file)
cli::cli_alert_success("Wrote {.file {CFG$out_file}} ({nrow(suitability)} rows)")

# ---- 7. performance report -----------------------------------------------

report_performance <- function(timings, gtimes, phases, workers,
                               bench_secs, bench_rasters,
                               par_secs, par_rasters, total_secs) {
  
  cli::cli_h1("Performance report")
  
  ph <- data.frame(phase = names(phases), secs = unlist(as.list(phases)),
                   row.names = NULL)
  ph <- ph[ph$secs > 0, ]; ph <- ph[order(-ph$secs), ]
  cli::cli_h2("Where the wall clock went")
  cli::cli_verbatim(paste(sprintf("  %-18s %8s  %6s", ph$phase,
                                  vapply(ph$secs, fmt_s, ""), vapply(ph$secs, pct_of, "", b = total_secs)),
                          collapse = "\n"))
  cli::cli_verbatim(sprintf("  %-18s %8s", "TOTAL", fmt_s(total_secs)))
  
  gp <- gtimes[gtimes$exec == "parallel", ]
  tp <- timings[timings$exec == "parallel" & timings$ok, ]
  ts <- timings[timings$exec == "sequential" & timings$ok, ]
  
  conc <- NA_real_
  if (nrow(gp) == 0L) {
    cli::cli_alert_warning("No groups ran in parallel; nothing to report.")
    return(invisible(NULL))
  }
  
  busy <- sum(gp$secs_group)
  span <- max(gp$t_end) - min(gp$t_start)
  conc <- if (isTRUE(span > 0)) busy / span else NA_real_
  
  cli::cli_h2("Parallel behaviour")
  cli::cli_dl(c(
    "Groups run in parallel" = "{nrow(gp)}",
    "Workers requested"      = "{workers}",
    "Distinct worker PIDs"   = "{length(unique(gp$pid))}",
    "Worker busy time"       = "{fmt_s(busy)}",
    "Parallel span"          = "{fmt_s(span)}",
    "Effective concurrency"  = "{ratio(busy, span)}",
    "Parallel efficiency"    = "{pct_of(conc, workers)}"))
  if (isTRUE(nrow(gp) < workers))
    cli::cli_alert_warning("Fewer groups than workers -- concurrency is capped at {nrow(gp)}x by the task count, not by the hardware.")
  
  per_w <- tapply(gp$secs_group, gp$pid, sum)
  ng    <- table(gp$pid)[names(per_w)]
  cli::cli_h2("Load balance")
  cli::cli_verbatim(paste(sprintf("  pid %-8s %8s  (%d group%s)", names(per_w),
                                  vapply(as.numeric(per_w), fmt_s, ""), as.integer(ng),
                                  ifelse(as.integer(ng) == 1L, "", "s")), collapse = "\n"))
  sr <- max(per_w) / stats::median(per_w)
  cli::cli_alert_info("Straggler ratio (max/median worker time): {sprintf('%.2f', sr)}")
  if (isTRUE(sr > 1.4))
    cli::cli_alert_warning("Uneven load -- the slowest worker's tail is wasted wall clock.")
  
  if (nrow(tp)) {
    rd <- sum(tp$secs_read, na.rm = TRUE); ct <- sum(tp$secs_count, na.rm = TRUE)
    mb <- sum(tp$bytes, na.rm = TRUE) / 2^20
    cli::cli_h2("Work profile")
    cli::cli_dl(c(
      "Read / decompress" = "{fmt_s(rd)} ({pct_of(rd, rd + ct)})",
      "Counting"          = "{fmt_s(ct)} ({pct_of(ct, rd + ct)})",
      "Throughput"        = "{sprintf('%.2f rasters/s', par_rasters / par_secs)}",
      "Read bandwidth"    = "{sprintf('%.0f MB/s aggregate (compressed bytes)', mb / par_secs)}",
      "Cell rate"         = "{sprintf('%.1f Mcell/s', sum(as.numeric(tp$ncell)) / par_secs / 1e6)}"))
    if (isTRUE(rd / (rd + ct) > 0.7))
      cli::cli_alert_warning("I/O dominated -- workers beyond disk bandwidth will raise concurrency but not throughput.")
  }
  
  cli::cli_h2("Speedup")
  if (isTRUE(bench_rasters > 0L && bench_secs > 0)) {
    seq_rate <- bench_rasters / bench_secs
    par_rate <- par_rasters / par_secs
    contention <- if (nrow(ts) && nrow(tp))
      mean(tp$secs_total) / mean(ts$secs_total) else NA_real_
    cli::cli_dl(c(
      "Sequential baseline" = "{sprintf('%.2f rasters/s', seq_rate)} (n = {bench_rasters}, 1 thread)",
      "Parallel throughput" = "{sprintf('%.2f rasters/s', par_rate)} (n = {par_rasters})",
      "Measured speedup"    = "{cli::col_green(ratio(par_rate, seq_rate))}",
      "Contention factor"   = "{if (is.finite(contention)) sprintf('%.2fx slower per raster under load', contention) else 'n/a'}",
      "Projected serial run"= "{fmt_s(nrow(timings) / seq_rate)}"))
    if (isTRUE(contention > 1.15))
      cli::cli_alert_warning("Workers are slowing each other by {pct_of(contention - 1, 1)} per raster -- memory bandwidth or disk queue is saturated.")
    cli::cli_alert_info("Caveat: the baseline runs first, so OS page cache state differs between the two measurements.")
  } else {
    cli::cli_alert_info("No sequential baseline; effective concurrency ({ratio(busy, span)}) is the best available estimate.")
  }
  
  serial_secs <- total_secs - par_secs
  est_serial  <- serial_secs + busy
  cli::cli_h2("Amdahl ceiling")
  cli::cli_dl(c(
    "Non-parallel portion" = "{fmt_s(serial_secs)} ({pct_of(serial_secs, total_secs)} of this run)",
    "Overall speedup"      = "{ratio(est_serial, total_secs)} vs an all-serial run",
    "Ceiling"              = "{ratio(est_serial, serial_secs)} with infinite workers"))
  if (isTRUE(serial_secs / total_secs > 0.25))
    cli::cli_alert_warning("Over a quarter of the run is serial setup -- optimise that before adding cores.")
  
  invisible(list(phases = ph, timings = timings, groups = gtimes, concurrency = conc))
}

total_secs <- as.numeric(difftime(Sys.time(), t_script_start, units = "secs"))
perf <- report_performance(timings, gtimes, PHASE, CFG$workers,
                           bench_secs = PHASE[["benchmark"]], bench_rasters = bench_rasters,
                           par_secs = PHASE[["parallel_sweep"]], par_rasters = n_par_rasters,
                           total_secs = total_secs)
saveRDS(perf, CFG$perf_file)
cli::cli_alert_success("Timing detail -> {.file {CFG$perf_file}}")

# ---- 8. optional: guard against tiny zones -------------------------------
# suitability |> mutate(pct_suitb = if_else(total < 500L, NA_real_, pct_suitb))