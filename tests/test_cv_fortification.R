suppressPackageStartupMessages({
  library(dplyr)
  library(sf)
  library(terra)
})

source(file.path("src", "helper_functions.R"))

failures <- character(0)
tests_run <- 0L

assert_true <- function(condition, message) {
  if (!isTRUE(condition)) {
    stop(message, call. = FALSE)
  }
}

expect_error <- function(expression, pattern = NULL) {
  error <- tryCatch(
    {
      force(expression)
      NULL
    },
    error = function(e) e
  )
  assert_true(inherits(error, "error"), "Expected an error, but none was raised.")
  if (!is.null(pattern)) {
    assert_true(
      grepl(pattern, conditionMessage(error), ignore.case = TRUE),
      paste0("Error did not match '", pattern, "': ", conditionMessage(error))
    )
  }
  invisible(error)
}

run_test <- function(name, expression) {
  tests_run <<- tests_run + 1L
  tryCatch(
    {
      force(expression)
      message("PASS: ", name)
    },
    error = function(e) {
      failures <<- c(failures, paste0(name, ": ", conditionMessage(e)))
      message("FAIL: ", name, " -- ", conditionMessage(e))
    }
  )
}

make_spatial_records <- function(n_sites,
                                 spacing = 200000,
                                 contexts = "europe") {
  grid_width <- ceiling(sqrt(n_sites))
  coordinates <- expand.grid(
    x = seq(0, by = spacing, length.out = grid_width),
    y = seq(0, by = spacing, length.out = grid_width)
  )[seq_len(n_sites), , drop = FALSE]
  records <- do.call(
    rbind,
    lapply(contexts, function(context) {
      do.call(
        rbind,
        lapply(0:1, function(species_value) {
          data.frame(
            coordinates,
            species = species_value,
            cv_context = context,
            cv_stratum = paste(context, species_value, sep = "::"),
            cv_group = paste(coordinates$x, coordinates$y, sep = "::"),
            stringsAsFactors = FALSE
          )
        })
      )
    })
  )
  sf::st_as_sf(records, coords = c("x", "y"), crs = 3035)
}

run_test("valid 100 km spatial CV", {
  records <- make_spatial_records(36L)
  plan <- make_spatial_cv_plan(
    records,
    requested_folds = 5L,
    block_sizes_m = 100000,
    min_train = 2L,
    min_test = 2L,
    seed = 41L,
    iteration = 50L
  )
  assert_true(plan$valid, plan$fallback_reason)
  assert_true(identical(plan$method, "spatial_block"), "Spatial plan was not selected.")
  assert_true(identical(as.numeric(plan$block_size_m), 100000), "Wrong block size selected.")
  assert_true(!anyNA(plan$fold_ids), "Some spatial records were left unassigned.")
})

run_test("spatial block-size adaptation", {
  records <- make_spatial_records(36L)
  plan <- make_spatial_cv_plan(
    records,
    requested_folds = 5L,
    block_sizes_m = c(5000000, 100000),
    min_train = 2L,
    min_test = 2L,
    seed = 42L,
    iteration = 50L
  )
  assert_true(plan$valid, plan$fallback_reason)
  assert_true(identical(as.numeric(plan$block_size_m), 100000),
              "The smaller valid spatial block size was not selected.")
  assert_true(any(!plan$diagnostics$valid), "The invalid larger-block attempt was not recorded.")
})

run_test("spatial fold-count adaptation", {
  records <- make_spatial_records(8L)
  plan <- make_spatial_cv_plan(
    records,
    requested_folds = 5L,
    block_sizes_m = 100000,
    min_train = 2L,
    min_test = 2L,
    seed = 43L,
    iteration = 100L
  )
  assert_true(plan$valid, plan$fallback_reason)
  assert_true(identical(plan$k, 4L), "The requested fold count was not reduced to four.")
})

run_test("spatial failure falls back to grouped stratified k-fold", {
  records <- make_spatial_records(24L)
  plan <- make_preferred_cv_plan(
    records,
    requested_folds = 5L,
    block_sizes_m = 1e9,
    enable_kfold_fallback = TRUE,
    min_train = 2L,
    min_test = 2L,
    seed = 44L,
    iteration = 50L
  )
  assert_true(plan$valid, plan$fallback_reason)
  assert_true(identical(plan$method, "stratified_kfold"), "K-fold fallback was not selected.")
  assert_true(!is.na(plan$fallback_reason), "The spatial fallback reason was not retained.")
})

run_test("disabled fallback reports not evaluable", {
  records <- make_spatial_records(24L)
  plan <- make_preferred_cv_plan(
    records,
    requested_folds = 5L,
    block_sizes_m = 1e9,
    enable_kfold_fallback = FALSE,
    min_train = 2L,
    min_test = 2L,
    seed = 45L,
    iteration = 20L
  )
  assert_true(!plan$valid, "An impossible spatial plan was marked valid.")
  assert_true(identical(plan$method, "not_evaluable"), "Invalid plan was not labelled not_evaluable.")
})

run_test("classes too small for two folds are not evaluable", {
  records <- data.frame(
    species = rep(0:1, each = 3L),
    cv_stratum = rep(c("europe::0", "europe::1"), each = 3L),
    cv_group = seq_len(6L)
  )
  plan <- make_stratified_kfold_plan(
    records,
    requested_folds = 5L,
    min_train = 2L,
    min_test = 2L,
    iteration = 20L
  )
  assert_true(!plan$valid, "A class with three records supported two-fold CV unexpectedly.")
  assert_true(identical(plan$method, "not_evaluable"), "Invalid k-fold plan has the wrong label.")
})

run_test("missing context classes are rejected before partitioning", {
  records <- data.frame(
    species = c(rep(0:1, each = 10L), rep(1L, 10L)),
    cv_context = c(rep("global", 20L), rep("europe", 10L)),
    cv_stratum = c(
      rep(c("global_0", "global_1"), each = 10L),
      rep("europe_1", 10L)
    ),
    cv_group = seq_len(30L)
  )
  plan <- make_stratified_kfold_plan(
    records,
    requested_folds = 5L,
    min_train = 2L,
    min_test = 2L,
    iteration = 20L
  )
  assert_true(!plan$valid, "A context lacking absences was accepted.")
  assert_true(grepl("europe class 0", plan$fallback_reason, fixed = TRUE),
              "The missing context/class was not diagnosed.")
})

run_test("grouped duplicate records remain in one fold", {
  records <- as.data.frame(make_spatial_records(24L))
  records <- rbind(records, records[c(1L, 25L), , drop = FALSE])
  plan <- make_stratified_kfold_plan(
    records,
    requested_folds = 5L,
    min_train = 2L,
    min_test = 2L,
    seed = 46L,
    iteration = 50L
  )
  assert_true(plan$valid, plan$fallback_reason)
  folds_by_group <- split(plan$fold_ids, records$cv_group)
  assert_true(all(vapply(folds_by_group, function(x) length(unique(x)) == 1L, logical(1))),
              "A duplicate/coincident group crossed train/test folds.")
})

run_test("raster NA filtering reports class-specific dropped records", {
  raster <- terra::rast(nrows = 2, ncols = 2, xmin = 0, xmax = 2,
                        ymin = 0, ymax = 2, crs = "EPSG:3857")
  terra::values(raster) <- c(NA, 2, 3, 4)
  names(raster) <- "predictor"
  xy <- terra::xyFromCell(raster, seq_len(4L))
  points <- sf::st_as_sf(
    data.frame(ID = paste0("P", 1:4), species = c(1L, 0L, 1L, 0L), xy),
    coords = c("x", "y"),
    crs = 3857
  )
  extracted <- extract_env(points, raster)
  assert_true(nrow(extracted$complete) == 3L, "The incomplete raster row was not removed.")
  assert_true(nrow(extracted$dropped) == 1L, "Dropped-record diagnostics are incomplete.")
  assert_true(extracted$dropped$species[[1]] == 1L, "The dropped class was not retained.")
})

run_test("rare European contexts remain independently plannable", {
  for (n_presence in c(21L, 25L, 30L)) {
    global <- data.frame(
      species = rep(0:1, each = 40L),
      cv_context = "global_climate",
      cv_group = paste0("global_", rep(seq_len(40L), 2L))
    )
    europe <- data.frame(
      species = rep(0:1, each = n_presence),
      cv_context = "europe_climate",
      cv_group = paste0("europe_", rep(seq_len(n_presence), 2L))
    )
    climate_records <- rbind(global, europe)
    climate_records$cv_stratum <- paste(
      climate_records$cv_context,
      climate_records$species,
      sep = "::"
    )
    climate_plan <- make_stratified_kfold_plan(
      climate_records,
      requested_folds = 5L,
      min_train = 2L,
      min_test = 2L,
      seed = n_presence,
      iteration = 50L
    )

    european_records <- climate_records
    european_records$cv_context[european_records$cv_context == "europe_climate"] <-
      "europe_components"
    european_records$cv_stratum <- paste(
      european_records$cv_context,
      european_records$species,
      sep = "::"
    )
    component_plan <- make_stratified_kfold_plan(
      european_records,
      requested_folds = 5L,
      min_train = 2L,
      min_test = 2L,
      seed = n_presence,
      iteration = 50L
    )

    assert_true(climate_plan$valid, paste("Climate plan failed for", n_presence, "presences."))
    assert_true(component_plan$valid, paste("Habitat/combined plan failed for", n_presence, "presences."))
    europe_rows <- european_records$cv_context == "europe_components"
    assert_true(
      identical(
        climate_plan$fold_ids[climate_records$cv_context == "europe_climate"],
        component_plan$fold_ids[europe_rows]
      ),
      "European habitat and combined assignments are not reusable."
    )
  }
})

predict.fake_cv_model <- function(object, newdata, method, ...) {
  n <- nrow(newdata)
  switch(
    method,
    good_1 = data.frame(probability = seq(0.1, 0.9, length.out = n)),
    good_2 = data.frame(probability = seq(0.2, 0.8, length.out = n)),
    good_3 = data.frame(probability = seq(0.3, 0.7, length.out = n)),
    error = stop("synthetic algorithm failure", call. = FALSE),
    empty = data.frame(probability = numeric(0)),
    wrong_length = data.frame(probability = rep(0.5, max(n - 1L, 0L))),
    non_finite = data.frame(probability = c(rep(0.5, max(n - 1L, 0L)), NA_real_)),
    stop("unknown method", call. = FALSE)
  )
}

run_test("prediction helper tolerates two failed algorithms", {
  model <- structure(list(), class = "fake_cv_model")
  datasets <- list(
    presence = data.frame(ID = paste0("P", 1:5), predictor = 1:5),
    absence = data.frame(ID = paste0("A", 1:5), predictor = 6:10)
  )
  result <- compute_median_favourability_safe(
    model,
    datasets,
    top5_methods = c("good_1", "good_2", "good_3", "error", "empty"),
    prev_ratio = 1,
    min_successful_methods = 3L
  )
  assert_true(result$valid, "Three successful methods did not produce an ensemble prediction.")
  assert_true(result$success_count == 3L, "The successful-method count is wrong.")
  assert_true(setequal(result$failed_methods, c("error", "empty")), "Failed methods were not recorded.")
  assert_true(any(grepl("synthetic algorithm failure", result$method_diagnostics$error,
                        fixed = TRUE), na.rm = TRUE), "The original prediction error was lost.")
})

run_test("prediction helper rejects wrong-length and non-finite outputs", {
  model <- structure(list(), class = "fake_cv_model")
  datasets <- list(
    presence = data.frame(ID = paste0("P", 1:5), predictor = 1:5)
  )
  result <- compute_median_favourability_safe(
    model,
    datasets,
    top5_methods = c("good_1", "good_2", "wrong_length", "non_finite", "error"),
    prev_ratio = 1,
    min_successful_methods = 3L
  )
  assert_true(!result$valid, "Only two valid algorithms incorrectly met the threshold.")
  diagnostics <- result$method_diagnostics$error
  assert_true(any(grepl("prediction length", diagnostics), na.rm = TRUE),
              "Wrong-length prediction was not diagnosed.")
  assert_true(any(grepl("non-finite", diagnostics), na.rm = TRUE),
              "Non-finite prediction was not diagnosed.")
})

run_test("validation metrics always return the complete schema", {
  metrics <- compute_validation_metrics(
    species = "Synthetic species",
    type = "Climate",
    region = "Europe",
    fold = 1L,
    all_suit_vals = c(0.1, 0.2),
    occ_suit_vals = numeric(0),
    abs_suit_vals = c(0.1, 0.2)
  )
  expected <- c(
    "Species", "Type", "Region", "test_fold", "n_pres", "n_abs",
    "auc", "boyce", "tss", "sens", "spec"
  )
  assert_true(identical(names(metrics), expected), "Metric schema changed for an unevaluable fold.")
  assert_true(nrow(metrics) == 1L && is.na(metrics$auc), "Unevaluable metrics were not returned as NA.")
})

run_test("validation summaries retain explicit CV metadata", {
  metrics <- dplyr::bind_rows(lapply(1:2, function(fold) {
    data.frame(
      Species = "Synthetic species",
      Type = "Climate",
      Region = "Global",
      test_fold = fold,
      n_pres = 2,
      n_abs = 2,
      auc = 0.75,
      boyce = NA_real_,
      tss = 0.5,
      sens = 0.5,
      spec = 1,
      cv_method = "spatial_block",
      requested_folds = 5L,
      effective_folds = 2L,
      block_size_m = 75000,
      fallback_reason = NA_character_
    )
  }))
  summary <- summarise_validation(metrics)
  required <- c(
    "cv_method", "requested_folds", "effective_folds", "block_size_m",
    "fallback_reason", "validation"
  )
  assert_true(all(required %in% names(summary)), "CV metadata was dropped from the summary.")
  assert_true(summary$n_folds[[1]] == 2L, "The effective evaluated fold count is wrong.")
  assert_true(identical(summary$validation[[1]], "spatial_block"),
              "The validation label is ambiguous.")
})

run_test("combined validation summary preserves the historical interface", {
  extended <- data.frame(
    Species = c("Historical_species", "Fallback_species", "Tiny_species"),
    Type = c("Climate", "Habitat", "Ensemble"),
    Region = c("Global", "Europe", "Europe"),
    cv_method = c("spatial_block", "stratified_kfold", "not_evaluable"),
    requested_folds = c(5L, 5L, 5L),
    effective_folds = c(5L, 3L, NA_integer_),
    block_size_m = c(100000, NA_real_, NA_real_),
    fallback_reason = c(NA_character_, "spatial plan invalid", "too few records"),
    n_folds = c(5L, 3L, 0L),
    mean_auc = c(0.9, 0.8, NA_real_),
    sd_auc = c(0.01, 0.02, NA_real_),
    mean_boyce = c(0.8, 0.7, NA_real_),
    sd_boyce = c(0.03, 0.04, NA_real_),
    mean_tss = c(0.7, 0.6, NA_real_),
    sd_tss = c(0.05, 0.06, NA_real_),
    mean_sens = c(0.8, 0.7, NA_real_),
    sd_sens = c(0.07, 0.08, NA_real_),
    mean_spec = c(0.9, 0.8, NA_real_),
    sd_spec = c(0.09, 0.1, NA_real_),
    validation = c("spatial_block", "stratified_kfold", "not_evaluable")
  )
  legacy <- as_legacy_validation_summary(extended)
  historical_columns <- c(
    "Species", "Type", "Region", "n_folds", "mean_auc", "sd_auc",
    "mean_boyce", "sd_boyce", "mean_tss", "sd_tss", "mean_sens",
    "sd_sens", "mean_spec", "sd_spec", "validation"
  )
  assert_true(identical(names(legacy), historical_columns),
              "The combined summary no longer has the historical columns in order.")
  assert_true(
    identical(
      legacy$validation,
      c("Cross-validation", "Cross-validation", "Not evaluable")
    ),
    "Legacy validation labels are incompatible."
  )
  assert_true(is.na(legacy$n_folds[[3]]),
              "An unevaluable row should not report zero completed folds in the legacy file.")
})

run_test("ensemble matching uses stable IDs", {
  habitat <- data.frame(ID = c("P1", "A1", "orphan_h"), median_favourability = c(0.4, 0.2, 0.8))
  climate <- data.frame(ID = c("A1", "P1", "orphan_c"), median_favourability = c(0.8, 0.9, 0.1))
  result <- suppressMessages(ensemble_geom_mean(habitat, climate, "test", return_data = TRUE))
  assert_true(identical(result$ID, c("P1", "A1")), "Predictions were not matched by stable ID.")
  assert_true(nrow(result) == 2L && all(is.finite(result$ensemble_favourability)),
              "Matched ensemble predictions are invalid.")
})

run_test("portable current-manifest resolution and signatures", {
  directory <- tempfile("cv-manifest-")
  dir.create(directory, recursive = TRUE)
  on.exit(unlink(directory, recursive = TRUE), add = TRUE)

  raster_1 <- terra::rast(nrows = 3, ncols = 4, xmin = 0, xmax = 4,
                          ymin = 0, ymax = 3, crs = "EPSG:3857")
  raster_2 <- raster_1
  terra::values(raster_1) <- 1:12
  terra::values(raster_2) <- 12:1
  file_1 <- file.path(directory, "bio1.tif")
  file_2 <- file.path(directory, "bio2.tif")
  terra::writeRaster(raster_1, file_1, overwrite = TRUE)
  terra::writeRaster(raster_2, file_2, overwrite = TRUE)

  manifest_path <- file.path(directory, "manifest.csv")
  utils::write.csv(
    data.frame(
      period = c("current", "current", "future"),
      scenario = c("current", "current", "ssp999"),
      var_name = c("bio1", "bio2", "bio1"),
      file_path = c("bio1.tif", "bio2.tif", "missing_future.tif")
    ),
    manifest_path,
    row.names = FALSE
  )

  manifest <- load_user_specific_current_manifest(
    manifest_path,
    "test_manifest",
    "climate"
  )
  assert_true(nrow(manifest$current_rows) == 2L,
              "Unavailable future projections blocked or altered current validation.")
  stack <- load_named_raster_stack(manifest$current_rows)
  signature <- build_predictor_signature(stack, c("bio1", "bio2"))

  saved_result <- resolve_model_current_predictors(
    saved_manifest_path = manifest_path,
    configured_manifest_path = file.path(directory, "does-not-exist.csv"),
    config_name = "test_manifest",
    data_label = "climate",
    selected_predictors = c("bio1", "bio2"),
    saved_crs = terra::crs(stack),
    saved_signature = signature
  )
  assert_true(identical(saved_result$source, "saved"), "A valid saved manifest was not preferred.")

  relocated_result <- suppressWarnings(resolve_model_current_predictors(
    saved_manifest_path = file.path(directory, "old-machine.csv"),
    configured_manifest_path = manifest_path,
    config_name = "test_manifest",
    data_label = "climate",
    selected_predictors = c("bio1", "bio2"),
    saved_crs = terra::crs(stack),
    saved_signature = NULL
  ))
  assert_true(identical(relocated_result$source, "configured"),
              "A compatible relocated manifest was not selected.")
  assert_true(relocated_result$compatibility$legacy,
              "Legacy predictor validation was not identified as weaker.")

  bad_path_manifest <- file.path(directory, "missing-raster.csv")
  utils::write.csv(
    data.frame(period = "current", scenario = "current", var_name = "bio1",
               file_path = "absent.tif"),
    bad_path_manifest,
    row.names = FALSE
  )
  expect_error(
    load_user_specific_current_manifest(bad_path_manifest, "test_manifest", "climate"),
    "do not exist"
  )
  expect_error(
    resolve_model_current_predictors(
      file.path(directory, "old-machine.csv"),
      manifest_path,
      "test_manifest",
      "climate",
      selected_predictors = "not_a_predictor",
      saved_crs = terra::crs(stack),
      saved_signature = NULL
    ),
    "missing selected predictor"
  )

  changed_grid <- terra::rast(nrows = 3, ncols = 5, xmin = 0, xmax = 5,
                              ymin = 0, ymax = 3, crs = "EPSG:3857")
  terra::values(changed_grid) <- seq_len(terra::ncell(changed_grid))
  changed_grid_path <- file.path(directory, "bio1-grid-change.tif")
  terra::writeRaster(changed_grid, changed_grid_path, overwrite = TRUE)
  changed_grid_manifest <- file.path(directory, "grid-change.csv")
  utils::write.csv(
    data.frame(period = "current", scenario = "current", var_name = "bio1",
               file_path = "bio1-grid-change.tif"),
    changed_grid_manifest,
    row.names = FALSE
  )
  single_signature <- build_predictor_signature(stack, "bio1")
  expect_error(
    resolve_model_current_predictors(
      file.path(directory, "old-machine.csv"),
      changed_grid_manifest,
      "test_manifest",
      "climate",
      selected_predictors = "bio1",
      saved_crs = terra::crs(stack),
      saved_signature = single_signature
    ),
    "signature mismatch"
  )

  changed_crs <- terra::rast(nrows = 3, ncols = 4, xmin = 0, xmax = 4,
                             ymin = 0, ymax = 3, crs = "EPSG:4326")
  terra::values(changed_crs) <- seq_len(terra::ncell(changed_crs))
  changed_crs_path <- file.path(directory, "bio1-crs-change.tif")
  terra::writeRaster(changed_crs, changed_crs_path, overwrite = TRUE)
  changed_crs_manifest <- file.path(directory, "crs-change.csv")
  utils::write.csv(
    data.frame(period = "current", scenario = "current", var_name = "bio1",
               file_path = "bio1-crs-change.tif"),
    changed_crs_manifest,
    row.names = FALSE
  )
  expect_error(
    resolve_model_current_predictors(
      file.path(directory, "old-machine.csv"),
      changed_crs_manifest,
      "test_manifest",
      "climate",
      selected_predictors = "bio1",
      saved_crs = terra::crs(stack),
      saved_signature = single_signature
    ),
    "signature mismatch"
  )
})

if (length(failures) > 0L) {
  stop(
    length(failures), " of ", tests_run, " CV fortification tests failed:\n- ",
    paste(failures, collapse = "\n- "),
    call. = FALSE
  )
}

message("All ", tests_run, " CV fortification tests passed.")
