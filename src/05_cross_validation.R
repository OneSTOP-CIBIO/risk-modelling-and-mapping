#-------------------------------------------------------------------------------
# Spatial cross-validation with a grouped stratified non-spatial fallback
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# Load packages
#-------------------------------------------------------------------------------
packages <- c(
  "dplyr", "qs", "digest", "terra", "tidyterra", "sf", "here", "matrixStats",
  "ggplot2", "dismo", "sdm", "purrr", "ecospat", "blockCV"
)

installed <- rownames(installed.packages())
for (package in packages) {
  print(package)
  if (!package %in% installed) install.packages(package)
  library(package, character.only = TRUE)
}

suppressWarnings(try(sdm::installAll(), silent = TRUE))


#-------------------------------------------------------------------------------
# Runtime and project configuration
#-------------------------------------------------------------------------------
options("rgdal_show_exportToProj4_warnings" = "none")
terra::setGDALconfig("GDAL_PAM_ENABLED", "FALSE")
terra::terraOptions(
  memfrac = 0.4,
  tempdir = file.path(tempdir()),
  todisk = TRUE
)

source(here::here("src", "helper_functions.R"))
source(here::here("src", "00_configurations.R"))

if (length(cv_folds) != 1L || is.na(cv_folds) ||
    cv_folds != as.integer(cv_folds) || cv_folds < 2L) {
  stop("cv_folds must be a single integer of at least 2.", call. = FALSE)
}
if (!is.numeric(cv_spatial_block_sizes_m) ||
    length(cv_spatial_block_sizes_m) == 0L ||
    any(!is.finite(cv_spatial_block_sizes_m)) ||
    any(cv_spatial_block_sizes_m <= 0) ||
    any(diff(cv_spatial_block_sizes_m) > 0)) {
  stop(
    "cv_spatial_block_sizes_m must contain positive metre values ordered from largest to smallest.",
    call. = FALSE
  )
}
if (!is.logical(enable_kfold_cv_fallback) ||
    length(enable_kfold_cv_fallback) != 1L ||
    is.na(enable_kfold_cv_fallback)) {
  stop("enable_kfold_cv_fallback must be TRUE or FALSE.", call. = FALSE)
}

minimum_test_class_records <- 2L
minimum_successful_methods <- 3L
cv_seed <- 123L
cv_iterations <- 200L


#-------------------------------------------------------------------------------
# Input and output paths
#-------------------------------------------------------------------------------
climate_path <- file.path(
  "data", "external", "climate", "chelsa_current", "processed",
  "globalclimpreds.tif"
)
habitat_path <- file.path(
  "data", "external", "habitat", "processed", "habitat_stack.tif"
)
biome_path <- file.path(
  "data", "external", "GIS", "official", "newRealms.shp"
)

taxa_info <- read.csv2(file.path(
  "data", "projects", project, paste0(project, "_taxa_info.csv")
))
accepted_taxonkeys <- unique(taxa_info$acceptedTaxonKey)

validation_dir <- file.path("data", "projects", project, "Model_validation")
dir.create(validation_dir, recursive = TRUE, showWarnings = FALSE)
validate_output_directory(validation_dir)
validation_summary <- list()
validation_diagnostics <- list()


#-------------------------------------------------------------------------------
# Stage-local helpers
#-------------------------------------------------------------------------------
sample_available_cells <- function(mask, requested_size, seed) {
  available <- terra::global(!is.na(mask), "sum", na.rm = TRUE)[[1]]
  if (!is.finite(available) || available < 1L) {
    stop("No non-NA cells are available for validation background sampling.",
         call. = FALSE)
  }
  sample_size <- as.integer(min(as.numeric(requested_size), available))
  set.seed(seed)
  terra::spatSample(
    mask,
    size = sample_size,
    method = "random",
    na.rm = TRUE,
    as.points = TRUE
  )
}


extract_background <- function(points, predictor_stack, ids) {
  point_vector <- as_spatvector_safe(points)
  if (!isTRUE(terra::same.crs(point_vector, predictor_stack))) {
    point_vector <- terra::project(point_vector, predictor_stack)
  }
  values <- terra::extract(
    predictor_stack,
    point_vector,
    ID = FALSE,
    xy = FALSE
  )
  if (nrow(values) != length(ids) || ncol(values) == 0L) {
    stop(
      "Background extraction returned ", nrow(values), " row(s) and ",
      ncol(values), " predictor column(s) for ", length(ids), " point(s).",
      call. = FALSE
    )
  }
  keep <- stats::complete.cases(values)
  numeric_columns <- vapply(values, is.numeric, logical(1))
  if (any(numeric_columns) && nrow(values) > 0L) {
    keep <- keep & apply(
      values[, numeric_columns, drop = FALSE],
      1,
      function(x) all(is.finite(x))
    )
  }
  values <- values[keep, , drop = FALSE]
  values$ID <- ids[keep]
  values
}


filter_records_for_stack <- function(records, predictor_stack) {
  extracted <- extract_env(records, predictor_stack)
  usable_ids <- extracted$complete$ID
  list(
    records = records[records$ID %in% usable_ids, , drop = FALSE],
    extracted = extracted,
    dropped_n = nrow(records) - length(usable_ids)
  )
}


make_spatial_group <- function(records, reference_raster) {
  point_vector <- as_spatvector_safe(records)
  if (!isTRUE(terra::same.crs(point_vector, reference_raster))) {
    point_vector <- terra::project(point_vector, reference_raster)
  }
  coordinates <- terra::crds(point_vector)
  cells <- terra::cellFromXY(reference_raster[[1]], coordinates)
  groups <- ifelse(
    is.na(cells),
    paste0(
      "xy_",
      format(round(coordinates[, 1], 6), scientific = FALSE, trim = TRUE),
      "_",
      format(round(coordinates[, 2], 6), scientific = FALSE, trim = TRUE)
    ),
    paste0("cell_", cells)
  )
  as.character(groups)
}


make_plan_context <- function(records, context, reference_raster) {
  output <- records[, c("species", "ID"), drop = FALSE]
  output$source_ID <- as.character(records$ID)
  output$cv_context <- context
  output$cv_stratum <- paste0(context, "_", records$species)
  output$cv_group <- make_spatial_group(records, reference_raster)
  output
}


predictor_filter_diagnostics <- function(before_records,
                                         after_records,
                                         species_name,
                                         context,
                                         component) {
  classes <- sort(unique(c(
    "0",
    "1",
    as.character(before_records$species),
    as.character(after_records$species)
  )))
  dplyr::bind_rows(lapply(classes, function(class_value) {
    before_ids <- as.character(
      before_records$ID[as.character(before_records$species) == class_value]
    )
    after_ids <- as.character(
      after_records$ID[as.character(after_records$species) == class_value]
    )
    dropped_ids <- setdiff(before_ids, after_ids)
    data.frame(
      diagnostic_type = "predictor_filter",
      cv_method = NA_character_,
      requested_folds = as.integer(cv_folds),
      effective_folds = NA_integer_,
      block_size_m = NA_real_,
      valid = length(after_ids) >= 2L,
      selected = NA,
      min_train_count = NA_integer_,
      min_test_count = NA_integer_,
      unassigned_records = NA_integer_,
      reason = if (length(dropped_ids) == 0L) {
        "no records removed by predictor filtering"
      } else {
        paste0(length(dropped_ids), " record(s) removed by predictor filtering")
      },
      fallback_reason = NA_character_,
      Species = species_name,
      validation_context = context,
      model_component = component,
      class_value = class_value,
      records_before = length(before_ids),
      records_after = length(after_ids),
      dropped_records = length(dropped_ids),
      dropped_record_ids = paste(dropped_ids, collapse = ";"),
      test_fold = NA_integer_,
      method = NA_character_,
      dataset = NA_character_,
      success = length(after_ids) >= 2L,
      error = if (length(after_ids) < 2L) {
        "fewer than two usable records remain in this class"
      } else {
        NA_character_
      },
      stringsAsFactors = FALSE
    )
  }))
}


plan_metadata <- function(metrics, plan) {
  metrics$cv_method <- plan$method
  metrics$requested_folds <- as.integer(plan$requested_folds)
  metrics$effective_folds <- if (is.null(plan$k)) NA_integer_ else as.integer(plan$k)
  metrics$block_size_m <- if (is.null(plan$block_size_m)) {
    NA_real_
  } else {
    as.numeric(plan$block_size_m)
  }
  metrics$fallback_reason <- if (is.null(plan$fallback_reason)) {
    NA_character_
  } else {
    as.character(plan$fallback_reason)
  }
  metrics
}


not_evaluable_plan <- function(reason) {
  list(
    valid = FALSE,
    method = "not_evaluable",
    requested_folds = as.integer(cv_folds),
    k = NA_integer_,
    block_size_m = NA_real_,
    fold_ids = integer(0),
    fallback_reason = reason,
    diagnostics = data.frame()
  )
}


not_evaluable_metrics <- function(species_name, type, region, plan) {
  plan_metadata(
    compute_validation_metrics(
      species = species_name,
      type = type,
      region = region,
      fold = NA_character_,
      all_suit_vals = numeric(0),
      occ_suit_vals = numeric(0),
      abs_suit_vals = numeric(0)
    ),
    plan
  )
}


decorate_partition_diagnostics <- function(plan, species_name, context) {
  diagnostics <- plan$diagnostics
  if (is.null(diagnostics) || nrow(diagnostics) == 0L) {
    return(data.frame())
  }
  diagnostics$Species <- species_name
  diagnostics$validation_context <- context
  diagnostics$model_component <- NA_character_
  diagnostics$requested_folds <- as.integer(plan$requested_folds)
  diagnostics$fallback_reason <- if (is.null(plan$fallback_reason)) {
    NA_character_
  } else {
    as.character(plan$fallback_reason)
  }
  diagnostics$selected <- diagnostics$valid &
    diagnostics$cv_method == plan$method &
    diagnostics$effective_folds == plan$k &
    (
      is.na(diagnostics$block_size_m) & is.na(plan$block_size_m) |
        diagnostics$block_size_m == plan$block_size_m
    )
  diagnostics$test_fold <- NA_integer_
  diagnostics$method <- NA_character_
  diagnostics$dataset <- NA_character_
  diagnostics$success <- diagnostics$valid
  diagnostics$error <- ifelse(diagnostics$valid, NA_character_, diagnostics$reason)
  diagnostics
}


decorate_prediction_diagnostics <- function(prediction,
                                            species_name,
                                            context,
                                            component,
                                            fold,
                                            plan) {
  diagnostics <- prediction$method_diagnostics
  if (is.null(diagnostics) || nrow(diagnostics) == 0L) {
    return(data.frame())
  }
  diagnostics$diagnostic_type <- "prediction"
  diagnostics$Species <- species_name
  diagnostics$validation_context <- context
  diagnostics$model_component <- component
  diagnostics$test_fold <- as.integer(fold)
  diagnostics$cv_method <- plan$method
  diagnostics$effective_folds <- as.integer(plan$k)
  diagnostics$block_size_m <- as.numeric(plan$block_size_m)
  diagnostics$requested_folds <- as.integer(plan$requested_folds)
  diagnostics$fallback_reason <- if (is.null(plan$fallback_reason)) {
    NA_character_
  } else {
    as.character(plan$fallback_reason)
  }
  diagnostics$selected <- TRUE
  diagnostics$valid <- diagnostics$success
  diagnostics$min_train_count <- NA_integer_
  diagnostics$min_test_count <- NA_integer_
  diagnostics$unassigned_records <- NA_integer_
  diagnostics$reason <- diagnostics$error
  diagnostics
}


runtime_diagnostic <- function(species_name,
                               context,
                               component,
                               fold,
                               plan,
                               reason) {
  data.frame(
    diagnostic_type = "runtime_failure",
    cv_method = plan$method,
    requested_folds = as.integer(plan$requested_folds),
    effective_folds = as.integer(plan$k),
    block_size_m = as.numeric(plan$block_size_m),
    valid = FALSE,
    selected = TRUE,
    min_train_count = NA_integer_,
    min_test_count = NA_integer_,
    unassigned_records = NA_integer_,
    reason = reason,
    fallback_reason = if (is.null(plan$fallback_reason)) {
      NA_character_
    } else {
      as.character(plan$fallback_reason)
    },
    Species = species_name,
    validation_context = context,
    model_component = component,
    test_fold = as.integer(fold),
    method = NA_character_,
    dataset = NA_character_,
    success = FALSE,
    error = reason,
    stringsAsFactors = FALSE
  )
}


fit_sdm_fold <- function(train_points, predictor_stack, methods) {
  tryCatch({
    train_points <- sf::st_transform(
      train_points[, "species", drop = FALSE],
      sf::st_crs(terra::crs(predictor_stack))
    )
    sdm_data <- sdm::sdmData(
      species ~ .,
      train = as_spatvector_safe(train_points),
      predictors = predictor_stack
    )
    sdm::sdm(species ~ ., data = sdm_data, methods = methods)
  }, error = function(e) e)
}


context_fold_map <- function(plan_records, plan, context) {
  keep <- plan_records$cv_context == context
  source_ids <- plan_records$source_ID[keep]
  folds <- plan$fold_ids[keep]
  if (length(source_ids) == 0L || anyDuplicated(source_ids) || anyNA(folds)) {
    stop("The selected CV context does not have a complete unique fold map.",
         call. = FALSE)
  }
  stats::setNames(as.integer(folds), source_ids)
}


run_climate_plan <- function(plan,
                             plan_records,
                             global_records,
                             eu_record_ids,
                             climate_selection,
                             methods,
                             global_background,
                             eu_background,
                             species_name) {
  diagnostics <- list()
  diagnostic_index <- 0L
  global_metrics <- list()
  eu_metrics <- list()

  if (!isTRUE(plan$valid)) {
    return(list(
      valid = FALSE,
      failure_reason = plan$fallback_reason,
      global_metrics = data.frame(),
      eu_metrics = data.frame(),
      diagnostics = data.frame()
    ))
  }

  global_map <- tryCatch(
    context_fold_map(plan_records, plan, "global"),
    error = function(e) e
  )
  if (inherits(global_map, "error")) {
    return(list(
      valid = FALSE,
      failure_reason = conditionMessage(global_map),
      global_metrics = data.frame(),
      eu_metrics = data.frame(),
      diagnostics = data.frame()
    ))
  }
  global_records$folds <- unname(global_map[as.character(global_records$ID)])
  include_eu <- length(eu_record_ids) > 0L

  if (include_eu) {
    eu_map <- tryCatch(
      context_fold_map(plan_records, plan, "eu_climate"),
      error = function(e) e
    )
    if (inherits(eu_map, "error") ||
        any(unname(eu_map[eu_record_ids]) != unname(global_map[eu_record_ids]))) {
      reason <- if (inherits(eu_map, "error")) {
        conditionMessage(eu_map)
      } else {
        "European climate records received inconsistent folds across contexts"
      }
      return(list(
        valid = FALSE,
        failure_reason = reason,
        global_metrics = data.frame(),
        eu_metrics = data.frame(),
        diagnostics = data.frame()
      ))
    }
  }

  for (fold in seq_len(plan$k)) {
    message(
      "Climate validation fold ", fold, "/", plan$k,
      " [", plan$method, "]"
    )
    train_data <- global_records[global_records$folds != fold, , drop = FALSE]
    test_data <- global_records[global_records$folds == fold, , drop = FALSE]
    pres_train <- sum(train_data$species == 1)
    abs_train <- sum(train_data$species == 0)
    if (pres_train < minimum_test_class_records ||
        abs_train < minimum_test_class_records) {
      reason <- "climate training split does not retain both classes"
      diagnostic_index <- diagnostic_index + 1L
      diagnostics[[diagnostic_index]] <- runtime_diagnostic(
        species_name, "climate", "climate", fold, plan, reason
      )
      return(list(valid = FALSE, failure_reason = reason,
                  global_metrics = data.frame(), eu_metrics = data.frame(),
                  diagnostics = dplyr::bind_rows(diagnostics)))
    }

    model <- fit_sdm_fold(train_data, climate_selection, methods)
    if (inherits(model, "error")) {
      reason <- paste0("climate model fitting failed: ", conditionMessage(model))
      diagnostic_index <- diagnostic_index + 1L
      diagnostics[[diagnostic_index]] <- runtime_diagnostic(
        species_name, "climate", "climate", fold, plan, reason
      )
      return(list(valid = FALSE, failure_reason = reason,
                  global_metrics = data.frame(), eu_metrics = data.frame(),
                  diagnostics = dplyr::bind_rows(diagnostics)))
    }

    test_environment <- extract_env(test_data, climate_selection)
    datasets <- list(
      global_points = global_background,
      occ_env = test_environment$presences,
      abs_env = test_environment$absences
    )
    if (include_eu) {
      eu_test <- test_data[test_data$ID %in% eu_record_ids, , drop = FALSE]
      eu_environment <- extract_env(eu_test, climate_selection)
      datasets$eu_points <- eu_background
      datasets$eu_occ_env <- eu_environment$presences
      datasets$eu_abs_env <- eu_environment$absences
    }

    prediction <- compute_median_favourability_safe(
      model = model,
      datasets = datasets,
      top5_methods = methods,
      prev_ratio = pres_train / abs_train,
      min_successful_methods = minimum_successful_methods
    )
    diagnostic_index <- diagnostic_index + 1L
    diagnostics[[diagnostic_index]] <- decorate_prediction_diagnostics(
      prediction, species_name, "climate", "climate", fold, plan
    )
    if (!isTRUE(prediction$valid)) {
      reason <- paste0(
        "only ", prediction$success_count, " of ", length(methods),
        " climate methods predicted successfully"
      )
      diagnostic_index <- diagnostic_index + 1L
      diagnostics[[diagnostic_index]] <- runtime_diagnostic(
        species_name, "climate", "climate", fold, plan, reason
      )
      return(list(valid = FALSE, failure_reason = reason,
                  global_metrics = data.frame(), eu_metrics = data.frame(),
                  diagnostics = dplyr::bind_rows(diagnostics)))
    }

    favourability <- prediction$median_favourability
    global_metrics[[fold]] <- compute_validation_metrics(
      species = species_name,
      type = "Climate",
      region = "Global",
      fold = fold,
      all_suit_vals = favourability$global_points$median_favourability,
      occ_suit_vals = favourability$occ_env$median_favourability,
      abs_suit_vals = favourability$abs_env$median_favourability
    )
    if (global_metrics[[fold]]$n_pres < minimum_test_class_records ||
        global_metrics[[fold]]$n_abs < minimum_test_class_records) {
      reason <- "climate test predictions do not retain enough records per class"
      return(list(valid = FALSE, failure_reason = reason,
                  global_metrics = data.frame(), eu_metrics = data.frame(),
                  diagnostics = dplyr::bind_rows(diagnostics)))
    }

    if (include_eu) {
      eu_metrics[[fold]] <- compute_validation_metrics(
        species = species_name,
        type = "Climate",
        region = "Europe",
        fold = fold,
        all_suit_vals = favourability$eu_points$median_favourability,
        occ_suit_vals = favourability$eu_occ_env$median_favourability,
        abs_suit_vals = favourability$eu_abs_env$median_favourability
      )
      if (eu_metrics[[fold]]$n_pres < minimum_test_class_records ||
          eu_metrics[[fold]]$n_abs < minimum_test_class_records) {
        reason <- "European climate test predictions do not retain enough records per class"
        return(list(valid = FALSE, failure_reason = reason,
                    global_metrics = data.frame(), eu_metrics = data.frame(),
                    diagnostics = dplyr::bind_rows(diagnostics)))
      }
    }
    terra::tmpFiles(remove = TRUE)
  }

  list(
    valid = TRUE,
    failure_reason = NA_character_,
    global_metrics = dplyr::bind_rows(global_metrics),
    eu_metrics = dplyr::bind_rows(eu_metrics),
    diagnostics = dplyr::bind_rows(diagnostics)
  )
}


run_european_plan <- function(plan,
                              plan_records,
                              global_records,
                              eu_records,
                              climate_selection,
                              habitat_selection,
                              climate_methods,
                              habitat_methods,
                              climate_background,
                              habitat_background,
                              species_name) {
  diagnostics <- list()
  diagnostic_index <- 0L
  habitat_metrics <- list()
  ensemble_metrics <- list()

  if (!isTRUE(plan$valid)) {
    return(list(
      valid = FALSE,
      failure_reason = plan$fallback_reason,
      habitat_metrics = data.frame(),
      ensemble_metrics = data.frame(),
      diagnostics = data.frame()
    ))
  }

  global_map <- tryCatch(context_fold_map(plan_records, plan, "global"),
                         error = function(e) e)
  eu_map <- tryCatch(context_fold_map(plan_records, plan, "eu_component"),
                     error = function(e) e)
  if (inherits(global_map, "error") || inherits(eu_map, "error")) {
    reason <- paste(
      if (inherits(global_map, "error")) conditionMessage(global_map) else NULL,
      if (inherits(eu_map, "error")) conditionMessage(eu_map) else NULL,
      collapse = " | "
    )
    return(list(valid = FALSE, failure_reason = reason,
                habitat_metrics = data.frame(), ensemble_metrics = data.frame(),
                diagnostics = data.frame()))
  }
  global_records$folds <- unname(global_map[as.character(global_records$ID)])
  eu_records$folds <- unname(eu_map[as.character(eu_records$ID)])

  for (fold in seq_len(plan$k)) {
    message(
      "European habitat/combined validation fold ", fold, "/", plan$k,
      " [", plan$method, "]"
    )
    global_train <- global_records[global_records$folds != fold, , drop = FALSE]
    eu_train <- eu_records[eu_records$folds != fold, , drop = FALSE]
    eu_test <- eu_records[eu_records$folds == fold, , drop = FALSE]
    global_pres <- sum(global_train$species == 1)
    global_abs <- sum(global_train$species == 0)
    eu_pres <- sum(eu_train$species == 1)
    eu_abs <- sum(eu_train$species == 0)
    if (min(global_pres, global_abs, eu_pres, eu_abs) < minimum_test_class_records) {
      reason <- "European component training split does not retain both classes"
      diagnostic_index <- diagnostic_index + 1L
      diagnostics[[diagnostic_index]] <- runtime_diagnostic(
        species_name, "european_component", "both", fold, plan, reason
      )
      return(list(valid = FALSE, failure_reason = reason,
                  habitat_metrics = data.frame(), ensemble_metrics = data.frame(),
                  diagnostics = dplyr::bind_rows(diagnostics)))
    }

    climate_model <- fit_sdm_fold(global_train, climate_selection, climate_methods)
    habitat_model <- fit_sdm_fold(eu_train, habitat_selection, habitat_methods)
    if (inherits(climate_model, "error") || inherits(habitat_model, "error")) {
      reason <- paste(
        if (inherits(climate_model, "error")) {
          paste0("climate fitting failed: ", conditionMessage(climate_model))
        } else NULL,
        if (inherits(habitat_model, "error")) {
          paste0("habitat fitting failed: ", conditionMessage(habitat_model))
        } else NULL,
        collapse = " | "
      )
      diagnostic_index <- diagnostic_index + 1L
      diagnostics[[diagnostic_index]] <- runtime_diagnostic(
        species_name, "european_component", "both", fold, plan, reason
      )
      return(list(valid = FALSE, failure_reason = reason,
                  habitat_metrics = data.frame(), ensemble_metrics = data.frame(),
                  diagnostics = dplyr::bind_rows(diagnostics)))
    }

    climate_environment <- extract_env(eu_test, climate_selection)
    habitat_environment <- extract_env(eu_test, habitat_selection)
    climate_prediction <- compute_median_favourability_safe(
      model = climate_model,
      datasets = list(
        eu_points = climate_background,
        ens_occ_env = climate_environment$presences,
        ens_abs_env = climate_environment$absences
      ),
      top5_methods = climate_methods,
      prev_ratio = global_pres / global_abs,
      min_successful_methods = minimum_successful_methods
    )
    habitat_prediction <- compute_median_favourability_safe(
      model = habitat_model,
      datasets = list(
        eu_habitat_points = habitat_background,
        occ_hab = habitat_environment$presences,
        abs_hab = habitat_environment$absences
      ),
      top5_methods = habitat_methods,
      prev_ratio = eu_pres / eu_abs,
      min_successful_methods = minimum_successful_methods
    )
    diagnostic_index <- diagnostic_index + 1L
    diagnostics[[diagnostic_index]] <- decorate_prediction_diagnostics(
      climate_prediction, species_name, "european_component", "climate",
      fold, plan
    )
    diagnostic_index <- diagnostic_index + 1L
    diagnostics[[diagnostic_index]] <- decorate_prediction_diagnostics(
      habitat_prediction, species_name, "european_component", "habitat",
      fold, plan
    )

    if (!isTRUE(climate_prediction$valid) || !isTRUE(habitat_prediction$valid)) {
      reason <- paste0(
        "insufficient successful component methods (climate=",
        climate_prediction$success_count, ", habitat=",
        habitat_prediction$success_count, ")"
      )
      diagnostic_index <- diagnostic_index + 1L
      diagnostics[[diagnostic_index]] <- runtime_diagnostic(
        species_name, "european_component", "both", fold, plan, reason
      )
      return(list(valid = FALSE, failure_reason = reason,
                  habitat_metrics = data.frame(), ensemble_metrics = data.frame(),
                  diagnostics = dplyr::bind_rows(diagnostics)))
    }

    climate_fav <- climate_prediction$median_favourability
    habitat_fav <- habitat_prediction$median_favourability
    habitat_metrics[[fold]] <- compute_validation_metrics(
      species = species_name,
      type = "Habitat",
      region = "Europe",
      fold = fold,
      all_suit_vals = habitat_fav$eu_habitat_points$median_favourability,
      occ_suit_vals = habitat_fav$occ_hab$median_favourability,
      abs_suit_vals = habitat_fav$abs_hab$median_favourability
    )

    ensemble_background <- ensemble_geom_mean(
      habitat_fav$eu_habitat_points,
      climate_fav$eu_points,
      type = "background",
      return_data = TRUE
    )
    ensemble_occurrences <- ensemble_geom_mean(
      habitat_fav$occ_hab,
      climate_fav$ens_occ_env,
      type = "occurrence",
      return_data = TRUE
    )
    ensemble_absences <- ensemble_geom_mean(
      habitat_fav$abs_hab,
      climate_fav$ens_abs_env,
      type = "absence",
      return_data = TRUE
    )
    ensemble_metrics[[fold]] <- compute_validation_metrics(
      species = species_name,
      type = "Ensemble",
      region = "Europe",
      fold = fold,
      all_suit_vals = ensemble_background$ensemble_favourability,
      occ_suit_vals = ensemble_occurrences$ensemble_favourability,
      abs_suit_vals = ensemble_absences$ensemble_favourability
    )

    if (habitat_metrics[[fold]]$n_pres < minimum_test_class_records ||
        habitat_metrics[[fold]]$n_abs < minimum_test_class_records ||
        ensemble_metrics[[fold]]$n_pres < minimum_test_class_records ||
        ensemble_metrics[[fold]]$n_abs < minimum_test_class_records) {
      reason <- "component predictions do not retain enough matched records per class"
      return(list(valid = FALSE, failure_reason = reason,
                  habitat_metrics = data.frame(), ensemble_metrics = data.frame(),
                  diagnostics = dplyr::bind_rows(diagnostics)))
    }
    terra::tmpFiles(remove = TRUE)
  }

  list(
    valid = TRUE,
    failure_reason = NA_character_,
    habitat_metrics = dplyr::bind_rows(habitat_metrics),
    ensemble_metrics = dplyr::bind_rows(ensemble_metrics),
    diagnostics = dplyr::bind_rows(diagnostics)
  )
}


run_with_runtime_fallback <- function(initial_plan,
                                      plan_records,
                                      runner,
                                      species_name,
                                      context) {
  diagnostics <- decorate_partition_diagnostics(
    initial_plan, species_name, context
  )
  if (!isTRUE(initial_plan$valid)) {
    return(list(
      plan = initial_plan,
      result = list(valid = FALSE, failure_reason = initial_plan$fallback_reason),
      diagnostics = diagnostics
    ))
  }

  result <- runner(initial_plan)
  diagnostics <- dplyr::bind_rows(diagnostics, result$diagnostics)
  if (isTRUE(result$valid) ||
      !identical(initial_plan$method, "spatial_block") ||
      !isTRUE(enable_kfold_cv_fallback)) {
    return(list(plan = initial_plan, result = result, diagnostics = diagnostics))
  }

  fallback_reason <- paste0(
    "Spatial CV passed structural checks but failed during model evaluation: ",
    result$failure_reason
  )
  message(fallback_reason, ". Retrying with grouped stratified k-fold CV.")
  kfold_plan <- make_stratified_kfold_plan(
    records = plan_records,
    requested_folds = cv_folds,
    stratum_col = "cv_stratum",
    group_col = "cv_group",
    min_train = minimum_test_class_records,
    min_test = minimum_test_class_records,
    seed = cv_seed,
    iteration = cv_iterations,
    fallback_reason = fallback_reason
  )
  diagnostics <- dplyr::bind_rows(
    diagnostics,
    decorate_partition_diagnostics(kfold_plan, species_name, context)
  )
  if (!isTRUE(kfold_plan$valid)) {
    return(list(
      plan = kfold_plan,
      result = list(valid = FALSE, failure_reason = fallback_reason),
      diagnostics = diagnostics
    ))
  }

  fallback_result <- runner(kfold_plan)
  diagnostics <- dplyr::bind_rows(diagnostics, fallback_result$diagnostics)
  list(plan = kfold_plan, result = fallback_result, diagnostics = diagnostics)
}


write_metric_outputs <- function(metrics,
                                 output_dir,
                                 per_fold_filename,
                                 summary_filename) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  validate_output_directory(output_dir)
  readr::write_csv(
    metrics,
    fortify_output_path(file.path(output_dir, per_fold_filename))
  )
  summary <- summarise_validation(metrics)
  readr::write_csv(
    summary,
    fortify_output_path(file.path(output_dir, summary_filename))
  )
  summary
}


#-------------------------------------------------------------------------------
# Species loop
#-------------------------------------------------------------------------------
for (i in seq_along(accepted_taxonkeys)) {
  species_validation_summary <- data.frame()
  species_diagnostics <- list()
  diagnostic_index <- 0L

  species <- taxa_info$acceptedScientificName[i]
  taxonkey <- taxa_info$acceptedTaxonKey[i]
  speciesName <- species_output_stem(species)
  message(
    "\n", strrep("=", 72),
    "\nSPECIES: ", speciesName, "  [taxonkey: ", taxonkey, "]",
    "\n", strrep("=", 72)
  )

  base_dir <- file.path(
    "data", "projects", project, paste0(speciesName, "_", taxonkey)
  )
  climate_qs_file <- file.path(
    base_dir, "Climate",
    paste0("Climate_model_", speciesName, "_", taxonkey, ".qs")
  )
  habitat_qs_file <- file.path(
    base_dir, "Habitat",
    paste0("Habitat_model_", speciesName, "_", taxonkey, ".qs")
  )
  if (!file.exists(climate_qs_file)) {
    warning("No climate model was found for ", species, "; skipping validation.")
    next
  }

  climate_model_artifact <- qs::qread(climate_qs_file)
  climate_methods <- climate_model_artifact$top5_models
  climate_predictors <- climate_model_artifact$selected_predictors
  global_presabs <- climate_model_artifact$global_presabs
  climate_input_mode <- climate_model_artifact$climate_input_mode
  if (is.null(climate_input_mode) || !nzchar(climate_input_mode)) {
    climate_input_mode <- "chelsa"
  }
  if (!climate_input_mode %in% c("chelsa", "user_specific")) {
    stop("Unsupported saved climate input mode: ", climate_input_mode,
         call. = FALSE)
  }

  if (identical(climate_input_mode, "user_specific")) {
    climate_resolution <- resolve_model_current_predictors(
      saved_manifest_path = climate_model_artifact$climate_manifest_path,
      configured_manifest_path = user_specific_climate_data,
      config_name = "user_specific_climate_data",
      data_label = "climate",
      selected_predictors = climate_predictors,
      saved_crs = climate_model_artifact$predictor_crs,
      saved_signature = climate_model_artifact$current_predictor_signature,
      materialize = function(manifest) {
        materialize_user_specific_current_stack(manifest)
      }
    )
    climate_stack_active <- climate_resolution$predictor_stack
  } else {
    climate_stack_active <- terra::rast(climate_path)
    compatibility <- validate_model_predictor_stack(
      climate_stack_active,
      climate_predictors,
      saved_crs = climate_model_artifact$predictor_crs,
      saved_signature = climate_model_artifact$current_predictor_signature
    )
    if (!isTRUE(compatibility$valid)) {
      stop("Default climate predictors are incompatible with the saved model: ",
           compatibility$reason, call. = FALSE)
    }
  }
  climate_selection <- terra::subset(climate_stack_active, climate_predictors)
  if (terra::nlyr(climate_selection) != length(climate_predictors) ||
      !identical(names(climate_selection), climate_predictors)) {
    stop("The active climate stack does not contain the saved predictors in order.",
         call. = FALSE)
  }

  has_habitat <- file.exists(habitat_qs_file)
  habitat_model_artifact <- NULL
  habitat_selection <- NULL
  habitat_methods <- character(0)
  if (has_habitat) {
    habitat_model_artifact <- qs::qread(habitat_qs_file)
    habitat_predictors <- habitat_model_artifact$selected_predictors
    habitat_methods <- habitat_model_artifact$top5_models
    landcover_input_mode <- habitat_model_artifact$landcover_input_mode
    if (is.null(landcover_input_mode) || !nzchar(landcover_input_mode)) {
      landcover_input_mode <- "default"
    }
    if (!landcover_input_mode %in% c("default", "user_specific")) {
      stop("Unsupported saved land-cover input mode: ", landcover_input_mode,
           call. = FALSE)
    }

    if (identical(landcover_input_mode, "user_specific")) {
      habitat_resolution <- resolve_model_current_predictors(
        saved_manifest_path = habitat_model_artifact$landcover_manifest_path,
        configured_manifest_path = user_specific_landcover_data,
        config_name = "user_specific_landcover_data",
        data_label = "land-cover",
        selected_predictors = habitat_predictors,
        saved_crs = habitat_model_artifact$landcover_predictor_crs,
        saved_signature = habitat_model_artifact$current_predictor_signature,
        materialize = function(manifest) {
          materialize_user_specific_landcover_stack(
            stack_rows = manifest$current_rows,
            period = "current",
            scenario = "current"
          )
        }
      )
      habitat_stack_active <- habitat_resolution$predictor_stack
    } else {
      habitat_stack_active <- terra::rast(habitat_path)
      compatibility <- validate_model_predictor_stack(
        habitat_stack_active,
        habitat_predictors,
        saved_crs = habitat_model_artifact$landcover_predictor_crs,
        saved_signature = habitat_model_artifact$current_predictor_signature
      )
      if (!isTRUE(compatibility$valid)) {
        stop("Default habitat predictors are incompatible with the saved model: ",
             compatibility$reason, call. = FALSE)
      }
    }
    habitat_selection <- terra::subset(habitat_stack_active, habitat_predictors)
    if (terra::nlyr(habitat_selection) != length(habitat_predictors) ||
        !identical(names(habitat_selection), habitat_predictors)) {
      stop("The active habitat stack does not contain the saved predictors in order.",
           call. = FALSE)
    }
  }

  # Stable IDs and predictor-complete model records.
  global_presabs$ID <- paste0("GLOBAL_", seq_len(nrow(global_presabs)))
  global_filtered <- filter_records_for_stack(global_presabs, climate_selection)
  global_records <- global_filtered$records[, c("species", "ID"), drop = FALSE]
  diagnostic_index <- diagnostic_index + 1L
  species_diagnostics[[diagnostic_index]] <- predictor_filter_diagnostics(
    global_presabs,
    global_records,
    speciesName,
    "global",
    "climate"
  )
  if (global_filtered$dropped_n > 0L) {
    warning(global_filtered$dropped_n,
            " global record(s) were removed because climate predictors were incomplete.")
  }
  global_class_counts <- table(
    factor(global_records$species, levels = c(0, 1))
  )
  if (any(global_class_counts < 2L)) {
    stop("Fewer than two usable global records remain in one class.", call. = FALSE)
  }

  euboundary_climate <- load_eu_boundary(
    custom_path = custom_eu_boundary_path,
    reference = climate_selection[[1]]
  )
  boundary_for_global <- sf::st_transform(
    euboundary_climate,
    sf::st_crs(global_records)
  )
  is_european_global <- lengths(sf::st_intersects(global_records, boundary_for_global)) > 0L
  eu_climate_records <- global_records[is_european_global, , drop = FALSE]
  boundary_for_original_global <- sf::st_transform(
    euboundary_climate,
    sf::st_crs(global_presabs)
  )
  original_european_global <- global_presabs[
    lengths(sf::st_intersects(global_presabs, boundary_for_original_global)) > 0L,
    ,
    drop = FALSE
  ]
  diagnostic_index <- diagnostic_index + 1L
  species_diagnostics[[diagnostic_index]] <- predictor_filter_diagnostics(
    original_european_global,
    eu_climate_records,
    speciesName,
    "eu_climate",
    "climate"
  )
  eu_climate_counts <- table(factor(eu_climate_records$species, levels = c(0, 1)))
  validate_climate_europe <- all(
    eu_climate_counts >= 2L * minimum_test_class_records
  )

  # Global Boyce background restricted to occupied biomes.
  biome <- sf::st_read(biome_path, quiet = TRUE) %>%
    sf::st_transform(sf::st_crs(terra::crs(climate_selection)))
  global_presences <- global_records[global_records$species == 1, , drop = FALSE] %>%
    sf::st_transform(sf::st_crs(biome))
  sf::sf_use_s2(FALSE)
  occupied_biomes <- lengths(sf::st_intersects(biome, global_presences)) > 0L
  sf::sf_use_s2(TRUE)
  biome <- biome[occupied_biomes, , drop = FALSE]
  if (nrow(biome) == 0L) {
    stop("No biome intersects the usable global occurrence records.", call. = FALSE)
  }
  biome_vector <- as_spatvector_safe(biome)
  global_mask <- terra::crop(climate_selection[[1]], biome_vector) %>%
    terra::mask(biome_vector)
  global_sample <- sample_available_cells(global_mask, boyce_background_size, 728L)
  global_background <- extract_background(
    global_sample,
    climate_selection,
    paste0("GLOBAL_BG_", seq_len(nrow(global_sample)))
  )

  # European background points are sampled once so climate and habitat IDs align.
  climate_eu_mask <- terra::crop(
    climate_selection[[1]], as_spatvector_safe(euboundary_climate)
  ) %>% terra::mask(as_spatvector_safe(euboundary_climate))
  if (has_habitat) {
    euboundary_habitat <- load_eu_boundary(
      custom_path = custom_eu_boundary_path,
      reference = habitat_selection[[1]]
    )
    habitat_eu_mask <- terra::crop(
      habitat_selection[[1]], as_spatvector_safe(euboundary_habitat)
    ) %>% terra::mask(as_spatvector_safe(euboundary_habitat))
    climate_on_habitat <- terra::project(climate_eu_mask, habitat_eu_mask)
    eu_sampling_mask <- terra::mask(habitat_eu_mask, climate_on_habitat)
  } else {
    eu_sampling_mask <- climate_eu_mask
  }
  eu_sample <- sample_available_cells(eu_sampling_mask, boyce_background_size, 728L)
  eu_ids <- paste0("EU_BG_", seq_len(nrow(eu_sample)))
  eu_climate_background <- extract_background(
    eu_sample, climate_selection, eu_ids
  )
  eu_habitat_background <- NULL
  if (has_habitat) {
    eu_habitat_background <- extract_background(
      eu_sample, habitat_selection, eu_ids
    )
    common_background_ids <- intersect(
      eu_climate_background$ID,
      eu_habitat_background$ID
    )
    eu_component_climate_background <- eu_climate_background[
      eu_climate_background$ID %in% common_background_ids, , drop = FALSE
    ]
    eu_component_habitat_background <- eu_habitat_background[
      eu_habitat_background$ID %in% common_background_ids, , drop = FALSE
    ]
  }

  # Climate plan: always global, and European climate metrics when eligible.
  climate_plan_records <- make_plan_context(
    global_records, "global", climate_selection[[1]]
  )
  eu_climate_ids <- character(0)
  if (validate_climate_europe) {
    eu_climate_ids <- as.character(eu_climate_records$ID)
    climate_plan_records <- rbind(
      climate_plan_records,
      make_plan_context(
        eu_climate_records, "eu_climate", climate_selection[[1]]
      )
    )
  }
  climate_plan <- make_preferred_cv_plan(
    records = climate_plan_records,
    requested_folds = cv_folds,
    block_sizes_m = cv_spatial_block_sizes_m,
    enable_kfold_fallback = enable_kfold_cv_fallback,
    min_train = minimum_test_class_records,
    min_test = minimum_test_class_records,
    seed = cv_seed,
    iteration = cv_iterations
  )
  climate_run <- run_with_runtime_fallback(
    initial_plan = climate_plan,
    plan_records = climate_plan_records,
    runner = function(plan) {
      run_climate_plan(
        plan = plan,
        plan_records = climate_plan_records,
        global_records = global_records,
        eu_record_ids = eu_climate_ids,
        climate_selection = climate_selection,
        methods = climate_methods,
        global_background = global_background,
        eu_background = eu_climate_background,
        species_name = speciesName
      )
    },
    species_name = speciesName,
    context = "climate"
  )
  diagnostic_index <- diagnostic_index + 1L
  species_diagnostics[[diagnostic_index]] <- climate_run$diagnostics

  if (isTRUE(climate_run$result$valid)) {
    global_climate_metrics <- plan_metadata(
      climate_run$result$global_metrics, climate_run$plan
    )
    if (validate_climate_europe) {
      eu_climate_metrics <- plan_metadata(
        climate_run$result$eu_metrics, climate_run$plan
      )
    } else {
      reason <- paste0(
        "European climate validation requires enough records of both classes for at least two folds. ",
        "Usable counts: absence=", eu_climate_counts[["0"]],
        ", presence=", eu_climate_counts[["1"]], "."
      )
      eu_climate_metrics <- not_evaluable_metrics(
        speciesName, "Climate", "Europe", not_evaluable_plan(reason)
      )
    }
  } else {
    failed_plan <- not_evaluable_plan(climate_run$result$failure_reason)
    global_climate_metrics <- not_evaluable_metrics(
      speciesName, "Climate", "Global", failed_plan
    )
    eu_climate_metrics <- not_evaluable_metrics(
      speciesName, "Climate", "Europe", failed_plan
    )
  }

  climate_validation_dir <- file.path(
    base_dir, "Climate", "Current", "Diagnostics", "Model_validation"
  )
  global_climate_summary <- write_metric_outputs(
    global_climate_metrics,
    climate_validation_dir,
    paste0(speciesName, "_global_climate_validation_per_fold.csv"),
    paste0(speciesName, "_global_climate_validation_summary.csv")
  )
  eu_climate_summary <- write_metric_outputs(
    eu_climate_metrics,
    climate_validation_dir,
    paste0(speciesName, "_eu_climate_validation_per_fold.csv"),
    paste0(speciesName, "_eu_climate_validation_summary.csv")
  )
  species_validation_summary <- dplyr::bind_rows(
    species_validation_summary,
    global_climate_summary,
    eu_climate_summary
  )

  # Habitat and combined validation use an independent European plan.
  if (has_habitat) {
    eu_presabs <- habitat_model_artifact$eu_presabs
    eu_presabs$ID <- paste0("EU_MODEL_", seq_len(nrow(eu_presabs)))
    eu_habitat_filtered <- filter_records_for_stack(
      eu_presabs, habitat_selection
    )
    eu_climate_filtered <- filter_records_for_stack(
      eu_presabs, climate_selection
    )
    common_eu_ids <- intersect(
      eu_habitat_filtered$records$ID,
      eu_climate_filtered$records$ID
    )
    eu_records <- eu_presabs[
      eu_presabs$ID %in% common_eu_ids,
      c("species", "ID"),
      drop = FALSE
    ]
    diagnostic_index <- diagnostic_index + 1L
    species_diagnostics[[diagnostic_index]] <- predictor_filter_diagnostics(
      eu_presabs,
      eu_habitat_filtered$records,
      speciesName,
      "eu_component",
      "habitat"
    )
    diagnostic_index <- diagnostic_index + 1L
    species_diagnostics[[diagnostic_index]] <- predictor_filter_diagnostics(
      eu_presabs,
      eu_climate_filtered$records,
      speciesName,
      "eu_component",
      "climate"
    )
    diagnostic_index <- diagnostic_index + 1L
    species_diagnostics[[diagnostic_index]] <- predictor_filter_diagnostics(
      eu_presabs,
      eu_records,
      speciesName,
      "eu_component",
      "matched_climate_habitat"
    )

    european_plan_records <- rbind(
      make_plan_context(global_records, "global", climate_selection[[1]]),
      make_plan_context(
        sf::st_transform(eu_records, sf::st_crs(global_records)),
        "eu_component",
        climate_selection[[1]]
      )
    )
    european_plan <- make_preferred_cv_plan(
      records = european_plan_records,
      requested_folds = cv_folds,
      block_sizes_m = cv_spatial_block_sizes_m,
      enable_kfold_fallback = enable_kfold_cv_fallback,
      min_train = minimum_test_class_records,
      min_test = minimum_test_class_records,
      seed = cv_seed,
      iteration = cv_iterations
    )
    european_run <- run_with_runtime_fallback(
      initial_plan = european_plan,
      plan_records = european_plan_records,
      runner = function(plan) {
        run_european_plan(
          plan = plan,
          plan_records = european_plan_records,
          global_records = global_records,
          eu_records = eu_records,
          climate_selection = climate_selection,
          habitat_selection = habitat_selection,
          climate_methods = climate_methods,
          habitat_methods = habitat_methods,
          climate_background = eu_component_climate_background,
          habitat_background = eu_component_habitat_background,
          species_name = speciesName
        )
      },
      species_name = speciesName,
      context = "european_component"
    )
    diagnostic_index <- diagnostic_index + 1L
    species_diagnostics[[diagnostic_index]] <- european_run$diagnostics

    if (isTRUE(european_run$result$valid)) {
      habitat_metrics <- plan_metadata(
        european_run$result$habitat_metrics, european_run$plan
      )
      ensemble_metrics <- plan_metadata(
        european_run$result$ensemble_metrics, european_run$plan
      )
    } else {
      failed_plan <- not_evaluable_plan(european_run$result$failure_reason)
      habitat_metrics <- not_evaluable_metrics(
        speciesName, "Habitat", "Europe", failed_plan
      )
      ensemble_metrics <- not_evaluable_metrics(
        speciesName, "Ensemble", "Europe", failed_plan
      )
    }

    habitat_validation_dir <- file.path(
      base_dir, "Habitat", "Current", "Diagnostics", "Model_validation"
    )
    habitat_summary <- write_metric_outputs(
      habitat_metrics,
      habitat_validation_dir,
      paste0(speciesName, "_habitat_validation_per_fold.csv"),
      paste0(speciesName, "_habitat_validation_summary.csv")
    )
    ensemble_validation_dir <- file.path(
      base_dir, "Combined", "Current", "Diagnostics", "Model_validation"
    )
    ensemble_summary <- write_metric_outputs(
      ensemble_metrics,
      ensemble_validation_dir,
      paste0(speciesName, "_combined_validation_per_fold.csv"),
      paste0(speciesName, "_combined_validation_summary.csv")
    )
    species_validation_summary <- dplyr::bind_rows(
      species_validation_summary,
      habitat_summary,
      ensemble_summary
    )
  }

  validation_summary[[speciesName]] <- species_validation_summary
  species_diagnostics_df <- dplyr::bind_rows(species_diagnostics)
  validation_diagnostics[[speciesName]] <- species_diagnostics_df
  readr::write_csv(
    species_diagnostics_df,
    fortify_output_path(file.path(
      validation_dir,
      paste0(speciesName, "_cross_validation_diagnostics.csv")
    ))
  )

  rm(climate_model_artifact, habitat_model_artifact)
  terra::tmpFiles(remove = TRUE)
  gc()
}


#-------------------------------------------------------------------------------
# Combined project outputs
#-------------------------------------------------------------------------------
final_validation_extended <- dplyr::bind_rows(validation_summary)
readr::write_csv(
  final_validation_extended,
  fortify_output_path(file.path(validation_dir, "Validation_summary_extended.csv"))
)
final_validation <- as_legacy_validation_summary(final_validation_extended)
readr::write_csv(
  final_validation,
  fortify_output_path(file.path(validation_dir, "Validation_summary.csv"))
)
final_diagnostics <- dplyr::bind_rows(validation_diagnostics)
readr::write_csv(
  final_diagnostics,
  fortify_output_path(file.path(validation_dir, "Cross_validation_diagnostics.csv"))
)
