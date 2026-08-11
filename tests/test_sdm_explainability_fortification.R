suppressPackageStartupMessages(library(dplyr))

source(file.path("src", "helper_functions.R"))

failures <- character(0)
tests_run <- 0L

assert_true <- function(condition, message) {
  if (!isTRUE(condition)) stop(message, call. = FALSE)
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

methods::setClass(
  "fakeSdmResponseCurve",
  slots = c(response = "list")
)
methods::setClass(
  "fakeSdmVarImportance",
  slots = c(varImportance = "data.frame")
)

fake_response <- function(values, responses = values) {
  methods::new(
    "fakeSdmResponseCurve",
    response = list(
      predictor_1 = data.frame(
        predictor_1 = values,
        response = responses
      )
    )
  )
}

fake_importance <- function(cor_test, auc_test = cor_test) {
  methods::new(
    "fakeSdmVarImportance",
    varImportance = data.frame(
      variables = "predictor_1",
      corTest = cor_test,
      AUCtest = auc_test,
      stringsAsFactors = FALSE
    )
  )
}

run_test("diagnostic failures are isolated by model ID and metric", {
  model_info <- data.frame(
    modelID = 1:4,
    method = c("glm", "rf", "cart", "mars"),
    success = TRUE,
    stringsAsFactors = FALSE
  )

  response_getter <- function(model, id) {
    switch(
      as.character(id),
      `1` = fake_response(1:3, c(0.2, 0.4, 0.6)),
      `2` = fake_response(1:3, c(0.3, 0.5, 0.7)),
      `3` = stop("synthetic response failure", call. = FALSE),
      `4` = fake_response(1:3, rep(NA_real_, 3))
    )
  }
  importance_getter <- function(model, id) {
    switch(
      as.character(id),
      `1` = NULL,
      `2` = fake_importance(0.4, 0.2),
      `3` = stop("synthetic importance failure", call. = FALSE),
      `4` = fake_importance(0.6, 0.3)
    )
  }

  result <- collect_sdm_explainability(
    model = list(),
    methods = model_info$method,
    model_info = model_info,
    response_getter = response_getter,
    importance_getter = importance_getter
  )

  assert_true(
    setequal(unique(result$response_df$Algorithm), c("glm", "rf")),
    "Valid response curves were lost or invalid curves were retained."
  )
  assert_true(
    setequal(unique(result$varimp_df$Algorithm), c("rf", "mars")),
    "Valid importance rows were lost or invalid rows were retained."
  )
  assert_true(nrow(result$diagnostics) == 8L, "The per-metric audit is incomplete.")
  assert_true(
    any(grepl("returned NULL", result$diagnostics$Reason, fixed = TRUE), na.rm = TRUE),
    "The NULL importance result was not diagnosed."
  )
  assert_true(
    any(grepl("synthetic response failure", result$diagnostics$Reason, fixed = TRUE), na.rm = TRUE),
    "The response error was not retained."
  )
  assert_true(
    any(grepl("no finite", result$diagnostics$Reason, fixed = TRUE), na.rm = TRUE),
    "The all-NA response curve was not rejected."
  )
})

run_test("multiple successful model IDs are averaged by algorithm", {
  model_info <- data.frame(
    modelID = c(11L, 12L),
    method = c("glm", "glm"),
    success = TRUE,
    stringsAsFactors = FALSE
  )

  result <- collect_sdm_explainability(
    model = list(),
    methods = "glm",
    model_info = model_info,
    response_getter = function(model, id) {
      fake_response(c(1, 2), rep(if (id == 11L) 0.2 else 0.4, 2))
    },
    importance_getter = function(model, id) {
      fake_importance(if (id == 11L) 0.2 else 0.4)
    }
  )

  assert_true(
    all(abs(result$response_df$Response - 0.3) < 1e-12),
    "Response curves were not averaged across model IDs."
  )
  assert_true(
    nrow(result$varimp_df) == 1L && abs(result$varimp_df$corTest - 0.3) < 1e-12,
    "Variable importance was not averaged across model IDs."
  )
  assert_true(
    nrow(result$diagnostics) == 4L && all(result$diagnostics$Success),
    "Successful replicated diagnostics were not audited independently."
  )
})

run_test("total diagnostic failure returns typed empty tables", {
  model_info <- data.frame(
    modelID = 1L,
    method = "glm",
    success = TRUE,
    stringsAsFactors = FALSE
  )

  result <- collect_sdm_explainability(
    model = list(),
    methods = "glm",
    model_info = model_info,
    response_getter = function(model, id) stop("response unavailable", call. = FALSE),
    importance_getter = function(model, id) NULL
  )

  assert_true(
    identical(
      names(result$response_df),
      c("Algorithm", "Predictor", "Predictor_value", "Response")
    ),
    "The empty response schema changed."
  )
  assert_true(
    identical(
      names(result$varimp_df),
      c("Algorithm", "Predictor", "corTest", "AUCtest")
    ),
    "The empty importance schema changed."
  )
  assert_true(
    is.numeric(result$response_df$Response) && is.numeric(result$varimp_df$corTest),
    "Empty metric columns have the wrong types."
  )
  assert_true(
    nrow(result$diagnostics) == 2L && !any(result$diagnostics$Success),
    "Total failure did not produce two failed audit rows."
  )
})

if (length(failures) > 0L) {
  stop(
    paste0(
      length(failures), " of ", tests_run, " SDM explainability test(s) failed:\n- ",
      paste(failures, collapse = "\n- ")
    ),
    call. = FALSE
  )
}

message("All ", tests_run, " SDM explainability fortification tests passed.")
