suppressPackageStartupMessages(library(terra))

source(file.path("src", "helper_functions.R"))

failures <- character(0)
tests_run <- 0L

assert_true <- function(condition, message) {
  if (!isTRUE(condition)) stop(message, call. = FALSE)
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

absolute_chars <- function(path) {
  nchar(normalizePath(path, winslash = "/", mustWork = FALSE))
}

run_test("ordinary species names retain the legacy stem", {
  assert_true(
    identical(
      species_output_stem("Ailanthus altissima (Mill.) Swingle"),
      "Ailanthus_altissima"
    ),
    "An ordinary binomial did not retain its historical output stem."
  )
})

run_test("hybrid markers and author citations are removed from filesystem stems", {
  names_to_test <- c(
    "Reynoutria × bohemica Chrtek & Chrtková",
    "Reynoutria ×bohemica Chrtek & Chrtková",
    "Reynoutria x bohemica Chrtek & Chrtková"
  )
  stems <- vapply(names_to_test, species_output_stem, character(1))
  assert_true(
    all(stems == "Reynoutria_bohemica"),
    "Hybrid spellings did not resolve to Reynoutria_bohemica."
  )

  parts <- scientific_name_output_parts(names_to_test[[1]])
  assert_true(
    identical(parts$taxon_label, "Reynoutria × bohemica") &&
      identical(parts$authorship, "Chrtek & Chrtková"),
    "The authoritative hybrid label or author citation was not preserved."
  )

  canonical_path <- file.path(
    "data", "projects", "onestop_03_reynoutria_bohemica",
    "Reynoutria_bohemica_4038485", "Climate", "2041-2070", "ssp126",
    "Diagnostics", "Confidence_maps", "Rasters",
    "Reynoutria_bohemica_Climate_2041-2070_ssp126_ensemble_SD.tif"
  )
  assert_true(
    absolute_chars(canonical_path) <= 240L &&
      identical(fortify_output_path(canonical_path), canonical_path),
    "The corrected reported path is not portable or was changed unexpectedly."
  )
})

run_test("taxon tokens are transliterated to portable ASCII", {
  assert_true(
    identical(species_output_stem("Caféa spécies Autôr"), "Cafea_species"),
    "Accented taxon tokens were not transliterated."
  )
})

run_test("invalid scientific names fail early", {
  expect_error(species_output_stem(""), "non-empty")
  expect_error(species_output_stem("Reynoutria"), "genus-and-species")
})

run_test("long species stems use stable collision-resistant digests", {
  first <- species_output_stem(
    paste0(strrep("Genus", 10), " ", strrep("species", 10), " Author"),
    max_chars = 40L
  )
  first_again <- species_output_stem(
    paste0(strrep("Genus", 10), " ", strrep("species", 10), " Author"),
    max_chars = 40L
  )
  same_taxon_other_author <- species_output_stem(
    paste0(strrep("Genus", 10), " ", strrep("species", 10), " DifferentAuthor"),
    max_chars = 40L
  )
  second <- species_output_stem(
    paste0(strrep("Genus", 10), " ", strrep("species", 9), "x Author"),
    max_chars = 40L
  )
  assert_true(nchar(first) == 40L, "The long stem was not trimmed to max_chars.")
  assert_true(identical(first, first_again), "Stem truncation was not deterministic.")
  assert_true(
    identical(first, same_taxon_other_author),
    "Author citations changed a truncated taxon stem."
  )
  assert_true(!identical(first, second), "Distinct long names produced the same stem.")
})

run_test("short safe paths remain byte-for-byte unchanged", {
  path <- file.path(tempdir(), "Ailanthus_altissima_Climate_current_ensemble.tif")
  assert_true(
    identical(fortify_output_path(path), path),
    "A short legacy-compatible path was changed."
  )
})

run_test("special and overlong leaf names are fortified deterministically", {
  directory <- file.path(
    tempdir(),
    paste(rep("segment1234", 8), collapse = .Platform$file.sep)
  )
  logical <- file.path(
    directory,
    paste0(
      paste(rep("Reynoutria × bohemica Chrtek & Chrtková", 5), collapse = "_"),
      "_Climate_2041-2070_ssp126_ensemble_SD.tif"
    )
  )
  resolved <- fortify_output_path(logical)
  resolved_again <- fortify_output_path(logical)

  assert_true(absolute_chars(logical) > 260L, "The logical test path is not long enough.")
  assert_true(absolute_chars(resolved) <= 240L, "The resolved path exceeds 240 characters.")
  assert_true(identical(resolved, resolved_again), "Path fortification was not deterministic.")
  assert_true(grepl("ensemble_SD\\.tif$", resolved), "The artifact suffix was not retained.")
  assert_true(
    !grepl("[×&á[:space:]]", basename(resolved)),
    "The resolved filename still contains unsafe characters."
  )
})

run_test("directories without a safe filename budget fail early", {
  directory <- file.path(tempdir(), paste(rep("segment1234", 20), collapse = "/"))
  expect_error(validate_output_directory(directory), "insufficient room")
})

run_test("Windows-incompatible directory components fail early", {
  expect_error(
    validate_output_directory(file.path(tempdir(), "bad?directory")),
    "Windows-incompatible"
  )
})

run_test("safe raster writing handles long paths and verified overwrites", {
  directory <- file.path(
    tempdir(),
    paste(rep("rasterseg123", 8), collapse = .Platform$file.sep)
  )
  dir.create(directory, recursive = TRUE, showWarnings = FALSE)
  logical <- file.path(
    directory,
    paste0(strrep("LongSpeciesArtifact_", 9), "ensemble_SD.tif")
  )
  target <- fortify_output_path(logical)
  raster <- terra::rast(nrows = 3, ncols = 4, vals = 1:12)

  written <- write_raster_safely(raster, logical, overwrite = TRUE)
  reopened <- terra::rast(written)
  assert_true(identical(written, target), "The writer returned an unexpected target path.")
  assert_true(file.exists(written) && file.info(written)$size > 0, "Raster was not written.")
  assert_true(
    terra::nlyr(reopened) == terra::nlyr(raster) &&
      isTRUE(terra::compareGeom(reopened, raster, lyrs = FALSE, stopOnError = FALSE)),
    "The reopened raster does not match the source geometry."
  )

  terra::values(raster) <- 12:1
  write_raster_safely(raster, logical, overwrite = TRUE)
  assert_true(
    identical(as.numeric(terra::values(terra::rast(target))), as.numeric(12:1)),
    "The verified overwrite did not replace raster values."
  )
  temporary_files <- list.files(directory, pattern = "^\\.w[rb]_", all.files = TRUE)
  assert_true(length(temporary_files) == 0L, "Staging or backup files were not cleaned up.")
})

run_test("unrecoverable writes report diagnostics and preserve existing output", {
  directory <- file.path(tempdir(), "raster_failure_diagnostics")
  dir.create(directory, recursive = TRUE, showWarnings = FALSE)
  target <- file.path(directory, "existing.tif")
  raster <- terra::rast(nrows = 2, ncols = 2, vals = 1:4)
  write_raster_safely(raster, target, overwrite = TRUE)
  original_values <- as.numeric(terra::values(terra::rast(target)))

  error <- expect_error(
    write_raster_safely("not-a-raster", target, overwrite = TRUE),
    "underlying write error"
  )
  message_text <- conditionMessage(error)
  assert_true(
    grepl("Logical path:", message_text, fixed = TRUE) &&
      grepl("Resolved path:", message_text, fixed = TRUE) &&
      grepl("chars/", message_text, fixed = TRUE),
    "The write failure omitted path diagnostics."
  )
  assert_true(
    identical(as.numeric(terra::values(terra::rast(target))), original_values),
    "A failed overwrite damaged the previous valid raster."
  )
  temporary_files <- list.files(directory, pattern = "^\\.w[rb]_", all.files = TRUE)
  assert_true(length(temporary_files) == 0L, "Failed-write temporary files were not cleaned up.")
})

if (length(failures) > 0L) {
  stop(
    paste0(
      length(failures), " of ", tests_run, " raster-output test(s) failed:\n- ",
      paste(failures, collapse = "\n- ")
    ),
    call. = FALSE
  )
}

message("All ", tests_run, " raster-output fortification tests passed.")
