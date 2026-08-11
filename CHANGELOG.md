# Changelog

All notable changes to the OneSTOP adaptation of wiSDM are documented here.
The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/).
The upstream release lineage is TrIAS wiSDM v2.0.0, and the exact upstream
development base is commit `2c7bc8712304e46bebd1bbcf2f90771c4814f865`.

## [Unreleased]

No changes yet.

## [2.1.0] - Forthcoming

Planned Git tag: `onestop-v2.1.0`.

### Added

- Manifest-driven support for user-supplied climate and land-cover raster
  inputs, including templates and validation of paths, predictor names,
  periods, scenarios, coordinate reference systems, and raster structure.
- Support for custom European and country boundaries and alignment of custom
  rasters to the active modelling grid.
- A functional cross-validation stage for climate, habitat, and combined
  suitability models, including per-fold and summary outputs.
- A chunked multi-species runner with preflight checks, locking, retries,
  recovery, stage reuse, dry-run support, and explicit full-model or
  climate-only completion states.
- A Shiny monitoring dashboard with progress and log inspection, stale and
  partial-run detection, and a no-overwrite workflow for importing and merging
  results produced across machines.
- Regression tests for cross-validation fallbacks, safe raster output,
  climate-only continuation, and model-explainability robustness.
- Post-modelling utilities for model checks, occurrence summaries, scenario
  statistics, raster discovery, validation summaries, and map production.
- Configuration options for custom inputs, custom European boundaries, and
  near-zero-variance habitat-predictor filtering.

### Changed

- Extended future projections to support custom climate inputs and either
  static or scenario-specific custom land-cover predictors.
- Strengthened raster cropping, masking, alignment, naming, writing, and
  overwrite verification for portable and deterministic outputs.
- Expanded saved model metadata so downstream stages can verify input modes,
  manifests, coordinate reference systems, and predictor names.
- Improved chunk execution and dashboard control files for clearer audit,
  recovery, and distributed-run handling.
- Updated the README with OneSTOP Task 5.1 scope, upstream provenance,
  operational implications, funding acknowledgement, and contributor details.
- Updated `.gitignore` to exclude generated dashboard error and output logs,
  local tooling files, custom data manifests, and downloaded boundary data.

### Fixed

- Made cross-validation adapt spatial block sizes and fold counts, keep
  duplicate locations together, and optionally fall back to grouped,
  stratified non-spatial folds.
- Prevented redundant `terra::vect()` conversion of background samples that
  are already `SpatVector` objects.
- Allowed stage 04 to finish cleanly when no European occurrences are
  available and continue through climate-only validation and retry workflows.
- Isolated algorithm-specific prediction, variable-importance, and response-
  curve failures so usable model results are retained.
- Added fallbacks for infeasible k-means thinning, unavailable taxon-specific
  bias grids, and individual model prediction failures.
- Fortified species-derived filenames and output paths against authorship,
  accents, hybrid markers, unsafe characters, and excessive path lengths.

### Removed

- Deprecated, inactive, duplicate, backup, and `_old` scripts that are not
  part of the OneSTOP v2.1.0 release workflow.
- Generated dashboard `.err.log` and `.out.log` files from version control.

### Security

- Removed embedded Zenodo access tokens from the active scripts.
- Read the Zenodo token from the `ZENODO_TOKEN` environment variable and stop
  with a clear error when it is not configured.

[Unreleased]: https://github.com/OneSTOP-CIBIO/risk-modelling-and-mapping/compare/onestop-v2.1.0...HEAD
[2.1.0]: https://github.com/OneSTOP-CIBIO/risk-modelling-and-mapping/compare/2c7bc8712304e46bebd1bbcf2f90771c4814f865...onestop-v2.1.0
