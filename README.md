
# wiSDM OneSTOP v2.1.0 (Task 5.1)  — Alien species risk modelling and mapping

This repository contains the framework and R code for predicting the distribution of alien species in Europe at 1 km<sup>2</sup> (climate model) and 1 km<sup>2</sup> resolution (habitat/land cover model) as part of the TrIAS project and **later adapted to the OneSTOP Task 5.1 modelling requirements**.    

More specifically, this is the OneSTOP downstream adaptation of TrIAS Alien species risk modelling and mapping, wiSDM v2.0.0:

- Upstream repository: https://github.com/trias-project/risk-modelling-and-mapping
- Upstream release lineage: v2.0.0
- Exact upstream development base: 2c7bc8712304e46bebd1bbcf2f90771c4814f865
- OneSTOP release: onestop-v2.1.0


## Repo structure

```
├── README.md              : Description of this repository
├── LICENSE                : Repository license
├── risk-modelling-and-mapping.Rproj : RStudio project file
├── .gitignore             : Files and directories to be ignored by git
│
├── data
│   ├── external          : external files required to run the model. The majority of these files will be downloaded and stored in the right folders by running script 01_prepare_files_and_folders.R.
│   
│
└── src                    : R Code
```


## Requirements to run this workflow
1.  **RStudio** installed on your local computer.
2.   Have an active **GBIF account**. Before running the workflow for the first time, store your GBIF username, password, and email address in your `~/.Renviron` file as instructed in the script 00_configurations.R.
3. **Clone this repository** to your local computer.
4. Set your specific configurations for the workflow in the **00_configurations** script (see below). 
<br>

## Prior to running the workflow

Before you execute this workflow, specify the following configurations in the **00_configurations script**, which is stored in the `src` folder:

* **project**: The name of your project. A folder with this name will be created automatically in the `./data/projects` folder, and all workflow outputs will be stored there.

* **species_to_model**: Specify the species you want to model as a character string (e.g., `"Vespa velutina"`). Multiple species can be provided as a character vector (e.g., `c("Vespa velutina", "Aedes albopictus")`). Use Latin binomials only (no authorship or year).

* **occurrence_thinning_method**: When a species has more than 10000 occurrences, records are thinned to 10000 to match the number of pseudoabsences. Options:<ul>
<li><code>"random"</code>: randomly samples 10000 occurrences.</li>
<li><code>"kmeans_clustering"</code>: performs k-means clustering in environmental space and selects 10000 cluster centroids, ensuring the thinned occurrences represent the broadest environmental variation.</li>
</ul>

* **mtp_probabilities**: Defines the minimum training presence (MTP) thresholds used to convert continuous favorability predictions into binary presence/absence maps. For each value in `mtp_probabilities`, the workflow removes the lowest *x*% of occurrence probabilities and uses the next-lowest value as the threshold. Example: `mtp_probabilities = c(0.01, 0.05)` will produce binarized maps where the threshold corresponds to the lowest favorability that remains after the 1% and 5% lowest-favorability occurrences are removed.

* **country_of_interest**: Either "Europe" or a European country. If a specific European country is provided, all output maps (PNG, PDF, and raster format) will be generated for that country only.

* **update_files**: Whether or not all needed raster input files should be downloaded again. Options are *"ask"*: the user will be prompted if a specific file should be downloaded again (any missing files are downloaded by default), *"yes"*: all files are downloaded again and *"no"*: only missing files are downloaded. When you select "ask", a pop-up window will appear for each layer asking whether it should be downloaded again. The workflow will pause and cannot continue automatically until you respond to each prompt. Note that these pop-ups may sometimes open behind the RStudio window, so you may need to minimize RStudio to locate them. 

## Executing the workflow

To execute this workflow, run the **06_run_wiSDM.R** script, stored in the `src` folder. This script automatically runs the following scripts in the designated order:

0. **Script 00_configurations**: Specifies the workflow's configurations (e.g., project name, species to model, thinning method, MTP threshold settings). IMPORTANT: users must actively set these fields themselves prior to running the workflow.
1. **Script 01_prepare_files_and_folders**: Sets up the folder structure and downloads the files (climate rasters, habitat predictors, spatial boundaries,...) necessary to run the workflow.
2. **Script 02_global_occurrence_download.R**: Retrieves occurrence data for the species of interest defined in script 00 from the Global Biodiversity Information Facility (GBIF). 
3. **Script 03_fit_climate_model.R**: Builds a global-scale climate-only species distribution model (SDM) for each species of interest, at a resolution of 1 km<sup>2</sup>. The results of this model are presented at the level of the country of interest that was specified in the configurations script and they can be found in the folder `./data/projects/<your project>/species/Climate`.
4. **Script 04_fit_habitat_model.R**: Generates a European-scale habitat-only species distribution model for the specified species at a resolution of 1 km<sup>2</sup> and integrates these predictions with the 1 km² climate-only predictions from script 03. The two prediction layers are integrated by using the geometric mean to generate a final suitability map at 1 km<sup>2</sup> resolution that reflects both climate suitability and habitat (land cover) suitability. Final predictions are generated for both current conditions and for two future periods, 2041-2070 and 2071-2100, under different climate change scenarios (SSP1-2.6, SSP3-7.0, and SSP5-8.5). The results of the habitat model can be found in the folder `./data/projects/<your project>/species/Habitat`, while the final predictions, combining both the habitat and the climate suitability, can be found in the `./data/projects/<your project>/species/Combined` folder.
<br>
 
## What does the TrIAS/wiSDM-2.0 modelling workflow do?
1.	**Generates habitat suitability maps using machine learning.**
The workflow requires only a species name and then fits an ensemble of ten machine-learning algorithms to estimate both climate suitability (using global occurrences) and habitat suitability (using European occurrences). For each type (habitat and climate), the favourability maps from all ten algorithms are first stacked and analysed using Principal Component Analysis (PCA). For each algorithm, the spatial variance of its predictions projected onto the first principal component (PC1) is calculated. The five algorithms with the highest variance along PC1, representing those that contribute most strongly to the dominant pattern of variation across the study area, are then used to generate the final suitability map. These two maps (one for the climate model and one for the habitat model) are then combined using the geometric mean in order to produce a final overall suitability prediction. All maps are generated automatically for current conditions. Climate suitability maps, along with the final combined maps, are also produced for three standard Shared Socioeconomic Pathways (SSP1-2.6, SSP3-7.0, and SSP5-8.5) for the periods 2041–2070 and 2071–2100.
2.	**Implements best practices for pseudoabsence placement:** 
    * **Climate model**: Pseudoabsences are sampled within the same biomes as species presences but excluded from presence grid cells. A taxonomic sampling effort grid (bias grid) captures the sampling intensity of the higher taxon and is used to weight grid cells, assigning greater weight to well-sampled areas.
    * **Habitat model**: Pseudoabsences are sampled across all of Europe, excluding presence grid cells. In ecoregions with presences, grid cells are weighted using the bias grid, while outside these ecoregions all cells are assigned the minimum weight (1).
3.	**Detects and removes highly correlated predictors.**<br>
Highly correlated predictors can have undesirable effects and confuse the interpretation of variable importance.
4.	**Automatic generation of confidence maps** for each suitability map. These maps illustrate prediction uncertainty across the study area by calculating the population standard deviation of the predictions produced by the five algorithms for both the climate and habitat models. The standard deviation of the combined prediction is then computed using the appropriate error-propagation formula for the geometric mean. Note that the standard deviation is always calculated from the mean of the five algorithms, whereas the final suitability predictions are based on their median. Consequently, the confidence maps provide a relative measure of model uncertainty but should not be interpreted as the uncertainty of the median prediction itself.
<br>

## Contributors

[List of contributors](https://github.com/trias-project/risk-modelling-and-mapping/contributors)
<br>

## References
Davis AJS, Groom Q, Adriaens T, Vanderhoeven S, De Troch R, Oldoni D, Desmet P, Reyserhove L, Lens L and Strubbe D (2024) Reproducible WiSDM: a workflow for reproducible invasive alien species risk maps under climate change scenarios using standardized open data. Front. Ecol. Evol. 12:1148895. doi: [10.3389/fevo.2024.1148895](https://doi.org/10.3389/fevo.2024.1148895)

**Important note**: The wiSDM version on the current main branch (v2.0.0) represents a complete restructuring and substantial update of the workflow and does not exactly reproduce the implementation described in Davis et al. (2024). The published workflow is preserved in release v1.1.0.


----
    
![](img/OneSTOP_logo_resize_w200.png)     

    
# Adaptations to the wiSDM 2.0 framework under OneSTOP Task 5.1

## Summary

The adapted framework keeps the core wiSDM 2.0 modelling logic: climate suitability, habitat suitability, and combined suitability are still produced for invasive alien species using the established staged workflow. The main changes make the workflow more flexible for project-specific inputs, more robust for large batches of species, and easier to monitor, validate, retry, and merge across multiple runs or machines.

The most relevant additions are:

1. Support for user-supplied climate and land-cover raster data.
2. Support for custom Europe and country boundaries.
3. A functional cross-validation stage integrated into the main workflow.
4. A chunked multi-species execution framework with logging, recovery, and a Shiny dashboard.
5. Several robustness improvements for raster alignment, pseudoabsence selection, failed model predictions, and batch reruns.

## Main alterations

### 1. User-specific environmental input data

The stock workflow was mainly wired to the default downloadable CHELSA climate layers and default land-cover predictors. The adapted version adds optional manifest-driven inputs:

- `user_specific_climate_data` allows the climate model to use a CSV manifest of user-provided climate rasters instead of the default CHELSA inputs.
- `user_specific_landcover_data` allows the habitat model to use a CSV manifest of user-provided land-cover or habitat rasters instead of the default downloaded land-cover stack.
- New template files document the required manifest structure:
  - `src/user_specific_climate_data_template.csv`
  - `src/user_specific_climate_data_template.md`
  - `src/user_specific_landcover_data_template.csv`
  - `src/user_specific_landcover_data_template.md`
- The framework validates period/scenario combinations, predictor names, file paths, CRS consistency, raster alignment, and single-layer raster assumptions.
- User-specific climate inputs support current and future combinations for the standard periods and scenarios.
- User-specific land-cover inputs support either current-only habitat predictors or dynamic future habitat predictors for all future period/scenario combinations.

This is one of the most important adaptations because it decouples the framework from the stock input data and enables project-specific IAS distribution modelling with bespoke environmental layers.

### 2. Custom spatial boundaries and raster alignment

The stock workflow assumed the default Europe/country spatial setup. The adapted version adds:

- `custom_eu_boundary_path` for a user-defined European modelling boundary.
- Continued support for `custom_country_boundary_path`, with broader use across the climate, habitat, and validation scripts.
- Helper functions to load, transform, crop, mask, and align boundary layers to the active raster CRS.
- Logic to handle user-specific rasters whose CRS, extent, and resolution may differ from the stock rasters.

These changes are important for project deployments that need a different geographic definition of Europe, a filtered set of countries, or input rasters in a non-stock projection.

### 3. Functional cross-validation integrated into the workflow

`src/06_run_wiSDM.R` now sources `src/05_cross_validation.R`, making validation part of the main staged workflow.

The cross-validation stage was expanded to:

- Validate climate-only, habitat-only, and combined climate-habitat suitability outputs.
- Work with both stock and user-specific climate and land-cover inputs.
- Work with custom Europe boundaries.
- Produce per-fold and summary validation outputs.
- Use robust Boyce-index and threshold-independent validation calculations.
- Fall back to climate-only validation where habitat or combined validation is not possible.
- Reuse saved model metadata, including input mode, manifest path, predictor CRS, and selected predictors.

This materially improves the internal evidence base for IAS distribution outputs, especially when modelling many species or comparing stock versus custom input data.

### 4. Chunked multi-species execution framework

A new chunk runner was added in `src/_run_wisdm_by_chunks.R` to support large species batches.

Key features include:

- Species-list driven execution from an Excel file, with optional filtering.
- Manual chunk assignment for running multiple RStudio background jobs in parallel.
- Automatic per-species project naming with sanitised and shortened names.
- Shared stage 01 setup with locking, so multiple chunk workers do not repeat setup unnecessarily.
- Registry, event, completed-species, and log files under `data/projects/<project_prefix>_chunk_control/`.
- Retry controls for failed or skipped species.
- Optional verification that completed outputs still exist before skipping reruns.
- Optional reuse of successful previous stages on retry attempts.
- Classification of completed models as full models or climate-only outputs.
- Dry-run and preflight modes for checking configuration before running.

This is a major operational adaptation. It turns the stock single/batch workflow into a more resilient production-style system for many IAS models.

### 5. Monitoring dashboard and import/merge workflow

Two new dashboard-related scripts were added:

- `src/_wisdm_chunk_dashboard.R`
- `src/_launch_wisdm_chunk_dashboard.R`

The dashboard reads the chunk-control files and provides:

- Progress monitoring by chunk, species, stage, status, and model outcome.
- Log inspection from the browser.
- Detection of active, completed, failed, skipped, stale, and partial climate-only runs.
- A detached launcher so monitoring does not interfere with RStudio background jobs.
- An Imports / Merge workflow for combining run outputs from different machines or transfer batches.
- A no-overwrite merge strategy that merges control CSVs/logs into a master control folder and adopts large species project folders by moving them into `data/projects/`.

This is especially relevant for internal delivery because it supports distributed or long-running IAS modelling campaigns and creates a clearer audit trail.

### 6. Model robustness improvements

Several targeted changes make the workflow less brittle:

- K-means occurrence and pseudoabsence thinning now uses a fallback that reduces cluster centers when the requested number is not feasible.
- Climate and habitat prediction code is more defensive around failed algorithms and raster prediction blocks.
- Cross-validation can use a safe median-favourability function that tolerates individual method prediction failures, while requiring a minimum number of successful methods.
- Habitat modelling can optionally drop near-zero-variance predictors before fitting, reducing instability from low-information moving-window land-cover variables.
- Bias grids are better aligned to the active climate/habitat raster grid and CRS.
- When a taxon-specific bias grid is unavailable, the workflow can fall back to a neutral bias-grid approach.
- Additional metadata is saved with model objects so downstream stages can check input modes, manifests, CRS, and predictor names.

These changes mainly reduce batch failures and improve reproducibility when species have sparse records, unusual predictor values, or custom raster inputs.

### 7. Future projections with custom inputs

The habitat and combined prediction logic was extended so future projections can use:

- Stock climate projections with static habitat predictors.
- User-specific climate projections.
- User-specific current-only land-cover predictors, treated as static for future projections.
- User-specific dynamic future land-cover predictors, when all required future combinations are provided.

The framework aligns habitat and climate rasters before combining them with the existing geometric-mean approach, preserving the original wiSDM concept while making the input pipeline more flexible.

### 8. Configuration and supporting files

The configuration script now includes new user-facing settings:

- `custom_eu_boundary_path`
- `user_specific_climate_data`
- `user_specific_landcover_data`
- `habitat_filter_near_zero_variance_predictors`

Additional supporting changes include:

- A species-list workbook added under `data/external/Species_list_v5.xlsx` for chunked runs.
- `.gitignore` updates for local Posit/Codex files, custom raster manifests, and downloaded GADM data.
- More helper functions in `src/helper_functions.R` for manifests, raster stack naming, raster alignment, custom boundaries, safe favourability prediction, validation metrics, and geometric-mean ensemble validation.

## Latest corrections and improvements

The `new-features/misc-adjusts-bug-corrections` branch adds the following targeted changes on top of `new-features/merge-chunks-dashboard`:

- Cross-validation is now more defensive: it adapts spatial block size and fold count, can fall back to grouped stratified k-fold validation, keeps duplicate locations in the same fold, records why data or partitions are not evaluable, and tolerates individual algorithm prediction failures. New configuration options control the number of folds, candidate block sizes, and whether the non-spatial fallback is enabled.
- A `terra::vect()` conversion bug in cross-validation was corrected so background samples already represented as `SpatVector` objects are not converted a second time.
- Raster outputs now use portable, deterministic species stems and fortified file paths. Hybrid markers, authorship, accents, unsafe characters, and excessive Windows path lengths are handled consistently, while writes and overwrites are verified and failures preserve existing outputs.
- Stage 04 now exits cleanly when a species has no European occurrences: the absence of a habitat model is recorded as an intentional climate-only outcome, stage 05 continues with climate validation, and chunk retries can reuse that decision.
- Climate and habitat variable-importance and response-curve generation now isolates failures by model and diagnostic type. Valid results are retained and averaged by algorithm, while unavailable diagnostics produce explicit audit information instead of aborting the modelling stage.
- Chunk execution accepts a directly defined, normalised species vector, uses the same safe species naming as model outputs, and validates the new cross-validation settings during preflight. The custom-climate manifest filename and related ignored local documentation entries were also adjusted.
- Regression tests were added for cross-validation partitioning and fallbacks, safe raster naming/writing, climate-only continuation, and robust model-explainability outputs.

## Practical implications

Compared with stock wiSDM 2.0, the adapted branch is better suited to internal IAS distribution modelling where:

- input rasters are supplied by the project instead of downloaded from stock sources;
- the modelling region differs from the default wiSDM Europe definition;
- many species need to be processed in parallel or across machines;
- interrupted or partial runs must be retried without losing previous successful stages;
- validation outputs are required as part of the modelling deliverable;
- the modelling team needs a readable execution audit trail through registries, logs, and dashboard summaries.

The core modelling structure remains recognisably wiSDM 2.0, but the surrounding data-ingestion, execution, validation, and monitoring layers have been substantially strengthened for operational IAS modelling.


----

## OneSTOP Task 5.1 contributors

[List of contributors](https://github.com/trias-project/risk-modelling-and-mapping/contributors)



## License

[MIT License](https://github.com/OneSTOP-CIBIO/risk-modelling-and-mapping/blob/master/LICENSE)

© 2026 — Developed under the OneSTOP Project

![](img/eu_funded_en.png)

OneSTOP receives funding from the European Union Horizon Europe Research and Innovation Programme (ID No 101180559). Views and opinions expressed are those of the author(s) only and do not necessarily reflect those of the European Union or the European Research Executive Agency (REA). Neither the EU nor REA can be held responsible for them.
