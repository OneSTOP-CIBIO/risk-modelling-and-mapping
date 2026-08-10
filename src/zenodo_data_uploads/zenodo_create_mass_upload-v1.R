###############################################
# Automate CHELSA v2.1 uploads to Zenodo
# One draft deposition per scenario (no publish)
###############################################

# install.packages("zen4R")   # if not yet installed
library(zen4R)

## 1. AUTH + MANAGER ----------------------------------------------

# Put your personal access token here, or rely on zenodo_pat()
token <- "4Dqs063GQX1lFAzfPRA7NAKbGtBOpT8veuqJ46UVuNaVyEJxx3s7cgmSLDRX"

zen <- ZenodoManager$new(
  token  = token,
  logger = "INFO"   # "DEBUG" if you want very verbose logs
)

## 2. COMMON METADATA ---------------------------------------------

version_string <- "v.2.1-osv.0.1-062025"
created_date   <- "2025-06-19"
project_url    <- "https://onestop-project.eu/"

## FUNDING BLOCK ------------------------------------------------------
# Zenodo expects:
# award_number, award_title, and funder
funding_award_number <- "101180559"
funding_award_title  <- "OneSTOP — OneBiosecurity Systems and Technology for People, Places and Pathways"
funding_funder       <- "European Commission"

common_description <- 
  "
This dataset provides bioclimatic and environmental predictor variables used in Species Distribution Models (SDMs) developed within the OneSTOP Project (Task 5.1). The original climatic data are based on the CHELSA v2.1 BIOCLIM+ dataset, which provides high-resolution (~1 km) bioclimatic variables derived from downscaled and bias-corrected climate data.
<br />
For future projections, CHELSA v2.1 climate layers were obtained for all five available CMIP6 Global Circulation Models (GCMs): GFDL-ESM4, UKESM1-0-LL, MPI-ESM1-2-HR, IPSL-CM6A-LR, and MRI-ESM2-0) and corresponding to the relevant SSP scenario. To produce a consistent climatic baseline that matches land-cover projections in the Chen et al. (2022) dataset, all available GCMs were averaged to generate a single ensemble mean for each time period and scenario. This ensemble approach reduces individual model biases and aims to provide a robust representation of mid- and late-century climatic conditions for SDMs. Raster data has been internally scaled and reprojected (bilinear method) in the terra R package.
<br />
All raster layers are provided as GeoTIFF (float) files in the coordinate reference system EPSG:6933, WGS 1984/NSIDC EASE-Grid 2.0 Global (Cylindrical Equal Area projection). The datasets correspond to one of the following temporal windows and SSP scenarios:
<br />
- Historical Baseline: 1981–2010 (used for model training)<br />
- Future Mid-Century (2041–2070): SSP1–2.6, SSP3–7.0, SSP5–8.5 (used for model projection)<br />
- Future Late-Century (2071–2100): SSP1–2.6, SSP3–7.0, SSP5–8.5 (used for model projection)<br />
<br />
These predictor/projection datasets are intended for use in SDM workflows to assess climatic suitability and potential species distributions under changing environmental conditions. They were prepared to support the OneSTOP project’s modeling of Invasive Alien Species (IAS) and integration with land-cover projections (Chen et al., 2022).
<br />
Climatic variables were obtained from CHELSA v2.1 (https://www.chelsa-climate.org/datasets), including the full set of bioclimatic variables describing temperature and precipitation regimes. In addition to the standard BIOCLIM variables, the predictor set includes BIOCLIM+ metrics such as growing-season length, growing-season precipitation, growing-season mean temperature, growing-degree days above 0 °C, 5 °C, and 10 °C, and net primary productivity, providing an expanded representation of climate conditions relevant to species’ ecological requirements.
<br />
---
<br />
Chen, G., Li, X., & Liu, X. (2022). Global land projection based on plant functional types with a 1-km resolution under socio-climatic scenarios. Scientific Data, 9, 125. https://doi.org/10.1038/s41597-022-01208-6
<br />
"

common_keywords <- c(
  "CHELSA",
  "bioclimatic variables", 
  "climate predictors",
  "environmental predictors", 
  "species distribution models",
  "SDM", 
  "climate change", 
  "CMIP6", 
  "SSP scenarios",
  "invasive alien species", "IAS",
  "OneSTOP project",
  "high-resolution climate data",
  "biodiversity modelling"
)

## 3. SCENARIO DEFINITIONS ----------------------------------------

base_dir <- "D:/invasions/projects/OneSTOP/CHELSA/v2.1_avg_gcm_scaled"

scenarios <- list(
  # list(
  #   title  = "CHELSA Bioclimatic and Environmental Predictors for Species Distribution Models (OneSTOP Project – Task 5.1): Historical Baseline (1981–2010)",
  #   folder = "CHELSA_1981-2010_v2.1"
  # ),
  # list(
  #   title  = "CHELSA Bioclimatic and Environmental Predictors for Species Distribution Models (OneSTOP Project – Task 5.1): SSP1-2.6 Scenario (2041–2070)",
  #   folder = "CHELSA_2041-2070_SSP126_avgGCM_v2.1"
  # ),
  list(
    title  = "CHELSA Bioclimatic and Environmental Predictors for Species Distribution Models (OneSTOP Project – Task 5.1): SSP3-7.0 Scenario (2041–2070)",
    folder = "CHELSA_2041-2070_SSP370_avgGCM_v2.1"
  ),
  list(
    title  = "CHELSA Bioclimatic and Environmental Predictors for Species Distribution Models (OneSTOP Project – Task 5.1): SSP5-8.5 Scenario (2041–2070)",
    folder = "CHELSA_2041-2070_SSP585_avgGCM_v2.1"
  ),
  list(
    title  = "CHELSA Bioclimatic and Environmental Predictors for Species Distribution Models (OneSTOP Project – Task 5.1): SSP1-2.6 Scenario (2071–2100)",
    folder = "CHELSA_2071-2100_SSP126_avgGCM_v2.1"
  ),
  list(
    title  = "CHELSA Bioclimatic and Environmental Predictors for Species Distribution Models (OneSTOP Project – Task 5.1): SSP3-7.0 Scenario (2071–2100)",
    folder = "CHELSA_2071-2100_SSP370_avgGCM_v2.1"
  ),
  list(
    title  = "CHELSA Bioclimatic and Environmental Predictors for Species Distribution Models (OneSTOP Project – Task 5.1): SSP5-8.5 Scenario (2071–2100)",
    folder = "CHELSA_2071-2100_SSP585_avgGCM_v2.1"
  )
)




## 4. LOOP: CREATE DRAFT + UPLOAD FILES ----------------------------

created <- list()

for (sc in scenarios) {
  message("\n--------------------------------------------------")
  message("Processing scenario: ", sc$title)
  
  # Build full path to folder
  folder_path <- file.path(base_dir, sc$folder)
  if (!dir.exists(folder_path)) {
    warning("Folder does not exist: ", folder_path, " – skipping.")
    next
  }
  
  # List all files to upload (here: everything; you can filter by pattern = \"\\.tif$\")
  files <- list.files(folder_path, full.names = TRUE)
  if (length(files) == 0) {
    warning("No files found in ", folder_path, " – skipping.")
    next
  }
  
  
  existing_deps <- zen$getDepositions(
    q            = sc$title,
    size         = 50,
    all_versions = TRUE,
    exact        = TRUE,
    quiet        = TRUE
  )
  
  dep <- NULL
  
  if (length(existing_deps) > 0) {
    existing_titles <- vapply(
      existing_deps,
      function(x) x$metadata$title,
      character(1)
    )
    
    idx <- which(existing_titles == sc$title)
    
    if (length(idx) > 0) {
      dep <- existing_deps[[idx[1]]]
      message("  -> Found existing deposition with same title (ID: ", dep$id, "). Reusing it.")
    }
  }
  
  if (is.null(dep)) {
    message("  -> No existing deposition with this title. Creating a new one.")
    
    
    ## 4.1 Create a local ZenodoRecord with metadata
    rec <- ZenodoRecord$new()
    rec$setTitle(sc$title)
    rec$setDescription(common_description)
    rec$setResourceType("dataset")
    
    rec$setPublisher("Zenodo")
    
    
    #rec$setAccessPolicyRecord("open")
    
    rec$setPublicationDate(as.character(Sys.Date()))
    
    rec$setLicense("cc0-1.0")     # change if you prefer another license
    rec$setSubjects(common_keywords)
    rec$setLanguage("eng")
    
    ## VERSION
    rec$setVersion(version_string)
    
    ## CREATED DATE
    rec$addDate(created_date, type = "created")
    
    ## IDENTIFIER (URL)
    rec$addRelatedIdentifier(
      identifier = project_url,
      scheme = "url",
      relation_type = "references"
    )
    
    ## ADD CREATOR (edit as needed)
    rec$addCreator(
      firstname   = "João",
      lastname    = "Gonçalves",
      affiliation = "Associação BIOPOLIS - Rede de Investigação em Biodiversidade e Biologia Evolutiva, Centro de Investigação em Biodiversidade e Recursos Genéticos, Universidade do Porto",
      orcid       = "0000-0002-6615-0218"
    )
    
    ## (v) FUNDING / AWARDS
    # rec$addGrant(
    #   grant = funding_award_number
    #   # award_title  = funding_award_title,
    #   # funder       = funding_funder
    # )
    
    ## 4.2 Deposit the record as a *draft* (publish = FALSE)
    dep <- zen$depositRecord(
      record      = rec,
      reserveDOI  = TRUE,
      publish     = FALSE   # IMPORTANT: keep draft only
    )
    
    message("  -> Created draft deposition with ID: ", dep$id)
    message("     DOI (concept or provisional): ", dep$getConceptDOI())
    
  } else {
    message("  -> Skipping creation: using existing draft deposition ID: ", dep$id)
  }
  
  ## 4.3 Upload all files for this scenario
  
  #files <- files[-c(1:8)]
  
  # Get list of files currently in the deposition
  existing_files <- zen$getFiles(dep$id)
  existing_fnames <- vapply(
    existing_files,
    function(x) {
      if (!is.null(x$key)) x$key else NA_character_
    },
    character(1)
  )
  
  # optionally drop NAs (in case any element has no key yet)
  existing_fnames <- existing_fnames[!is.na(existing_fnames)] 
  
  
  for (f in files) {
    fname <- basename(f)
    
    if (fname %in% existing_fnames) {
      message("    Skipping (already exists): ", fname)
    } else {
      message("    Uploading: ", fname)
      zen$uploadFile(
        path   = f,
        record = dep
      )
    }
  }
  
  # Store info for later use
  created[[sc$folder]] <- list(
    id   = dep$id,
    doi  = dep$getConceptDOI(),
    path = folder_path
  )
}

## 5. Inspect what was created ------------------------------------

created
# This list will contain IDs and DOIs of all draft depositions.
# You can copy-paste those into a small table for your records.
