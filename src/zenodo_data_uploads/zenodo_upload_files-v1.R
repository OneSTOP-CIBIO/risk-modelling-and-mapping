
library(zen4R)

token <- Sys.getenv("ZENODO_TOKEN")

if (!nzchar(token)) {
  stop("ZENODO_TOKEN is not configured.")
}


zen <- ZenodoManager$new(token = token, logger = "INFO")

#dep <- zen$getDepositionByDOI("10.5281/zenodo.17644815")
dep <- zen$getDepositionByDOI("10.5281/zenodo.17655743")

# Vector of files you want to upload
# files <- list.files("D:/invasions/projects/OneSTOP/CHELSA/v2.1_avg_gcm_scaled/CHELSA_1981-2010_v2.1", 
#                     full.names = TRUE)

files <- list.files("D:/invasions/projects/OneSTOP/Global_LU_DataScenarios/Chen2022_Global PFT-based land projection", 
                    full.names = TRUE, pattern="\\.tif$", recursive = TRUE)

#files <- files[grepl("_recl",files)]
files <- c(
  files[2],
  files[grepl("(SSP1_RCP26|SSP3_RCP70|SSP5_RCP85).*_recl\\.tif$", files)])

# Upload each file
for (f in files) {
  zen$uploadFile(
    record = dep,
    path = f
  )
}

