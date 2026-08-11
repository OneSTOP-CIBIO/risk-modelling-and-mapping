

library(sf)
library(terra)
library(tidyverse)
library(readr)
library(readxl)


extract_pft_metadata <- function(paths, strict = TRUE) {
  
  files <- basename(paths)
  
  # Expected future filename structure:
  # global_PFT_SSP1_RCP19_2020_recl_eur.tif
  scenario_pattern <- paste0(
    "^global_PFT_",
    "(SSP[1-5])_",       # SSP
    "(RCP[0-9]{2})_",    # RCP
    "([0-9]{4})_",       # Year
    "recl_eur\\.tif$"
  )
  
  # Expected baseline filename:
  # global_PFT_2015_recl_eur.tif
  baseline_pattern <- "^global_PFT_(2015)_recl_eur\\.tif$"
  
  output <- data.frame(
    path     = paths,
    filename = files,
    year     = NA_integer_,
    ssp      = NA_character_,
    rcp      = NA_character_,
    scenario = NA_character_,
    stringsAsFactors = FALSE
  )
  
  ## Extract 2015 baseline
  is_baseline <- grepl(baseline_pattern, files)
  
  output$year[is_baseline] <- 2015L
  output$ssp[is_baseline] <- NA_character_
  output$rcp[is_baseline] <- NA_character_
  output$scenario[is_baseline] <- "baseline"
  
  ## Extract future SSP–RCP scenarios
  is_scenario <- grepl(scenario_pattern, files)
  
  matches <- regexec(scenario_pattern, files[is_scenario])
  extracted <- regmatches(files[is_scenario], matches)
  
  if (length(extracted) > 0L) {
    extracted <- do.call(rbind, extracted)
    
    output$ssp[is_scenario] <- extracted[, 2]
    output$rcp[is_scenario] <- extracted[, 3]
    output$year[is_scenario] <- as.integer(extracted[, 4])
    
    output$scenario[is_scenario] <- paste(
      output$ssp[is_scenario],
      output$rcp[is_scenario],
      sep = "_"
    )
  }
  
  ## Identify filenames that do not follow either convention
  invalid <- !(is_baseline | is_scenario)
  
  if (any(invalid)) {
    message_text <- paste0(
      "The following filenames do not match the expected naming conventions:\n",
      paste0("  - ", files[invalid], collapse = "\n")
    )
    
    if (strict) {
      stop(message_text, call. = FALSE)
    } else {
      warning(message_text, call. = FALSE)
    }
  }
  
  output
}

##_________________________________________________________________________##

tmpFiles(remove=TRUE)

eu_bounds <- 
  vect("./data/external/gadm/europe_selected_countries_wgs84cea_v3-1.gpkg")

eu_bounds$ID <- 1:nrow(eu_bounds)

countries_list <- as.data.frame(eu_bounds)

colnames(countries_list) <- c("country","ID")

reclass_chen <- read_excel(
  "./data/external/Chen2022_LULC_Projections_Reclass-v1.xlsx",
                           sheet = 2)

##_________________________________________________________________________##

#paths <- read_csv("./data/external/file_paths_custom_data.csv")
paths_lc <- read_csv("./data/external/file_paths_custom_landcover_data.csv")

paths <- list.files(
  "E:/DATA/Global_LU_DataScenarios/Chen2022_Global_PFT_Projections_eur",
                    pattern="\\.tif$", full.names=TRUE)

lc_meta <- extract_pft_metadata(paths) |> 
  filter(year %in% c(2015, 2055, 2085)) |> 
  filter(scenario %in% c("baseline","SSP1_RCP26",
                         "SSP3_RCP70","SSP5_RCP85"))


##_________________________________________________________________________##

ftabs <- list()

r_eu_bounds <- rasterize(eu_bounds, r_lcover, field="ID")


##_________________________________________________________________________##



pb <- txtProgressBar(1,nrow(lc_meta),style=3)

for(i in 1:nrow(lc_meta)){
  
  period_ <- lc_meta$year[i]
  scenario <- lc_meta$scenario[i]
  
  if(period_ =="2015")
    period_n <- "current"
  if(period_==2055)
    period_n <- "2041-2070"
  if(period_==2085)
    period_n <- "2071-2100"
  
  file_path <- lc_meta$path[i]
  
  r_lcover <- rast(file_path)
  
  
  
  freq_tab <- freq(r_lcover, zones=r_eu_bounds)
  
  ftabs[[i]] <- freq_tab |> 
    select(-layer) |> 
    left_join(reclass_chen, by = c("value"="code")) |> 
    left_join(countries_list, by = c("zone"="ID")) |> 
    mutate(period = period_,
           period_name = period_n,
           scenario = scenario)
  
  
  setTxtProgressBar(pb, i)
}

##_________________________________________________________________________##


freq_tables_all <- bind_rows(ftabs)

##_________________________________________________________________________##


output_dir <- "./data/post_model_outputs/landcover_stats"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

write_csv(
  freq_tables_all,
  file.path(output_dir, "landcover_freq_tables_all-v1.csv")
)
write_rds(
  freq_tables_all,
  file.path(output_dir, "landcover_freq_tables_all-v1.rds")
)

##_________________________________________________________________________##



