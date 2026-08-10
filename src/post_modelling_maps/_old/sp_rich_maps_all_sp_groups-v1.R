
library(terra)
library(dplyr)
library(tidyverse) 
library(purrr)
library(future)
library(furrr)
library(progressr)
library(cli)

# Load the raster file paths crawler/retriever
source("./src/post_modelling_maps/rast_files_crawler-v1.R")

clean_for_path <- function(x) {
  x |>
    iconv(from = "UTF-8", to = "ASCII//TRANSLIT") |>  # remove accents/en dashes
    tolower() |>
    gsub("[^a-z0-9]+", "_", x = _) |>
    gsub("_+", "_", x = _) |>
    gsub("^_|_$", "", x = _)
}

##__________________________________________________________________________##

# Load the full paths list by map type, period, scenario
#
if(file.exists("./data/processed/best_mod_set_paths.rds")){
  best_mod_set_paths <- readRDS("./data/processed/best_mod_set_paths.rds")
}else{
  source("./src/post_modelling_analyses/check_best_model-v1.R")
}

# Load species list with taxa groups
sp_df <- readxl::read_excel(
  "./data/sel_species-standardised_organism_types_v2.xlsx")

tb <- table(sp_df$Organism_Type) |> sort()

# Select taxa groups with 3 or more species
selected_sp_groups <- c("all",names(tb[tb>=3]))

# Make a table with all scenarios
per_scn <- tibble(
  period   = c(rep(c("mid_2041_2070", "late_2071_2100"), each = 3), "hist"),
  scenario = c(rep(c("ssp126", "ssp370", "ssp585"), times = 2), "hist")
) |> mutate(ps_name = paste(period, scenario, sep = "_"))

##__________________________________________________________________________##
## Analysis options


map_type <- "binary"

out_dir <- "./data/post_model_maps/sp_rich_maps/"

##__________________________________________________________________________##


bounds_gpkg <- "./data/external/gadm/europe_selected_countries_wgs84cea_v3-1.gpkg"

bounds_vec <- vect(bounds_gpkg)

rst_templ <- rast(best_mod_set_paths[[1]]$paths_bin[1])
rst_mask <- rasterize(bounds_vec, rst_templ, field=1)


# pb <- cli::cli_progress_bar(
#   name = "Calculating maps",
#   total = total_jobs,
#   format = paste(
#     "{cli::pb_spin}",
#     "{cli::pb_name}",
#     "{cli::pb_current}/{cli::pb_total}",
#     "{cli::pb_bar}",
#     "{cli::pb_percent}",
#     #"| {cli::pb_status}",
#     "| ETA {cli::pb_eta}"
#   ),
#   clear = TRUE,
#   current = FALSE,
#   auto_terminate = FALSE
# )
# 
# ii <- 0
# for(ps in seq_len(nrow(per_scn))){
#   
#   per <- per_scn[ps,]$period
#   scn <- per_scn[ps,]$scenario
#   
#   for(spg in selected_sp_groups){
#     
#     
#     if(spg!= "all"){
#       # Changed to column Sp_name which matches the species names in file names
#       # retrieved by wiSDM-2.0
#       sp_names <- sp_df |> 
#         filter(Organism_Type == spg) |> 
#         pull(`Sp_name`)
#     }else{
#       sp_names <- "all"
#     }
#       
#     #cli_progress_update(id = pb, status = "Getting raster paths", force = TRUE)
#     
#     # Get the raster files paths for all or species sets
#     paths_df <- crawl_raster_paths(species = sp_names, 
#                                    map_type = "binary", 
#                                    period = per, 
#                                    scenario = scn,
#                                    fuzzy = TRUE,
#                                    max_fuzzy_distance = 0.1)
#     
#     #cli_progress_update(id = pb, status = "Loading rasters", force = TRUE)
#     
#     r <- rast(paths_df$path)
#     
#     out_path <- paste0(out_dir,"sprich_",clean_for_path(spg),"_",per,"_",
#                        scn,".tif")
#     
#     # cli_progress_update(id = pb, status = "Aggregating/summing rasters", 
#     #                     force = TRUE)
#     
#     if(!file.exists(out_path)){
#       r_sum <- app(r, 
#                    fun="sum", 
#                    na.rm = TRUE, 
#                    filename = out_path, 
#                    overwrite = TRUE)      
#     }
#   
#     # cli_progress_update(id = pb, status = paste("Finished species group:", spg),
#     #                     force = TRUE)
#     ii <- ii + 1L
#     
#     stopifnot(ii <= total_jobs)
#     
#     cli::cli_progress_update(
#       id = pb,
#       set = ii,
#       force = TRUE
#     )
# 
#   }
#   
# }
# stopifnot(ii == total_jobs)
# cli::cli_progress_done(id = pb)


#########################################
#total_jobs <- nrow(per_scn) * length(selected_sp_groups)


jobs <- tidyr::crossing(
  ps = seq_len(nrow(per_scn)),
  spg = selected_sp_groups
)

pb <- cli::cli_progress_bar(
  name = "Calculating maps",
  total = nrow(jobs),
  format = paste(
    "{cli::pb_spin}",
    "{cli::pb_name}",
    "{cli::pb_current}/{cli::pb_total}",
    "{cli::pb_bar}",
    "{cli::pb_percent}",
    "| {cli::pb_status}",
    "| ETA {cli::pb_eta}"
  ),
  current = FALSE,
  auto_terminate = FALSE,
  clear = FALSE
)

for (i in seq_len(nrow(jobs))) {
  
  ps  <- jobs$ps[[i]]
  spg <- jobs$spg[[i]]
  
  per <- per_scn$period[[ps]]
  scn <- per_scn$scenario[[ps]]
  
  cli::cli_progress_update(
    id = pb,
    set = i - 1L,
    inc = 0,
    status = paste("Processing", spg),
    force = TRUE
  )
  
  if(spg!= "all"){
    # Changed to column Sp_name which matches the species names in file names
    # retrieved by wiSDM-2.0
    sp_names <- sp_df |> 
      filter(Organism_Type == spg) |> 
      pull(`Sp_name`)
  }else{
    sp_names <- "all"
  }
  
  # Get the raster files paths for all or species sets
  paths_df <- crawl_raster_paths(species = sp_names, 
                                 map_type = "binary", 
                                 period = per, 
                                 scenario = scn,
                                 fuzzy = TRUE,
                                 max_fuzzy_distance = 0.1)
  
  r <- rast(paths_df$path)
  
  out_path <- paste0(out_dir,"sprich_",clean_for_path(spg),"_",per,"_",
                     scn,".tif")
  
  cli_progress_update(id = pb, status = "Aggregating/summing rasters",
                      force = TRUE, inc = 0)
  
  if(!file.exists(out_path)){

    if(scn != "hist"){
      
      r_sum <- app(r, 
                   fun="sum", 
                   na.rm = TRUE)
      
      raster::mask(r_sum, rst_mask,
                   filename = out_path, 
                   overwrite = TRUE)
    }else{
      r_sum <- app(r, 
                   fun="sum", 
                   na.rm = TRUE, 
                   filename = out_path, 
                   overwrite = TRUE)
    }
    
  }
  
  cli::cli_progress_update(
    id = pb,
    set = i - 1L,
    inc = 0,
    status = paste("Finished", spg),
    force = TRUE
  )
}

cli::cli_progress_done(id = pb)

