

library(sf)
library(terra)
library(tidyverse)
library(readr)

source("./src/post_modelling_analyses/check_best_model-v1.R")


##-------------------------------------------------------------------------##

tmpFiles(remove=TRUE)

eu_bounds <- 
  vect("./data/external/gadm/europe_selected_countries_wgs84cea_v3-1.gpkg")

eu_bounds$ID <- 1:nrow(eu_bounds)

countries_list <- as.data.frame(eu_bounds)

##-------------------------------------------------------------------------##

sample_path <- best_mod_set_paths[[1]]$paths_bin[1]

sample_rast <- rast(sample_path)

eu_bounds_rst <- rasterize(eu_bounds, sample_rast, field="ID")

##-------------------------------------------------------------------------##


# Create list of periods and scenarios to invoke them in the file list
periods <- c(rep(c("mid_2041_2070","late_2071_2100"), each=3),"hist")
scenarios <- c(rep(c("ssp126","ssp370","ssp585"), 2), "hist")

per_scn <- data.frame(
  period = periods,
  scenario = scenarios,
  ps_name = paste(periods, scenarios, sep="_")
)



pb <- cli::cli_progress_bar(
  name = "Sweeping species and predictions",
  total = 120*7,
  format = "{cli::pb_name} {cli::pb_current}/{cli::pb_total} {cli::pb_percent}",
  clear = FALSE
)

suitb_tables <- list()

k<-0
for(i in seq_along(best_mod_set_paths)){
  
  sp_name <- best_mod_set_paths[[i]]$sp_name
  sp_key <- best_mod_set_paths[[i]]$tkey
  map_used <- "binary10pct"
  
  
  for(j in 1:7){
    
    p <- per_scn[j,"period"]
    s <- per_scn[j,"scenario"]
    nm <- per_scn[j,"ps_name"]
    k<-k+1
    
    fpath <- best_mod_set_paths[[i]]$paths_list_bin[[p]][[s]]
    
    active_rast <- rast(fpath)
    
    comp <- compareGeom(active_rast, eu_bounds_rst, 
                        ext         = TRUE,
                        rowcol      = FALSE,
                        crs         = FALSE,
                        res         = FALSE,
                        stopOnError = FALSE)
    
    if(!comp){
      warning("Modifying the zonal raster support")
      eu_bounds_rst <- rasterize(eu_bounds, active_rast, field="ID")
    }
    
    
    r <- c(eu_bounds_rst, active_rast)
    names(r) <- c("zone", "class")          # names must be unique
    
    ct <- crosstab(r, long = TRUE, useNA = FALSE)
    
    tot <- tapply(ct$n, ct$zone, sum)
    n1  <- tapply(ct$n * (ct$class == 1), ct$zone, sum)
    
    suitability <- data.frame(
      sp_name   = sp_name,
      ps_name   = nm,
      id        = as.integer(names(tot)),
      n1        = as.vector(n1),
      total     = as.vector(tot),
      pct_suitb = round(100 * as.vector(n1) / as.vector(tot), 3)
    )
    
    suitb_tb_tmp <- data.frame(id=1:nrow(eu_bounds),
               countries_list) |> 
      left_join(suitability, by = "id")
    
    suitb_tables[[k]] <- suitb_tb_tmp
    
    cli::cli_progress_update(id = pb, set = k)
    
  }

}





