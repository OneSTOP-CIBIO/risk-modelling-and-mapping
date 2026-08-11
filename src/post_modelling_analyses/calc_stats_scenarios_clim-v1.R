
library(sf)
library(terra)
library(tidyverse)
library(readr)

##_________________________________________________________________________##

tmpFiles(remove=TRUE)

eu_bounds <- 
  vect("./data/external/gadm/europe_selected_countries_wgs84cea_v3-1.gpkg")

eu_bounds$ID <- 1:nrow(eu_bounds)

countries_list <- as.data.frame(eu_bounds)

##_________________________________________________________________________##

#paths <- read_csv("./data/external/file_paths_custom_data.csv")
paths <- read_csv("./data/external/file_paths_custom_landcover_data.csv")

funs <- c("mean","sd","median","mad")

n_total <- nrow(paths) * length(funs)

##_________________________________________________________________________##

r0 <- rast(paths$file_path[1])

eu_bounds_rst <- rasterize(eu_bounds, r0, field="ID")

##_________________________________________________________________________##

pb <- txtProgressBar(1,n_total, style=3)
j<-0
for(fun_to_use in funs){
  for(i in 1:nrow(paths)){
    
    j <- j +1
    file_path <- paths$file_path[i]
    period_ <- paths$period[i]
    scenario <- paths$scenario[i]
    var_name <- paths$var_name[i]
    
    r <- rast(file_path)
    
    comp <- compareGeom(r, eu_bounds_rst, 
                        ext         = TRUE,
                        rowcol      = FALSE,
                        crs         = FALSE,
                        res         = FALSE,
                        stopOnError = FALSE)
    
    if(!comp){
      warning("Modifying the zonal raster support")
      eu_bounds_rst <- rasterize(eu_bounds, r, field="ID")
    }
    
    zon_analysis <- zonal(r, eu_bounds_rst, fun = fun_to_use, na.rm=TRUE)
    colnames(zon_analysis) <- c("ID", "value")
    
    zonal_df <-  cbind(countries_list,
                  data.frame(
                           summary_fun = fun_to_use,
                           period      = period_,
                           scenario    = scenario,
                           var_name    = var_name
                           )) |> left_join(
                             zon_analysis, by = "ID"
                           )
    
    if(j==1){
      zonal_full_df <- zonal_df
    }else{
      zonal_full_df <- rbind(zonal_full_df, zonal_df)
    }
    
    cat("[",j,"/",n_total,"]","Finished:",
        var_name,"|",
        period_,"|",
        scenario,"|",
        fun_to_use,"\n")
    setTxtProgressBar(pb, j)
    cat("\n")
  }
  # write_csv(zonal_full_df,
  #           "./data/post_model_outputs/clim_stats/clim_stats_zonal-country_full_df-v1.csv")
}


# write_csv(zonal_full_df,
#           "./data/post_model_outputs/clim_stats/clim_stats_zonal-country_full_df-v1.csv")

write_csv(zonal_full_df,
          "./data/post_model_outputs/landcover_stats/landcover_stats_zonal-country_full_df-v1.csv")


