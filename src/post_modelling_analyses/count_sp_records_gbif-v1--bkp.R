

library(tidyverse)
library(qs)
library(sf)
library(terra)

##___________________________________________________________________________##


paths <- read_csv("./data/external/file_paths_custom_data.csv")
file_path <- paths$file_path[1]
r <- rast(file_path)
#values(r)<-0
#plot(r)

##___________________________________________________________________________##

sp_data_files <- list.files("./data/projects", 
                            pattern="_processed_occurrences.qs",
                            recursive = TRUE,
                            full.names = TRUE) |> sort()

custom_eu_boundary_path <- 
  read_sf("./data/external/gadm/europe_selected_countries_wgs84cea_v3-1.gpkg")

eu_bounds_vec <- 
  vect("./data/external/gadm/europe_selected_countries_wgs84cea_v3-1.gpkg")

countries <- custom_eu_boundary_path$NAME_0
countries <- gsub("\\ +","_",countries)

eu_idx <- c(
  4, 6, 8, 9, 10, 11, 12, 13, 14, 15,
  17, 18, 19, 21, 22, 24, 26, 27, 28,
  32, 35, 36, 37, 40, 41, 42, 43
)

eu27_countries <- countries[eu_idx]

##___________________________________________________________________________##


pb <- txtProgressBar(1,length(sp_data_files),style=3)
i=0
for(sp_data_path in sp_data_files){
  
  i = i+1
  sp_data <- qread(sp_data_path)$cleaned_1km
  sp_data_sf <- st_as_sf(sp_data,coords = c("decimalLongitude",
                                            "decimalLatitude"),
                         crs="EPSG:4326")
  
  sp_name <- sp_data[1,"species"]

  sp_data_sf_transf <- st_transform(sp_data_sf, 
                                    crs = st_crs(custom_eu_boundary_path))
  
  sp_data_transf_vect <- as(sp_data_sf_transf, "SpatVector")
  
  # Get intersection
  sp_int <- st_intersects(sp_data_sf_transf, custom_eu_boundary_path, 
                          sparse=FALSE)
  
  sp_int_vect <- terra::intersect(sp_data_transf_vect, eu_bounds_vec)
  
  # By country stats/count
  colnames(sp_int)<-countries
  by_country_count <- apply(sp_int,2,sum)
  
  sp_int_vec <- rowSums(sp_int) > 0
  
  sp_data_global_vect <- as(sp_data_sf_transf,"SpatVector")
  
  sp_data_eur <- sp_data_sf_transf[sp_int_vec, ]
  sp_data_eur_vect <- as(sp_data_eur,"SpatVector")
  
  extr_global <- terra::extract(r, sp_data_global_vect,cells=TRUE) |> 
    na.omit()
  
  extr_eur <- terra::extract(r, sp_data_eur_vect, cells=TRUE) |> 
    na.omit()
  
  by_country_counts_unique <- cbind(sp_int_vect, extr_eur) |> 
    as.data.frame() |> 
    group_by(NAME_0) |> 
    summarize(n_unique_by_country = n_distinct(cell)) |> 
    as.data.frame()
    
  
  n_total_global <- nrow(sp_data_global_vect)
  n_total_eur <- nrow(sp_data_eur_vect)
  
  n_unique_global <- length(unique(extr_global$cell))
  n_unique_eur <- length(unique(extr_eur$cell))
  
  tmp_counts <- data.frame(sp_name = sp_name, 
                           n_total_global = n_total_global,
                           n_total_eur = n_total_eur,
                           
                           n_unique_global = n_unique_global,
                           n_unique_eur = n_unique_eur,
                           
                           r_uni_global = n_unique_global/n_total_global,
                           r_uni_eur = n_unique_eur/n_total_eur,
                           r_conc_eur = n_unique_eur / n_unique_global)
  
  tmp_by_country_counts <- data.frame(sp_name = sp_name,
                                     country = names(by_country_count),
                                     n_records = by_country_count,
                                     n_records_unique = 0,
                                     row.names = names(by_country_count))
  
  tmp_by_country_counts[by_country_counts_unique$NAME_0,
                        "n_records_unique"] <- 
    by_country_counts_unique$n_unique_by_country
  
  if(i==1){
    sp_data_all <- sp_data
    sp_counts <- tmp_counts
  }else{# accumulate data
    sp_data_all <- rbind(sp_data_all,sp_data)
    sp_counts <- rbind(sp_counts, tmp_counts)
  }
  
  setTxtProgressBar(pb,i)
}

##___________________________________________________________________________##


write_csv(sp_data_all, 
          "./data/post_model_outputs/gbif_records_stats/sp_data_all.csv")
write_rds(sp_data_all, 
          "./data/post_model_outputs/gbif_records_stats/sp_data_all.rds")

##___________________________________________________________________________##

write_csv(sp_counts, 
          "./data/post_model_outputs/gbif_records_stats/sp_counts.csv")
write_rds(sp_counts, 
          "./data/post_model_outputs/gbif_records_stats/sp_counts.rds")






