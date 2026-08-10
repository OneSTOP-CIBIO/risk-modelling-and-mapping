

project_folder <- "./data/projects/onestop_v02/"

dirs <- list.dirs(project_folder,recursive = FALSE)

check_data <- 
data.frame(
  
  species = basename(dirs),
  
  global_model=NA,
  pres_occ=NA,
  qs_file=NA,
  
  ## Hist
  proj_hist=NA,
  proj_hist_bin=NA,
  
  ## 2070
  proj_ssp1_2070=NA,
  proj_ssp1_2070_bin=NA,
  
  proj_ssp3_2070=NA,
  proj_ssp3_2070_bin=NA,
  
  proj_ssp5_2070=NA,
  proj_ssp5_2070_bin=NA,
  
  ## 2100
  proj_ssp1_2100=NA,
  proj_ssp1_2100_bin=NA,
  
  proj_ssp3_2100=NA,
  proj_ssp3_2100_bin=NA,
  
  proj_ssp5_2100=NA,
  proj_ssp5_2100_bin=NA,
  
  ## Maps
  map_hist=NA,
  map_ssp1_2070=NA,
  map_ssp3_2070=NA,
  map_ssp5_2070=NA,
  map_ssp1_2100=NA,
  map_ssp3_2100=NA,
  map_ssp5_2100=NA
  
)

for(i in 1:length(dirs)){
  
  active_dir <- dirs[i]
  sp <- basename(active_dir)
  fl <- list.files(active_dir, full.names = TRUE, recursive = TRUE)
  
  # Acacia_dealbata_2979474.tif
  check_data[i,"global_model"] <- any(grepl(paste0("Global_model_",sp,".tif"), fl))
  check_data[i,"pres_occ"] <- any(grepl(paste0("Global_model_",sp,".tif"), fl))
  check_data[i,"qs_file"] <- any(grepl(paste0("Global_model_",sp,".tif"), fl))
  
  check_data[i,"proj_hist"] <- any(grepl(paste0(sp,"_hist_eur.tif"), fl))
  check_data[i,"proj_hist_bin"] <- any(grepl(paste0(sp,"_hist_bin_eur.tif"), fl))
  
  check_data[i,"proj_ssp1_2070"] <- any(grepl(paste0(sp,"_2041_2070_ssp126_eur.tif"), fl))
  check_data[i,"proj_ssp1_2070_bin"] <- any(grepl(paste0(sp,"_2041_2070_ssp126_bin_eur.tif"), fl))
  
  check_data[i,"proj_ssp3_2070"] <- any(grepl(paste0(sp,"_2041_2070_ssp370_eur.tif"), fl))
  check_data[i,"proj_ssp3_2070_bin"] <- any(grepl(paste0(sp,"_2041_2070_ssp370_bin_eur.tif"), fl))
  
  check_data[i,"proj_ssp5_2070"] <- any(grepl(paste0(sp,"_2041_2070_ssp585_eur.tif"), fl))
  check_data[i,"proj_ssp5_2070_bin"] <- any(grepl(paste0(sp,"_2041_2070_ssp585_bin_eur.tif"), fl))
  
  check_data[i,"proj_ssp1_2100"] <- any(grepl(paste0(sp,"_2071_2100_ssp126_eur.tif"), fl))
  check_data[i,"proj_ssp1_2100_bin"] <- any(grepl(paste0(sp,"_2071_2100_ssp126_bin_eur.tif"), fl))
  
  check_data[i,"proj_ssp3_2100"] <- any(grepl(paste0(sp,"_2071_2100_ssp370_eur.tif"), fl))
  check_data[i,"proj_ssp3_2100_bin"] <- any(grepl(paste0(sp,"_2071_2100_ssp370_bin_eur.tif"), fl))
  
  check_data[i,"proj_ssp5_2100"] <- any(grepl(paste0(sp,"_2071_2100_ssp585_eur.tif"), fl))
  check_data[i,"proj_ssp5_2100_bin"] <- any(grepl(paste0(sp,"_2071_2100_ssp585_bin_eur.tif"), fl))
  
  check_data[i,"map_hist"] <- any(grepl(paste0("__hist_bin_.png"), fl))
  
  check_data[i,"map_ssp1_2070"] <- any(grepl(paste0("__2041_2070_ssp126_bin_.png"), fl))
  check_data[i,"map_ssp3_2070"] <- any(grepl(paste0("__2041_2070_ssp370_bin_.png"), fl))
  check_data[i,"map_ssp5_2070"] <- any(grepl(paste0("__2041_2070_ssp585_bin_.png"), fl))
  
  check_data[i,"map_ssp1_2100"] <- any(grepl(paste0("__2071_2100_ssp126_bin_.png"), fl))
  check_data[i,"map_ssp3_2100"] <- any(grepl(paste0("__2071_2100_ssp370_bin_.png"), fl))
  check_data[i,"map_ssp5_2100"] <- any(grepl(paste0("__2071_2100_ssp585_bin_.png"), fl))
  
}

View(check_data)


