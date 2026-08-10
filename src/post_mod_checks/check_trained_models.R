

options("rgdal_show_exportToProj4_warnings"="none")

packages <- c( "dplyr", "stringr", "here", "qs","CoordinateCleaner","terra", 
               "raster", "rnaturalearth", "rnaturalearthdata", 
               "ggplot2","tidyterra","mapview", "dismo", "sdm", "caret", 
               "viridisLite", "kableExtra","future", "future.apply",
               "randomForest","earth", "progressr", "sf", "gbm", 
               "PresenceAbsence","tictoc")

for(package in packages) {
  print(package)
  if( ! package %in% rownames(installed.packages()) ) { install.packages( package ) }
  library(package, character.only = TRUE)
}


#--------------------------------------------

source("./src/aux_funs.R")


#--------------------------------------------

projectname <- "onestop_v01"

#--------------------------------------------


#global<-qs::qread( paste0("./data/projects/",projectname,"/",projectname,"_occurrences.qs"))
taxa_info<-read.csv2(paste0("./data/projects/",projectname,"/",projectname,"_taxa_info.csv"))

taxa_info <- taxa_info %>% 
  arrange(acceptedScientificName) %>% 
  mutate(has_global_model = FALSE) %>%
  mutate(has_global_mod_map = FALSE) %>%
  mutate(model_date = as.Date(NA))

accepted_taxonkeys<-taxa_info%>%
  dplyr::pull(speciesKey)%>%
  unique()


#--------------------------------------------
# Check and list full model runs
#--------------------------------------------


for(i in 1:nrow(taxa_info)){ 
  

  species <- taxa_info$acceptedScientificName[i]
  first_two_words <- sub("^(\\w+)\\s+(\\w+).*", "\\1_\\2", species)  # Extract first two words of species name
  
  taxonkey <- accepted_taxonkeys[i]
  
  spcode <- paste(first_two_words, taxonkey, sep="_")
  
  qs_mod_file <- paste0("./data/projects/",projectname,"/",spcode)
  tif_mod_file <- paste0("./data/projects/",projectname,"/",spcode,"/Rasters/Global")
    
  has_model <- length(list.files(qs_mod_file, pattern=".qs$")) > 0
  has_preds <- length(list.files(tif_mod_file, pattern=".tif$")) > 0
  
    if(has_model){
    taxa_info[i,"model_date"] <- get_creation_date(qs_mod_file)
  } 
  
  taxa_info[i,"has_global_model"] <- has_model
  taxa_info[i,"has_global_mod_map"] <- has_preds
  
}

View(taxa_info)


#--------------------------------------------
# Delete incomplete model runs
#--------------------------------------------

for(i in 1:nrow(taxa_info)){
  
  if(taxa_info$has_global_model[i]){
    next
  } else{
    species <- taxa_info$acceptedScientificName[i]
    first_two_words <- sub("^(\\w+)\\s+(\\w+).*", "\\1_\\2", species)  # Extract first two words of species name
    
    taxonkey <- taxa_info$speciesKey[i]
    spcode <- paste(first_two_words, taxonkey, sep="_")
    
    path <- paste0("./data/projects/",projectname,"/",spcode)
    status <- unlink(path, recursive = TRUE, force = FALSE)
    print(status)
    message("Deleting: ", path)
  }

}



