library(tidyverse)
library(qs)
library(terra)

##___________________________________________________________________________##


paths <- read_csv("./data/external/file_paths_custom_data.csv")
file_path <- paths$file_path[1]
r <- terra::rast(file_path)
# values(r) <- 0
# plot(r)

##___________________________________________________________________________##

sp_data_files <- list.files(
  "E:/wiSDM_v2/risk-modelling-and-mapping",
  pattern = "_processed_occurrences.qs",
  recursive = TRUE,
  full.names = TRUE
) |>
  sort()


i=101

sp_data_path <- sp_data_files[[i]]
sp_data <- qread(sp_data_path)$cleaned_1km

length(unique(paste(sp_data$decimalLatitude,
                    sp_data$decimalLongitude,
                    sp_data$coordinateUncertaintyInMeters)))
