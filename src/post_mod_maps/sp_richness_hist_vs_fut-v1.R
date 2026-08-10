
library(terra)
library(sf)
library(tidyverse)

model_out_dir <- "D:/DATA/OneSTOP/onestop_v02/"

fl <- list.files(model_out_dir, pattern="\\.tif$", full.names = TRUE, recursive = TRUE)


fl <- fl[!grepl("Global_model_|Biasgrid_", fl)]


fl_bin <- fl[grepl("_bin_", fl)]
fl_hsu <- fl[!grepl("_bin_", fl)]

#sp <- sub(".*/([^/]+)/Rasters/.*", "\\1", fl)


## --- Continuous habitat suitability 
##
fl_hsu_hist <- fl_hsu[grepl("_hist_", fl_hsu)]

fl_hsu_2070_ssp1 <- fl_hsu[grepl("_2041_2070_ssp126_", fl_hsu)]
fl_hsu_2070_ssp3 <- fl_hsu[grepl("_2041_2070_ssp370_", fl_hsu)]
fl_hsu_2070_ssp5 <- fl_hsu[grepl("_2041_2070_ssp585_", fl_hsu)]

fl_hsu_2100_ssp1 <- fl_hsu[grepl("_2071_2100_ssp126_", fl_hsu)]
fl_hsu_2100_ssp3 <- fl_hsu[grepl("_2071_2100_ssp370_", fl_hsu)]
fl_hsu_2100_ssp5 <- fl_hsu[grepl("_2071_2100_ssp585_", fl_hsu)]

## --- Binary suitability
##
fl_bin_hist <- fl_bin[grepl("_hist_", fl_bin)]

fl_bin_2070_ssp1 <- fl_bin[grepl("_2041_2070_ssp126_", fl_bin)]
fl_bin_2070_ssp3 <- fl_bin[grepl("_2041_2070_ssp370_", fl_bin)]
fl_bin_2070_ssp5 <- fl_bin[grepl("_2041_2070_ssp585_", fl_bin)]

fl_bin_2100_ssp1 <- fl_bin[grepl("_2071_2100_ssp126_", fl_bin)]
fl_bin_2100_ssp3 <- fl_bin[grepl("_2071_2100_ssp370_", fl_bin)]
fl_bin_2100_ssp5 <- fl_bin[grepl("_2071_2100_ssp585_", fl_bin)]


## ------------------------------------------------------------- ##

out_path_spr <- "D:/DATA/OneSTOP/onestop_v02_post/spr"

# Ensure output directory exists
dir.create(out_path_spr, recursive = TRUE, showWarnings = FALSE)

# Define scenarios and their file vectors
fl_by_scn <- list(
  hist         = fl_bin_hist,
  `2070_ssp126` = fl_bin_2070_ssp1,
  `2070_ssp370` = fl_bin_2070_ssp3,
  `2070_ssp585` = fl_bin_2070_ssp5,
  `2100_ssp126` = fl_bin_2100_ssp1,
  `2100_ssp370` = fl_bin_2100_ssp3,
  `2100_ssp585` = fl_bin_2100_ssp5
)

# --- Loop through each scenario ---
for (label in names(fl_by_scn)) {
  
  fl_vec <- fl_by_scn[[label]]
  
  if (length(fl_vec) == 0) {
    message(sprintf("Skipping '%s': no files found.", label))
    next
  }
  
  message(sprintf("Processing scenario: %s", label))
  
  # Load binary rasters into a SpatRaster stack
  r_stack <- rast(fl_vec)
  
  # Calculate species richness (sum of binary layers)
  spr <- sum(r_stack, na.rm = TRUE)
  
  # Output filename
  outfile <- file.path(out_path_spr, paste0("spr_", label, "_bin_v02.tif"))
  
  # Write to disk with compression and integer format
  writeRaster(
    spr, filename = outfile, overwrite = TRUE,
    wopt = list(datatype = "INT2U", gdal = c("COMPRESS=LZW"))
  )
  
  message(sprintf(" -> Wrote species richness map: %s", outfile))
}

message("✅ All species richness maps created successfully.")


