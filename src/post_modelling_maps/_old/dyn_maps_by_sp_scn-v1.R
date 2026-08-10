

library(terra)


model_out_dir <- "D:/invasions/projects/OneSTOP/risk-modelling-and-mapping/data/projects/onestop_v02"

fl <- list.files(model_out_dir, pattern="\\.tif$", full.names = TRUE, recursive = TRUE)
fl <- fl[!grepl("Global_model_|Biasgrid_", fl)]
fl_bin <- fl[grepl("_bin_", fl)]
## --- Binary suitability
##
fl_bin_hist <- fl_bin[grepl("_hist_", fl_bin)]

fl_bin_2070_ssp1 <- fl_bin[grepl("_2041_2070_ssp126_", fl_bin)]
fl_bin_2070_ssp3 <- fl_bin[grepl("_2041_2070_ssp370_", fl_bin)]
fl_bin_2070_ssp5 <- fl_bin[grepl("_2041_2070_ssp585_", fl_bin)]

fl_bin_2100_ssp1 <- fl_bin[grepl("_2071_2100_ssp126_", fl_bin)]
fl_bin_2100_ssp3 <- fl_bin[grepl("_2071_2100_ssp370_", fl_bin)]
fl_bin_2100_ssp5 <- fl_bin[grepl("_2071_2100_ssp585_", fl_bin)]

# -----------------------------
# Inputs (your vectors already prepared)
# -----------------------------
# fl_bin_hist, fl_bin_2070_ssp1, fl_bin_2070_ssp3, fl_bin_2070_ssp5,
# fl_bin_2100_ssp1, fl_bin_2100_ssp3, fl_bin_2100_ssp5

# Output folder (parameter)
#out_dir <- "D:/DATA/OneSTOP/onestop_v02_post/dynamics_bin"
out_dir <- "D:/invasions/projects/OneSTOP/risk-modelling-and-mapping/data/projects/onestop_v02_post/dynamics_bin"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# -----------------------------
# Helper: extract species ID from full path
# Example: .../Acacia_dealbata_2979474/Rasters/Acacia_dealbata_2979474_...tif
# -----------------------------
get_species_id <- function(paths) sub(".*/([^/]+)/Rasters/.*", "\\1", paths)

# Build named vectors (names = species IDs)
name_by_species <- function(paths) {
  sp <- get_species_id(paths)
  names(paths) <- sp
  paths
}

fl_hist_named <- name_by_species(fl_bin_hist)

scenarios <- list(
  `2070_ssp126` = name_by_species(fl_bin_2070_ssp1),
  `2070_ssp370` = name_by_species(fl_bin_2070_ssp3),
  `2070_ssp585` = name_by_species(fl_bin_2070_ssp5),
  `2100_ssp126` = name_by_species(fl_bin_2100_ssp1),
  `2100_ssp370` = name_by_species(fl_bin_2100_ssp3),
  `2100_ssp585` = name_by_species(fl_bin_2100_ssp5)
)

# -----------------------------
# RAT / levels table
# -----------------------------
dyn_levels <- data.frame(
  id      = 1:4,
  dynamics = c("Stable unsuitable (0->0)",
               "Gain (0->1)",
               "Loss (1->0)",
               "Stable suitable (1->1)")
)

# -----------------------------
# Per-scenario, per-species loop
# -----------------------------
for (lab in names(scenarios)) {
  fut_named <- scenarios[[lab]]
  
  # species present in both historical and this scenario
  spp <- intersect(names(fl_hist_named), names(fut_named))
  if (length(spp) == 0) {
    message(sprintf("No matching species for scenario %s; skipping.", lab))
    next
  }
  
  message(sprintf("Scenario %s: %d species", lab, length(spp)))
  
  for (sp in spp) {
    hist_path <- fl_hist_named[[sp]]
    fut_path  <- fut_named[[sp]]
    
    # Load aligned 0/1 rasters
    rH <- rast(hist_path)
    rF <- rast(fut_path)
    
    # Boolean dynamics encoded as 1..4:
    # code = hist*2 + fut + 1
    # 0,0 -> 1 ; 0,1 -> 2 ; 1,0 -> 3 ; 1,1 -> 4
    dyn_code <- rH * 2 + rF + 1
    
    # Keep NA where either input is NA
    # (The arithmetic above already yields NA if any NA is present.)
    
    # Make categorical and attach RAT
    dyn_cat <- as.factor(dyn_code)
    levels(dyn_cat) <- dyn_levels
    
    # Output filename: dyn_<scenario>_<species>_bin_v02.tif
    outfile <- file.path(out_dir, paste0("dyn_", lab, "_", sp, "_bin_v02.tif"))
    
    writeRaster(
      dyn_cat, outfile, overwrite = TRUE,
      wopt = list(datatype = "INT1U", gdal = c("COMPRESS=LZW"))
    )
    
    message(sprintf("  -> wrote %s", outfile))
  }
}

message("✅ Done: per-species binary dynamics with RAT written.")
