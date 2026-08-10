library(terra)

# -----------------------------
# User params
# -----------------------------
in_dir  <- "D:/DATA/OneSTOP/onestop_v02_post/spr"
out_dir <- "D:/DATA/OneSTOP/onestop_v02_post/analysis"
labels_future <- c("2070_ssp126", "2070_ssp370", "2070_ssp585",
                   "2100_ssp126", "2100_ssp370", "2100_ssp585")

window_km <- 10   # e.g., 10 => ~10x10 km window
min_pairs <- 5    # min valid pairs required to compute rho

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# Load historical richness (aligned with futures)
r_hist <- rast(file.path(in_dir, "spr_hist_bin_v02.tif"))

# Window size (cells) from km (works for non-square pixels)
res_m <- res(r_hist)                 # meters per cell (x, y)
wx <- max(3, 2 * floor((window_km * 1000) / res_m[1] / 2) + 1)
wy <- max(3, 2 * floor((window_km * 1000) / res_m[2] / 2) + 1)

# Custom focal function: gets a vector of length (wx*wy*2)
# (first wx*wy belong to hist, next wx*wy to fut)
spearman_fun <- function(v, ...) {
  # Split vector into two columns (hist, fut)
  m <- matrix(v, ncol = 2)
  ok <- is.finite(m[,1]) & is.finite(m[,2])
  if (sum(ok) < min_pairs) return(NA_real_)
  suppressWarnings(stats::cor(m[ok,1], m[ok,2], method = "spearman"))
}

# Run for each future scenario
for (lab in labels_future) {
  message(sprintf("Local Spearman for %s using %dx%d cells (~%dkm):",
                  lab, wx, wy, window_km))
  
  r_fut <- rast(file.path(in_dir, paste0("spr_", lab, "_bin_v02.tif")))
  
  # Stack hist + fut => 2-layer SpatRaster
  S <- c(r_hist, r_fut)
  
  # Local Spearman via terra::focal (no raster pkg needed)
  loc_cor <- focal(
    S, w = c(wy, wx),
    fun = spearman_fun,
    na.policy = "all",   # evaluate at all centers; fun handles NAs inside
    fillvalue = NA
  )
  
  outfile <- file.path(out_dir, paste0("spr_", lab, "_loc_spearman_", window_km, "km.tif"))
  writeRaster(loc_cor, outfile, overwrite = TRUE,
              wopt = list(datatype = "FLT4S", gdal = c("COMPRESS=LZW")))
  message(sprintf("  -> wrote %s", outfile))
}

message("✅ Done.")
