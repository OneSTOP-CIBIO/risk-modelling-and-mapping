library(terra)

# -----------------------------
# User parameters
# -----------------------------
in_dir  <- "D:/invasions/projects/OneSTOP/risk-modelling-and-mapping/data/projects/onestop_v02_post/spr"
out_dir <- "D:/invasions/projects/OneSTOP/risk-modelling-and-mapping/data/projects/onestop_v02_post/analysis"

labels_future <- c("2070_ssp126", "2070_ssp370", "2070_ssp585",
                   "2100_ssp126", "2100_ssp370", "2100_ssp585")

window_km <- 10    # window size in km (e.g., 10 => ~10x10 km)
min_pairs <- 5     # min valid pairs in the window to compute correlation
method    <- "spearman"   # "pearson" (fast) or "spearman" (slow)

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# Historical richness (assumed aligned with futures)
H <- rast(file.path(in_dir, "spr_hist_bin_v02.tif"))

# Window size (cells) from km (handles non-square pixels)
res_m <- res(H)                       # meters per cell (x, y)
wx <- max(3, 2 * floor((window_km * 1000) / res_m[1] / 2) + 1)
wy <- max(3, 2 * floor((window_km * 1000) / res_m[2] / 2) + 1)
w  <- c(wy, wx)                       # rows, cols

# -----------------------------
# Fast local Pearson correlation via focal sums
# corr = cov(H,F)/sqrt(var(H)*var(F)), where
# cov(H,F) = E(HF) - E(H)E(F)
# E(·) computed with focal sums / N
# -----------------------------
local_corr_pearson <- function(H, F, w, min_pairs) {
  # valid pairs mask: 1 where both finite, NA elsewhere
  valid <- ifel(is.finite(H) & is.finite(F), 1, NA)
  
  # count valid pairs in each window
  N  <- focal(valid, w = w, fun = "sum", na.rm = TRUE)
  
  # replace N < min_pairs with NA later
  # focal sums (na.rm=TRUE ignores NA cells automatically)
  SH  <- focal(H,         w = w, fun = "sum", na.rm = TRUE)
  SF  <- focal(F,         w = w, fun = "sum", na.rm = TRUE)
  SH2 <- focal(H * H,     w = w, fun = "sum", na.rm = TRUE)
  SF2 <- focal(F * F,     w = w, fun = "sum", na.rm = TRUE)
  SHF <- focal(H * F,     w = w, fun = "sum", na.rm = TRUE)
  
  # means
  EH <- SH / N
  EF <- SF / N
  
  # covariance and variances
  covHF <- SHF - N * EH * EF
  varH  <- SH2 - N * EH * EH
  varF  <- SF2 - N * EF * EF
  
  # correlation
  rho <- covHF / sqrt(varH * varF)
  
  # enforce requirements: enough pairs, positive variance
  rho <- ifel(N >= min_pairs & varH > 0 & varF > 0, rho, NA)
  
  rho
}

# -----------------------------
# Spearman fallback (slower): uses R-level function inside focal
# -----------------------------
local_corr_spearman <- function(H, F, w, min_pairs) {
  S <- c(H, F)  # 2-layer SpatRaster
  spearman_fun <- function(v, ...) {
    m <- matrix(v, ncol = 2)
    ok <- is.finite(m[,1]) & is.finite(m[,2])
    if (sum(ok) < min_pairs) return(NA_real_)
    suppressWarnings(stats::cor(m[ok,1], m[ok,2], method = "spearman"))
  }
  focal(S, w = w, fun = spearman_fun, na.policy = "all", fillvalue = NA)
}

# -----------------------------
# Run for each future scenario
# -----------------------------
for (lab in labels_future) {
  message(sprintf("Local %s correlation for %s using %dx%d cells (~%dkm):",
                  method, lab, wx, wy, window_km))
  
  F <- rast(file.path(in_dir, paste0("spr_", lab, "_bin_v02.tif")))
  
  if (method == "pearson") {
    loc_cor <- local_corr_pearson(H, F, w, min_pairs)
  } else {
    loc_cor <- local_corr_spearman(H, F, w, min_pairs)
  }
  
  outfile <- file.path(out_dir, paste0("spr_", lab, "_loc_", method, "_", window_km, "km.tif"))
  writeRaster(loc_cor, outfile, overwrite = TRUE,
              wopt = list(datatype = "FLT4S", gdal = c("COMPRESS=LZW")))
  message(sprintf("  -> wrote %s", outfile))
}

message("✅ Done.")
