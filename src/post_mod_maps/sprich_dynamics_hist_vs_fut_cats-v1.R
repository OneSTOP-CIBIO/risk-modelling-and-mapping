

library(terra)

# -----------------------------
# User parameters
# -----------------------------
in_dir  <- "D:/DATA/OneSTOP/onestop_v02_post/spr"
out_dir <- "D:/DATA/OneSTOP/onestop_v02_post/analysis"

labels_future <- c("2070_ssp126", "2070_ssp370", "2070_ssp585",
                   "2100_ssp126", "2100_ssp370", "2100_ssp585")

window_km <- 25  # focal window size in kilometers (e.g., 10 => ~10x10 km)

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# -----------------------------
# Load historical richness (assumed aligned with futures)
# -----------------------------
r_hist <- rast(file.path(in_dir, "spr_hist_bin_v02.tif"))

# -----------------------------
# 1) Compute thresholds (terciles) from the entire historical map (exclude zeros)
# -----------------------------
v_hist <- values(r_hist, mat = FALSE)
v_hist <- v_hist[is.finite(v_hist) & v_hist > 0]
if (length(v_hist) < 10) stop("Not enough non-zero historical values to compute terciles.")

qs <- quantile(v_hist, probs = c(1/3, 2/3), na.rm = TRUE)
q1 <- as.numeric(qs[1])
q2 <- as.numeric(qs[2])

# Reclassification matrix for 3 classes (1=Low, 2=Moderate, 3=High)
# Intervals are (from, to] in terra::classify
rcl <- matrix(c(-Inf, q1, 1,
                q1, q2, 2,
                q2,  Inf, 3),
              ncol = 3, byrow = TRUE)

# -----------------------------
# 2) Focal window in cells from km
# -----------------------------
res_m <- res(r_hist)                     # meters per cell (x, y)
wx <- max(3, 2 * floor((window_km * 1000) / res_m[1] / 2) + 1)
wy <- max(3, 2 * floor((window_km * 1000) / res_m[2] / 2) + 1)
w  <- c(wy, wx)                          # rows, cols (must be odd)

message(sprintf("Using focal window %dx%d cells (~%dkm)", wx, wy, window_km))

# -----------------------------
# 3) Focal summary for historical (use median; na.rm=TRUE)
#     Then classify with historical thresholds
# -----------------------------
hist_focal <- focal(r_hist, w = w, fun = "median", na.rm = TRUE)
hist_cls   <- classify(hist_focal, rcl)

# -----------------------------
# 4) Loop futures: focal, classify (with same thresholds), build dynamics, write RAT
# -----------------------------
# Dynamics coding: code = (hist_cls - 1) * 3 + fut_cls  ∈ {1..9}
dyn_labels <- c(
  "Low -> Low",       # 1
  "Low -> Moderate",  # 2
  "Low -> High",      # 3
  "Moderate -> Low",  # 4
  "Moderate -> Moderate", # 5
  "Moderate -> High", # 6
  "High -> Low",      # 7
  "High -> Moderate", # 8
  "High -> High"      # 9
)

# Encoder goes as:
# hist_cls	fut_cls	(hist_cls-1)*3	dyn_code	Transition class
# 1	1	0	1	Low->Low
# 1	2	0	2	Low->Moderate
# 1	3	0	3	Low->High
# 2	1	3	4	Moderate->Low
# 2	2	3	5	Moderate->Moderate
# 2	3	3	6	Moderate->High
# 3	1	6	7	High->Low
# 3	2	6	8	High->Moderate
# 3	3	6	9	High->High

# Decode as: 
# hist = ((code - 1) %/% 3) + 1
# fut = ((code - 1) %% 3) + 1)


for (lab in labels_future) {
  message(sprintf("Processing dynamics for %s ...", lab))
  
  r_fut     <- rast(file.path(in_dir, paste0("spr_", lab, "_bin_v02.tif")))
  fut_focal <- focal(r_fut, w = w, fun = "median", na.rm = TRUE)
  fut_cls   <- classify(fut_focal, rcl)
  
  # Combine into 1..9 dynamics codes
  dyn_code <- (hist_cls - 1) * 3 + fut_cls
  
  # Set cells to NA where either class is NA (no valid focal neighborhood)
  dyn_code <- mask(dyn_code, is.na(hist_cls) | is.na(fut_cls), maskvalues = TRUE, 
                   updatevalue = NA)
  
  # Make it categorical and attach RAT/levels
  dyn_cat <- as.factor(dyn_code)
  lev <- data.frame(id = 1:9, dynamics = dyn_labels)
  levels(dyn_cat) <- lev
  
  # Write to disk (RAT stored in levels)
  outfile <- file.path(out_dir, paste0("spr_", lab, "_focal_dynamics_", window_km, "km.tif"))
  writeRaster(dyn_cat, outfile, overwrite = TRUE,
              wopt = list(datatype = "INT1U", gdal = c("COMPRESS=LZW")))
  message(sprintf("   ->  wrote %s", outfile))
}

message("✅ Done: focal dynamics maps and RAT written.")
