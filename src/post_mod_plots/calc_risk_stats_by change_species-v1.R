

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
# Inputs: file lists you already have
# -----------------------------
# fl_bin_hist, fl_bin_2070_ssp1, fl_bin_2070_ssp3, fl_bin_2070_ssp5,
# fl_bin_2100_ssp1, fl_bin_2100_ssp3, fl_bin_2100_ssp5

# -----------------------------
# User parameters
# -----------------------------
out_dir_tables <- "D:/invasions/projects/OneSTOP/risk-modelling-and-mapping/data/projects/onestop_v02_post/risk_tables"
dir.create(out_dir_tables, recursive = TRUE, showWarnings = FALSE)

# weights for the risk index (tunable)
w_s <- 1   # weight for stable suitable (1->1)
w_g <- 2   # weight for gain (0->1)
w_l <- 1   # weight for loss (1->0)
eps <- 1e-9

# Build named vectors: names = species IDs
get_species_id <- function(paths) sub(".*/([^/]+)/Rasters/.*", "\\1", paths)
name_by_species <- function(paths) { sp <- get_species_id(paths); names(paths) <- sp; paths }

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
# Helper: fast area stats for one species (hist vs fut)
# Assumes aligned 0/1 rasters in an equal-area CRS
# -----------------------------
species_stats <- function(hist_path, fut_path, species, scenario_label) {
  rH <- rast(hist_path)
  rF <- rast(fut_path)
  
  # Coerce to clean 0/1 just in case
  rH <- clamp(round(rH), 0, 1)
  rF <- clamp(round(rF), 0, 1)
  
  # Pixel area in km² (handles your 1km grid but keeps it general)
  res_m <- res(rH)
  pix_km2 <- (res_m[1] * res_m[2]) / 1e6
  
  # Dynamics (logical -> 0/1 via as.integer)
  stable_suit <- ifel(rH == 1 & rF == 1, 1, 0)
  gain        <- ifel(rH == 0 & rF == 1, 1, 0)
  loss        <- ifel(rH == 1 & rF == 0, 1, 0)
  
  # Areas (km²); NAs ignored
  A_stable <- as.numeric(global(stable_suit, "sum", na.rm = TRUE)) * pix_km2
  A_gain   <- as.numeric(global(gain,        "sum", na.rm = TRUE)) * pix_km2
  A_loss   <- as.numeric(global(loss,        "sum", na.rm = TRUE)) * pix_km2
  
  A_hist   <- as.numeric(global(rH, "sum", na.rm = TRUE)) * pix_km2
  A_future <- as.numeric(global(rF, "sum", na.rm = TRUE)) * pix_km2
  
  net_change_km2 <- A_future - A_hist
  pct_change     <- if (A_hist > 0) 100 * net_change_km2 / A_hist else NA_real_
  stability_rate <- if (A_hist > 0) A_stable / A_hist else NA_real_
  expansion_rate <- if (A_hist > 0) A_gain   / A_hist else NA_real_
  
  # Risk index (higher = more expansion while retaining area)
  risk_index <- (w_s * A_stable + w_g * A_gain - w_l * A_loss) / (A_hist + eps)
  
  data.frame(
    scenario      = scenario_label,
    species       = species,
    hist_area_km2 = A_hist,
    future_area_km2 = A_future,
    stable_km2    = A_stable,
    gain_km2      = A_gain,
    loss_km2      = A_loss,
    net_change_km2 = net_change_km2,
    pct_change     = pct_change,
    stability_rate = stability_rate,
    expansion_rate = expansion_rate,
    risk_index     = risk_index,
    stringsAsFactors = FALSE
  )
}

# -----------------------------
# Run: compute table per scenario, then save ranked CSVs
# -----------------------------
all_results <- list()

for (lab in names(scenarios)) {
  fut_named <- scenarios[[lab]]
  spp <- intersect(names(fl_hist_named), names(fut_named))
  if (length(spp) == 0) {
    message(sprintf("No matching species in scenario %s; skipping.", lab))
    next
  }
  
  message(sprintf("Scenario %s: %d species", lab, length(spp)))
  
  res_list <- vector("list", length(spp))
  i <- 1
  for (sp in spp) {
    res_list[[i]] <- species_stats(fl_hist_named[[sp]], fut_named[[sp]], sp, lab)
    i <- i + 1
  }
  df <- do.call(rbind, res_list)
  
  # Rank by risk_index (desc), with tie-breakers: net change (desc), future area (desc)
  df$rank <- rank(-df$risk_index, ties.method = "min")
  df <- df[order(-df$risk_index, -df$net_change_km2, -df$future_area_km2), ]
  
  # Save CSV per scenario
  csv_path <- file.path(out_dir_tables, paste0("risk_table_", lab, "_w", 
                                               w_s, "_g", w_g, "_l", w_l, ".csv"))
  write.csv(df, csv_path, row.names = FALSE)
  message(sprintf("  -> wrote %s", csv_path))
  
  all_results[[lab]] <- df
}

# Optionally, bind all scenarios together into one big table:
risk_all <- do.call(rbind, all_results)
csv_all  <- file.path(out_dir_tables, paste0("risk_table_ALL_w", w_s, "_g", 
                                             w_g, "_l", w_l, ".csv"))
write.csv(risk_all, csv_all, row.names = FALSE)
message(sprintf("✅ Wrote combined table: %s", csv_all))

# Print a preview (top 10 across all scenarios)
print(head(risk_all[order(risk_all$scenario, -risk_all$risk_index), ], 10))


## ------------------------------------------------------------------------------- ##

library(dplyr)
library(ggplot2)
library(forcats)
library(stringr)


risk_all_avg <- risk_all %>%
  #select(-scenario) %>%
  group_by(species) %>%
  summarise(across(where(is.numeric), \(x) mean(x, na.rm = TRUE))) %>% 
  arrange(desc(risk_index))

write_csv(risk_all_avg, "./data/projects/onestop_v02_post/risk_tables/risk_all_avg_v02.csv")


# ---- 1) Compute mean, sd, and n (for SE) ----
risk_stats <- risk_all %>%
  group_by(species) %>%
  summarise(
    across(
      where(is.numeric),
      list(
        avg = \(x) mean(x, na.rm = TRUE),
        std = \(x) sd(x, na.rm = TRUE),
        n   = \(x) sum(!is.na(x))
      ),
      .names = "{.col}_{.fn}"
    )
  ) %>%
  arrange(desc(risk_index_avg))

write_csv(risk_stats, "./data/projects/onestop_v02_post/risk_tables/risk_all_stats_v02.csv")


# ---- 2) Plot top-k species by risk_index_avg ----
k <- 30  # number of species to show

plot_df <- risk_stats %>%
  slice_head(n = k) %>%
  mutate(
    species_clean = str_replace(species, "_\\d+$", ""),       # remove numeric code
    species_clean = str_replace_all(species_clean, "_", " "), # replace underscores with spaces
    species_clean = fct_reorder(species_clean, risk_index_avg, .desc = FALSE),
    risk_index_se = risk_index_std / sqrt(risk_index_n)       # standard error
  )

g<-ggplot(plot_df, aes(x = species_clean, y = risk_index_avg, fill = risk_index_avg)) +
  geom_col() +
  geom_errorbar(
    aes(ymin = risk_index_avg - risk_index_se,
        ymax = risk_index_avg + risk_index_se),
    width = 0.25
  ) +
  scale_fill_gradient(
    name = "Risk index",
    low = "#FFBFBF",      # light brick-red
    high = "#7a0000"      # dark brick-red
  ) +
  labs(
    x = "Species",
    y = "Risk index (avg ± SE)",
    title = paste0("Top ", k, " species by risk index"),
    subtitle = "Average across all scenarios/periods"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "right",
    axis.text.y = element_text(angle = 0, hjust = 1, face = "italic"), # italic species names
    axis.title.x = element_text(face = "bold"),
    axis.title.y = element_text(face = "bold"),
    text = element_text(size = 20)
  ) +
  coord_flip()
plot(g)

ggsave(plot = g,
       filename = 
         sprintf("./data/projects/onestop_v02_post/risk_tables/risk_top_%d_barplot.png", k), 
       width = 11, height = 9, dpi = 300)



