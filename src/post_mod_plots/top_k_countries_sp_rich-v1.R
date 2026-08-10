
# Packages
library(readr)       # read_csv / read_tsv
library(dplyr)
library(forcats)
library(ggplot2)
library(ggtext)      # to show <img> in axis text
library(countrycode) # to get ISO2 codes for flags

# ---------- 1) Read the table ----------
# Set your file path
csv_path <- "D:/DATA/OneSTOP/onestop_v02_post/risk_tables/sprich_hist_zonal_by_country_v2.csv"  # adjust if needed

# If your file is comma-separated, use read_csv();
# if it's tab-separated, use read_tsv().
# Uncomment the one that applies:

dat <- read_csv(csv_path, show_col_types = FALSE)
# dat <- read_tsv(csv_path, show_col_types = FALSE)

# ---------- 2) Prepare stats (SE from STD and COUNT) ----------
dat2 <- dat %>%
  mutate(
    COUNT = as.numeric(COUNT),
    MEAN  = as.numeric(MEAN),
    STD   = as.numeric(STD),
    se    = STD / sqrt(pmax(COUNT, 1))   # guard against COUNT==0
  ) %>%
  filter(COUNT > 0, is.finite(MEAN))

# ---------- 3) Pick top-k by MEAN ----------
k <- 20  # change as needed
topk <- dat2 %>%
  filter(!(NAME_0 %in% c("Monaco", "San Marino", "Liechtenstein", "Vatican City"))) %>% 
  arrange(desc(MEAN)) %>%
  slice_head(n = k)

# ---------- 4) Build flag labels for axis ----------
# Map country names to ISO2 codes (for flagcdn URLs)
topk <- topk %>%
  mutate(
    iso2 = countrycode(NAME_0, origin = "country.name", destination = "iso2c", warn = FALSE),
    # flagcdn uses lowercase ISO2 in URLs; some special cases might be NA
    flag_url = ifelse(!is.na(iso2),
                      paste0("https://flagcdn.com/w20/", tolower(iso2), ".png"),
                      NA_character_),
    # Build an HTML label with the flag image (if available) + country name
    name_label = ifelse(
      is.na(flag_url),
      NAME_0,
      paste0("<img src='", flag_url, "' width='16'/> ", NAME_0)
    ),
    # Order for plotting (ascending so coord_flip shows highest at top)
    name_label = fct_reorder(name_label, MEAN, .desc = FALSE)
  )

# ---------- 5) Plot (horizontal bars + horizontal error bars via coord_flip) ----------
g <- ggplot(topk, aes(x = name_label, y = MEAN, fill = MEAN)) +
  geom_col() +
  # error bars: mean ± SE
  geom_errorbar(aes(ymin = MEAN - se, ymax = MEAN + se), width = 0.25) +
  scale_fill_gradient(
    name = "Mean richness",
    low = "#FFBFBF",   # light brick-red
    high = "#7a0000"   # dark brick-red
  ) +
  labs(
    x = NULL,
    y = "Species richness (mean ± SE)",
    title = paste0("Top ", k, " countries by avg. sp. richness"),
    subtitle = "Hist. 1971–2024 | Zonal statistics by country"
  ) +
  theme_minimal(base_size = 21) +
  theme(
    legend.position = "right",
    axis.text.y = element_markdown(hjust = 0),  # enable flag images in labels
    axis.title.y = element_text(face = "bold"),
    axis.title.x = element_text(face = "bold"),
    plot.title = element_text(face = "bold")
  ) +
  coord_flip()

plot(g)

ggsave(plot = g,
       filename = 
         sprintf("D:/DATA/OneSTOP/onestop_v02_post/risk_tables/sprich_avg_zonal_by_country_top_%d_barplot.png", k), 
       width = 11, height = 9, dpi = 300, bg = "white")


