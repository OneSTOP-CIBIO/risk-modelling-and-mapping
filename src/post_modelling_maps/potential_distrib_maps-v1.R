
# ==============================================================================
# Automated Species Distribution Mapping for Binary SDM Outputs
# ==============================================================================
# This script creates publication-quality maps from binary species distribution
# model outputs with proper cartographic elements (legend, north arrow, scale bar)

# ==============================================================================

# Required Libraries
# ------------------------------------------------------------------------------
library(terra)          # Modern raster handling
library(sf)             # Spatial features
library(ggplot2)        # Plotting
library(ggspatial)      # North arrow and scale bar
library(rnaturalearth)  # Country boundaries
library(rnaturalearthdata)
library(cowplot)        # For enhanced layouts
library(scales)         # For better color scales

# ==============================================================================
# CONFIGURATION
# ==============================================================================

# Base directory to read species folders
#base_dir <- "D:/OneSTOP/onestop_v02"
base_dir <- "./data/projects/onestop_v02/"

# Reference grid to use
grid_ref <- read_sf("./data/external/grids/Grid_20km_Eur_OneSTOP_EPSG6933_v1.shp")


# Map settings
map_dpi <- 300                 # DPI for print quality
map_width <- 11                # Width in inches
map_height <- map_width*0.536 # Height in inches #9
color_suitable <- "#2E7D32"    # Dark green for suitable habitat
color_unsuitable <- "#EEEEEE"  # Light gray for unsuitable areas
color_border <- "#757575"      # Gray for country borders

# Projection scenarios data
proj_scenarios <- data.frame(
  
  proj_codes = c(
  "_hist_bin_",
  "_2041_2070_ssp126_bin_",
  "_2041_2070_ssp370_bin_",
  "_2041_2070_ssp585_bin_",
  "_2071_2100_ssp126_bin_",
  "_2071_2100_ssp370_bin_",
  "_2071_2100_ssp585_bin_"),
  proj_names = c(
    "Historical (1971-2024)",
    "SSP1/2.6 (2070)",
    "SSP3/7.0 (2070)",
    "SSP5/8.5 (2070)",
    "SSP1/2.6 (2100)",
    "SSP3/7.0 (2100)",
    "SSP5/8.5 (2100)"
  )
)
# ==============================================================================
# FUNCTION DEFINITIONS
# ==============================================================================

#' Load and prepare binary SDM raster
#'
#' @param file_path Path to the raster file
#' @return terra SpatRaster object
load_sdm_raster <- function(file_path) {
  cat("Loading raster data...\n")
  r <- rast(file_path)
  
  # Ensure binary values
  # if (any(!is.na(values(r)) & !values(r) %in% c(0, 1))) {
  #   warning("Raster contains non-binary values. Converting to binary (threshold = 0.5)")
  #   r <- ifel(r >= 0.5, 1, 0)
  # }
  
  cat("Raster loaded successfully\n")
  cat(sprintf("  - Dimensions: %d x %d\n", nrow(r), ncol(r)))
  cat(sprintf("  - Extent: xmin=%.2f, xmax=%.2f, ymin=%.2f, ymax=%.2f\n", 
              ext(r)[1], ext(r)[2], ext(r)[3], ext(r)[4]))
  #cat(sprintf("  - CRS: %s\n", crs(r)))
  
  return(r)
}

#' Convert raster to data frame for ggplot
#'
#' @param r terra SpatRaster
#' @return data frame with x, y, and value columns
raster_to_df <- function(r) {
  cat("Converting raster to data frame...\n")
  
  # Convert to data frame
  df <- as.data.frame(r, xy = TRUE, na.rm = FALSE)
  colnames(df) <- c("x", "y", "suitability")
  
  # Convert to factor for discrete colors
  df$suitability <- factor(df$suitability, 
                           levels = c(0, 1),
                           labels = c("Unsuitable", "Suitable"))
  
  return(df)
}

#' Load and prepare European country boundaries
#'
#' @param target_crs CRS to transform boundaries to
#' @return sf object with country boundaries
load_europe_boundaries <- function(target_crs) {
  cat("Loading European boundaries...\n")
  
  # Get world countries
  world <- ne_countries(scale = "medium", returnclass = "sf")
  
  # Define European countries (including nearby regions)
  europe_countries <- c(
    "Albania", "Andorra", "Armenia", "Austria", "Azerbaijan", "Belarus", "Belgium",
    "Bosnia and Herzegovina", "Bulgaria", "Croatia", "Cyprus", "Czechia",
    "Denmark", "Estonia", "Finland", "France", "Georgia", "Germany", "Greece",
    "Hungary", "Iceland", "Ireland", "Italy", "Latvia", "Liechtenstein",
    "Lithuania", "Luxembourg", "Malta", "Moldova", "Monaco", "Montenegro",
    "Netherlands", "North Macedonia", "Norway", "Poland", "Portugal", "Romania",
    "San Marino", "Serbia", "Slovakia", "Slovenia", "Spain", "Sweden",
    "Switzerland", "Turkey", "Ukraine", "United Kingdom", "Vatican City", "Kosovo"
  )
  
  # Filter European countries
  #europe <- world[world$name %in% europe_countries | world$continent == "Europe", ]
  europe <- world[world$name %in% europe_countries, ]
  
  # Transform to target CRS
  europe_transformed <- st_transform(europe, target_crs)
  
  cat(sprintf("  - Loaded %d country boundaries\n", nrow(europe_transformed)))
  
  return(europe_transformed)
}

#' Create species distribution map
#'
#' @param raster_df Data frame from raster_to_df
#' @param boundaries sf object with country boundaries
#' @param title Map title
#' @param extent_buffer Buffer factor for map extent (default 1.05 for 5% buffer)
#' @return ggplot object

create_sdm_map <- function(raster_df, boundaries, grid, title, 
                           subtitle, extent_buffer = 1.025) {
  cat("Creating map...\n")
  
  # Calculate extent with buffer
  x_range <- range(raster_df$x, na.rm = TRUE)
  y_range <- range(raster_df$y, na.rm = TRUE)
  
  x_center <- mean(x_range)
  y_center <- mean(y_range)
  x_span <- diff(x_range) * extent_buffer
  y_span <- diff(y_range) * extent_buffer
  
  xlim <- c(x_center - x_span/2, x_center + x_span/2)
  ylim <- c(y_center - y_span/2, y_center + y_span/2)
  
  # Create the map
  p <- ggplot() +
    # Raster layer
    geom_raster(data = raster_df, aes(x = x, y = y, fill = suitability)) +
    
    # Country boundaries
    geom_sf(data = boundaries, fill = NA, color = color_border, 
            linewidth = 0.1, alpha = 0.6) +
    
    # ADDED: intersected grid (blue, thin)
    geom_sf(data = grid, fill = NA, color = "#3C64AD",
            linewidth = 0.15, alpha = 0.8) +
    
    # Color scale
    scale_fill_manual(
      values = c("Unsuitable" = color_unsuitable, 
                 "Suitable" = color_suitable),
      na.value = "white",
      name = "Habitat Suitability",
      drop = FALSE
    ) +
    
    # Scale bar
    annotation_scale(
      location = "br",
      width_hint = 0.25,
      pad_x = unit(0.2, "in"),
      pad_y = unit(0.2, "in"),
      style = "ticks",
      line_col = "grey40",
      text_col = "grey40"
    ) +
    
    
    # North arrow
    # annotation_north_arrow(
    #   location = "tr",
    #   which_north = "true",
    #   pad_x = unit(0.1, "in"),
    #   pad_y = unit(0.1, "in"),
    #   style = north_arrow_fancy_orienteering(
    #     fill = c("grey40", "white"),
    #     line_col = "grey20"
    #   )
    # ) +
    
    # Coordinate system
    coord_sf(xlim = xlim, ylim = ylim, expand = FALSE) +
    
    # Theme
    theme_minimal(base_size = 12) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "italic", size = 16),
      plot.subtitle = element_text(hjust = 0.5, size = 10, color = "grey40"),
      legend.position = "none",            # CHANGED: hide legend
      legend.title = element_text(face = "bold", size = 11),
      legend.text = element_text(size = 10),
      legend.key.height = unit(0.8, "cm"),
      legend.key.width = unit(0.5, "cm"),
      panel.grid.major = element_line(color = "grey90", linewidth = 0.2),
      panel.grid.minor = element_blank(),
      panel.background = element_rect(fill = "white", color = NA),
      plot.background = element_rect(fill = "white", color = NA),
      axis.title = element_text(size = 10),
      plot.margin = margin(10, 10, 10, 10)
    ) +
    
    # Labels
    labs(
      title = title,
      subtitle = subtitle,
      x = "Longitude",
      y = "Latitude"
    )
  
  return(p)
}




#' Save map to multiple formats
#'
#' @param plot ggplot object
#' @param filename Base filename (without extension)
#' @param output_dir Output directory
#' @param width Width in inches
#' @param height Height in inches
#' @param dpi DPI for raster output
save_map <- function(plot, filename, output_dir, width, height, dpi) {
  cat("Saving maps...\n")
  
  # PNG output
  png_path <- file.path(output_dir, paste0(filename, ".png"))
  ggsave(
    filename = png_path,
    plot = plot,
    width = width,
    height = height,
    dpi = dpi,
    bg = "white"
  )
  cat(sprintf("  - PNG saved: %s\n", png_path))
  
  # PDF output
  # pdf_path <- file.path(output_dir, paste0(filename, ".pdf"))
  # ggsave(
  #   filename = pdf_path,
  #   plot = plot,
  #   width = width,
  #   height = height,
  #   device = cairo_pdf,
  #   bg = "white"
  # )
  # cat(sprintf("  - PDF saved: %s\n", pdf_path))
  
  #invisible(list(png = png_path, pdf = pdf_path))
  invisible(list(png = png_path))
}

#' Calculate and print summary statistics
#'
#' @param r terra SpatRaster
print_summary_stats <- function(r,sp_name=NA,proj_code=NA) {
  
  cat("\n=== Summary Statistics ===\n")
  
  vals <- values(r)
  total_cells <- length(vals)
  na_cells <- sum(is.na(vals))
  valid_cells <- total_cells - na_cells
  suitable_cells <- sum(vals == 1, na.rm = TRUE)
  unsuitable_cells <- sum(vals == 0, na.rm = TRUE)
  
  cat(sprintf("Total cells: %d\n", total_cells))
  cat(sprintf("Valid cells: %d (%.1f%%)\n", valid_cells, 
              100 * valid_cells / total_cells))
  cat(sprintf("Suitable habitat: %d cells (%.1f%% of valid)\n", 
              suitable_cells, 100 * suitable_cells / valid_cells))
  cat(sprintf("Unsuitable habitat: %d cells (%.1f%% of valid)\n", 
              unsuitable_cells, 100 * unsuitable_cells / valid_cells))
  
  # Calculate area if cell size is known
  cell_area <- prod(res(r))  # Area in map units squared
  cat(sprintf("\nSuitable area: %.2f km² (assuming meters)\n", 
              suitable_cells * cell_area / 1e6))
  
  return(
    data.frame(sp_name     = sp_name,
               proj_code   = proj_code,
               total_cells = total_cells,
               na_cells    = na_cells,
               valid_cells = valid_cells,
               suitable_cells   = suitable_cells,
               unsuitable_cells = unsuitable_cells))
  
}




# ==============================================================================
# MAIN EXECUTION
# ==============================================================================


cat("\n")
cat("==============================================================================\n")
cat("Species Distribution Mapping Pipeline\n")
cat("==============================================================================\n\n")


flist <- list.dirs(base_dir, recursive = FALSE)

if(file.exists(paste0(base_dir,"/summary_stats_by_scn.csv"))){
  all_sum_stats <- read.csv(paste0(base_dir,"/summary_stats_by_scn.csv"))
}

k <- 0
for(fname in flist){
  
  species_name <- paste(unlist(strsplit(basename(fname),"_"))[1:2],collapse=" ")
  
  fl <- list.files(paste(fname,"Rasters",sep="/"), full.names = TRUE)
  
  output_dir <- fname
  
  shp_path <- list.files(paste(output_dir,"Rasters/Global",sep="/"), 
                         pattern="\\.shp$", full.names = TRUE)
  
  # Load presence points used for model training
  sp_pts <- read_sf(shp_path)
  
  if(is.na(st_crs(sp_pts)))
    st_crs(sp_pts) <- st_crs(grid_ref)
  
  # Get intersected cells based on presence records
  grid_intersected <- st_filter(grid_ref, sp_pts)
  

  # Make one map per future scenario
  for(i in 1:nrow(proj_scenarios)){
    
    k<-k+1
    
    proj_code <- proj_scenarios[i,"proj_codes"]
    proj_name <- proj_scenarios[i,"proj_names"]
    
    input_raster_path <- fl[grepl(sprintf("%s(?!.*(aux|xml))", proj_code), fl, perl = TRUE)]
    
    # Check if input file exists
    if (!file.exists(input_raster_path)) {
      stop(sprintf("Input raster file not found: %s\n", input_raster_path))
    }
    
    output_filename <- gsub(" ", "_", tolower(species_name))
    output_filename <- paste0("sdm_map_", output_filename, "_", proj_code)
    output_path <- paste0(output_dir, "/", output_filename,".png")
    
    if(file.exists(output_path)){
      cat(":: Skipping file:",output_path,"... already exists\n\n")
      next
    }
    

    # Load raster with model projections
    sdm_raster <- try(load_sdm_raster(input_raster_path))
    
    if(inherits(sdm_raster,"try-error")){
      warning("Could not load file: ", input_raster_path)
      next
    }
    
    # Get CRS
    target_crs <- crs(sdm_raster)
    
    # Load boundaries
    europe_boundaries <- load_europe_boundaries(target_crs)
    
    # Convert raster to dataframe
    raster_df <- raster_to_df(sdm_raster)
    
    # Create map
    #map_title <- sprintf("Species distribution: %s", species_name)
    map_title <- species_name
    sdm_map <- create_sdm_map(raster_df = raster_df, 
                              boundaries = europe_boundaries, 
                              grid = grid_intersected, 
                              title = map_title, 
                              subtitle = proj_name)
    
    # Save map
    save_map(
      plot = sdm_map,
      filename = output_filename,
      output_dir = output_dir,
      width = map_width,
      height = map_height,
      dpi = map_dpi
    )
    
    # Print summary statistics
    sum_stats <- print_summary_stats(sdm_raster, species_name, proj_name)
    
    if(k==1){
      all_sum_stats <- sum_stats
    }else{
      all_sum_stats <- rbind(all_sum_stats,sum_stats)
    }
    
    write.csv(all_sum_stats, paste0(base_dir,"/summary_stats_by_scn.csv"),
              row.names = FALSE)
    
    cat("\n==============================================================================\n")
    cat("Mapping completed successfully!\n")
    cat(output_filename,"\n")
    cat("==============================================================================\n\n")
    
    
  }
  
}

