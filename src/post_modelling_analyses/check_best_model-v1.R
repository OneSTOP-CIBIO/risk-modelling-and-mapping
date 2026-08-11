
library(terra)
library(tidyverse)
library(cli)
library(utils)
library(tools)

match_raster_geometry <- function(
    x,
    template,
    method = "near",
    fill = NA,
    filename = "",
    overwrite = FALSE
) {
  
  if (!inherits(x, "SpatRaster")) {
    stop("'x' must be a terra SpatRaster.")
  }
  
  if (!inherits(template, "SpatRaster")) {
    stop("'template' must be a terra SpatRaster.")
  }
  
  # Different CRS: project directly onto the template geometry
  if (!terra::same.crs(x, template)) {
    
    out <- terra::project(
      x,
      template,
      method = method
    )
    
  } else {
    
    same_resolution <- isTRUE(
      all.equal(
        terra::res(x),
        terra::res(template),
        tolerance = 1e-10
      )
    )
    
    same_origin <- isTRUE(
      all.equal(
        terra::origin(x),
        terra::origin(template),
        tolerance = 1e-10
      )
    )
    
    if (same_resolution && same_origin) {
      
      # Remove cells outside the template
      out <- terra::crop(
        x,
        template,
        snap = "near"
      )
      
      # Add cells where x does not cover the complete template
      out <- terra::extend(
        out,
        template,
        snap = "near",
        fill = fill
      )
      
    } else {
      
      # Different resolution/origin: transfer values to template grid
      out <- terra::resample(
        x,
        template,
        method = method
      )
    }
  }
  
  # Confirm that the resulting geometry exactly matches the template
  geometry_matches <- terra::compareGeom(
    out,
    template,
    lyrs = FALSE,
    stopOnError = FALSE,
    messages = TRUE
  )
  
  if (!geometry_matches) {
    stop("The resulting raster does not match the template geometry.")
  }
  
  if (nzchar(filename)) {
    out <- terra::writeRaster(
      out,
      filename,
      overwrite = overwrite
    )
  }
  
  out
}


mod_dir <- "E:/wiSDM_v2/risk-modelling-and-mapping/data/projects"

sp_mod_dirs <- list.dirs(mod_dir,recursive = FALSE)[1:120]

# Start variables
sp_names <- c()
tkeys <- c()
full_mod_set <- c()
best_mod_set_paths <- list()

n_dirs <- length(sp_mod_dirs)

pb <- cli::cli_progress_bar(
  name = "Sweeping folders",
  total = n_dirs,
  format = "{cli::pb_name} {cli::pb_current}/{cli::pb_total} {cli::pb_percent}",
  clear = FALSE
)

tryCatch(
  {
    for (i in seq_along(sp_mod_dirs)) {
      
      csv_taxa_info_path <- list.files(
        sp_mod_dirs[i],
        pattern = "\\.csv$",
        full.names = TRUE
      )
      
      # Prevent unexpected multiple-file input and resulting warnings
      if (length(csv_taxa_info_path) != 1L) {
        stop(
          "Expected exactly one CSV in: ", sp_mod_dirs[i],
          "; found ", length(csv_taxa_info_path)
        )
      }
      
      # Species info in the individual project folder 
      tinfo <- readr::read_delim(
        csv_taxa_info_path,
        delim = ";",
        show_col_types = FALSE,
        progress = FALSE       # Important: disable the inner progress bar
      )
      
      # Check for only one species info
      if (nrow(tinfo) != 1L) {
        stop(
          "Expected one taxon-information row in: ",
          csv_taxa_info_path,
          "; found ", nrow(tinfo)
        )
      }
      
      # Get species key and name
      sp_names[i] <- tinfo$acceptedScientificName[[1]]
      tkeys[i]    <- tinfo$acceptedTaxonKey[[1]]
      
      # List files
      mod_fl_bin <- list.files(
        sp_mod_dirs[i],
        pattern = "_binary10pct\\.tif$",
        recursive = TRUE,
        full.names = TRUE
      )
      
      mod_fl_fav <- list.files(
        sp_mod_dirs[i],
        pattern = "_ensemble\\.tif$",
        recursive = TRUE,
        full.names = TRUE
      )
      
      has_full_mod_set <- length(mod_fl_bin) == 21L
      full_mod_set[i] <- as.integer(has_full_mod_set)
      
      if (has_full_mod_set) {
        mod_bin_paths <- mod_fl_bin[8:14]
        mod_fav_paths <- mod_fl_fav[8:14]
      } else {
        mod_bin_paths <- mod_fl_bin
        mod_fav_paths <- mod_fl_fav
      }
      
      ## Change raster maps extent for overlapping 
      ## current vs future projections
      ##
      ##
      r_bin_templ <- rast(mod_bin_paths[1])
      rc_bin <-  rast(mod_bin_paths[7])
      
      r_fav_templ <- rast(mod_fav_paths[1])
      rc_fav <-  rast(mod_fav_paths[7])
      
      out_bin_rmatch_path <- paste0(
        file_path_sans_ext(mod_bin_paths[7]),"_matched_ext.tif")
      out_fav_rmatch_path <- paste0(
        file_path_sans_ext(mod_fav_paths[7]),"_matched_ext.tif")
      
      # Make a new file matching the extent of future projections
      # Binary maps
      if(!file.exists(out_bin_rmatch_path)){
        rc_bin_matched <- match_raster_geometry(
          x=rc_bin,
          template=r_bin_templ,
          method = "near",
          fill = NA,
          filename = out_bin_rmatch_path,
          overwrite = TRUE
        )
      }
      # Continuous favourability maps
      if(!file.exists(out_fav_rmatch_path)){
        rc_fav_matched <- match_raster_geometry(
          x=rc_fav,
          template=r_fav_templ,
          method = "near",
          fill = NA,
          filename = out_fav_rmatch_path,
          overwrite = TRUE
        )
      }

      # Replace the matched extent file
      mod_bin_paths[7] <- out_bin_rmatch_path
      mod_fav_paths[7] <- out_fav_rmatch_path
      
      mod_bin_path_list <- list(
        mid_2041_2070 = list(
          ssp126 = mod_bin_paths[1],
          ssp370 = mod_bin_paths[2],
          ssp585 = mod_bin_paths[3]
        ),
        late_2071_2100 = list(
          ssp126 = mod_bin_paths[4],
          ssp370 = mod_bin_paths[5],
          ssp585 = mod_bin_paths[6]
        ),
        hist = list(
          hist = mod_bin_paths[7]
        )
      )
      
      mod_fav_path_list <- list(
        mid_2041_2070 = list(
          ssp126 = mod_fav_paths[1],
          ssp370 = mod_fav_paths[2],
          ssp585 = mod_fav_paths[3]
        ),
        late_2071_2100 = list(
          ssp126 = mod_fav_paths[4],
          ssp370 = mod_fav_paths[5],
          ssp585 = mod_fav_paths[6]
        ),
        hist = list(
          hist = mod_fav_paths[7]
        )
      )
      
      best_mod_set_paths[[i]] <- list(
        sp_name        = tinfo$acceptedScientificName[[1]],
        tkey           = tinfo$acceptedTaxonKey[[1]],
        paths_bin      = mod_bin_paths,
        paths_list_bin = mod_bin_path_list,
        paths_fav      = mod_fav_paths,
        paths_list_fav = mod_fav_path_list
      )
      
      # Explicit ID prevents another progress bar becoming "current"
      cli::cli_progress_update(id = pb, set = i)
    }
  },
  finally = {
    cli::cli_progress_done(id = pb)
  }
)

write_rds(best_mod_set_paths,"./data/processed/best_mod_set_paths.rds")

check_df <- data.frame(tkeys        = tkeys,
                       sp_names     = sp_names,
                       full_mod_set = full_mod_set,
                       sp_mod_dirs  = sp_mod_dirs)


#check_df |> filter(full_mod_set==0) |> select(sp_names)
