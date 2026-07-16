
#-----------------------------------------------------------------------------------
#This function calculates the number of decimal places in any given numeric value 
# eg., 15.21 has 2 decimal places, 15.2569 has 4 decimal places, 15.25690 also has 4, as 0 in the end doesn't count
#-----------------------------------------------------------------------------------
decimalplaces <- function(x) {
  if (abs(x - round(x)) > .Machine$double.eps^0.5) {
    # Remove trailing zeros and split at the decimal point
    split_result <- strsplit(sub('0+$', '', as.character(x)), ".", fixed = TRUE)[[1]]
    # Check if there are any decimals
    if (length(split_result) > 1) {
      nchar(split_result[[2]]) # Count characters in the decimal part
    } else {
      return(0) # No decimal part
    }
  } else {
    return(0) # No decimals for whole numbers
  }
}

#-----------------------------------------------------------------------------------
#Divide a numerical value by 10
#-----------------------------------------------------------------------------------
divide10<-function(x){
  value<-x/10
  return(value)
}


#-----------------------------------------------------------------------------------
# Only keep one occurrence point per grid cell of a spatRaster
#-----------------------------------------------------------------------------------
remove_duplicates <- function(occurrences, rast_template){
  
  #Initial dataset
  initial_occurrences<-nrow(occurrences)
  
  #Indicate for each occurrence point in which cell of the raster it falls
  occurrences$cell <- terra::cellFromXY( rast_template, occurrences) 
  
  #Remove occurrences that don't fall in any cell of the raster and duplicate occurrences
  occurrences <- occurrences %>%
    dplyr::filter(!is.na(cell)) %>% 
    dplyr::distinct(cell, .keep_all=TRUE) %>% 
    dplyr::select(1:2)
  
  #Print how many occurrences were removed
  print(paste(initial_occurrences - nrow(occurrences), "duplicate occurrence records removed."))
  
  return(occurrences)
  
}


#-----------------------------------------------------------------------------------
# Remove occurrences that fall in grid cells with NA values
#-----------------------------------------------------------------------------------
remove_nodata_occurrences <- function(occurrences, rast_template, crs){
  
  #Store number of initial occurrences
  initial_occurrences<-nrow(occurrences)
  
  #Remove occurrences in NA cells and convert to sf
  env<- terra::extract(rast_template, occurrences, xy = F, ID = F)
  
  occurrences <- cbind(env, occurrences)%>%
    dplyr::filter(!is.na(.[[1]])) %>%   # keep only rows where first column is not NA
    dplyr::select(c(decimalLongitude, decimalLatitude)) %>%
    sf::st_as_sf(coords = c("decimalLongitude", "decimalLatitude"), crs = crs, remove = FALSE)
  
  #Print how many occurrences were removed
  print(paste(initial_occurrences - nrow(occurrences), "occurrence records in grid cells with NAs removed."))
  
  return(occurrences)
}


#-----------------------------------------------------------------------------------
# Run k-means with fallback to fewer cluster centers
#-----------------------------------------------------------------------------------
kmeans_with_center_fallback <- function(data, center_number, step = 500,
                                        iter.max = 10, nstart = 1) {
  current_centers <- center_number
  last_error <- NULL
  
  if (length(current_centers) != 1 || is.na(current_centers) || current_centers <= 0) {
    stop("K-means clustering failed: center_number must be greater than 0.")
  }
  
  if (length(step) != 1 || is.na(step) || step <= 0) {
    stop("K-means clustering failed: step must be greater than 0.")
  }
  
  while (current_centers > 0) {
    result <- tryCatch(
      kmeans(data, centers = current_centers, iter.max = iter.max, nstart = nstart),
      error = function(e) {
        last_error <<- e
        NULL
      }
    )
    
    if (!is.null(result)) {
      if (current_centers < center_number) {
        message(
          "K-means clustering succeeded after reducing centers from ",
          center_number, " to ", current_centers, "."
        )
      }
      return(list(cluster = result$cluster, centers = current_centers))
    }
    
    next_centers <- current_centers - step
    current_centers <- if (next_centers <= 0) {
      if (current_centers > 1) 1 else 0
    } else {
      next_centers
    }
  }
  
  last_error_message <- if (!is.null(last_error)) conditionMessage(last_error) else "No k-means attempts were made."
  stop(
    "K-means clustering failed after reducing centers from ", center_number,
    " to 0. Last error: ", last_error_message
  )
}


#-----------------------------------------------------------------------------------
#Divide occurrence column with either y=0 (absences) or y=1 (presences)
#-----------------------------------------------------------------------------------
add.occ<-function(x,y){
  occ<-rep(y,nrow(x))
  cbind(x,occ)
}


#-----------------------------------------------------------------------------------
#Define the favourability transformation function
#-----------------------------------------------------------------------------------
favourability_from_prob <- function(prob_raster, prev_ratio) {
  odds <- prob_raster / (1 - prob_raster)
  fav <- odds / (prev_ratio + odds)
  fav[is.infinite(fav)] <- NA
  fav[fav < 0] <- 0
  fav[fav > 1] <- 1
  return(fav)
}


#-----------------------------------------------------------------------------------
#Function to return threshold where sens=spec from caret results 
#-----------------------------------------------------------------------------------
findThresh<-function(df){
  df<-df[c("rowIndex","obs","present")]
  df<-df %>%
    dplyr::mutate(observed= ifelse(obs == "present",1,0)) %>%
    dplyr::select(rowIndex,observed,predicted=present)
  result<-PresenceAbsence::optimal.thresholds(df,opt.methods = 2)
  return(result)
}


#-----------------------------------------------------------------------------------
#Recalculate accuracy for a given model with the threshold that has been optimized
#-----------------------------------------------------------------------------------
accuracyStats<-function(df,y){
  df<-df[c("rowIndex","obs","present")]
  df<-df %>%
    dplyr::mutate(observed= ifelse(obs == "present",1,0)) %>%
    dplyr::select(rowIndex,observed,predicted=present)
  result<-PresenceAbsence::presence.absence.accuracy(df,threshold = y,st.dev=FALSE)
  return(result)
}


#-----------------------------------------------------------------------------------
# Model predictions for a large raster in a more efficient way using parallellization 
#-----------------------------------------------------------------------------------
predict_large_raster<-function(rasterstack, model, type) {
  
  # Ensure that connections are closed even in case of an error
  on.exit({
    plan(strategy = "sequential")  # Ensure that the parallel plan is returned to sequential
    gc()  # Trigger garbage collection
    closeAllConnections()  # Close any open file connections
  }, add = TRUE)
  
  gc() #Free up memory
  
  ncores<-min(4, availableCores()-2)  #Set up number of cores
  
  if(class(rasterstack)!="SpatRaster"){
    raster_terra<-terra::rast(rasterstack)  #Convert raster to terra raster format if not already
  }else{
    raster_terra<-rasterstack
  }
  
  chunk_size <- ceiling(nrow(raster_terra) / ncores)   # Define chunk size
  
  # Create a list of row indices for each chunk
  chunk_indices <- split(seq_len(nrow(raster_terra)), ceiling(seq_along(seq_len(nrow(raster_terra))) / chunk_size))
  
  # Extract raster chunks and put raster chunks in list
  r_list<- vector(mode = "list", length = ncores)
  for (core in 1:ncores) {
    r_list[[core]]<- terra::wrap(raster_terra[min(chunk_indices[[core]]):max(chunk_indices[[core]]), ,drop=FALSE])
  } #SpatRasters need to be wrapped before sending out to different cores
  
  # Save model to disk if it’s large
  saveRDS(model, "model.rds")
  options(future.globals.maxSize = 4.5 * 1024^3)
  plan(strategy = "multisession", workers=ncores) #Set up parallel
  
  out_list <- future_lapply(r_list,  function(chunk) {
    model <- readRDS("model.rds")  # Load model from disk
    unwrapped_raster <- terra::unwrap(chunk)  # Unwrap raster for processing
    predicted_raster <- terra::predict(unwrapped_raster, model, type = type, na.rm = TRUE)
    rm(unwrapped_raster)
    terra::wrap(predicted_raster)  # Wrap the raster again
  }, future.seed = TRUE)
  
  
  plan(strategy = "sequential")   #Close parallel processing
  file.remove("model.rds")
  rm(r_list) #Remove large objects we don't need anymore
  out_list<- lapply(out_list, terra::unwrap) #unwrap chunks
  gc() # Clean up memory after processing
  model_parallel<- do.call(terra::merge, out_list)  # Merge the chunks 
  rm(out_list) #Remove large objects we don't need anymore
  gc()  #Final garbage collect
  options(future.globals.maxSize = 500 * 1024^2)  # Reset to 500 MB
  return(model_parallel)
}


#-----------------------------------------------------------------------------------
# Export PNG function
#-----------------------------------------------------------------------------------
exportPNG<-function(rst,taxonkey,taxonName,nameextension,is.diff="FALSE"){
  filename=file.path(pdfOutput,paste("be_",taxonkey, "_",nameextension,sep=""))
  png(file=filename)
  par(bty="n")#to turn off box around plot
  ifelse(is.diff=="TRUE", brks<-seq(-1, 1, by=0.25), brks <- seq(0, 1, by=0.1)) 
  nb <- length(brks)-1 
  pal <- grDevices::colorRampPalette(rev(brewer.pal(11, 'Spectral')))
  cols<-pal(nb)
  maintitle<-paste(taxonName,taxonkey,"_",nameextension, sep= " ")
  plot(rst, breaks=brks, col=cols,main=maintitle, lab.breaks=brks,axes=FALSE)
  dev.off() 
} 


#-----------------------------------------------------------------------------------
# Generate pseudoabsences
#-----------------------------------------------------------------------------------
generate_pseudoabs <- function(index = NULL,mask, alternative_mask, n, p) {
  tryf_values <- c(50,100, 150)  # tryf values to attempt in each stage
  current_raster <- mask  # Start with the initial raster layer
  
  # Attempt to generate points
  for (tryf in tryf_values) {
    # Generate random points
    suppressWarnings(pseudoabs <- as.data.frame(
      dismo::randomPoints(
        current_raster, 
        n, 
        p, 
        ext = NULL, 
        extf = 1.1, 
        excludep = TRUE, 
        prob = TRUE, 
        cellnumbers = FALSE, 
        tryf = tryf, 
        warn = 2, 
        lonlatCorrection = TRUE
      )
    )
    )
    # Check if the number of pseudoabsences reaches required amount
    if (nrow(pseudoabs) == n) {
      # If index is provided, include it in the message (only for lists)
      if (!is.null(index)) {
        message(paste0(n, " out of ", n, " pseudoabsences generated while accounting for observer bias in set ", index))
      } else {
        message(paste0(n, " out of ", n, " pseudoabsences generated while accounting for observer bias."))
      }
      return(pseudoabs)  # Return dataset if the required amount of pseudoabsences are generated
    }
  }
  
  # If unsuccessful with biasgrid ecoregions raster, switch to the full ecoregions raster and retry
  current_raster <- alternative_mask
  
  for (tryf in tryf_values) {
    pseudoabs <- as.data.frame(
      dismo::randomPoints(
        current_raster, 
        n, 
        p, 
        ext = NULL, 
        extf = 1.1, 
        excludep = TRUE, 
        prob = TRUE, 
        cellnumbers = FALSE, 
        tryf = tryf, 
        warn = 2, 
        lonlatCorrection = TRUE
      )
    )
    
    # Check if the number of rows meets the desired count
    if (nrow(pseudoabs) == n) {
      # If index is provided, include it in the warning (only for lists)
      if (!is.null(index)) {
        warning(paste0(n, " out of ", n, " pseudoabsences generated without accounting for observer bias in set ", index))
      } else {
        warning(paste0(n, " out of ", n, " pseudoabsences generated without accounting for observer bias."))
      }
      return(pseudoabs)  # Return dataset if enough pseudoabsences were generated
    }
  }
  
  # If all attempts fail, return the last generated dataframe with fewer pseudoabsences than requested
  # If index is provided, include it in the warning
  if (!is.null(index)) {
    warning(paste0("Could not generate the required number of pseudoabsences: ", n, " out of ", n, " pseudoabsences generated without accounting for observer bias in set ", index))
  } else {
    warning(paste0("Could not generate the required number of pseudoabsences: ", n, " out of ", n, " pseudoabsences generated without accounting for observer bias."))
  }
  
  return(pseudoabs)  # Return the pseudoabs data, even if incomplete
}


#-----------------------------------------------------------------------------------
# Recode factor levels to absent (0) and present(1), and set present as the reference level
#-----------------------------------------------------------------------------------
factorVars<-function(df,var){
  df[,c(var)]<-as.factor(df[,c(var)])
  levels(df[,c(var)])<-c("absent","present")
  df[,c(var)]<-relevel(df[,c(var)], ref = "present")
  return(df)
}


#-----------------------------------------------------------------------------------
#----------------Create folders when they don't exist yet---------------------------
#-----------------------------------------------------------------------------------
create_folder <- function(path, name) {
  if (!dir.exists(path)) {
    dir.create(path, recursive = TRUE)
    message(paste0("Folder '", name, "' created at path: '", path, "' 🎉"))
  }
}


#-----------------------------------------------------------------------------------
#-------------------------Export PDF new function-----------------------------------
#-----------------------------------------------------------------------------------
exportPDF <- function(predictions=NULL, period=NULL, scenario, occ_data=NULL, dataType, returnPredictions=FALSE,returnPNG=FALSE, providedPNG=NULL, exportPNG=FALSE, LabelValue=NULL, LabelName=NULL, Label2Value=NULL, Label2Name=NULL, PDF_title, PNG_folder=NULL, PDF_folder, filename){
  
  #Set scenario to "" if period is Current
  if(period=="Current") scenario<-""
  
  #Define scenario title
  scenarioTitle<- switch(paste0(period,scenario),
                         "Current" = "Current",
                         "all" = "all",
                         "2041-2070ssp126" = "2041-2070: SSP1-2.6",
                         "2041-2070ssp370" = "2041-2070: SSP3-7.0",
                         "2041-2070ssp585" = "2041-2070: SSP5-8.5",
                         "2071-2100ssp126" = "2071-2100: SSP1-2.6",
                         "2071-2100ssp370" = "2071-2100: SSP3-7.0",
                         "2071-2100ssp585" = "2071-2100: SSP5-8.5"
  )
  
  PNG_filename <- paste0(filename, ".png")
  PDF_filename <- paste0(filename, ".pdf")
  
  # Define file paths
  if(!is.null(PNG_folder)){
    plot_png_path <- file.path(PNG_folder, PNG_filename)# If not EU or Global
  }
  plot_pdf_path <- file.path(PDF_folder, PDF_filename)
  
  
  #If png is not provided, create a PNG based on the input predictions
  if(is.null(providedPNG)){
    
    #Get extent and padding in map units so projected rasters plot correctly
    exten<-as.vector(terra::ext(predictions))
    x_span <- exten[2] - exten[1]
    y_span <- exten[4] - exten[3]
    x_span <- ifelse(is.finite(x_span) && x_span > 0, x_span, max(abs(exten[1:2]), na.rm = TRUE))
    y_span <- ifelse(is.finite(y_span) && y_span > 0, y_span, max(abs(exten[3:4]), na.rm = TRUE))
    x_pad_left <- x_span * 0.12
    x_pad_right <- x_span * 0.05
    y_pad_top <- y_span * 0.04
    x_limits <- c(exten[1] - x_pad_left, exten[2] + x_pad_right)
    y_limits <- c(exten[3], exten[4] + y_pad_top)
    x_value <- x_limits[1] + 0.02 * diff(x_limits)
    
    #Settings for plot
    if (dataType == "Diff") {
      brks <- seq(-1, 1, by = 0.25)
    } else if (dataType %in% c("Suit", "Conf", "Masked_Suit", "Stdev")) {
      brks <- seq(0, 1, by = 0.2)
    }
    
    if (dataType != "Binary"){
      nb <- length(brks) - 1
      viridis_palette <- viridis::viridis(nb)
      
      #Create dummy raster with all values
      template <- predictions
      values(template) <- rep(brks, length.out = ncell(template))
      template <- mask(template, predictions)
      
      #Create plot
      suppressMessages(
        country_plot <- ggplot() + 
          geom_spatraster(data=template)+
          geom_spatraster(data = predictions) +
          scale_fill_gradientn(colors = viridis_palette, 
                               breaks = brks, 
                               labels = brks, 
                               na.value = "transparent") +
          theme_bw() +
          theme(axis.title = element_blank())+
          theme(plot.margin = unit(c(0.2,0.2,0.2,0.2), "cm"))+
          coord_sf(xlim = x_limits,
                   ylim = y_limits)
      )
      
    }else{
      #Create plot
      suppressMessages(
        country_plot <- ggplot() + 
          geom_spatraster(data = predictions) +
          scale_fill_manual(values = c("Absent" = "lightgrey", "Present" = "#085099"),
                            na.value = "transparent",
                            na.translate=FALSE)+
          theme_bw() +
          theme(axis.title = element_blank())+
          theme(plot.margin = unit(c(0.2,0.2,0.2,0.2), "cm"))+
          coord_sf(xlim = x_limits,
                   ylim = y_limits)
      )
    }
    # Define text label, fill label, and hjust based on dataType
    text_label <- ifelse(dataType=="Diff",paste(scenarioTitle, "- current"), scenarioTitle)
    
    fill_label <- switch(dataType,
                         "Suit" = "Suitability",
                         "Diff" = "Suitability difference",
                         "Conf" = "Confidence",
                         "Masked_Suit" = "Suitability",
                         "Binary" = "Suitability",
                         "Stdev" = "Standard deviation")
    
    # Update the plot
    country_plot <- country_plot +
      labs(fill = fill_label) +
      annotate("text",
               x = x_value, y = Inf,       # Position at top-right
               label = text_label,     # Text to display
               hjust = 0,
               vjust = 2.4,            # Adjust text alignment to the right and above
               size = 4.8,
               color = "#636363",
               fontface = "bold")+
      theme(aspect.ratio=NULL)
    
    if(!is.null(occ_data)){
      crs_value<-st_crs(occ_data)
      #Only show occurrences that fall within raster cells
      occ_data<- terra::extract(predictions, occ_data, xy = T, ID=F)%>%
        dplyr::filter(!is.na(.[, 1]))%>% #Keep rows that do not have any NA values in value column
        dplyr::select(c(x,y))%>%
        dplyr::rename(decimalLongitude=x,
                      decimalLatitude=y)%>%
        st_as_sf(coords=c("decimalLongitude", "decimalLatitude"), crs=crs_value)
      
      suppressMessages(
        country_plot<-country_plot +
          geom_sf(data = occ_data, color = "black", fill = "red", 
                  size = 1.5, shape = 21)+
          coord_sf(xlim = x_limits,
                   ylim = y_limits)
      )
    }
    
    if(!is.null(LabelValue)){
      
      assertthat::assert_that(!is.null(LabelName), msg = "LabelValue is provided but LabelName is not.")
      
      country_plot<-country_plot +
        annotate("text",
                 x = x_value, y = Inf,       # Position at top-right
                 label = paste(LabelName,"=",LabelValue),     # Text to display
                 hjust = 0,
                 vjust = 4.5,            # Adjust text alignment to the right and above
                 size = 4.4,
                 color = "#919191",
                 fontface = "bold")
    }
    
    if(!is.null(Label2Value)){
      
      assertthat::assert_that(!is.null(Label2Name), msg = "Label2Value is provided but Label2Name is not.")
      
      country_plot<-country_plot +
        annotate("text",
                 x = x_value, y = Inf,       # Position at top-right
                 label = paste(Label2Name,"=",Label2Value),     # Text to display
                 hjust = 0,
                 vjust = 6.5,            # Adjust text alignment to the right and above
                 size = 4.4,
                 color = "#919191",
                 fontface = "bold")
    }
    
    #Create an empty plot to fill PDF
    empty_plot <- ggplot() + 
      theme_void() + 
      theme(plot.background = element_blank()) 
    
    #Create final plot
    plot_final<-country_plot 
    
    # Save plot temporarily as a PNG file
    ggplot2::ggsave(filename = PNG_filename, plot = plot_final, 
                    device = "png", width =7.7 , height = 6.94, units = "in", dpi= 300, path= PDF_folder)
    
  }else{
    # Save plot temporarily as a PNG file
    ggplot2::ggsave(filename = PNG_filename, plot =providedPNG, 
                    device = "png", width =7.7 , height = 6.94, units = "in", dpi= 300, path= PDF_folder)
  }
  
  # Read the PNG image back in
  img <- magick::image_read(here::here(PDF_folder, PNG_filename))
  
  # Start a PDF device for output (A4 portrait in inches)
  pdf(plot_pdf_path, width = 8.27, height = 11.69)
  
  grid.newpage()
  
  # Add title at the top of the PDF
  grid.text(
    label = PDF_title,
    x = 0.5, y = 0.95, just = "center", gp = gpar(fontsize = 14, fontface = "bold")
  )
  
  # Insert the PNG centered on the page
  # - x = 0.5 centers horizontally
  # - y = 0.5 centers vertically
  # - just = "center" anchors the image by its center point
  # - width/height use npc units so the image scales to the page
  grid::grid.raster(
    img,
    x = unit(0.5, "npc"),
    y = unit(0.627, "npc"),  # Put plot under title
    width = unit(0.95, "npc"),
    height = unit(0.6, "npc"),  # roughly matches PNG aspect ratio
    just = "center"
  )
  
  # Close the PDF device (single close)
  dev.off()
  
  # Print confirmation
  print(paste(PDF_filename," has been created.", sep=""))
  
  # Remove the temporary PNG file
  file.remove(here::here(PDF_folder, PNG_filename))
  gc()
  
  #Store PNG file in PNG folder if exportPNG is TRUE
  if(exportPNG){
    ggplot2::ggsave(filename = PNG_filename, plot = country_plot, 
                    device = "png", width = 7.7, height = 6.94, units = "in", dpi = 300,
                    path = PNG_folder)
    print(paste(PNG_filename," has been created.", sep="")) 
  }
  
  
  #Store plots or models
  if (returnPredictions & returnPNG) {
    
    return(list("model" = predictions,
                "png"=country_plot,
                "scenario"=scenario))
  }
  
  if (returnPredictions & !returnPNG) {
    # return(setNames(list("model" = predictions), scenario))
    return(list("model" = predictions,
                "scenario"=scenario))
  }
  
  if (!returnPredictions & returnPNG) {
    # return(setNames(list("model" = predictions), scenario))
    return(list("png" = country_plot,
                "scenario"=scenario))
  }
}


#-----------------------------------------------------------------------------------
#------------------------- Standardize residuals -----------------------------------
#-----------------------------------------------------------------------------------
stdres<-function(obs.numeric, yhat){
  num<-obs.numeric-yhat #Obtain residuals (the difference between observed and predicted values)
  denom<-sqrt(yhat * (1 - yhat) + 1e-10) #Approximates the residual variance in logistic regression,  + 1e-10 is added in case predicted values are 0
  return(num/denom)#Standardize
}


#-----------------------------------------------------------------------------------
#--------Return number of elements that are equal or less than a threshold ---------
#-----------------------------------------------------------------------------------

GetLength <- function(x, y) {
  sum(x <= y)
}

#-----------------------------------------------------------------------------------
#--------- classify based on probabilities compared to a confidence level ----------
#-----------------------------------------------------------------------------------
CPconf<-function(pA,pB,confidence){
  if(pA > confidence && pB< confidence){
    predClass<-"classA"
  }else if(pA < confidence && pB> confidence){
    predClass<-"classB"
  }else if(pA< confidence && pB< confidence){
    predClass<-"noClass"
  }else{
    predClass<-"bothClasses"
  }
  return(predClass)
}


#-----------------------------------------------------------------------------------
#---------------------- calculate confidence of each prediction --------------------
#-----------------------------------------------------------------------------------
get.confidence<-function(pvalA,pvalB){
  secondHighest<-ifelse(pvalA>pvalB,pvalB,pvalA)
  conf<-(1-secondHighest)
  return(conf)
}



#-----------------------------------------------------------------------------------
#-------------- Return presence/absence based on values a and b --------------------
#-----------------------------------------------------------------------------------
forcedCp<-function(pvalA,pvalB){
  ifelse(pvalA>pvalB,"presence","absence")
}


#-----------------------------------------------------------------------------------
#--- Extract probability of presence and absence from prediction raster ------------
#-----------------------------------------------------------------------------------
extractVals<-function(predras){
  vals <-  as.numeric(terra::values(predras))
  vals[is.nan(vals)] <- NA
  coord <-  terra::xyFromCell(predras,1:terra::ncell(predras))
  raster_fitted <- cbind(coord,vals)
  raster_fitted.df<-as.data.frame(raster_fitted)
  raster_fitted.df1<-na.omit(raster_fitted.df)
  raster_fitted.df1$presence<-raster_fitted.df1$vals
  raster_fitted.df1$absence<- (1-raster_fitted.df1$presence)
  return(raster_fitted.df1)
}


#-----------------------------------------------------------------------------------
#-------------- Class conformal prediction --------------------
#-----------------------------------------------------------------------------------
classConformalPrediction<-function(x,y){
  #Extract model results
  ens_results <- x
  ens_calib<-ens_results$ens_model$pred
  
  # Filter and extract calibration data for presence and absence
  calibPresence<-ens_calib %>%
    dplyr::filter(obs=='present')%>%
    dplyr::select(present)
  calibPresence<-unname(unlist(calibPresence[c("present")]))
  
  calibAbsence<-ens_calib %>%
    dplyr::filter(obs=='absent')%>%
    dplyr::select(absent)
  calibAbsence<-unname(unlist(calibAbsence[c("absent")]))
  
  #Extract predicted values
  predicted.values<-extractVals(y)
  testPresence<-predicted.values$presence
  testAbsence<-predicted.values$absence
  
  #derive p.Values for class A
  smallrA<-lapply(testPresence,function(x) GetLength(calibPresence,x))#For each value in testPresence, you calculate the number of values in calibPresence that are smaller or equal to the testPresence value
  smallrA_1<- unlist (smallrA)+1 #Create a vector of resulting values and add 1 to each value
  nCalibSet<-length(calibPresence)+1
  pvalA<-smallrA_1+1/nCalibSet
  
  # derive p.Values for Class B
  smallrB<-lapply(testAbsence,function(x) GetLength(calibAbsence,x))
  smallrB_1<- unlist (smallrB)+1
  nCalibSetB<-length(calibAbsence)
  pvalB<-smallrB_1/nCalibSetB
  
  pvalsdf<-as.data.frame(cbind(pvalA,pvalB,0.20))
  #raster_cp_20<-mapply(CPconf,pvalsdf$pvalA,pvalsdf$pvalB,pvalsdf[3])
  #table(raster_cp_20)
  
  pvalsdf$conf<-get.confidence(pvalsdf$pvalA,pvalsdf$pvalB)
  pvalsdf_1<-cbind(pvalsdf,predicted.values)
  
  return(pvalsdf_1)  
}


#-----------------------------------------------------------------------------------
#--------------------------- Create confidence maps --------------------------------
#-----------------------------------------------------------------------------------
confidenceMaps<-function(x,original_raster,taxonName, taxonNameTitle, nameExtension, taxonKey ,scenario, regionName, scenarioTitle, dataType, folder, GlobalModel=FALSE, resampling_rast=NULL, country_sf=NULL){
  # Create a SpatVector from the data.xyz
  data.xyz <- x[c("x","y","conf")]
  points <- terra::vect(data.xyz, geom = c("x", "y"), crs = terra::crs(original_raster))
  
  # Rasterize the points using the original SpatRaster as a template
  rst <- terra::rasterize(points, original_raster, field = "conf")
  
  #If global model is used, resample map
  if(GlobalModel){
    rst_to_export<- rst%>%
      terra::project(terra::crs(country_sf))%>%
      terra::resample(resampling_rast, method="bilinear") 
  }else{
    
    rst_to_export<-rst
  }
  
  #Export raster
  raster_file<-paste(taxonName, "_", taxonKey, "_", scenario, "_confidence_", regionName, ".tif", sep="")
  terra::writeRaster(rst_to_export,
                     filename=file.path(folder, raster_file),
                     overwrite=TRUE)
  #Print
  print(paste(raster_file," has been created.", sep=""))
  
  exportPDF(predictions=rst_to_export,
            taxonName,
            nameExtension, 
            taxonNameTitle,
            taxonKey=taxonKey, 
            scenario, 
            regionName,
            returnPredictions=FALSE,
            returnPNG=FALSE, 
            dataType="Conf")
  return(rst)
  
}


#----------------------------------------------------------
#---------------- Assess response curves-------------------
#----------------------------------------------------------
#evaluate predictions while varying only the selected variable (x) and keeping all other variables at their observed values
partial_gbm<-function(x){
  m.gbm<-pdp::partial(bestModel$models$gbm$finalModel,pred.var=paste(x),train = bestModel.train,type="classification",
                      prob=TRUE,n.trees= bestModel$models$gbm$finalModel$n.trees, which.class = 1,grid.resolution=nrow(bestModel.train))
}

partial_glm<-function(x){
  m.glm<-pdp::partial(bestModel$models$glm$finalModel,pred.var=paste(x),train = bestModel.train,type="classification",
                      prob=TRUE,which.class = 1,grid.resolution=nrow(bestModel.train))
}

partial_mars<-function(x){
  m.mars<-pdp::partial(bestModel$models$earth$finalModel,pred.var=paste(x),train = bestModel.train,type="classification",
                       prob=TRUE,which.class = 2,grid.resolution=nrow(bestModel.train)) # class=2 because in earth pkg, absense is the first class
}

partial_rf<-function(x){
  pdp::partial(bestModel$models$rf$finalModel,pred.var=paste(x),train = bestModel.train,type="classification",
               prob=TRUE,which.class = 1,grid.resolution=nrow(bestModel.train))
}

#----------------------------------------------------------
#---------------- Plot response curves-------------------
#----------------------------------------------------------

responseCurves<-function(x,y) {
  colors <- c("GLM" = "gray", "GBM"="red","RF"="blueviolet","MARS"= "hotpink") 
  ggplot(all_dfs,(aes(x=.data[[x]],y=.data[[y]]))) +
    geom_line(aes(color = data), size =1.2, position=position_dodge(width=0.2))+
    theme_bw()+
    labs(y="Partial probability", x= gsub("//..*","",x),color="Legend") +
    scale_color_manual(values = colors)
}  



#----------------------------------------------------------
#------------- Evaluate model predictions------------------
#----------------------------------------------------------
eu_eval<-function (ras,y){
  indep.bil<-terra::extract(ras,y,method="bilinear")
  indep.bil.df<-as.data.frame(indep.bil)
  indep.bil.df<-indep.bil.df %>%
    dplyr::mutate(predicted= ifelse(indep.bil >= 0.5,"present","absent")) 
  indep.bil.df$observed<-rep("present",nrow(indep.bil.df))
  indep.bil.df$predicted<-as.factor(indep.bil.df$predicted)
  indep.bil.df$observed<-as.factor(indep.bil.df$observed)
  xtab<-table(indep.bil.df$predicted,indep.bil.df$observed)
  return(xtab)
}

#----------------------------------------------------------
#-------------   update_files_logic   ----------------------
#----------------------------------------------------------
#' @param dest_file output file to test existance of
#' @param update_files whether to ask, or ignore existance of input environmental layers
update_files_logic <- function(dest_file,
                               dest_folder = NULL,
                               update_files) {
  
  # normalize update_files
  if (is.factor(update_files)) update_files <- as.character(update_files)
  update_files <- trimws(tolower(update_files))  # remove spaces, make lowercase
  
  if (!update_files %in% c("yes","no","ask")) {
    stop("update_files must be 'yes', 'no', or 'ask'")
  }
  
  # ---------- SINGLE FILE ----------
  if (!is.data.frame(dest_file)) {
    
    if (update_files == "yes") {
      
      update_files_final <- TRUE
      
    } else if (update_files == "no") {
      
      update_files_final <- !file.exists(dest_file)
      
    } else {  # "ask"
      
      if (file.exists(dest_file)) {
        update_files_final <- askYesNo(
          paste("Download\n", basename(dest_file), "\n again?")
        )
      } else {
        update_files_final <- TRUE
      }
    }
    
    # ---------- MULTIPLE FILES ----------
  } else {
    
    if (update_files == "yes") {
      
      dest_file$update_file <- TRUE
      
    } else if (update_files == "no") {
      
      dest_file$update_file <- !file.exists(
        file.path(dest_folder, dest_file$file)
      )
      
    } else {  # "ask"
      
      dest_file$update_file <- FALSE
      
      for (i in seq_len(nrow(dest_file))) {
        f <- file.path(dest_folder, dest_file$file[i])
        
        if (!file.exists(f)) {
          dest_file$update_file[i] <- TRUE
        } else {
          dest_file$update_file[i] <- askYesNo(
            paste("Download\n", basename(f), "\n again?")
          )
        }
      }
    }
    
    update_files_final <- dplyr::filter(dest_file, update_file)
  }
  
  return(update_files_final)
}

#----------------------------------------------------------
#-------------   safe_download_zenodo   -------------------
#----------------------------------------------------------

#' a wrapper for zen4R::download_zenodo to delete the failed file and thus trigger 
#' a redownload with update_files == FALSE

safe_download_zenodo <- function(doi, path, files, timeout = 600, quiet = FALSE) {
  tryCatch(
    {
      zen4R::download_zenodo(
        doi    = doi,
        path   = path,
        files  = files,
        timeout = timeout,
        quiet  = quiet
      )
    },
    error = function(e) {
      msg <- conditionMessage(e)
      
      # Try to extract the dest_file name from the error message
      # Adapt the regex to match the actual message format
      m <- regexpr("dest_file ['\"]?([^'\" ]+)['\"]?", msg)
      if (m[1] != -1) {
        fname <- regmatches(msg, m)
        # fname now contains something like \"dest_file 'myfile.ext'\"
        # extract just the filename part
        fname_only <- sub(".*dest_file ['\"]?([^'\" ]+)['\"]?.*", "\\1", fname)
        
        file_to_remove <- file.path(path, fname_only)
        if (file.exists(file_to_remove)) {
          unlink(file_to_remove)
        }
      }
      
      # Re-display the original error
      stop(e)
    }
  )
}


#----------------------------------------------------------
#---- Check, and if necessary, redownloaded tif files -----
#----------------------------------------------------------
read_or_redownload <- function(file, folder, doi, max_attempts = 3) {
  
  file_path <- file.path(folder, file)
  attempt <- 1
  
  while (attempt <= max_attempts) {
    
    r <- tryCatch({
      r <- terra::rast(file_path)
      terra::ncell(r)
      
    }, error = function(e) {
      NULL
    })
    
    if (!is.null(r)) {
      return(r)  # success
    }
    
    message(paste("Corrupt file detected:", file, "- redownloading (attempt", attempt, ")"))
    
    # Remove corrupt file if it exists
    if (file.exists(file_path)) {
      file.remove(file_path)
    }
    
    # Redownload
    zen4R::download_zenodo(
      doi = doi,
      path = folder,
      files = file,
      timeout = 600,
      quiet = FALSE
    )
    
    attempt <- attempt + 1
  }
  
  stop(paste("Failed to obtain valid raster after", max_attempts, "attempts:", file))
}



#-----------------------------------------------------------------
#--Resolve a user-provided path against a base directory-----------
#-----------------------------------------------------------------
resolve_input_path <- function(path, base_dir = getwd()) {
  if (is.na(path)) {
    return(NA_character_)
  }
  
  path <- trimws(path)
  
  if (!nzchar(path)) {
    return(path)
  }
  
  is_absolute <- grepl("^(?:[A-Za-z]:[\\\\/]|/|\\\\\\\\)", path)
  resolved <- if (is_absolute) path else file.path(base_dir, path)
  
  normalizePath(resolved, winslash = "/", mustWork = FALSE)
}


#-----------------------------------------------------------------
#--Return an sf CRS from a supported reference object-------------
#-----------------------------------------------------------------
get_reference_crs <- function(reference) {
  if (inherits(reference, "SpatRaster") || inherits(reference, "SpatVector")) {
    reference_crs <- terra::crs(reference)
    if (!nzchar(reference_crs)) {
      stop("The reference raster/vector does not have a defined CRS.", call. = FALSE)
    }
    return(sf::st_crs(reference_crs))
  }
  
  if (inherits(reference, "sf") || inherits(reference, "sfc")) {
    reference_crs <- sf::st_crs(reference)
    if (is.na(reference_crs)) {
      stop("The reference vector does not have a defined CRS.", call. = FALSE)
    }
    return(reference_crs)
  }
  
  if (inherits(reference, "crs")) {
    if (is.na(reference)) {
      stop("The reference CRS is not defined.", call. = FALSE)
    }
    return(reference)
  }
  
  if ((is.character(reference) || is.numeric(reference)) && length(reference) == 1) {
    reference_crs <- sf::st_crs(reference)
    if (is.na(reference_crs)) {
      stop("The reference CRS could not be parsed.", call. = FALSE)
    }
    return(reference_crs)
  }
  
  stop("Unsupported reference object supplied for CRS matching.", call. = FALSE)
}


#-----------------------------------------------------------------
#--Create the default EU boundary from the habitat raster---------
#-----------------------------------------------------------------
create_default_eu_boundary <- function(habitat_boundary_raster = file.path("data", "external", "habitat", "Agriculture.tif")) {
  habitat_boundary_raster <- resolve_input_path(habitat_boundary_raster)
  
  if (!file.exists(habitat_boundary_raster)) {
    stop("The habitat raster used to derive the default EU boundary does not exist: ", habitat_boundary_raster, call. = FALSE)
  }
  
  euboundary <- terra::rast(habitat_boundary_raster)
  euboundary <- (euboundary * 0) + 1
  euboundary <- terra::as.polygons(euboundary, dissolve = TRUE)
  euboundary <- sf::st_as_sf(euboundary)
  
  if (!all(sf::st_is_valid(euboundary))) {
    euboundary <- suppressWarnings(sf::st_make_valid(euboundary))
  }
  
  euboundary
}


#-----------------------------------------------------------------
#--Load the climate masking layer while preserving legacy default--
#-----------------------------------------------------------------
load_climate_eu_boundary <- function(custom_path = NULL,
                                     reference,
                                     habitat_boundary_raster = file.path("data", "external", "habitat", "Agriculture.tif"),
                                     legacy_extent = c(-38, 50, 24.29152732065, 72.66652712715)) {
  if (missing(reference) || is.null(reference)) {
    stop("A reference raster must be supplied to derive the climate EU boundary.", call. = FALSE)
  }
  
  if (is.null(custom_path)) {
    habitat_boundary_raster <- resolve_input_path(habitat_boundary_raster)
    
    if (!file.exists(habitat_boundary_raster)) {
      stop("The habitat raster used to derive the default climate EU boundary does not exist: ", habitat_boundary_raster, call. = FALSE)
    }
    
    return(
      terra::rast(habitat_boundary_raster) %>%
        terra::project(reference) %>%
        terra::crop(terra::ext(legacy_extent[1], legacy_extent[2], legacy_extent[3], legacy_extent[4]))
    )
  }
  
  load_eu_boundary(
    custom_path = custom_path,
    reference = reference,
    habitat_boundary_raster = habitat_boundary_raster
  ) %>%
    terra::vect()
}


#-----------------------------------------------------------------
#--Load the active EU boundary and align it to a reference CRS----
#-----------------------------------------------------------------
load_eu_boundary <- function(custom_path = NULL,
                             reference = NULL,
                             default_path = file.path("data", "external", "GIS", "Europe", "EUboundary.shp"),
                             habitat_boundary_raster = file.path("data", "external", "habitat", "Agriculture.tif")) {
  if (is.null(custom_path)) {
    default_path <- resolve_input_path(default_path)
    
    if (file.exists(default_path)) {
      euboundary <- sf::st_read(default_path, quiet = TRUE)
    } else {
      euboundary <- create_default_eu_boundary(habitat_boundary_raster)
    }
  } else {
    custom_path <- resolve_input_path(custom_path)
    
    if (!file.exists(custom_path)) {
      stop("The file provided in 'custom_eu_boundary_path' does not exist: ", custom_path, call. = FALSE)
    }
    
    euboundary <- sf::st_read(custom_path, quiet = TRUE)
  }
  
  if (nrow(euboundary) == 0) {
    stop("The EU boundary layer does not contain any features.", call. = FALSE)
  }
  
  empty_features <- sf::st_is_empty(euboundary)
  if (any(empty_features)) {
    euboundary <- euboundary[!empty_features, , drop = FALSE]
  }
  
  if (nrow(euboundary) == 0) {
    stop("The EU boundary layer only contains empty geometries.", call. = FALSE)
  }
  
  if (!all(sf::st_is_valid(euboundary))) {
    euboundary <- suppressWarnings(sf::st_make_valid(euboundary))
  }
  
  boundary_crs <- sf::st_crs(euboundary)
  if (is.na(boundary_crs)) {
    stop("The EU boundary layer must have a defined CRS.", call. = FALSE)
  }
  
  if (!is.null(reference)) {
    reference_crs <- get_reference_crs(reference)
    
    if (!isTRUE(boundary_crs == reference_crs)) {
      euboundary <- sf::st_transform(euboundary, reference_crs)
    }
    
    if (!isTRUE(sf::st_crs(euboundary) == reference_crs)) {
      stop("The EU boundary layer CRS could not be aligned to the reference raster/vector CRS.", call. = FALSE)
    }
  }
  
  euboundary
}


#-----------------------------------------------------------------
#--Load a named raster stack from manifest rows--------------------
#-----------------------------------------------------------------
load_named_raster_stack <- function(stack_rows) {
  required_cols <- c("var_name", "file_path")
  missing_cols <- setdiff(required_cols, names(stack_rows))
  
  if (length(missing_cols) > 0) {
    stop(
      "The raster stack rows are missing required columns: ",
      paste(missing_cols, collapse = ", "),
      call. = FALSE
    )
  }
  
  if (nrow(stack_rows) == 0) {
    stop("No raster rows were supplied to build the climate stack.", call. = FALSE)
  }
  
  if (anyDuplicated(stack_rows$var_name)) {
    stop("The raster stack rows contain duplicated 'var_name' values.", call. = FALSE)
  }
  
  stack <- terra::rast(stack_rows$file_path)
  
  if (terra::nlyr(stack) != nrow(stack_rows)) {
    stop(
      "The number of raster layers loaded does not match the number of manifest rows provided.",
      call. = FALSE
    )
  }
  
  expected_names <- as.character(stack_rows$var_name)
  names(stack) <- expected_names
  
  if (!identical(names(stack), expected_names)) {
    stop("The raster stack could not be named exactly as specified in 'var_name'.", call. = FALSE)
  }
  
  stack
}


#-----------------------------------------------------------------
#--Build a signature for a normalized raster stack----------------
#-----------------------------------------------------------------
build_manifest_stack_signature <- function(stack_rows,
                                           signature_scope = "manifest_stack_v1") {
  if (nrow(stack_rows) == 0) {
    stop("No raster rows were supplied to build a stack signature.", call. = FALSE)
  }
  
  raster_info <- file.info(stack_rows$file_path)
  if (any(is.na(raster_info$size)) || any(is.na(raster_info$mtime))) {
    stop("One or more raster files could not be inspected for signature generation.", call. = FALSE)
  }
  
  key_lines <- c(
    paste("signature_scope", signature_scope, sep = "="),
    paste("n_layers", nrow(stack_rows), sep = "="),
    paste(
      stack_rows$var_name,
      stack_rows$file_path,
      raster_info$size,
      format(raster_info$mtime, tz = "UTC", usetz = TRUE),
      sep = "|"
    )
  )
  
  key_file <- tempfile(pattern = "manifest_stack_signature_", fileext = ".txt")
  on.exit(unlink(key_file), add = TRUE)
  writeLines(key_lines, key_file, useBytes = TRUE)
  unname(tools::md5sum(key_file))
}


#-----------------------------------------------------------------
#--Write a deterministic processed raster stack when needed-------
#-----------------------------------------------------------------
materialize_processed_manifest_stack <- function(stack_rows,
                                                 output_file,
                                                 signature_scope = "manifest_stack_v1",
                                                 signature_file = paste0(output_file, ".signature.txt"),
                                                 apply_common_na_mask = TRUE) {
  output_file <- resolve_input_path(output_file)
  signature_file <- resolve_input_path(signature_file)
  
  if (!dir.exists(dirname(output_file))) {
    dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)
  }
  
  expected_names <- as.character(stack_rows$var_name)
  expected_signature <- build_manifest_stack_signature(stack_rows, signature_scope = signature_scope)
  
  if (file.exists(output_file) && file.exists(signature_file)) {
    current_signature <- paste(readLines(signature_file, warn = FALSE), collapse = "\n")
    cached_stack <- tryCatch(terra::rast(output_file), error = function(e) NULL)
    
    if (!is.null(cached_stack) &&
        identical(trimws(current_signature), expected_signature) &&
        terra::nlyr(cached_stack) == length(expected_names) &&
        identical(names(cached_stack), expected_names)) {
      return(cached_stack)
    }
  }
  
  stack <- load_named_raster_stack(stack_rows)
  if (apply_common_na_mask) {
    stack <- terra::mask(stack, anyNA(stack), maskvalue = 1)
  }
  
  terra::writeRaster(
    stack,
    filename = output_file,
    overwrite = TRUE,
    wopt = list(gdal = c("COMPRESS=LZW"))
  )
  writeLines(expected_signature, signature_file, useBytes = TRUE)
  
  cached_stack <- terra::rast(output_file)
  if (terra::nlyr(cached_stack) != length(expected_names) ||
      !identical(names(cached_stack), expected_names)) {
    stop(
      "The processed raster stack does not match the expected predictor names.",
      call. = FALSE
    )
  }
  
  cached_stack
}


#-----------------------------------------------------------------
#--Build a stable cache key for current user-specific climate------
#-----------------------------------------------------------------
build_user_specific_climate_cache_key <- function(climate_manifest) {
  current_rows <- climate_manifest$current_rows
  if (is.null(current_rows) || nrow(current_rows) == 0) {
    stop("The user-specific climate manifest does not contain current/current rows.", call. = FALSE)
  }
  
  raster_info <- file.info(current_rows$file_path)
  if (any(is.na(raster_info$size)) || any(is.na(raster_info$mtime))) {
    stop("One or more current climate rasters could not be inspected for caching.", call. = FALSE)
  }
  
  key_lines <- c(
    "cache_scope=current_rows_v2",
    paste("n_current_layers", nrow(current_rows), sep = "="),
    paste(
      current_rows$var_name,
      current_rows$file_path,
      raster_info$size,
      format(raster_info$mtime, tz = "UTC", usetz = TRUE),
      sep = "|"
    )
  )
  
  key_file <- tempfile(pattern = "user_specific_climate_key_", fileext = ".txt")
  on.exit(unlink(key_file), add = TRUE)
  writeLines(key_lines, key_file, useBytes = TRUE)
  unname(tools::md5sum(key_file))
}


#-----------------------------------------------------------------
#--Materialize current user-specific climate stack to disk--------
#-----------------------------------------------------------------
materialize_user_specific_current_stack <- function(climate_manifest,
                                                    cache_dir = file.path("data", "external", "climate", "chelsa_current", "processed"),
                                                    cache_prefix = "user_specific_current_stack") {
  cache_dir <- resolve_input_path(cache_dir)
  if (!dir.exists(cache_dir)) {
    dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
  }
  
  cache_key <- build_user_specific_climate_cache_key(climate_manifest)
  cache_file <- file.path(cache_dir, paste0(cache_prefix, "_", cache_key, ".tif"))
  expected_names <- as.character(climate_manifest$current_rows$var_name)
  
  if (file.exists(cache_file)) {
    cached_stack <- tryCatch(terra::rast(cache_file), error = function(e) NULL)
    
    if (!is.null(cached_stack) &&
        terra::nlyr(cached_stack) == length(expected_names) &&
        identical(names(cached_stack), expected_names)) {
      return(cached_stack)
    }
  }
  
  stack <- load_named_raster_stack(climate_manifest$current_rows)
  stack <- terra::mask(stack, anyNA(stack), maskvalue = 1)
  
  terra::writeRaster(
    stack,
    filename = cache_file,
    overwrite = TRUE,
    wopt = list(gdal = c("COMPRESS=LZW"))
  )
  
  cached_stack <- terra::rast(cache_file)
  
  if (terra::nlyr(cached_stack) != length(expected_names) ||
      !identical(names(cached_stack), expected_names)) {
    stop(
      "The cached user-specific climate stack does not match the expected predictor names.",
      call. = FALSE
    )
  }
  
  cached_stack
}


#-----------------------------------------------------------------
#--Align a continuous raster to a template raster------------------
#-----------------------------------------------------------------
align_continuous_raster <- function(raster, template, method = "bilinear") {
  if (!inherits(raster, "SpatRaster")) {
    stop("'raster' must be a terra SpatRaster.", call. = FALSE)
  }
  
  if (!inherits(template, "SpatRaster")) {
    stop("'template' must be a terra SpatRaster.", call. = FALSE)
  }
  
  if (!nzchar(terra::crs(raster))) {
    stop("The raster to align does not have a defined CRS.", call. = FALSE)
  }
  
  if (!nzchar(terra::crs(template))) {
    stop("The template raster does not have a defined CRS.", call. = FALSE)
  }
  
  if (!isTRUE(terra::same.crs(raster, template))) {
    return(terra::project(raster, template, method = method))
  }
  
  if (!isTRUE(terra::compareGeom(raster, template, lyrs = FALSE, stopOnError = FALSE))) {
    return(terra::resample(raster, template, method = method))
  }
  
  raster
}


#-----------------------------------------------------------------
#--Create an approximately 5 km template for pseudoabsences-------
#-----------------------------------------------------------------
create_pseudoabsence_template <- function(raster_layer, target_resolution_m = 5000) {
  stopifnot(inherits(raster_layer, "SpatRaster"))
  
  x_res <- abs(terra::xres(raster_layer))
  y_res <- abs(terra::yres(raster_layer))
  
  target_resolution <- if (terra::is.lonlat(raster_layer)) {
    target_resolution_m / 111320
  } else {
    target_resolution_m
  }
  
  if (x_res >= target_resolution && y_res >= target_resolution) {
    return(raster_layer)
  }
  
  fact_x <- max(1, round(target_resolution / x_res))
  fact_y <- max(1, round(target_resolution / y_res))
  
  terra::aggregate(raster_layer, fact = c(fact_x, fact_y), fun = mean, na.rm = TRUE)
}


#-----------------------------------------------------------------
#--Supported current/future combinations for user raster manifests-
#-----------------------------------------------------------------
get_supported_user_raster_combos <- function() {
  c(
    "current__current",
    "2041-2070__ssp126",
    "2041-2070__ssp370",
    "2041-2070__ssp585",
    "2071-2100__ssp126",
    "2071-2100__ssp370",
    "2071-2100__ssp585"
  )
}


#-----------------------------------------------------------------
#--Validate and normalize a user raster manifest-------------------
#-----------------------------------------------------------------
load_user_specific_raster_manifest <- function(manifest_path,
                                               config_name,
                                               data_label,
                                               require_all_future = TRUE) {
  manifest_path <- resolve_input_path(manifest_path)
  
  if (!file.exists(manifest_path)) {
    stop("The file provided in '", config_name, "' does not exist: ", manifest_path, call. = FALSE)
  }
  
  manifest <- utils::read.csv(manifest_path, stringsAsFactors = FALSE, check.names = FALSE)
  required_cols <- c("period", "scenario", "var_name", "file_path")
  missing_cols <- setdiff(required_cols, names(manifest))
  
  if (length(missing_cols) > 0) {
    stop(
      "The ", data_label, " manifest is missing required columns: ",
      paste(missing_cols, collapse = ", "),
      call. = FALSE
    )
  }
  
  manifest <- manifest[, required_cols]
  manifest[] <- lapply(manifest, trimws)
  
  if (nrow(manifest) == 0) {
    stop("The ", data_label, " manifest is empty.", call. = FALSE)
  }
  
  manifest$period <- tolower(manifest$period)
  manifest$scenario <- tolower(manifest$scenario)
  manifest$combo_id <- paste(manifest$period, manifest$scenario, sep = "__")
  
  valid_combos <- get_supported_user_raster_combos()
  required_future_ids <- setdiff(valid_combos, "current__current")
  
  invalid_combo_ids <- setdiff(unique(manifest$combo_id), valid_combos)
  if (length(invalid_combo_ids) > 0) {
    stop(
      "The ", data_label, " manifest contains unsupported period/scenario combinations: ",
      paste(invalid_combo_ids, collapse = ", "),
      call. = FALSE
    )
  }
  
  if (any(is.na(manifest$period)) || any(is.na(manifest$scenario))) {
    stop("The ", data_label, " manifest contains missing 'period' or 'scenario' values.", call. = FALSE)
  }
  
  if (any(is.na(manifest$var_name)) || any(!nzchar(manifest$var_name))) {
    stop("The ", data_label, " manifest contains empty 'var_name' values.", call. = FALSE)
  }
  
  if (any(is.na(manifest$file_path)) || any(!nzchar(manifest$file_path))) {
    stop("The ", data_label, " manifest contains empty 'file_path' values.", call. = FALSE)
  }
  
  manifest_dir <- dirname(manifest_path)
  manifest$file_path <- vapply(
    manifest$file_path,
    resolve_input_path,
    character(1),
    base_dir = manifest_dir
  )
  
  missing_files <- unique(manifest$file_path[!file.exists(manifest$file_path)])
  if (length(missing_files) > 0) {
    stop(
      "The ", data_label, " manifest references files that do not exist: ",
      paste(missing_files, collapse = ", "),
      call. = FALSE
    )
  }
  
  if (!"current__current" %in% manifest$combo_id) {
    stop("The ", data_label, " manifest must contain rows with period='current' and scenario='current'.", call. = FALSE)
  }
  
  stack_rows <- split(manifest, manifest$combo_id)
  master_rows <- stack_rows[["current__current"]]
  
  if (anyDuplicated(master_rows$var_name)) {
    stop("The current/current ", data_label, " rows contain duplicated 'var_name' values.", call. = FALSE)
  }
  
  master_var_names <- master_rows$var_name
  master_reference <- NULL
  normalized_rows <- list()
  future_combo_ids_present <- intersect(required_future_ids, unique(manifest$combo_id))
  
  if (require_all_future) {
    missing_future_ids <- setdiff(required_future_ids, future_combo_ids_present)
    if (length(missing_future_ids) > 0) {
      stop(
        "The ", data_label, " manifest is missing required future period/scenario combinations: ",
        paste(missing_future_ids, collapse = ", "),
        call. = FALSE
      )
    }
  } else if (length(future_combo_ids_present) > 0 && !setequal(future_combo_ids_present, required_future_ids)) {
    missing_future_ids <- setdiff(required_future_ids, future_combo_ids_present)
    stop(
      "The ", data_label, " manifest must either omit all future period/scenario combinations or provide all of them. Missing: ",
      paste(missing_future_ids, collapse = ", "),
      call. = FALSE
    )
  }
  
  for (combo_id in valid_combos) {
    rows <- stack_rows[[combo_id]]
    
    if (is.null(rows)) {
      next
    }
    
    if (anyDuplicated(rows$var_name)) {
      stop("The ", data_label, " manifest contains duplicated 'var_name' values for ", combo_id, ".", call. = FALSE)
    }
    
    if (!setequal(rows$var_name, master_var_names)) {
      stop(
        "The predictor names for ",
        combo_id,
        " do not exactly match the current/current predictor names.",
        call. = FALSE
      )
    }
    
    rows <- rows[match(master_var_names, rows$var_name), , drop = FALSE]
    rasters <- lapply(rows$file_path, terra::rast)
    
    if (any(vapply(rasters, terra::nlyr, numeric(1)) != 1)) {
      stop("Each file in '", config_name, "' must contain exactly one raster layer.", call. = FALSE)
    }
    
    reference_raster <- rasters[[1]]
    if (!nzchar(terra::crs(reference_raster))) {
      stop("The ", data_label, " raster CRS is missing for ", rows$file_path[1], ".", call. = FALSE)
    }
    
    alignment_ok <- vapply(
      rasters[-1],
      function(x) terra::compareGeom(x, reference_raster, lyrs = FALSE, stopOnError = FALSE),
      logical(1)
    )
    
    if (length(alignment_ok) > 0 && !all(alignment_ok)) {
      stop(
        "All rasters within ",
        combo_id,
        " must share the same grid, extent, resolution, and CRS.",
        call. = FALSE
      )
    }
    
    if (is.null(master_reference)) {
      master_reference <- reference_raster
    } else if (!isTRUE(terra::same.crs(reference_raster, master_reference))) {
      stop(
        "All user-specific ", data_label, " rasters must share the same CRS in this workflow.",
        call. = FALSE
      )
    }
    
    normalized_rows[[combo_id]] <- rows
  }
  
  future_rows <- normalized_rows[required_future_ids[required_future_ids %in% names(normalized_rows)]]
  
  list(
    manifest_path = manifest_path,
    current_rows = normalized_rows[["current__current"]],
    future_rows = future_rows,
    master_var_names = master_var_names,
    has_future = length(future_rows) > 0,
    predictor_crs = terra::crs(master_reference)
  )
}


#-----------------------------------------------------------------
#--Validate and normalize a user climate manifest------------------
#-----------------------------------------------------------------
load_user_specific_climate_manifest <- function(manifest_path) {
  load_user_specific_raster_manifest(
    manifest_path = manifest_path,
    config_name = "user_specific_climate_data",
    data_label = "climate",
    require_all_future = TRUE
  )
}


#-----------------------------------------------------------------
#--Validate and normalize a user land-cover manifest--------------
#-----------------------------------------------------------------
load_user_specific_landcover_manifest <- function(manifest_path) {
  load_user_specific_raster_manifest(
    manifest_path = manifest_path,
    config_name = "user_specific_landcover_data",
    data_label = "land-cover",
    require_all_future = FALSE
  )
}


#-----------------------------------------------------------------
#--Validate only the current/current rows needed by validation----
#-----------------------------------------------------------------
load_user_specific_current_manifest <- function(manifest_path,
                                                config_name,
                                                data_label) {
  manifest_path <- resolve_input_path(manifest_path)

  if (!file.exists(manifest_path)) {
    stop("The file provided in '", config_name, "' does not exist: ",
         manifest_path, call. = FALSE)
  }

  manifest <- utils::read.csv(
    manifest_path,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  required_cols <- c("period", "scenario", "var_name", "file_path")
  missing_cols <- setdiff(required_cols, names(manifest))

  if (length(missing_cols) > 0L) {
    stop(
      "The ", data_label, " manifest is missing required columns: ",
      paste(missing_cols, collapse = ", "),
      call. = FALSE
    )
  }

  manifest <- manifest[, required_cols, drop = FALSE]
  manifest[] <- lapply(manifest, function(x) trimws(as.character(x)))
  manifest$period <- tolower(manifest$period)
  manifest$scenario <- tolower(manifest$scenario)
  current_rows <- manifest[
    manifest$period == "current" & manifest$scenario == "current",
    ,
    drop = FALSE
  ]

  if (nrow(current_rows) == 0L) {
    stop(
      "The ", data_label,
      " manifest must contain rows with period='current' and scenario='current'.",
      call. = FALSE
    )
  }
  if (any(is.na(current_rows$var_name)) || any(!nzchar(current_rows$var_name))) {
    stop("The current ", data_label, " rows contain empty predictor names.",
         call. = FALSE)
  }
  if (anyDuplicated(current_rows$var_name)) {
    stop("The current ", data_label, " rows contain duplicated predictor names.",
         call. = FALSE)
  }
  if (any(is.na(current_rows$file_path)) || any(!nzchar(current_rows$file_path))) {
    stop("The current ", data_label, " rows contain empty raster paths.",
         call. = FALSE)
  }

  manifest_dir <- dirname(manifest_path)
  current_rows$file_path <- vapply(
    current_rows$file_path,
    resolve_input_path,
    character(1),
    base_dir = manifest_dir
  )
  missing_files <- unique(current_rows$file_path[!file.exists(current_rows$file_path)])
  if (length(missing_files) > 0L) {
    stop(
      "The current ", data_label, " manifest rows reference files that do not exist: ",
      paste(missing_files, collapse = ", "),
      call. = FALSE
    )
  }

  rasters <- lapply(current_rows$file_path, terra::rast)
  if (any(vapply(rasters, terra::nlyr, numeric(1)) != 1L)) {
    stop("Each current raster in '", config_name,
         "' must contain exactly one layer.", call. = FALSE)
  }
  reference <- rasters[[1]]
  if (!nzchar(terra::crs(reference))) {
    stop("The current ", data_label, " raster CRS is missing.", call. = FALSE)
  }
  alignment_ok <- vapply(
    rasters[-1],
    function(x) terra::compareGeom(
      x,
      reference,
      lyrs = FALSE,
      stopOnError = FALSE
    ),
    logical(1)
  )
  if (length(alignment_ok) > 0L && !all(alignment_ok)) {
    stop(
      "All current ", data_label,
      " rasters must share the same grid, extent, resolution, and CRS.",
      call. = FALSE
    )
  }

  list(
    manifest_path = manifest_path,
    current_rows = current_rows,
    future_rows = list(),
    master_var_names = current_rows$var_name,
    has_future = FALSE,
    predictor_crs = terra::crs(reference)
  )
}


#-----------------------------------------------------------------
#--Build and compare a portable predictor-stack signature---------
#-----------------------------------------------------------------
build_predictor_signature <- function(predictor_stack,
                                      predictor_names = names(predictor_stack)) {
  missing_predictors <- setdiff(predictor_names, names(predictor_stack))
  if (length(missing_predictors) > 0L) {
    stop(
      "Cannot build predictor signature; missing layer(s): ",
      paste(missing_predictors, collapse = ", "),
      call. = FALSE
    )
  }
  selected <- terra::subset(predictor_stack, predictor_names)
  if (terra::nlyr(selected) == 0L) {
    stop("Cannot build a predictor signature for an empty stack.", call. = FALSE)
  }

  list(
    predictor_names = names(selected),
    crs = terra::crs(selected),
    nrow = as.integer(terra::nrow(selected)),
    ncol = as.integer(terra::ncol(selected)),
    resolution = as.numeric(terra::res(selected)),
    extent = c(
      xmin = terra::xmin(selected),
      xmax = terra::xmax(selected),
      ymin = terra::ymin(selected),
      ymax = terra::ymax(selected)
    )
  )
}


validate_model_predictor_stack <- function(predictor_stack,
                                           selected_predictors,
                                           saved_crs = NULL,
                                           saved_signature = NULL,
                                           tolerance = 1e-7) {
  missing_predictors <- setdiff(selected_predictors, names(predictor_stack))
  if (length(missing_predictors) > 0L) {
    return(list(
      valid = FALSE,
      legacy = is.null(saved_signature),
      reason = paste0(
        "missing selected predictor(s): ",
        paste(missing_predictors, collapse = ", ")
      ),
      signature = NULL
    ))
  }

  active_signature <- build_predictor_signature(
    predictor_stack,
    selected_predictors
  )

  if (!is.null(saved_signature)) {
    required_fields <- c(
      "predictor_names", "crs", "nrow", "ncol", "resolution", "extent"
    )
    missing_fields <- setdiff(required_fields, names(saved_signature))
    if (length(missing_fields) > 0L) {
      return(list(
        valid = FALSE,
        legacy = FALSE,
        reason = paste0(
          "saved predictor signature is missing field(s): ",
          paste(missing_fields, collapse = ", ")
        ),
        signature = active_signature
      ))
    }

    checks <- c(
      predictor_names = identical(
        as.character(saved_signature$predictor_names),
        as.character(active_signature$predictor_names)
      ),
      crs = isTRUE(
        sf::st_crs(saved_signature$crs) == sf::st_crs(active_signature$crs)
      ),
      dimensions = identical(
        as.integer(c(saved_signature$nrow, saved_signature$ncol)),
        as.integer(c(active_signature$nrow, active_signature$ncol))
      ),
      resolution = isTRUE(all.equal(
        as.numeric(saved_signature$resolution),
        active_signature$resolution,
        tolerance = tolerance,
        check.attributes = FALSE
      )),
      extent = isTRUE(all.equal(
        as.numeric(saved_signature$extent),
        active_signature$extent,
        tolerance = tolerance,
        check.attributes = FALSE
      ))
    )
    if (!all(checks)) {
      return(list(
        valid = FALSE,
        legacy = FALSE,
        reason = paste0(
          "predictor signature mismatch: ",
          paste(names(checks)[!checks], collapse = ", ")
        ),
        signature = active_signature
      ))
    }
  } else if (!is.null(saved_crs) && length(saved_crs) == 1L &&
             !is.na(saved_crs) && nzchar(saved_crs) &&
             !isTRUE(sf::st_crs(saved_crs) == sf::st_crs(active_signature$crs))) {
    return(list(
      valid = FALSE,
      legacy = TRUE,
      reason = "predictor CRS does not match the legacy model metadata",
      signature = active_signature
    ))
  }

  list(
    valid = TRUE,
    legacy = is.null(saved_signature),
    reason = if (is.null(saved_signature)) {
      "compatible by selected predictor names and CRS (legacy metadata)"
    } else {
      "compatible predictor signature"
    },
    signature = active_signature
  )
}


#-----------------------------------------------------------------
#--Resolve a saved manifest with a portable configured fallback---
#-----------------------------------------------------------------
resolve_model_current_predictors <- function(saved_manifest_path,
                                             configured_manifest_path,
                                             config_name,
                                             data_label,
                                             selected_predictors,
                                             saved_crs = NULL,
                                             saved_signature = NULL,
                                             materialize = function(x) {
                                               load_named_raster_stack(x$current_rows)
                                             }) {
  candidates <- c(
    saved = if (is.null(saved_manifest_path)) NA_character_ else saved_manifest_path,
    configured = if (is.null(configured_manifest_path)) NA_character_ else configured_manifest_path
  )
  candidates <- candidates[!is.na(candidates) & nzchar(candidates)]
  if (length(candidates) > 1L) {
    normalized <- vapply(candidates, resolve_input_path, character(1))
    candidates <- candidates[!duplicated(tolower(normalized))]
  }
  if (length(candidates) == 0L) {
    stop(
      "The saved ", data_label,
      " model requires user-specific predictors, but neither a saved nor configured manifest path is available.",
      call. = FALSE
    )
  }

  failures <- character(0)
  for (candidate_name in names(candidates)) {
    candidate_path <- candidates[[candidate_name]]
    attempt <- tryCatch({
      manifest <- load_user_specific_current_manifest(
        manifest_path = candidate_path,
        config_name = config_name,
        data_label = data_label
      )
      predictor_stack <- materialize(manifest)
      compatibility <- validate_model_predictor_stack(
        predictor_stack = predictor_stack,
        selected_predictors = selected_predictors,
        saved_crs = saved_crs,
        saved_signature = saved_signature
      )
      if (!isTRUE(compatibility$valid)) {
        stop(compatibility$reason, call. = FALSE)
      }
      list(
        manifest = manifest,
        predictor_stack = predictor_stack,
        compatibility = compatibility
      )
    }, error = function(e) e)

    if (inherits(attempt, "error")) {
      failures <- c(
        failures,
        paste0(candidate_name, " [", candidate_path, "]: ", conditionMessage(attempt))
      )
      next
    }

    if (identical(candidate_name, "configured") && length(failures) > 0L) {
      warning(
        "The saved ", data_label, " manifest could not be reused; using the configured manifest '",
        attempt$manifest$manifest_path, "'. Saved-path failure: ",
        paste(failures, collapse = " | "),
        call. = FALSE
      )
    }
    if (isTRUE(attempt$compatibility$legacy)) {
      warning(
        "The ", data_label,
        " model predates portable predictor signatures. Compatibility was checked using selected predictor names and CRS only.",
        call. = FALSE
      )
    }

    return(c(
      attempt,
      list(source = candidate_name, failures = failures)
    ))
  }

  stop(
    "No compatible current ", data_label, " manifest could be resolved. ",
    paste(failures, collapse = " | "),
    call. = FALSE
  )
}


#-----------------------------------------------------------------
#--Get deterministic processed paths for user land-cover stacks---
#-----------------------------------------------------------------
get_user_specific_landcover_processed_target <- function(period = "current",
                                                         scenario = "current",
                                                         processed_dir = file.path("data", "external", "habitat", "processed")) {
  processed_dir <- resolve_input_path(processed_dir)
  
  if (identical(period, "current") && identical(scenario, "current")) {
    raster_file <- file.path(processed_dir, "habitat_stack.tif")
  } else {
    raster_file <- file.path(
      processed_dir,
      "future",
      period,
      scenario,
      paste0("habitat_stack_", period, "_", scenario, ".tif")
    )
  }
  
  list(
    raster_file = raster_file,
    signature_file = paste0(raster_file, ".signature.txt")
  )
}


#-----------------------------------------------------------------
#--Materialize deterministic user land-cover stacks---------------
#-----------------------------------------------------------------
materialize_user_specific_landcover_stack <- function(stack_rows,
                                                      period = "current",
                                                      scenario = "current",
                                                      processed_dir = file.path("data", "external", "habitat", "processed")) {
  target <- get_user_specific_landcover_processed_target(
    period = period,
    scenario = scenario,
    processed_dir = processed_dir
  )
  
  materialize_processed_manifest_stack(
    stack_rows = stack_rows,
    output_file = target$raster_file,
    signature_file = target$signature_file,
    signature_scope = paste0("user_specific_landcover__", period, "__", scenario, "__v1"),
    apply_common_na_mask = TRUE
  )
}


#-----------------------------------------------------------------
#--Validate fold assignments across all required strata----------
#-----------------------------------------------------------------
validate_cv_context_classes <- function(records,
                                        context_col = "cv_context",
                                        class_col = "species") {
  if (!all(c(context_col, class_col) %in% names(records))) {
    return(list(valid = TRUE, reason = NA_character_))
  }
  contexts <- sort(unique(as.character(records[[context_col]])))
  counts <- table(
    factor(as.character(records[[context_col]]), levels = contexts),
    factor(records[[class_col]], levels = c(0, 1))
  )
  missing_rows <- which(counts == 0L, arr.ind = TRUE)
  if (nrow(missing_rows) == 0L) {
    return(list(valid = TRUE, reason = NA_character_))
  }
  missing_labels <- paste0(
    rownames(counts)[missing_rows[, "row"]],
    " class ",
    colnames(counts)[missing_rows[, "col"]]
  )
  list(
    valid = FALSE,
    reason = paste0(
      "required context/class strata contain no records: ",
      paste(missing_labels, collapse = ", ")
    )
  )
}


validate_cv_fold_assignments <- function(records,
                                         fold_ids,
                                         k,
                                         stratum_col = "cv_stratum",
                                         min_train = 2L,
                                         min_test = 2L) {
  if (!stratum_col %in% names(records)) {
    stop("The CV records are missing '", stratum_col, "'.", call. = FALSE)
  }
  if (length(fold_ids) != nrow(records)) {
    return(list(
      valid = FALSE,
      reason = paste0(
        "fold assignment length ", length(fold_ids),
        " differs from record count ", nrow(records)
      ),
      counts = data.frame(),
      min_train_count = 0L,
      min_test_count = 0L,
      unassigned_count = as.integer(nrow(records))
    ))
  }
  invalid_assignment <- is.na(fold_ids) | !fold_ids %in% seq_len(k)
  if (any(invalid_assignment)) {
    return(list(
      valid = FALSE,
      reason = "one or more records were not assigned to a valid fold",
      counts = data.frame(),
      min_train_count = 0L,
      min_test_count = 0L,
      unassigned_count = as.integer(sum(invalid_assignment))
    ))
  }

  strata <- sort(unique(as.character(records[[stratum_col]])))
  counts <- do.call(
    rbind,
    lapply(seq_len(k), function(fold) {
      do.call(
        rbind,
        lapply(strata, function(stratum) {
          in_stratum <- as.character(records[[stratum_col]]) == stratum
          data.frame(
            fold = fold,
            stratum = stratum,
            train_n = sum(in_stratum & fold_ids != fold),
            test_n = sum(in_stratum & fold_ids == fold),
            stringsAsFactors = FALSE
          )
        })
      )
    })
  )
  min_train_count <- if (nrow(counts) == 0L) 0L else min(counts$train_n)
  min_test_count <- if (nrow(counts) == 0L) 0L else min(counts$test_n)
  valid <- nrow(counts) > 0L &&
    min_train_count >= min_train &&
    min_test_count >= min_test

  list(
    valid = valid,
    reason = if (valid) {
      "all folds meet class-count requirements"
    } else {
      paste0(
        "minimum train/test stratum counts are ",
        min_train_count, "/", min_test_count,
        "; required ", min_train, "/", min_test
      )
    },
    counts = counts,
    min_train_count = min_train_count,
    min_test_count = min_test_count,
    unassigned_count = 0L
  )
}


cv_plan_attempt_row <- function(method,
                                k,
                                block_size_m,
                                validation,
                                reason = validation$reason) {
  data.frame(
    diagnostic_type = "partition_attempt",
    cv_method = method,
    effective_folds = as.integer(k),
    block_size_m = if (is.null(block_size_m)) NA_real_ else as.numeric(block_size_m),
    valid = isTRUE(validation$valid),
    min_train_count = as.integer(validation$min_train_count),
    min_test_count = as.integer(validation$min_test_count),
    unassigned_records = if (is.null(validation$unassigned_count)) {
      NA_integer_
    } else {
      as.integer(validation$unassigned_count)
    },
    reason = as.character(reason),
    stringsAsFactors = FALSE
  )
}


#-----------------------------------------------------------------
#--Build adaptive spatial block folds-----------------------------
#-----------------------------------------------------------------
make_spatial_cv_plan <- function(records,
                                 requested_folds,
                                 block_sizes_m,
                                 stratum_col = "cv_stratum",
                                 min_train = 2L,
                                 min_test = 2L,
                                 seed = 123L,
                                 iteration = 200L) {
  if (!inherits(records, "sf") && !inherits(records, "SpatVector")) {
    stop("Spatial CV records must be an sf or SpatVector object.", call. = FALSE)
  }
  requested_folds <- as.integer(requested_folds)
  block_sizes_m <- as.numeric(block_sizes_m)
  attempts <- list()
  attempt_index <- 0L
  context_validation <- validate_cv_context_classes(records)
  if (!isTRUE(context_validation$valid)) {
    validation <- list(
      valid = FALSE,
      reason = context_validation$reason,
      min_train_count = 0L,
      min_test_count = 0L,
      unassigned_count = 0L
    )
    diagnostics <- cv_plan_attempt_row(
      "spatial_block", requested_folds, NA_real_, validation
    )
    return(list(
      valid = FALSE,
      method = "spatial_block",
      requested_folds = requested_folds,
      k = NA_integer_,
      block_size_m = NA_real_,
      fold_ids = rep(NA_integer_, nrow(records)),
      fold_counts = data.frame(),
      fallback_reason = context_validation$reason,
      diagnostics = diagnostics,
      blockcv = NULL
    ))
  }

  for (k in seq.int(requested_folds, 2L, by = -1L)) {
    stratum_counts <- table(records[[stratum_col]])
    if (length(stratum_counts) == 0L || any(stratum_counts < k * min_test)) {
      attempt_index <- attempt_index + 1L
      validation <- list(
        valid = FALSE,
        reason = "at least one stratum is too small for this fold count",
        min_train_count = 0L,
        min_test_count = 0L
      )
      attempts[[attempt_index]] <- cv_plan_attempt_row(
        "spatial_block", k, NA_real_, validation
      )
      next
    }

    for (block_size_m in block_sizes_m) {
      spatial_result <- tryCatch(
        blockCV::cv_spatial(
          x = records,
          column = stratum_col,
          k = k,
          hexagon = TRUE,
          selection = "random",
          iteration = as.integer(iteration),
          size = block_size_m,
          seed = as.integer(seed),
          progress = FALSE,
          report = FALSE,
          plot = FALSE
        ),
        error = function(e) e
      )
      attempt_index <- attempt_index + 1L

      if (inherits(spatial_result, "error")) {
        validation <- list(
          valid = FALSE,
          reason = conditionMessage(spatial_result),
          min_train_count = 0L,
          min_test_count = 0L
        )
        attempts[[attempt_index]] <- cv_plan_attempt_row(
          "spatial_block", k, block_size_m, validation
        )
        next
      }

      fold_ids <- as.integer(spatial_result$folds_ids)
      validation <- validate_cv_fold_assignments(
        records = records,
        fold_ids = fold_ids,
        k = k,
        stratum_col = stratum_col,
        min_train = min_train,
        min_test = min_test
      )
      attempts[[attempt_index]] <- cv_plan_attempt_row(
        "spatial_block", k, block_size_m, validation
      )

      if (isTRUE(validation$valid)) {
        return(list(
          valid = TRUE,
          method = "spatial_block",
          requested_folds = requested_folds,
          k = k,
          block_size_m = block_size_m,
          fold_ids = fold_ids,
          fold_counts = validation$counts,
          fallback_reason = NA_character_,
          diagnostics = dplyr::bind_rows(attempts),
          blockcv = spatial_result
        ))
      }
    }
  }

  list(
    valid = FALSE,
    method = "spatial_block",
    requested_folds = requested_folds,
    k = NA_integer_,
    block_size_m = NA_real_,
    fold_ids = rep(NA_integer_, nrow(records)),
    fold_counts = data.frame(),
    fallback_reason = paste(unique(vapply(attempts, function(x) x$reason[[1]], character(1))), collapse = " | "),
    diagnostics = dplyr::bind_rows(attempts),
    blockcv = NULL
  )
}


#-----------------------------------------------------------------
#--Build grouped, stratified, deterministic non-spatial folds-----
#-----------------------------------------------------------------
make_stratified_kfold_plan <- function(records,
                                       requested_folds,
                                       stratum_col = "cv_stratum",
                                       group_col = "cv_group",
                                       min_train = 2L,
                                       min_test = 2L,
                                       seed = 123L,
                                       iteration = 200L,
                                       fallback_reason = NA_character_) {
  if (!stratum_col %in% names(records)) {
    stop("The CV records are missing '", stratum_col, "'.", call. = FALSE)
  }
  if (!group_col %in% names(records)) {
    records[[group_col]] <- seq_len(nrow(records))
  }

  requested_folds <- as.integer(requested_folds)
  context_validation <- validate_cv_context_classes(records)
  if (!isTRUE(context_validation$valid)) {
    validation <- list(
      valid = FALSE,
      reason = context_validation$reason,
      min_train_count = 0L,
      min_test_count = 0L,
      unassigned_count = 0L
    )
    diagnostics <- cv_plan_attempt_row(
      "stratified_kfold", requested_folds, NA_real_, validation
    )
    combined_reason <- paste(
      stats::na.omit(c(fallback_reason, context_validation$reason)),
      collapse = " | "
    )
    return(list(
      valid = FALSE,
      method = "not_evaluable",
      requested_folds = requested_folds,
      k = NA_integer_,
      block_size_m = NA_real_,
      fold_ids = rep(NA_integer_, nrow(records)),
      fold_counts = data.frame(),
      fallback_reason = combined_reason,
      diagnostics = diagnostics,
      blockcv = NULL
    ))
  }
  strata <- sort(unique(as.character(records[[stratum_col]])))
  group_ids <- unique(as.character(records[[group_col]]))
  group_strata <- table(
    factor(as.character(records[[group_col]]), levels = group_ids),
    factor(as.character(records[[stratum_col]]), levels = strata)
  )
  attempts <- list()
  attempt_index <- 0L

  for (k in seq.int(requested_folds, 2L, by = -1L)) {
    stratum_counts <- colSums(group_strata)
    if (any(stratum_counts < k * min_test)) {
      attempt_index <- attempt_index + 1L
      validation <- list(
        valid = FALSE,
        reason = "at least one stratum is too small for this fold count",
        min_train_count = 0L,
        min_test_count = 0L
      )
      attempts[[attempt_index]] <- cv_plan_attempt_row(
        "stratified_kfold", k, NA_real_, validation
      )
      next
    }

    best <- NULL
    best_score <- Inf
    for (iteration_id in seq_len(as.integer(iteration))) {
      set.seed(as.integer(seed) + iteration_id + k * 1000L)
      group_totals <- rowSums(group_strata)
      order_groups <- order(
        -apply(group_strata, 1, max),
        -group_totals,
        stats::runif(length(group_ids))
      )
      fold_counts <- matrix(0, nrow = k, ncol = length(strata))
      fold_groups <- integer(k)
      group_fold <- integer(length(group_ids))
      target <- matrix(
        rep(stratum_counts / k, each = k),
        nrow = k,
        ncol = length(strata)
      )

      for (group_index in order_groups) {
        candidate_scores <- vapply(seq_len(k), function(fold) {
          candidate_counts <- fold_counts
          candidate_counts[fold, ] <- candidate_counts[fold, ] + group_strata[group_index, ]
          imbalance <- sum(
            ((candidate_counts - target)^2) /
              matrix(rep(pmax(target[1, ], 1), each = k), nrow = k)
          )
          group_balance <- sum((replace(fold_groups, fold, fold_groups[fold] + 1L) -
                                  (sum(fold_groups) + 1) / k)^2)
          imbalance + group_balance * 1e-3
        }, numeric(1))
        chosen <- sample(which(candidate_scores == min(candidate_scores)), 1L)
        group_fold[group_index] <- chosen
        fold_counts[chosen, ] <- fold_counts[chosen, ] + group_strata[group_index, ]
        fold_groups[chosen] <- fold_groups[chosen] + 1L
      }

      fold_map <- stats::setNames(group_fold, group_ids)
      fold_ids <- as.integer(fold_map[as.character(records[[group_col]])])
      validation <- validate_cv_fold_assignments(
        records = records,
        fold_ids = fold_ids,
        k = k,
        stratum_col = stratum_col,
        min_train = min_train,
        min_test = min_test
      )
      score <- if (isTRUE(validation$valid)) {
        stats::sd(validation$counts$test_n)
      } else {
        Inf
      }
      if (score < best_score) {
        best <- list(fold_ids = fold_ids, validation = validation)
        best_score <- score
      }
      if (isTRUE(validation$valid) && isTRUE(all.equal(score, 0))) {
        break
      }
    }

    attempt_index <- attempt_index + 1L
    if (is.null(best)) {
      validation <- list(
        valid = FALSE,
        reason = "no grouped stratified assignment was generated",
        min_train_count = 0L,
        min_test_count = 0L
      )
    } else {
      validation <- best$validation
    }
    attempts[[attempt_index]] <- cv_plan_attempt_row(
      "stratified_kfold", k, NA_real_, validation
    )

    if (!is.null(best) && isTRUE(validation$valid)) {
      return(list(
        valid = TRUE,
        method = "stratified_kfold",
        requested_folds = requested_folds,
        k = k,
        block_size_m = NA_real_,
        fold_ids = best$fold_ids,
        fold_counts = validation$counts,
        fallback_reason = fallback_reason,
        diagnostics = dplyr::bind_rows(attempts),
        blockcv = NULL
      ))
    }
  }

  list(
    valid = FALSE,
    method = "not_evaluable",
    requested_folds = requested_folds,
    k = NA_integer_,
    block_size_m = NA_real_,
    fold_ids = rep(NA_integer_, nrow(records)),
    fold_counts = data.frame(),
    fallback_reason = fallback_reason,
    diagnostics = dplyr::bind_rows(attempts),
    blockcv = NULL
  )
}


make_preferred_cv_plan <- function(records,
                                   requested_folds,
                                   block_sizes_m,
                                   enable_kfold_fallback = TRUE,
                                   stratum_col = "cv_stratum",
                                   group_col = "cv_group",
                                   min_train = 2L,
                                   min_test = 2L,
                                   seed = 123L,
                                   iteration = 200L) {
  spatial_plan <- make_spatial_cv_plan(
    records = records,
    requested_folds = requested_folds,
    block_sizes_m = block_sizes_m,
    stratum_col = stratum_col,
    min_train = min_train,
    min_test = min_test,
    seed = seed,
    iteration = iteration
  )
  if (isTRUE(spatial_plan$valid)) {
    return(spatial_plan)
  }
  if (!isTRUE(enable_kfold_fallback)) {
    spatial_plan$method <- "not_evaluable"
    return(spatial_plan)
  }

  kfold_plan <- make_stratified_kfold_plan(
    records = records,
    requested_folds = requested_folds,
    stratum_col = stratum_col,
    group_col = group_col,
    min_train = min_train,
    min_test = min_test,
    seed = seed,
    iteration = iteration,
    fallback_reason = spatial_plan$fallback_reason
  )
  kfold_plan$diagnostics <- dplyr::bind_rows(
    spatial_plan$diagnostics,
    kfold_plan$diagnostics
  )
  kfold_plan
}


#----------------------------------------------------------------------
#- Make predictions per model algorithm and dataset and obtain median -
#----------------------------------------------------------------------

compute_median_favourability <- function(model,
                                         datasets,
                                         top5_methods,
                                         prev_ratio) {
  
  #---------------------------
  #---- Make predictions -----
  #---------------------------
  
  env_favourability <- list()
  for(modelmethod in top5_methods){
    
    message("Predicting for method: ", modelmethod,".")
    
    for(dataset_name in names(datasets)) {
      
      #Load datasets
      dataset <- datasets[[dataset_name]]
      IDs <-dataset$ID
      dataset<-dplyr::select(dataset, -ID)
      
      #Predict for dataset
      dataset_suit <- predict(model,
                              newdata = dataset,
                              method = modelmethod)
      
      #Convert suitability to favourability
      dataset_fav<- favourability_from_prob(dataset_suit[[1]], prev_ratio)
      
      #Store in list
      env_favourability[[modelmethod]][[dataset_name]] <- data.frame(ID = IDs,
                                                                     fav = dataset_fav)
      
      #Clean up
      rm(dataset_suit, dataset_fav, IDs, dataset)
      
    }
  }  
  
  
  #-----------------------------------------
  #---- Calculate median favourability  ----
  #-----------------------------------------
  median_favourability<-lapply(
    names(datasets),
    function(dataset_name) {
      
      fav_matrix <- do.call(
        cbind,
        lapply(env_favourability, function(x) x[[dataset_name]]$fav)
      )
      
      data.frame(ID = env_favourability[[1]][[dataset_name]]$ID,
                 median_favourability = matrixStats::rowMedians(fav_matrix,na.rm = TRUE))
    }
  )
  
  names(median_favourability) <- names(datasets)
  
  return(median_favourability)
}


#----------------------------------------------------------------------
#- Make predictions per model algorithm and dataset and obtain median -
#- while tolerating individual method prediction failures             -
#----------------------------------------------------------------------
compute_median_favourability_safe <- function(model,
                                              datasets,
                                              top5_methods,
                                              prev_ratio,
                                              min_successful_methods = 3L) {
  if (!is.list(datasets) || length(datasets) == 0L || is.null(names(datasets)) ||
      any(!nzchar(names(datasets)))) {
    stop("Prediction datasets must be a named, non-empty list.", call. = FALSE)
  }
  if (length(prev_ratio) != 1L || !is.finite(prev_ratio) || prev_ratio <= 0) {
    stop("The prevalence ratio must be a finite positive value.", call. = FALSE)
  }

  env_favourability <- list()
  successful_methods <- character(0)
  failed_methods <- character(0)
  failed_reasons <- list()
  method_diagnostics <- list()
  diagnostic_index <- 0L
  
  for (modelmethod in top5_methods) {
    
    message("Predicting for method: ", modelmethod, ".")
    
    method_predictions <- list()
    method_failed <- FALSE
    method_error <- NULL
    
    for (dataset_name in names(datasets)) {
      dataset_input <- datasets[[dataset_name]]
      prediction_result <- tryCatch({
        if (!is.data.frame(dataset_input)) {
          dataset_input <- as.data.frame(dataset_input)
        }
        if (!"ID" %in% names(dataset_input)) {
          stop("dataset is missing its stable 'ID' column", call. = FALSE)
        }
        if (nrow(dataset_input) == 0L) {
          stop("dataset contains no usable rows", call. = FALSE)
        }
        if (anyDuplicated(dataset_input$ID)) {
          stop("dataset contains duplicated stable IDs", call. = FALSE)
        }
        IDs <- dataset_input$ID
        dataset <- dplyr::select(dataset_input, -ID)
        if (ncol(dataset) == 0L) {
          stop("dataset contains no predictor columns", call. = FALSE)
        }
        dataset_suit <- predict(model,
                                newdata = dataset,
                                method = modelmethod)
        dataset_suit_df <- as.data.frame(dataset_suit)
        if (ncol(dataset_suit_df) == 0L) {
          stop("prediction returned no columns", call. = FALSE)
        }
        dataset_prob <- dataset_suit_df[[1]]
        if (length(dataset_prob) != nrow(dataset_input)) {
          stop(
            "prediction length ", length(dataset_prob),
            " differs from input row count ", nrow(dataset_input),
            call. = FALSE
          )
        }
        if (any(!is.finite(dataset_prob))) {
          stop("prediction contains non-finite values", call. = FALSE)
        }
        dataset_fav <- favourability_from_prob(dataset_prob, prev_ratio)
        if (length(dataset_fav) != length(IDs) || any(!is.finite(dataset_fav))) {
          stop("favourability transformation returned invalid values", call. = FALSE)
        }
        data.frame(ID = IDs, fav = dataset_fav)
      }, error = function(e) {
        e
      })

      diagnostic_index <- diagnostic_index + 1L
      if (inherits(prediction_result, "error")) {
        method_failed <- TRUE
        method_error <- paste0(dataset_name, ": ", conditionMessage(prediction_result))
        method_diagnostics[[diagnostic_index]] <- data.frame(
          method = modelmethod,
          dataset = dataset_name,
          success = FALSE,
          error = conditionMessage(prediction_result),
          stringsAsFactors = FALSE
        )
        break
      }
      method_diagnostics[[diagnostic_index]] <- data.frame(
        method = modelmethod,
        dataset = dataset_name,
        success = TRUE,
        error = NA_character_,
        stringsAsFactors = FALSE
      )
      method_predictions[[dataset_name]] <- prediction_result
    }
    
    if (method_failed) {
      failed_methods <- c(failed_methods, modelmethod)
      failed_reasons[[modelmethod]] <- method_error
      next
    }
    
    env_favourability[[modelmethod]] <- method_predictions
    successful_methods <- c(successful_methods, modelmethod)
  }
  
  if (length(successful_methods) < min_successful_methods) {
    return(list(
      valid = FALSE,
      median_favourability = NULL,
      successful_methods = successful_methods,
      failed_methods = failed_methods,
      failed_reasons = failed_reasons,
      success_count = length(successful_methods),
      method_diagnostics = dplyr::bind_rows(method_diagnostics)
    ))
  }
  
  median_favourability <- lapply(
    names(datasets),
    function(dataset_name) {
      fav_matrix <- do.call(
        cbind,
        lapply(successful_methods, function(method_name) {
          env_favourability[[method_name]][[dataset_name]]$fav
        })
      )
      
      data.frame(
        ID = env_favourability[[successful_methods[[1]]]][[dataset_name]]$ID,
        median_favourability = matrixStats::rowMedians(fav_matrix, na.rm = TRUE)
      )
    }
  )
  
  names(median_favourability) <- names(datasets)
  
  list(
    valid = TRUE,
    median_favourability = median_favourability,
    successful_methods = successful_methods,
    failed_methods = failed_methods,
    failed_reasons = failed_reasons,
    success_count = length(successful_methods),
    method_diagnostics = dplyr::bind_rows(method_diagnostics)
  )
}


#------------------------------------------
#----- Define Boyce helper functions -----
#------------------------------------------
#fit_vals are suitability or favourability values of a background sample or all available pixels in the target region
#obs_vals are suitability or favourability values at occurrence locations

compute_boyce_robust <- function(fit_vals, obs_vals) {
  
  #------------------------------------------
  #----------- Basic checks -----------------
  #------------------------------------------
  #Define conditions for number of occurrences and total (or background) points
  if (length(obs_vals) < 5 || length(fit_vals) < 200) return(NA_real_)
  if (length(unique(fit_vals)) < 3)                   return(NA_real_)
  
  
  #------------------------------------------
  #---- Try moving window boyce ---
  #-----------------------------------------
  boyce_result <- try(ecospat::ecospat.boyce(fit = fit_vals, 
                                             obs = obs_vals,
                                             nclass = 0, #moving window boyce index
                                             PEplot = FALSE), 
                      silent = TRUE)
  
  
  #------------------------------------------
  #--- Return boyce index if all went well -----
  #------------------------------------------
  if (!inherits(boyce_result, "try-error") && !is.null(boyce_result$cor) && is.finite(boyce_result$cor)){
    return(boyce_result$cor)
    
  }else{
    
    #------------------------------------------
    #------------ Try binned Boyce ------------
    #------------------------------------------
    for (nc in c(10L, 20L)) {
      boyce_result2 <- try(
        ecospat::ecospat.boyce(
          fit    = fit_vals,
          obs    = obs_vals,
          nclass = nc,
          PEplot = FALSE
        ),
        silent = TRUE
      )
      
      if (!inherits(boyce_result2, "try-error") && !is.null(boyce_result2$cor) && is.finite(boyce_result2$cor)){
        return(boyce_result2$cor)
      }
    }
  }  
  #------------------------------------------
  #--Return NA if binned boyce also fails ---
  #------------------------------------------
  return(NA_real_)
  
}


#------------------------------------------
#---Calculate model validation metrics ----
#------------------------------------------
compute_validation_metrics <- function(species,
                                       type,
                                       region,
                                       fold,
                                       all_suit_vals,
                                       occ_suit_vals,
                                       abs_suit_vals) {
  all_suit_vals <- as.numeric(all_suit_vals)
  occ_suit_vals <- as.numeric(occ_suit_vals)
  abs_suit_vals <- as.numeric(abs_suit_vals)
  all_suit_vals <- all_suit_vals[is.finite(all_suit_vals)]
  occ_suit_vals <- occ_suit_vals[is.finite(occ_suit_vals)]
  abs_suit_vals <- abs_suit_vals[is.finite(abs_suit_vals)]

  n_pres <- length(occ_suit_vals)
  n_abs <- length(abs_suit_vals)
  empty_metrics <- function() {
    data.frame(
      Species = species,
      Type = type,
      Region = region,
      test_fold = fold,
      n_pres = as.numeric(n_pres),
      n_abs = as.numeric(n_abs),
      auc = NA_real_,
      boyce = NA_real_,
      tss = NA_real_,
      sens = NA_real_,
      spec = NA_real_,
      stringsAsFactors = FALSE
    )
  }

  if (n_pres < 2L || n_abs < 2L || length(all_suit_vals) < 2L) {
    return(empty_metrics())
  }
  
  #---------------
  #----- AUC -----
  #---------------
  # AUC via Wilcoxon statistic (exact, no ties correction needed for our use)
  auc_val <- tryCatch({
    wt <- wilcox.test(occ_suit_vals, abs_suit_vals, alternative = "greater", exact = FALSE)
    as.numeric(wt$statistic) / (n_pres * n_abs)
  }, error = function(e) NA_real_)
  
  
  #---------------
  #- Boyce index -
  #---------------
  boyce_val<-compute_boyce_robust(all_suit_vals, occ_suit_vals)
  
  
  #---------------
  #----- tss -----
  #---------------
  # max TSS: sweep candidate thresholds = unique sorted predicted values
  tss_val <- tryCatch({
    all_preds   <- c(occ_suit_vals, abs_suit_vals)
    all_labels  <- c(rep(1L, n_pres), rep(0L, n_abs))
    ord         <- order(all_preds, decreasing = TRUE)
    preds_s     <- all_preds[ord]
    labels_s    <- all_labels[ord]
    
    # Cumulative TP and FP as threshold descends
    tp_cum <- cumsum(labels_s == 1L)
    fp_cum <- cumsum(labels_s == 0L)
    sens   <- tp_cum / n_pres      # sensitivity 
    spec   <- 1 - fp_cum / n_abs   # specificity
    tss_v  <- sens + spec - 1
    
    # Identify the optimal threshold
    best_idx <- which.max(tss_v)
    
    list(TSS        = tss_v[best_idx],
         sensitivity = sens[best_idx],
         specificity = spec[best_idx],
         threshold   = preds_s[best_idx]
    )
  }, error = function(e) {
    list(TSS = NA_real_,
         sensitivity = NA_real_,
         specificity = NA_real_,
         threshold = NA_real_
    )
  })
  
  
  #------------------
  #- Return metrics -
  #------------------
  data.frame(
    Species = species,
    Type = type,
    Region = region,
    test_fold = fold,
    n_pres = as.numeric(n_pres),
    n_abs = as.numeric(n_abs),
    auc = as.numeric(auc_val),
    boyce = as.numeric(boyce_val),
    tss = as.numeric(tss_val$TSS),
    sens = as.numeric(tss_val$sensitivity),
    spec = as.numeric(tss_val$specificity),
    stringsAsFactors = FALSE
  )
 
}


#-----------------------------------------------------------------
#--Normalize supported spatial point objects to a SpatVector------
#-----------------------------------------------------------------
as_spatvector_safe <- function(x) {
  if (inherits(x, "SpatVector")) {
    return(x)
  }
  terra::vect(x)
}


#-----------------------------------------------------------------
#--Extract raster values at presence/absence points and return----
#-----------------------------------------------------------------
extract_env <- function(pres_abs_points, raster) {
  if (terra::nlyr(raster) == 0L) {
    stop("Cannot extract environmental values from an empty raster stack.",
         call. = FALSE)
  }
  if (!"species" %in% names(pres_abs_points)) {
    stop("Presence/absence points are missing the 'species' column.",
         call. = FALSE)
  }
  if (!"ID" %in% names(pres_abs_points)) {
    pres_abs_points$ID <- seq_len(nrow(pres_abs_points))
  }
  point_vector <- as_spatvector_safe(pres_abs_points)
  if (!isTRUE(terra::same.crs(point_vector, raster))) {
    point_vector <- terra::project(point_vector, raster)
  }
  env_values <- terra::extract(raster,
                               point_vector,
                               ID = FALSE,
                               xy = FALSE)
  if (nrow(env_values) != nrow(pres_abs_points) || ncol(env_values) == 0L) {
    stop(
      "Environmental extraction returned ", nrow(env_values), " row(s) and ",
      ncol(env_values), " predictor column(s) for ", nrow(pres_abs_points),
      " input point(s).",
      call. = FALSE
    )
  }
  df <- data.frame(
    species = pres_abs_points$species,
    ID = pres_abs_points$ID,
    env_values,
    check.names = FALSE
  )
  usable <- stats::complete.cases(df)
  numeric_predictors <- vapply(df[, -(1:2), drop = FALSE], is.numeric, logical(1))
  if (any(numeric_predictors)) {
    usable <- usable & apply(
      as.data.frame(df[, -(1:2), drop = FALSE][, numeric_predictors, drop = FALSE]),
      1,
      function(x) all(is.finite(x))
    )
  }
  dropped <- df[!usable, c("species", "ID"), drop = FALSE]
  df <- df[usable, , drop = FALSE]

  list(
    presences = df[df$species == 1, -1, drop = FALSE],
    absences = df[df$species == 0, -1, drop = FALSE],
    complete = df,
    dropped = dropped
  )
}


#-----------------------------------------------------------------
#--Calculate the geometric mean for ensemble validation------------
#-----------------------------------------------------------------
ensemble_geom_mean <- function(hab_df,
                               clim_df,
                               type,
                               value_col_hab = "median_favourability",
                               value_col_clim = "median_favourability",
                               return_data = FALSE) {
  required_hab <- c("ID", value_col_hab)
  required_clim <- c("ID", value_col_clim)
  if (length(setdiff(required_hab, names(hab_df))) > 0L ||
      length(setdiff(required_clim, names(clim_df))) > 0L) {
    stop("Ensemble inputs are missing stable IDs or favourability values.",
         call. = FALSE)
  }
  if (anyDuplicated(hab_df$ID) || anyDuplicated(clim_df$ID)) {
    stop("Ensemble inputs contain duplicated stable IDs.", call. = FALSE)
  }
  merged <- dplyr::inner_join(hab_df, clim_df, by = "ID", suffix = c("_hab", "_clim"))
  geom_mean <- sqrt(merged[[paste0(value_col_hab, "_hab")]] *
                      merged[[paste0(value_col_clim, "_clim")]])
  keep <- is.finite(geom_mean)
  merged <- merged[keep, , drop = FALSE]
  geom_mean <- geom_mean[keep]

  #Create a message for printing
  n_all <- length(unique(c(hab_df$ID, clim_df$ID)))
  n_merged <- nrow(merged)
  perc_retained <- if (n_all == 0L) 0 else 100 * n_merged / n_all
  n_lost <- n_all - n_merged

  
  message(sprintf(paste("Ensemble validation:", "%d European",type ,"points retained (%.0f%%).",
    "%d point(s) excluded due to missing predictions in the climate or habitat model."
  ), n_merged, perc_retained, n_lost))
  
  if (isTRUE(return_data)) {
    return(data.frame(ID = merged$ID, ensemble_favourability = geom_mean))
  }
  geom_mean
}


#-----------------------------------------------------------------
#--Put NO CV validation data in right format ----------
#-----------------------------------------------------------------
summarise_validation <- function(df, validation = NULL) {
  #Check that all necessary columns are present
  required_cols <- c("Species", "Type", "Region", "test_fold","auc", "boyce", "tss", "sens", "spec")
  missing_cols <- setdiff(required_cols, names(df))
  
  if (length(missing_cols) > 0) {
    stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
  }
  
  if (!"cv_method" %in% names(df)) {
    inferred_method <- if (identical(validation, "Cross-validation")) {
      "spatial_block"
    } else if (identical(validation, "No cross-validation")) {
      "no_cv"
    } else if (!is.null(validation)) {
      as.character(validation)
    } else {
      "not_evaluable"
    }
    df$cv_method <- inferred_method
  }
  defaults <- list(
    requested_folds = NA_integer_,
    effective_folds = NA_integer_,
    block_size_m = NA_real_,
    fallback_reason = NA_character_
  )
  for (column in names(defaults)) {
    if (!column %in% names(df)) {
      df[[column]] <- defaults[[column]]
    }
  }

  safe_mean <- function(x) {
    x <- x[is.finite(x)]
    if (length(x) == 0L) NA_real_ else mean(x)
  }
  safe_sd <- function(x) {
    x <- x[is.finite(x)]
    if (length(x) < 2L) NA_real_ else stats::sd(x)
  }

  df %>%
    dplyr::group_by(
      Species, Type, Region, cv_method, requested_folds,
      effective_folds, block_size_m, fallback_reason
    ) %>%
    dplyr::summarise(
      n_folds = if (dplyr::first(cv_method) %in% c("spatial_block", "stratified_kfold")) {
        dplyr::n_distinct(test_fold[!is.na(test_fold)])
      } else {
        0L
      },
      mean_auc = safe_mean(auc),
      sd_auc = safe_sd(auc),
      mean_boyce = safe_mean(boyce),
      sd_boyce = safe_sd(boyce),
      mean_tss = safe_mean(tss),
      sd_tss = safe_sd(tss),
      mean_sens = safe_mean(sens),
      sd_sens = safe_sd(sens),
      mean_spec = safe_mean(spec),
      sd_spec = safe_sd(spec),
      .groups = "drop"
    ) %>%
    dplyr::mutate(validation = cv_method)
}


#-----------------------------------------------------------------
#--Convert fortified summaries to the historical project schema--
#-----------------------------------------------------------------
as_legacy_validation_summary <- function(df) {
  legacy_columns <- c(
    "Species", "Type", "Region", "n_folds",
    "mean_auc", "sd_auc", "mean_boyce", "sd_boyce",
    "mean_tss", "sd_tss", "mean_sens", "sd_sens",
    "mean_spec", "sd_spec", "validation"
  )
  metric_columns <- setdiff(legacy_columns, "validation")
  missing_columns <- setdiff(metric_columns, names(df))
  if (length(missing_columns) > 0L) {
    stop(
      "Cannot create the legacy validation summary; missing column(s): ",
      paste(missing_columns, collapse = ", "),
      call. = FALSE
    )
  }

  method <- if ("cv_method" %in% names(df)) {
    as.character(df$cv_method)
  } else if ("validation" %in% names(df)) {
    as.character(df$validation)
  } else {
    rep(NA_character_, nrow(df))
  }
  validation_label <- dplyr::case_when(
    method %in% c("spatial_block", "stratified_kfold", "Cross-validation") ~
      "Cross-validation",
    method %in% c("no_cv", "No cross-validation") ~
      "No cross-validation",
    method %in% c("not_evaluable", "Not evaluable") ~
      "Not evaluable",
    TRUE ~ method
  )
  legacy_n_folds <- as.numeric(df$n_folds)
  legacy_n_folds[validation_label %in% c("No cross-validation", "Not evaluable")] <-
    NA_real_

  output <- df[, metric_columns, drop = FALSE]
  output$n_folds <- legacy_n_folds
  output$validation <- validation_label
  output[, legacy_columns, drop = FALSE]
}

#-----------------------------------------------------------------
#--Make species list string for the chunk runner--
#-----------------------------------------------------------------

rand_strings <- function(count = 10, n = 12, chars = c(letters, LETTERS, 0:9)) {
  replicate(count, paste0(sample(chars, n, replace = TRUE), collapse = ""))
}

make_species_list <- function(sp_name, chunk_nr, n_chunks){
  out<-rand_strings(n_chunks)
  out[chunk_nr] <- sp_name
  return(out)
}


