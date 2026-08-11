
library(readr)
library(tidyverse)

fl <- list.files("./data/projects", 
           recursive = TRUE,
           pattern = "Validation_summary.csv",
           full.names = TRUE)
val_summs <- list()
i<-0
for(file_path in fl){
  i<-i+1
  val_summs[[i]] <- suppressMessages(
    suppressWarnings(
      read_csv(file_path)))
}

val_summs_all <- bind_rows(val_summs)

write_csv(val_summs_all,"./data/post_model_outputs/validation_summaries_all-v1.csv")
write_rds(val_summs_all,"./data/post_model_outputs/validation_summaries_all-v1.rds")
