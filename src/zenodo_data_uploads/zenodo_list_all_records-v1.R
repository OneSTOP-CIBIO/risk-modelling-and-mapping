

library(zen4R)

library(dplyr)
library(tidyverse)


## 1. Connect to Zenodo -------------------------------------------------
token <- "4Dqs063GQX1lFAzfPRA7NAKbGtBOpT8veuqJ46UVuNaVyEJxx3s7cgmSLDRX"

zen <- ZenodoManager$new(
  token  = token,
  logger = "INFO"   # or "DEBUG" for very verbose logs, or NULL for quiet
)

## 2. Get all your depositions (workspace records) ----------------------
# This returns a list of ZenodoRecord objects (drafts + your deposited records)
deps <- zen$getDepositions(
  q            = "",     # no filter: all your depositions
  size         = 200,    # increase if needed
  all_versions = TRUE,
  exact        = FALSE,
  quiet        = TRUE
)

## 3. Build a simplified data.frame -------------------------------------

if (length(deps) == 0) {
  dep_df <- data.frame(
    id    = character(0),
    title = character(0),
    url   = character(0),
    doi   = character(0),
    stringsAsFactors = FALSE
  )
} else {
  dep_df <- data.frame(
    id = vapply(deps, function(x) x$id, character(1)),
    
    title = vapply(
      deps,
      function(x) if (!is.null(x$metadata$title)) x$metadata$title else NA_character_,
      character(1)
    ),
    
    # Try html link first; fall back to self if html is missing
    url = vapply(
      deps,
      function(x) {
        if (!is.null(x$links$html)) {
          x$links$html
        } else if (!is.null(x$links$self)) {
          x$links$self
        } else {
          NA_character_
        }
      },
      character(1)
    ),
    
    doi = vapply(
      deps,
      function(x) if (!is.null(x$getConceptDOI())) x$getConceptDOI() else NA_character_,
      character(1)
    ),
    
    stringsAsFactors = FALSE
  )
}

dep_df <- dep_df %>% mutate(url = gsub("api/","", url))

## 4. Inspect the result -------------------------------------------------
sel_entries <- dep_df[1:8,]
# View(dep_df)  # in RStudio, if you want a spreadsheet-like view

readr::write_csv(sel_entries,"C:/Users/JG/Desktop/ZenodoRecords_OneSTOP-T5.1-data-v1.csv")
