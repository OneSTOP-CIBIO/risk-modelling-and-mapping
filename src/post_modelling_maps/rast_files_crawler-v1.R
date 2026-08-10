# Raster path crawler for the post-modelling SDM outputs.
#
# Examples:
# crawl_raster_paths(
#   "Acacia mearnsii De Wild.", "binary", "hist", "hist"
# )
# crawl_raster_paths(
#   c("Acacia mearnsii De Wild.", "Ailanthus altissima (Mill.) Swingle"),
#   "favourability", "mid_2041_2070", "ssp126"
# )
# crawl_raster_paths(
#   "all", "binary", "late_2071_2100", "ssp585"
# )
# crawl_raster_paths(
#   "Acacia mearnsi De Wild.", "binary", "hist", "hist", fuzzy = TRUE
# )

best_mod_set_paths <- readRDS(
  file.path("data", "processed", "best_mod_set_paths.rds")
)

raster_paths_mark_utf8 <- function(x) {
  x <- as.character(x)
  valid <- validUTF8(x)
  if (any(valid)) {
    Encoding(x[valid]) <- "UTF-8"
  }
  x
}

raster_paths_normalize_species <- function(x) {
  tolower(raster_paths_mark_utf8(trimws(x)))
}

raster_paths_validate_scalar_choice <- function(value, argument, choices) {
  if (!is.character(value) ||
      length(value) != 1L ||
      is.na(value) ||
      !nzchar(value) ||
      !value %in% choices) {
    stop(
      "`",
      argument,
      "` must be one of: ",
      paste(sprintf('"%s"', choices), collapse = ", "),
      ".",
      call. = FALSE
    )
  }
  invisible(value)
}

raster_paths_validate_catalog <- function(catalog) {
  required_fields <- c(
    "sp_name",
    "tkey",
    "paths_list_bin",
    "paths_list_fav"
  )

  if (!is.list(catalog) || length(catalog) == 0L) {
    stop("`catalog` must be a non-empty list.", call. = FALSE)
  }

  for (i in seq_along(catalog)) {
    entry <- catalog[[i]]
    if (!is.list(entry)) {
      stop("Catalog entry ", i, " must be a list.", call. = FALSE)
    }

    missing_fields <- setdiff(required_fields, names(entry))
    if (length(missing_fields) > 0L) {
      stop(
        "Catalog entry ",
        i,
        " is missing required field(s): ",
        paste(missing_fields, collapse = ", "),
        ".",
        call. = FALSE
      )
    }

    if (!is.character(entry$sp_name) ||
        length(entry$sp_name) != 1L ||
        is.na(entry$sp_name) ||
        !nzchar(trimws(entry$sp_name))) {
      stop(
        "Catalog entry ",
        i,
        " must contain one non-empty `sp_name`.",
        call. = FALSE
      )
    }
    if (length(entry$tkey) != 1L || is.na(entry$tkey)) {
      stop(
        "Catalog entry ",
        i,
        " must contain one non-missing `tkey`.",
        call. = FALSE
      )
    }
    if (!is.list(entry$paths_list_bin) ||
        !is.list(entry$paths_list_fav)) {
      stop(
        "Catalog entry ",
        i,
        " must contain list-valued path fields.",
        call. = FALSE
      )
    }
  }

  species_names <- vapply(
    catalog,
    function(entry) entry$sp_name,
    character(1)
  )
  normalized_names <- raster_paths_normalize_species(species_names)
  duplicated_names <- unique(normalized_names[duplicated(normalized_names)])
  if (length(duplicated_names) > 0L) {
    stop(
      "`catalog` contains duplicate species names after normalization.",
      call. = FALSE
    )
  }

  invisible(species_names)
}

raster_paths_format_suggestions <- function(names, distances, limit = 3L) {
  ranked <- order(distances)
  ranked <- ranked[seq_len(min(length(ranked), limit))]
  paste(
    sprintf("%s (distance %.3f)", names[ranked], distances[ranked]),
    collapse = ", "
  )
}

raster_paths_match_species <- function(species, catalog_names, fuzzy,
                                       max_fuzzy_distance) {
  normalized_catalog <- raster_paths_normalize_species(catalog_names)
  normalized_queries <- raster_paths_normalize_species(species)
  matched_indices <- integer(length(species))

  for (i in seq_along(species)) {
    query <- normalized_queries[[i]]
    exact <- which(normalized_catalog == query)
    if (length(exact) == 1L) {
      matched_indices[[i]] <- exact
      next
    }

    raw_distances <- as.numeric(utils::adist(
      query,
      normalized_catalog,
      partial = FALSE,
      ignore.case = FALSE
    ))
    denominators <- pmax(
      nchar(query, type = "chars"),
      nchar(normalized_catalog, type = "chars")
    )
    distances <- raw_distances / denominators

    if (!fuzzy) {
      stop(
        "No exact species match for \"",
        species[[i]],
        "\". Closest catalog names: ",
        raster_paths_format_suggestions(catalog_names, distances),
        ". Set `fuzzy = TRUE` to allow conservative fuzzy matching.",
        call. = FALSE
      )
    }

    candidates <- which(distances <= max_fuzzy_distance)
    if (length(candidates) == 0L) {
      stop(
        "No fuzzy species match for \"",
        species[[i]],
        "\" within `max_fuzzy_distance = ",
        format(max_fuzzy_distance),
        "`. Closest catalog names: ",
        raster_paths_format_suggestions(catalog_names, distances),
        ".",
        call. = FALSE
      )
    }
    if (length(candidates) > 1L) {
      candidate_order <- order(distances[candidates])
      candidates <- candidates[candidate_order]
      stop(
        "Ambiguous fuzzy species match for \"",
        species[[i]],
        "\": ",
        paste(
          sprintf(
            "%s (distance %.3f)",
            catalog_names[candidates],
            distances[candidates]
          ),
          collapse = ", "
        ),
        ". Use a fuller name or lower `max_fuzzy_distance`.",
        call. = FALSE
      )
    }

    matched_indices[[i]] <- candidates
  }

  matched_indices[!duplicated(matched_indices)]
}

#' Retrieve raster paths from the best-model catalog.
#'
#' @param species A character vector of complete species names, or "all".
#' @param map_type Either "binary" or "favourability".
#' @param period One of "hist", "mid_2041_2070", or "late_2071_2100".
#' @param scenario "hist" for historical data, otherwise one of "ssp126",
#'   "ssp370", or "ssp585".
#' @param fuzzy Whether to use conservative fuzzy species-name matching.
#' @param max_fuzzy_distance Maximum normalized Levenshtein distance in [0, 1].
#' @param catalog The species raster-path catalog.
#'
#' @return A data frame with one row per selected species and raster path.
crawl_raster_paths <- function(species,
                               map_type,
                               period,
                               scenario,
                               fuzzy = FALSE,
                               max_fuzzy_distance = 0.10,
                               catalog = best_mod_set_paths) {
  map_types <- c("binary", "favourability")
  periods <- c("hist", "mid_2041_2070", "late_2071_2100")
  future_scenarios <- c("ssp126", "ssp370", "ssp585")
  scenarios <- c("hist", future_scenarios)

  if (!is.character(species) ||
      length(species) == 0L ||
      anyNA(species) ||
      any(!nzchar(trimws(species)))) {
    stop(
      "`species` must be a non-empty character vector without missing or empty names.",
      call. = FALSE
    )
  }
  if (!is.logical(fuzzy) || length(fuzzy) != 1L || is.na(fuzzy)) {
    stop("`fuzzy` must be TRUE or FALSE.", call. = FALSE)
  }
  if (!is.numeric(max_fuzzy_distance) ||
      length(max_fuzzy_distance) != 1L ||
      is.na(max_fuzzy_distance) ||
      !is.finite(max_fuzzy_distance) ||
      max_fuzzy_distance < 0 ||
      max_fuzzy_distance > 1) {
    stop(
      "`max_fuzzy_distance` must be one finite number between 0 and 1.",
      call. = FALSE
    )
  }

  raster_paths_validate_scalar_choice(map_type, "map_type", map_types)
  raster_paths_validate_scalar_choice(period, "period", periods)
  raster_paths_validate_scalar_choice(scenario, "scenario", scenarios)

  if (identical(period, "hist") && !identical(scenario, "hist")) {
    stop(
      'Historical period "hist" requires scenario "hist".',
      call. = FALSE
    )
  }
  if (!identical(period, "hist") && identical(scenario, "hist")) {
    stop(
      'Future periods require one of: "ssp126", "ssp370", or "ssp585".',
      call. = FALSE
    )
  }

  catalog_names <- raster_paths_validate_catalog(catalog)
  normalized_species <- raster_paths_normalize_species(species)
  all_requested <- normalized_species == "all"
  if (any(all_requested) && length(species) != 1L) {
    stop(
      '`species = "all"` cannot be combined with explicit species names.',
      call. = FALSE
    )

  }

  selected_indices <- if (any(all_requested)) {
    seq_along(catalog)
  } else {
    raster_paths_match_species(
      species,
      catalog_names,
      fuzzy,
      max_fuzzy_distance
    )
  }

  path_field <- switch(
    map_type,
    binary = "paths_list_bin",
    favourability = "paths_list_fav"
  )

  rows <- lapply(selected_indices, function(index) {
    entry <- catalog[[index]]
    period_paths <- entry[[path_field]][[period]]
    if (!is.list(period_paths)) {
      stop(
        "Catalog entry for \"",
        entry$sp_name,
        "\" has no list-valued path data for period \"",
        period,
        "\".",
        call. = FALSE
      )
    }

    path <- period_paths[[scenario]]
    if (!is.character(path) ||
        length(path) != 1L ||
        is.na(path) ||
        !nzchar(path)) {
      stop(
        "Catalog entry for \"",
        entry$sp_name,
        "\" has no valid ",
        map_type,
        " path for ",
        period,
        "/",
        scenario,
        ".",
        call. = FALSE
      )
    }

    path <- raster_paths_mark_utf8(path)
    data.frame(
      sp_name = raster_paths_mark_utf8(entry$sp_name),
      tkey = entry$tkey,
      map_type = map_type,
      period = period,
      scenario = scenario,
      path = path,
      file_exists = file.exists(path),
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
  })

  result <- do.call(rbind, rows)
  rownames(result) <- NULL
  result
}
