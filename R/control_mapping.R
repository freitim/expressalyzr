controls_map_is_configured <- function(controls_map) {
  !is.null(controls_map) &&
    !identical(controls_map, FALSE) &&
    !identical(controls_map, NA) &&
    length(controls_map) > 0
}

control_design_column <- function(config) {
  column <- config$controls_design_column
  if (is.null(column) || is.na(column) || !nzchar(column)) {
    column <- "Channel"
  }
  column
}

control_population_column <- function(config, design = NULL) {
  column <- config$controls_population_column
  if (!is.null(column) && !is.na(column) && nzchar(column)) {
    return(column)
  }

  candidates <- c("gated_population", "GatedPopulation", "Population",
                  "population", "control_population", "ControlPopulation")
  if (!is.null(design)) {
    match <- candidates[candidates %in% names(design)]
    if (length(match) > 0L) {
      return(match[[1L]])
    }
  }

  "gated_population"
}

sample_design_candidates <- function(sample_names) {
  lapply(sample_names, function(sample_name) {
    base <- basename(sample_name)
    sans_ext <- tools::file_path_sans_ext(base)

    well <- sub("^.*-([A-Z]{1}\\d{1,2})\\.fcs$", "\\1", base)
    if (identical(well, base)) {
      well <- sub("^.*([A-Z]{1}\\d{1,2}).*$", "\\1", sans_ext)
    }

    run_well <- sub("^0*(\\d{1,})-.*-([A-Z]{1}\\d{1,2})\\.fcs$",
                    "\\1-\\2", base)
    if (identical(run_well, base)) {
      run_well <- character()
    }

    unique(c(sample_name, base, sans_ext, run_well, well))
  })
}

load_control_design <- function(data_path, config) {
  if (is.null(config$spec_file) || is.na(config$spec_file)) {
    return(NULL)
  }

  design_path <- file.path(data_path, config$spec_file)
  if (!file.exists(design_path)) {
    return(NULL)
  }

  design <- data.table::fread(design_path)
  sample_id_col <- config$merge_by
  if (is.null(sample_id_col) || is.na(sample_id_col) ||
      !sample_id_col %in% names(design)) {
    return(NULL)
  }

  channel_col <- control_design_column(config)
  if (!channel_col %in% names(design) &&
      identical(channel_col, "Channel") &&
      "control_channel" %in% names(design)) {
    channel_col <- "control_channel"
  }

  if (!channel_col %in% names(design)) {
    return(NULL)
  }

  list(
    design = design,
    sample_id_col = sample_id_col,
    channel_col = channel_col,
    population_col = control_population_column(config, design)
  )
}

control_design_values <- function(design, channel_col) {
  values <- trimws(as.character(design[[channel_col]]))
  values[is.na(values)] <- ""
  values
}

filter_control_design_for_population <- function(control_design,
                                                 population_label = NULL) {
  if (is.null(control_design) || is.null(population_label) ||
      is.na(population_label) || !nzchar(population_label)) {
    return(control_design)
  }

  population_col <- control_design$population_col
  if (is.null(population_col) || is.na(population_col) ||
      !nzchar(population_col) ||
      !population_col %in% names(control_design$design)) {
    return(control_design)
  }

  design <- data.table::as.data.table(data.table::copy(control_design$design))
  channel_col <- control_design$channel_col
  if (!channel_col %in% names(design)) {
    return(control_design)
  }

  population_values <- trimws(as.character(design[[population_col]]))
  population_values[is.na(population_values)] <- ""
  control_values <- control_design_values(design, channel_col)
  matching_rows <- tolower(population_values) == tolower(population_label)

  specific_controls <- matching_rows & nzchar(control_values)
  if (any(specific_controls)) {
    design <- design[matching_rows]
  } else {
    design <- design[!nzchar(population_values)]
  }

  out <- control_design
  out$design <- design
  out
}

match_design_ids_to_samples <- function(sample_names, design_ids) {
  candidates <- sample_design_candidates(sample_names)

  vapply(design_ids, function(design_id) {
    hits <- which(vapply(candidates,
                         function(candidate) design_id %in% candidate,
                         logical(1)))
    if (length(hits) == 0) {
      return(NA_integer_)
    }
    if (length(hits) > 1) {
      stop("Design ID '", design_id, "' matches more than one FCS sample.",
           call. = FALSE)
    }
    hits
  }, integer(1))
}

identify_control_samples_from_design <- function(sample_names, control_design,
                                                 control_index = NULL) {
  design <- data.table::copy(control_design$design)
  sample_id_col <- control_design$sample_id_col
  channel_col <- control_design$channel_col

  control_values <- control_design_values(design, channel_col)
  design[, control_channel := control_values]
  design[, sample_match_index := match_design_ids_to_samples(
    sample_names,
    as.character(design[[sample_id_col]])
  )]

  if (!is.null(control_index)) {
    control_sample_indices <- which(control_index)
    control_design_dt <- design[sample_match_index %in% control_sample_indices]

    if (nrow(control_design_dt) == 0) {
      return(NULL)
    }

    missing_design_samples <- setdiff(control_sample_indices,
                                      control_design_dt$sample_match_index)
    if (length(missing_design_samples) > 0) {
      stop("Control samples matching controls_pattern are missing from the design file: ",
           paste(sample_names[missing_design_samples], collapse = ", "),
           call. = FALSE)
    }

    if (!any(nzchar(control_design_dt$control_channel))) {
      return(NULL)
    }

    missing_channel <- control_design_dt[!nzchar(control_channel)]
    if (nrow(missing_channel) > 0) {
      stop("Control rows in the design file need a non-empty ",
           channel_col, " value: ",
           paste(missing_channel[[sample_id_col]], collapse = ", "),
           call. = FALSE)
    }
  } else {
    control_rows <- which(nzchar(control_values))

    if (length(control_rows) == 0) {
      return(NULL)
    }

    control_design_dt <- design[control_rows]
  }

  missing_rows <- control_design_dt[is.na(sample_match_index)]
  if (nrow(missing_rows) > 0) {
    stop("Control rows in the design file do not match FCS samples: ",
         paste(missing_rows[[sample_id_col]], collapse = ", "),
         call. = FALSE)
  }

  duplicated_samples <- control_design_dt[
    duplicated(sample_match_index),
    unique(sample_match_index)
  ]
  if (length(duplicated_samples) > 0) {
    stop("More than one control design row matches the same FCS sample: ",
         paste(sample_names[duplicated_samples], collapse = ", "),
         call. = FALSE)
  }

  control_index <- rep(FALSE, length(sample_names))
  control_index[control_design_dt$sample_match_index] <- TRUE

  list(
    control_index = control_index,
    design = control_design_dt
  )
}

normalize_control_label <- function(x) {
  toupper(gsub("[^A-Za-z0-9]", "", trimws(as.character(x))))
}

strip_area_suffix <- function(x) {
  sub("(A|H|W|WIDTH)$", "", x)
}

parameter_channel_metadata <- function(cs, channels) {
  if (length(cs) == 0) {
    return(data.table::data.table(channel = character(), desc = character()))
  }

  parameters_dt <- data.table::as.data.table(
    Biobase::pData(flowCore::parameters(cs[[1]]))
  )
  if (!"desc" %in% names(parameters_dt)) {
    parameters_dt[, desc := NA_character_]
  }

  parameters_dt <- parameters_dt[parameters_dt[["name"]] %in% channels]
  data.table::data.table(
    channel = parameters_dt$name,
    desc = as.character(parameters_dt$desc)
  )
}

channel_alias_table <- function(channel_metadata) {
  alias_rows <- lapply(seq_len(nrow(channel_metadata)), function(i) {
    channel <- channel_metadata$channel[i]
    desc <- channel_metadata$desc[i]
    desc_values <- if (is.na(desc) || !nzchar(desc)) {
      character()
    } else {
      unlist(strsplit(desc, "[[:space:],;/|()]+"))
    }

    raw_aliases <- unique(c(channel, desc, desc_values))
    normalized <- unique(normalize_control_label(raw_aliases))
    normalized <- normalized[nzchar(normalized)]

    data.table::data.table(
      alias = unique(c(normalized, strip_area_suffix(normalized))),
      channel = channel
    )
  })

  data.table::rbindlist(alias_rows, use.names = TRUE)
}

resolve_control_value_to_channel <- function(value, channels, channel_aliases) {
  normalized <- normalize_control_label(value)
  if (!nzchar(normalized)) {
    return(NA_character_)
  }

  channel_hits <- channels[normalize_control_label(channels) == normalized]
  if (length(channel_hits) == 1) {
    return(channel_hits)
  }

  alias_hits <- unique(channel_aliases[alias == normalized, channel])
  if (length(alias_hits) == 0) {
    alias_hits <- unique(channel_aliases[alias == strip_area_suffix(normalized),
                                         channel])
  }

  if (length(alias_hits) == 1) {
    return(alias_hits)
  }

  if (length(alias_hits) > 1) {
    stop("Design control value '", value,
         "' matches more than one FCS channel: ",
         paste(alias_hits, collapse = ", "),
         call. = FALSE)
  }

  stop("Design control value '", value,
       "' does not match a compensation channel or FCS fluorophore label.",
       call. = FALSE)
}

resolve_control_indices_from_design <- function(sample_names, channels,
                                                matched_design,
                                                sample_id_col,
                                                channel_metadata = NULL) {
  if (length(channels) == 0) {
    stop("No compensation channels matched channel_pattern for the design controls.",
         call. = FALSE)
  }

  design <- data.table::copy(matched_design)
  design[, sample_match_index := match_design_ids_to_samples(
    sample_names,
    as.character(design[[sample_id_col]])
  )]

  design <- design[!is.na(sample_match_index)]
  design[, control_channel_normalized := tolower(trimws(control_channel))]

  unstained_values <- c("unstained", "negative", "empty", "empty vector",
                        "empty_vector", "none", "no stain", "no_stain")
  unstained_rows <- design[control_channel_normalized %in% unstained_values]

  if (nrow(unstained_rows) != 1) {
    stop("Design controls must include exactly one unstained/negative control.",
         call. = FALSE)
  }

  if (is.null(channel_metadata)) {
    channel_metadata <- data.table::data.table(
      channel = channels,
      desc = NA_character_
    )
  }
  channel_aliases <- channel_alias_table(channel_metadata)

  single_stain_rows <- design[!control_channel_normalized %in% unstained_values]
  single_stain_rows[, inferred_channel := vapply(
    control_channel,
    resolve_control_value_to_channel,
    character(1),
    channels = channels,
    channel_aliases = channel_aliases
  )]

  find_channel_control <- function(channel) {
    rows <- single_stain_rows[inferred_channel == channel]
    if (nrow(rows) == 0) {
      stop("Design controls are missing compensation control for ", channel,
           ".", call. = FALSE)
    }
    if (nrow(rows) > 1) {
      stop("Design controls include more than one compensation control for ",
           channel, ".", call. = FALSE)
    }
    rows$sample_match_index
  }

  channel_indices <- vapply(channels, find_channel_control, integer(1))
  controls_index <- c(unstained_rows$sample_match_index, channel_indices)

  duplicated_indices <- unique(controls_index[duplicated(controls_index)])
  if (length(duplicated_indices) > 0) {
    stop("Design controls assign the same FCS sample more than once: ",
         paste(sample_names[duplicated_indices], collapse = ", "),
         call. = FALSE)
  }

  matched_control_values <- single_stain_rows$control_channel[
    match(channels, single_stain_rows$inferred_channel)
  ]

  mapping <- data.table::data.table(
    role = c("unstained", rep("single_stain", length(channels))),
    channel = c(NA_character_, channels),
    control_value = c(unstained_rows$control_channel, matched_control_values),
    control_position = unname(controls_index),
    sample_name = sample_names[unname(controls_index)]
  )

  list(
    controls_index = unname(controls_index),
    mapping = mapping
  )
}

compensation_channel_names <- function(channels, channel_pattern) {
  fluorescence_channel_names(channels, channel_pattern)
}

resolve_control_indices <- function(sample_names, channels, controls_map) {
  if (!controls_map_is_configured(controls_map)) {
    stop("controls_map is not configured.", call. = FALSE)
  }

  if (length(channels) == 0) {
    stop("No compensation channels matched channel_pattern for controls_map.",
         call. = FALSE)
  }

  if (is.null(controls_map$unstained)) {
    stop("controls_map must include an 'unstained' entry.", call. = FALSE)
  }

  channel_map <- controls_map$channels
  if (is.null(channel_map)) {
    channel_map <- controls_map$single_color
  }

  if (is.null(channel_map) || length(channel_map) == 0) {
    stop("controls_map must include a non-empty 'channels' map.",
         call. = FALSE)
  }

  channel_controls <- unlist(channel_map, recursive = FALSE, use.names = TRUE)
  missing_channels <- setdiff(channels, names(channel_controls))
  if (length(missing_channels) > 0) {
    stop(
      paste0(
        "controls_map is missing compensation controls for: ",
        paste(missing_channels, collapse = ", "),
        ". Add these channels under controls_map$channels or narrow channel_pattern."
      ),
      call. = FALSE
    )
  }

  lookup_control <- function(value, label) {
    if (length(value) != 1) {
      stop("Control mapping for ", label, " must resolve to one value.",
           call. = FALSE)
    }

    if (is.numeric(value) || is.integer(value)) {
      index <- as.integer(value)
      if (is.na(index) || index < 1 || index > length(sample_names)) {
        stop(
          "Control mapping for ", label, " points to position ", value,
          ", but there are only ", length(sample_names), " detected controls.",
          call. = FALSE
        )
      }
      return(index)
    }

    if (is.character(value)) {
      exact_hits <- which(sample_names == value)
      if (length(exact_hits) == 1) {
        return(exact_hits)
      }

      pattern_hits <- grep(value, sample_names, ignore.case = TRUE)
      if (length(pattern_hits) == 1) {
        return(pattern_hits)
      }

      if (length(pattern_hits) == 0) {
        stop("Control mapping for ", label, " did not match any detected controls: ",
             value, call. = FALSE)
      }

      stop("Control mapping for ", label, " matched more than one detected control: ",
           value, call. = FALSE)
    }

    stop("Control mapping for ", label,
         " must be a numeric position or character sample pattern.",
         call. = FALSE)
  }

  unstained_index <- lookup_control(controls_map$unstained, "unstained")
  channel_indices <- vapply(
    channels,
    function(channel) lookup_control(channel_controls[[channel]], channel),
    integer(1)
  )

  controls_index <- c(unstained = unstained_index, channel_indices)
  duplicated_indices <- unique(controls_index[duplicated(controls_index)])
  if (length(duplicated_indices) > 0) {
    stop(
      "controls_map assigns the same detected control position more than once: ",
      paste(duplicated_indices, collapse = ", "),
      call. = FALSE
    )
  }

  mapping <- data.table::data.table(
    role = c("unstained", rep("single_stain", length(channels))),
    channel = c(NA_character_, channels),
    control_value = as.character(c(controls_map$unstained,
                                   unname(channel_controls[channels]))),
    control_position = unname(controls_index),
    sample_name = sample_names[unname(controls_index)]
  )

  list(
    controls_index = unname(controls_index),
    mapping = mapping
  )
}
