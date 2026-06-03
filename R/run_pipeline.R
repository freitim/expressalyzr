#' Run the data processing pipeline
#'
#' Run the data processing pipeline that links together other functions in the
#' expressalyzr package.
#'
#' @param data_path Path to the experiment data directory.
#' @param view_config Whether to open the config file for viewing.
#' @param gating_output Legacy compatibility option. Use \code{NULL} for the
#'   default run or \code{"report"} as an alias for \code{report = TRUE}.
#'   Legacy \code{"inspect"} and \code{"set"} values fail clearly.
#' @param interactive Whether the pipeline may pause for user-facing setup or
#'   review steps. When \code{TRUE}, missing starter files can be created,
#'   \code{config.yml} can be opened for review, the gate designer can launch,
#'   and manual compensation/cutoff prompts may wait for input. When
#'   \code{FALSE}, no editors, browser apps, or console prompts are allowed;
#'   required setup files must already exist.
#' @param gate_designer If \code{TRUE}, launch
#'   \code{\link{design_gating_template}} before applying gates. This is
#'   additive and requires \code{interactive = TRUE}.
#' @param gating_templates Optional named character vector/list of OpenCyto
#'   template files to apply independently. Paths are resolved relative to
#'   \code{data_path}. In interactive runs, missing explicit template files are
#'   created from the starter \code{gt_samples.csv} template. When \code{NULL},
#'   the legacy \code{gt_samples.csv} workflow is used.
#' @param gated_populations Optional character vector of biological population
#'   labels to pair with \code{gating_templates}. When omitted, names on
#'   \code{gating_templates} are used, falling back to template file names.
#' @param report If \code{TRUE}, generate one comprehensive HTML report with QC,
#'   gating, channel, compensation, and summary sections.
#' @param include_gating_plots If \code{TRUE}, include per-sample gate scatter
#'   grids in the unified HTML report. These plots are useful for gate review
#'   but can make reports large for plate-scale experiments.
#' @param summaries If \code{TRUE}, write additive population-level and
#'   marker/channel-level summary CSVs. Legacy outputs and return values are
#'   unchanged.
#' @param save_gating_set If \code{TRUE}, save the applied
#'   \code{flowWorkspace::GatingSet}. In multi-template mode, one GatingSet is
#'   saved per template.
#' @param config_profile Starter config profile to copy when
#'   \code{config.yml} does not exist. Existing config files are loaded
#'   unchanged. One of \code{"default"}, \code{"legacy"}, or
#'   \code{"advanced"}.
#' @param ... Reserved. Deprecated options such as \code{qc},
#'   \code{qc_output}, \code{gating_report}, \code{channel_report},
#'   \code{provenance}, \code{plots}, and \code{use_gui} are no longer
#'   accepted.
#'
#' @return The event-level \code{data.table}. The same table is always written
#'   to \code{<experiment>.csv} in \code{data_path}. Report, summary, and
#'   GatingSet options add files on disk without changing the return type.
#' @export
run_pipeline <- function(data_path, view_config = TRUE, gating_output = NULL,
                         interactive = TRUE,
                         gate_designer = TRUE,
                         gating_templates = NULL,
                         gated_populations = NULL,
                         report = TRUE, include_gating_plots = TRUE,
                         summaries = FALSE,
                         save_gating_set = FALSE,
                         config_profile = c("default", "legacy", "advanced"),
                         ...) {
  data_path <- normalizePath(data_path, mustWork = FALSE)
  config_profile <- match.arg(config_profile)
  extra_args <- list(...)
  removed_args <- intersect(
    names(extra_args),
    c("qc", "qc_output", "gating_report", "channel_report",
      "provenance", "plots", "use_gui")
  )
  if (length(removed_args) > 0L) {
    stop(
      "run_pipeline() no longer accepts: ",
      paste(removed_args, collapse = ", "),
      ". Use report = TRUE for QC/report output, summaries = TRUE for ",
      "population/marker CSVs, save_gating_set = TRUE for GatingSets, ",
      "or gate_designer = TRUE for interactive gate editing.",
      call. = FALSE
    )
  }
  unknown_args <- setdiff(names(extra_args), removed_args)
  if (length(unknown_args) > 0L) {
    stop("Unknown run_pipeline() argument(s): ",
         paste(unknown_args, collapse = ", "), call. = FALSE)
  }

  qc <- isTRUE(report)

  if (!is.null(gating_output) &&
      !gating_output %in% c("inspect", "report", "set")) {
    stop("gating_output must be one of NULL, 'report', 'inspect', or 'set'. ",
         "Use report = TRUE for the current reporting workflow.",
         call. = FALSE)
  }

  if (identical(gating_output, "inspect") || identical(gating_output, "set")) {
    stop(
      "gating_output = '", gating_output, "' is no longer supported. ",
      "Use report = TRUE for reviewable gate plots, gate_designer = TRUE ",
      "to inspect/edit gates interactively, or save_gating_set = TRUE ",
      "to save reusable GatingSets.",
      call. = FALSE
    )
  }

  if (identical(gating_output, "report")) {
    report <- TRUE
    qc <- TRUE
  }

  if (!interactive && isTRUE(gate_designer)) {
    stop("gate_designer = TRUE requires interactive = TRUE.", call. = FALSE)
  }

  # initialize
  experiment_name <- basename(data_path)

  create_data_subdir(data_path)

  config_file_path <- file.path(data_path, "config.yml")
  config <- load_config(config_file_path, view_config,
                        interactive = interactive,
                        config_profile = config_profile)

  gt_file <- file.path(data_path, "gt_samples.csv")
  if (is.null(gating_templates)) {
    ensure_default_gating_template(gt_file, interactive = interactive)
  }
  template_specs <- resolve_pipeline_gating_templates(
    data_path,
    default_template = gt_file,
    gating_templates = gating_templates,
    gated_populations = gated_populations,
    interactive = interactive
  )
  if (isTRUE(config$adjust_gating)) {
    if (!interactive) {
      stop("config$adjust_gating = TRUE requires interactive = TRUE.",
           call. = FALSE)
    }

    utils::file.edit(gt_file)

    cat("\n")
    readline(prompt = "Press [Enter] to continue.")
  }

  if (isTRUE(gate_designer)) {
    template_labels <- ifelse(
      is.na(template_specs$gated_population),
      basename(template_specs$template_file),
      paste0(template_specs$gated_population, " (",
             template_specs$gating_template, ")")
    )
    message("Launching combined gate designer for ",
            nrow(template_specs), " template(s).")
    design_gating_template(
      data_path = data_path,
      output_file = template_specs$template_file[[1L]],
      template_file = template_specs$template_file[[1L]],
      template_files = template_specs$template_file,
      output_files = template_specs$template_file,
      template_labels = template_labels,
      launch = TRUE
    )
  }

  # start analysis
  cs <- load_fcs(data_path)

  channel_map <- config_channel_map(config)
  cs <- normalize_channels(cs, channel_map)
  chs <- flowWorkspace::colnames(cs)

  bead_index <- grepl(config$beads_pattern, flowCore::sampleNames(cs))

  if (sum(bead_index) == 1) {

    cs_beads <- cs[bead_index]
    cs <- cs[!bead_index]

    if (config$mefl_transform) {
      t_fun <- generate_transformation(cs_beads,
                                       channel_pattern = config$channel_pattern)
      trans <- TRUE
    } else {
      trans <- FALSE
    }
  } else if (sum(bead_index) > 1) {

    stop("More than one bead sample found in the dataset. Please restrict ",
         "beads_pattern to one sample.", call. = FALSE)

  } else {

    message("Bead sample not found. Values will not be transformed to MEFL.")
    trans <- FALSE
  }

  # gating
  openCyto::register_plugins(fun = density_gate, "density_gate", dep = NA, "gating")
  openCyto::register_plugins(fun = mixture_gate, "mixture_gate", dep = NA, "gating")
  openCyto::register_plugins(fun = polygon_gate, "polygon_gate", dep = NA, "gating")

  multi_template_mode <- isTRUE(any(template_specs$explicit))
  fresh_analysis_cs <- function() {
    fresh_cs <- normalize_channels(load_fcs(data_path), channel_map)
    fresh_bead_index <- grepl(config$beads_pattern,
                              flowCore::sampleNames(fresh_cs))
    if (sum(fresh_bead_index) == 1) {
      fresh_cs <- fresh_cs[!fresh_bead_index]
    } else if (sum(fresh_bead_index) > 1) {
      stop("More than one bead sample found in the dataset. Please restrict ",
           "beads_pattern to one sample.", call. = FALSE)
    }
    fresh_cs
  }

  gating_results <- vector("list", nrow(template_specs))
  for (i in seq_len(nrow(template_specs))) {
    cs_for_gating <- if (multi_template_mode) fresh_analysis_cs() else cs
    gating_results[[i]] <- run_gating_template(
      cs_for_gating,
      template_specs$template_file[[i]],
      extraction_population = template_specs$extraction_population[[i]]
    )
    gating_results[[i]]$label <- template_specs$gating_template[[i]]
    gating_results[[i]]$population_label <- template_specs$gated_population[[i]]
    gating_results[[i]]$extraction_population <- template_specs$extraction_population[[i]]
    gating_results[[i]]$template_file <- template_specs$template_file[[i]]
  }
  gs <- gating_results[[1L]]$gs

  gating_artifact_label <- function(result) make.names(result$label)

  gating_set_path <- function(result) {
    if (multi_template_mode) {
      file.path(data_path, paste0("gs_", gating_artifact_label(result)))
    } else {
      file.path(data_path, "gs")
    }
  }

  if (isTRUE(save_gating_set)) {
    for (result in gating_results) {
      flowWorkspace::save_gs(result$gs, path = gating_set_path(result))
    }
  }

  control_design <- load_control_design(data_path, config)
  process_extracted_population <- function(data_cs, template_label = NULL,
                                           population_label = NULL) {
    if (trans) {
      data_cs <- apply_transform(data_cs, t_fun,
                                 channel_pattern = config$channel_pattern)
    }

    cont_index <- grepl(config$controls_pattern, flowCore::sampleNames(data_cs))
    design_controls <- NULL
    population_control_design <- filter_control_design_for_population(
      control_design,
      population_label
    )
    if (!is.null(control_design) && any(cont_index)) {
      design_controls <- identify_control_samples_from_design(
        flowCore::sampleNames(data_cs),
        population_control_design,
        control_index = cont_index
      )
    }

    n_controls <- sum(cont_index)
    so_mat <- NULL
    comp_controls_index <- config$controls_index
    comp_cont_cs <- NULL
    control_mapping <- NULL
    cont_cs <- NULL
    mapping_suffix <- if (multi_template_mode && !is.na(template_label)) {
      paste0("_", make.names(template_label))
    } else {
      ""
    }

    if (n_controls > 0) {
      cont_cs <- data_cs[cont_index]
      cont_names <- flowCore::sampleNames(cont_cs)
      cont_order <- gtools::mixedorder(gsub("^.*([A-Z]{0,1}\\d{1,2}).fcs",
                                            "\\1", cont_names))
      cont_cs <- cont_cs[cont_order]

      comp_cont_cs <- cont_cs
      if (!is.null(design_controls)) {
        comp_channels <- compensation_channel_names(flowCore::colnames(cont_cs),
                                                    config$channel_pattern)
        channel_metadata <- parameter_channel_metadata(cont_cs, comp_channels)
        control_mapping <- resolve_control_indices_from_design(
          flowCore::sampleNames(cont_cs),
          comp_channels,
          design_controls$design,
          population_control_design$sample_id_col,
          channel_metadata = channel_metadata
        )
        comp_controls_index <- control_mapping$controls_index
        data.table::fwrite(
          control_mapping$mapping,
          file = file.path(data_path,
                           paste0(experiment_name, mapping_suffix,
                                  "_compensation_controls.csv"))
        )
        comp_cont_cs <- cont_cs[comp_controls_index]
        comp_controls_index <- seq_along(comp_controls_index)
      } else if (controls_map_is_configured(config$controls_map)) {
        comp_channels <- compensation_channel_names(flowCore::colnames(cont_cs),
                                                    config$channel_pattern)
        control_mapping <- resolve_control_indices(
          flowCore::sampleNames(cont_cs),
          comp_channels,
          config$controls_map
        )
        comp_controls_index <- control_mapping$controls_index
        data.table::fwrite(
          control_mapping$mapping,
          file = file.path(data_path,
                           paste0(experiment_name, mapping_suffix,
                                  "_compensation_controls.csv"))
        )
        comp_cont_cs <- cont_cs[comp_controls_index]
        comp_controls_index <- seq_along(comp_controls_index)
      }
    }

    if (n_controls > 2) {
      so_mat_path <- if (multi_template_mode && !is.na(template_label)) {
        file.path(data_path, paste0("so_mat_", make.names(template_label),
                                    ".RData"))
      } else {
        file.path(data_path, "so_mat.RData")
      }

      if (!file.exists(so_mat_path) || config$redo_comp) {
        if (isTRUE(config$manual_comp) && !interactive) {
          stop("config$manual_comp = TRUE requires interactive = TRUE.",
               call. = FALSE)
        }

        so_mat <- spillover_matrix(comp_cont_cs,
                                   comp_controls_index,
                                   config$channel_pattern,
                                   config$density_th,
                                   config$manual_comp,
                                   interactive = interactive)

        save(so_mat, file = so_mat_path)
      } else {
        load(so_mat_path)
      }

      data_cs <- flowWorkspace::compensate(data_cs, so_mat)
    } else {
      message("No control samples found. Proceeding without compensation.")
    }

    data_dt <- cs_to_dt(data_cs)

    chs <- fluorescence_channel_names(colnames(data_dt),
                                      config$channel_pattern)
    value_scale <- if (trans) {
      "mefl"
    } else if (n_controls > 2) {
      "compensated"
    } else {
      "raw"
    }

    if (length(chs) > 0L) {
      data_dt[, no_negative := rowSums(data_dt[, chs, with = FALSE] <= 0) == 0]
      data_dt[, positive := rowSums(data_dt[, chs, with = FALSE] < 0) == 0]
    } else {
      data_dt[, no_negative := TRUE]
      data_dt[, positive := TRUE]
    }

    if (n_controls > 0) {
      cont_dt <- cs_to_dt(cont_cs)
      negative_control_index <- if (is.null(control_mapping)) {
        config$controls_index[1]
      } else {
        control_mapping$controls_index[1]
      }
      neg_dt <- cs_to_dt(cont_cs[negative_control_index])

      if (isTRUE(config$manual_cutoff)) {
        if (!interactive) {
          stop("config$manual_cutoff = TRUE requires interactive = TRUE.",
               call. = FALSE)
        }

        cont_sub_dt <- cont_dt[, c("File", chs), with = FALSE]
        config$bg_cutoff <- adjust_threshold(cont_sub_dt, config$bg_cutoff)
        write_value <- paste("bg_cutoff:", config$bg_cutoff)
        write_config(config_file_path, "default", write_value)
      }

      neg_dt[, lapply(mget(chs), quantile, config$bg_cutoff),
             by = .(File)]

      bg_dt <- neg_dt[, lapply(mget(chs), quantile, config$bg_cutoff)]

      f_pos <- function(channel, values) values >= bg_dt[[channel]]

      pos_chs <- paste0(chs, "_pos")
      data_dt[, (pos_chs) := lapply(chs, function(ch) f_pos(ch, get(ch))),
              by = .(File)]
    }

    if (!is.null(config$bg_channels)) {
      data_dt[(positive), (paste0(config$bg_channels, "_bg")) :=
                lapply(mget(config$bg_channels),
                       assign_bg,
                       n_comp = NULL,
                       rm = 1,
                       inspect = FALSE),
              by = .(File)]
    }

    if (multi_template_mode && !is.na(template_label)) {
      data_dt[, gating_template := template_label]
      data_dt[, gated_population := population_label]
    }

    list(
      data_cs = data_cs,
      data_dt = data_dt,
      chs = chs,
      value_scale = value_scale,
      so_mat = so_mat,
      n_controls = n_controls
    )
  }

  processed_results <- lapply(gating_results, function(result) {
    process_extracted_population(result$data_cs, result$label,
                                 result$population_label)
  })
  data_cs <- processed_results[[1L]]$data_cs
  data_dt <- data.table::rbindlist(
    lapply(processed_results, `[[`, "data_dt"),
    use.names = TRUE,
    fill = TRUE
  )
  chs <- unique(unlist(lapply(processed_results, `[[`, "chs"),
                       use.names = FALSE))
  chs <- chs[chs %in% colnames(data_dt)]
  value_scales <- unique(vapply(processed_results, `[[`, character(1L),
                                "value_scale"))
  value_scale <- if (length(value_scales) == 1L) value_scales else "mixed"
  so_mat <- processed_results[[1L]]$so_mat

  qc_results <- NULL
  if (qc) {
    if (multi_template_mode) {
      qc_by_template <- Map(function(gating_result, processed_result) {
        collect_qc(gating_result$gs, processed_result$data_cs,
                   processed_result$data_dt, processed_result$chs)
      }, gating_results, processed_results)

      tag_qc <- function(tbl, result) {
        if (is.null(tbl)) {
          return(NULL)
        }
        tbl <- data.table::copy(tbl)
        tbl[, gating_template := result$label]
        tbl[, gated_population := result$population_label]
        tbl
      }

      combine_qc <- function(metric) {
        tables <- Map(function(qc_item, result) {
          tag_qc(qc_item[[metric]], result)
        }, qc_by_template, gating_results)
        tables <- Filter(Negate(is.null), tables)
        if (length(tables) == 0L) {
          return(NULL)
        }
        data.table::rbindlist(tables, use.names = TRUE, fill = TRUE)
      }

      qc_results <- list(
        event_counts  = combine_qc("event_counts"),
        scatter_stats = combine_qc("scatter_stats"),
        channel_stats = combine_qc("channel_stats"),
        time_check    = combine_qc("time_check")
      )
    } else {
      qc_results <- collect_qc(gs, data_cs, data_dt, chs)
    }
  }

  # assign experimental specifications
  s_file_path <- file.path(data_path, config$spec_file)

  if (file.exists(s_file_path)) {
    s_file <- data.table::fread(s_file_path)
    s_file[, (config$merge_by) := gsub("^0(\\d{1})(.*)$", "\\1\\2", get(config$merge_by))]
    data_dt[, (config$merge_by) := gsub("^0(\\d{1})-.*-([A-Z]{1}\\d{1,2})\\.fcs$", "\\1-\\2", File)]
    data_dt <- merge(data_dt, s_file, by = config$merge_by, all.x = TRUE)
  }

  data.table::fwrite(data_dt, file = file.path(data_path, paste0(experiment_name, ".csv")))

  population_table <- NULL
  marker_table <- NULL
  if (summaries || report) {
    if (multi_template_mode) {
      population_table <- data.table::rbindlist(
        Map(function(result) {
          table <- compute_population_table(result$gs, experiment = experiment_name)
          table[, gating_template := result$label]
          table[, gated_population := result$population_label]
          table
        }, gating_results),
        use.names = TRUE,
        fill = TRUE
      )
      marker_table <- data.table::rbindlist(
        Map(function(processed, result) {
          table <- compute_marker_table(
            processed$data_dt,
            processed$chs,
            experiment = experiment_name,
            population = if (multi_template_mode) {
              result$population_label
            } else {
              result$extraction_population
            },
            scale = processed$value_scale
          )
          table[, gating_template := result$label]
          table[, gated_population := result$population_label]
          table
        }, processed_results, gating_results),
        use.names = TRUE,
        fill = TRUE
      )
    } else {
      population_table <- compute_population_table(gs, experiment = experiment_name)
      marker_table <- compute_marker_table(
        data_dt,
        chs,
        experiment = experiment_name,
        population = gating_results[[1L]]$extraction_population,
        scale = value_scale
      )
    }

    if (summaries) {
      data.table::fwrite(
        population_table,
        file = file.path(data_path, paste0(experiment_name, "_populations.csv"))
      )
      data.table::fwrite(
        marker_table,
        file = file.path(data_path, paste0(experiment_name, "_markers.csv"))
      )
    }
  }

  # generate reports if requested
  if (report) {
    generate_qc_report(
      data_path,
      gs,
      qc_results,
      so_mat = so_mat,
      gating_results = gating_results,
      data_dt = data_dt,
      population_table = population_table,
      marker_table = marker_table,
      save_gating_set = save_gating_set,
      summaries_written = summaries,
      include_gating_plots = include_gating_plots
    )
  }

  return(data_dt)
}
