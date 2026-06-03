#' Run the data processing pipline
#'
#' Run the data processing pipline that links together other functions in the
#'  expressalyzr package.
#'
#' @return
#' @export
#'
run_pipeline <- function(data_path, view_config = TRUE, gating_output = NULL) {

  # initialize
  experiment_name <- basename(data_path)

  create_data_subdir(data_path)

  config_file_path <- file.path(data_path, "config.yml")
  config <- load_config(config_file_path, view_config)

  gt_file <- file.path(data_path, "gt_samples.csv")

  if (config$adjust_gating) {
    file.edit(gt_file)

    cat("\n")
    readline(prompt = "Press [Enter] to continue.")
  }

  # start analysis
  cs <- load_fcs(data_path)

  to_cytoflex_channels <- list(
    "445 [B]-A" = "FL6-A",
    "445 [B]-H" = "FL6-H",
    "445 [B]-W" = "FL6-Width",
    "488 [C]-A" = "FL1-A",
    "488 [C]-H" = "FL1-H",
    "488 [C]-W" = "FL1-Width",
    "640 [C]-A" = "FL3-A",
    "640 [C]-H" = "FL3-H",
    "640 [C]-W" = "FL3-Width"
  )

  chs <- flowWorkspace::colnames(cs)

  for (i in seq_along(chs)) {
    if (!is.null(to_cytoflex_channels[[chs[i]]])) {
      chs[i] <- to_cytoflex_channels[[chs[i]]]
    }
  }

  flowWorkspace::colnames(cs) <- chs

  bead_index <- grepl(config$beads_pattern, flowCore::sampleNames(cs))

  if (sum(bead_index) == 1) {

    cs_beads <- cs[bead_index]
    cs <- cs[!bead_index]

    if (config$mefl_transform) {
      t_fun <- generate_transformation(cs_beads)
      trans <- TRUE
    } else {
      trans <- FALSE
    }
  } else if (sum(bead_index) > 1) {

    stop("More than one bead sample found in the dataset. Please restricit the
         bead_pattern to one sample.")

  } else {

    message("Bead sample not found. Values will not be transformed to MEFL.")
    trans <- FALSE
  }

  # gating
  openCyto::register_plugins(fun = density_gate, "density_gate", dep = NA, "gating")
  openCyto::register_plugins(fun = mixture_gate, "mixture_gate", dep = NA, "gating")

  gt <- openCyto::gatingTemplate(gt_file)
  gs <- flowWorkspace::GatingSet(cs)

  openCyto::gt_gating(gt, gs)

  if (!is.null(gating_output)) {
    if (gating_output == "inspect") {
      for (i in 1:length(gs)) {
        print(autoplot(gs[[i]]))
        readline()
      }
    } else if (gating_output == "report") {
      rmarkdown::render(system.file("tools","generate_gating_report.R", package = "expressalyzr",
                                    mustWork = TRUE),
                        output_file = file.path(data_path, paste0(experiment_name, "_gating.html")))
    } else if (gating_output == "set") {
      flowWorkspace::save_gs(gs, path = "gs")    }
  }

  data_cs <- flowWorkspace::gs_pop_get_data(gs, y = "singlets")

  if (trans) {
    data_cs <- apply_transform(data_cs, t_fun)
  }

  # compensation
  cont_index <- grepl(config$controls_pattern, flowCore::sampleNames(data_cs))

  print(flowCore::sampleNames(data_cs)[cont_index])

  n_controls <- sum(cont_index)
  if (n_controls > 0) {
    cont_cs <- data_cs[cont_index]
    cont_names <- flowCore::sampleNames(cont_cs)
    cont_order <- gtools::mixedorder(gsub("^.*([A-Z]{0,1}\\d{1,2}).fcs", "\\1", cont_names))
    cont_cs <- cont_cs[cont_order]
  }

  if (n_controls > 2) {
    so_mat_path <- file.path(data_path, "so_mat.RData")

    if (!file.exists(so_mat_path) || config$redo_comp) {
      so_mat <- spillover_matrix(cont_cs,
                                 config$controls_index,
                                 config$channel_pattern,
                                 config$density_th,
                                 config$manual_comp)

      save(so_mat, file = so_mat_path)
    } else {
      load(so_mat_path)
    }

    data_cs <- flowWorkspace::compensate(data_cs, so_mat)
    # data_cs <- data_cs[!cont_index]
  } else {
    message("No control samples found. Proceeding with out compensation.")
  }

  # convert cytoset to data table
  data_dt <- cs_to_dt(data_cs)

  chs <- colnames(data_dt)[grepl("FL", colnames(data_dt))]

  # gate populations
  data_dt[, no_negative := rowSums(data_dt[, chs, with = FALSE] <= 0) == 0]
  data_dt[, positive := rowSums(data_dt[, chs, with = FALSE] < 0) == 0]

  if (n_controls > 0) {
    cont_dt <- cs_to_dt(cont_cs)
    neg_dt <- cs_to_dt(cont_cs[config$controls_index[1]])

    if (config$manual_cutoff) {
      select_chs <- chs[grepl(config$channel_pattern, chs)]
      cont_sub_dt <- cont_dt[, c("File", select_chs), with = FALSE]
      config$bg_cutoff <- adjust_threshold(cont_sub_dt, config$bg_cutoff)
      write_value <- paste("bg_cutoff:", config$bg_cutoff)
      write_config(config_file_path, "default", write_value)
    }

    neg_dt[, lapply(mget(chs), quantile, config$bg_cutoff),
           by = .(File)]

    bg_dt <- neg_dt[, lapply(mget(chs), quantile, config$bg_cutoff)]

    f_pos <- function(channel, values) values >= bg_dt[[channel]]

    pos_chs <- paste0(chs, "_pos")
    data_dt[, (pos_chs) := lapply(chs, function(ch) f_pos(ch, get(ch))), by = .(File)]
  }

  # new background removal
  if (!is.null(config$bg_channels)) {
    data_dt[(positive), (paste0(config$bg_channels, "_bg")) := lapply(mget(config$bg_channels),
                                                                      assign_bg,
                                                                      n_comp = NULL,
                                                                      rm = 1,
                                                                      inspect = FALSE),
            by = .(File)]
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
  return(data_dt)
}
