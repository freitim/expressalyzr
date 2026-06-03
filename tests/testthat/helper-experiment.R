setup_compat_experiment <- function(config_overrides = character(),
                                    include_config = TRUE) {
  test_file_path <- system.file("extdata", "example_fcs_files",
                                package = "expressalyzr", mustWork = TRUE)
  fcs_files <- list.files(file.path(test_file_path, "data"), full.names = TRUE)

  test_dir <- file.path(tempdir(), paste0("compat_pipeline_", Sys.getpid(),
                                          "_", sample.int(100000, 1)))
  test_data_dir <- file.path(test_dir, "data")
  dir.create(test_data_dir, recursive = TRUE, showWarnings = FALSE)
  file.copy(fcs_files, test_data_dir)

  gt_tmpl <- system.file("tools", "gt_samples.csv",
                         package = "expressalyzr", mustWork = TRUE)
  file.copy(gt_tmpl, file.path(test_dir, "gt_samples.csv"))

  if (include_config) {
    config_lines <- c(
      "default:",
      "  mefl_transform: FALSE",
      "  beads_pattern: \"Beads\"",
      "  redo_comp: FALSE",
      "  controls_pattern: \"NOMATCH\"",
      "  controls_index: !expr c(1)",
      "  controls_design_column: \"Channel\"",
      "  controls_map: NULL",
      "  channel_map: \"cytoflex\"",
      "  channel_pattern: \"-A\"",
      "  manual_comp: FALSE",
      "  density_th: 0.0002",
      "  bg_cutoff: 0.99",
      "  manual_cutoff: FALSE",
      "  spec_file: \"does_not_exist.csv\"",
      "  merge_by: \"ID\"",
      "  adjust_gating: FALSE",
      "  bg_channels: NULL"
    )

    if (length(config_overrides) > 0) {
      override_names <- sub("^  ([^:]+):.*$", "\\1", config_overrides)
      config_names <- sub("^  ([^:]+):.*$", "\\1", config_lines)
      keep <- !config_names %in% override_names
      config_lines <- c(config_lines[keep], config_overrides)
    }

    writeLines(config_lines, file.path(test_dir, "config.yml"))
  }

  test_dir
}

run_with_messages <- function(expr) {
  messages <- character()
  value <- withCallingHandlers(
    expr,
    message = function(m) {
      messages <<- c(messages, conditionMessage(m))
      invokeRestart("muffleMessage")
    }
  )
  list(value = value, messages = messages)
}
