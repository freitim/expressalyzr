#### helper: create a temporary experiment directory ####
setup_test_experiment <- function() {
  test_file_path <- system.file("extdata", "example_fcs_files",
                                package = "expressalyzr", mustWork = TRUE)
  fcs_files <- list.files(file.path(test_file_path, "data"), full.names = TRUE)

  test_dir <- file.path(tempdir(), paste0("pipeline_test_", Sys.getpid()))
  test_data_dir <- file.path(test_dir, "data")
  dir.create(test_data_dir, recursive = TRUE, showWarnings = FALSE)
  file.copy(fcs_files, test_data_dir)

  gt_tmpl <- system.file("tools", "gt_samples.csv",
                         package = "expressalyzr", mustWork = TRUE)
  file.copy(gt_tmpl, file.path(test_dir, "gt_samples.csv"))

  # non-interactive config that skips compensation & bg removal for the
  # example dataset (only 2 FL channels, no proper control set)
  config_text <- paste(
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
    "  bg_channels: NULL",
    sep = "\n"
  )
  writeLines(config_text, file.path(test_dir, "config.yml"))

  test_dir
}

#### test run_pipeline backward compatibility ####

test_that("run_pipeline can still return a legacy data.table", {
  test_dir <- setup_test_experiment()
  on.exit(unlink(test_dir, recursive = TRUE), add = TRUE)

  result <- run_pipeline(test_dir, view_config = FALSE,
                         gate_designer = FALSE, report = FALSE)

  expect_s3_class(result, "data.table")
  expect_true(nrow(result) > 0)
  expect_true(file.exists(file.path(test_dir, paste0(basename(test_dir), ".csv"))))
})

#### test run_pipeline with QC ####

test_that("run_pipeline with report = TRUE returns data table and writes report", {
  test_dir <- setup_test_experiment()
  on.exit(unlink(test_dir, recursive = TRUE), add = TRUE)

  result <- run_pipeline(test_dir, view_config = FALSE,
                         gate_designer = FALSE, report = TRUE)

  expect_s3_class(result, "data.table")
  expect_false(file.exists(file.path(test_dir, paste0(basename(test_dir), "_qc.csv"))))
  expect_true(file.exists(file.path(test_dir, paste0(basename(test_dir), "_report.html"))))
})

#### test run_pipeline report outputs ####

test_that("run_pipeline renders unified report without changing data output", {
  test_dir <- setup_test_experiment()
  on.exit(unlink(test_dir, recursive = TRUE), add = TRUE)

  result <- run_pipeline(test_dir, view_config = FALSE, interactive = FALSE,
                         gate_designer = FALSE,
                         report = TRUE)

  expect_s3_class(result, "data.table")
  expect_false(file.exists(file.path(test_dir, paste0(basename(test_dir), "_qc.csv"))))
  expect_true(file.exists(file.path(test_dir, paste0(basename(test_dir), "_report.html"))))
  expect_false(file.exists(file.path(test_dir, paste0(basename(test_dir), "_channels.html"))))
  expect_false(file.exists(file.path(test_dir, paste0(basename(test_dir), "_gating.html"))))
})
