test_that("run_pipeline preserves legacy return and output contract", {
  test_dir <- setup_compat_experiment()
  on.exit(unlink(test_dir, recursive = TRUE), add = TRUE)

  result <- run_pipeline(test_dir, view_config = FALSE, interactive = FALSE,
                         gate_designer = FALSE, report = FALSE)

  expect_s3_class(result, "data.table")
  expect_true(nrow(result) > 0)
  expect_true(file.exists(file.path(test_dir, paste0(basename(test_dir), ".csv"))))
  expect_true(all(c("positive", "no_negative") %in% names(result)))
  expect_false(file.exists(file.path(test_dir, paste0(basename(test_dir), "_manifest.csv"))))
  expect_false(file.exists(file.path(test_dir, paste0(basename(test_dir), "_populations.csv"))))
  expect_false(file.exists(file.path(test_dir, paste0(basename(test_dir), "_markers.csv"))))
})

test_that("run_pipeline keeps singlets as the default extraction population", {
  test_dir <- setup_compat_experiment()
  on.exit(unlink(test_dir, recursive = TRUE), add = TRUE)

  result <- run_pipeline(test_dir, view_config = FALSE, interactive = FALSE,
                         gate_designer = FALSE, report = TRUE)

  expect_s3_class(result, "data.table")
  expect_true(nrow(result) > 0)
})

test_that("run_pipeline can render QC report without standalone QC CSV", {
  test_dir <- setup_compat_experiment()
  on.exit(unlink(test_dir, recursive = TRUE), add = TRUE)

  result <- run_pipeline(
    test_dir,
    view_config = FALSE,
    interactive = FALSE,
    gate_designer = FALSE,
    report = TRUE
  )

  experiment_name <- basename(test_dir)
  expect_s3_class(result, "data.table")
  expect_true(file.exists(file.path(test_dir,
                                    paste0(experiment_name, "_report.html"))))
  expect_false(file.exists(file.path(test_dir,
                                     paste0(experiment_name, "_qc.csv"))))
})

test_that("run_pipeline preserves no-bead and no-control compatibility path", {
  test_dir <- setup_compat_experiment()
  on.exit(unlink(test_dir, recursive = TRUE), add = TRUE)

  run <- run_with_messages(
    run_pipeline(test_dir, view_config = FALSE, interactive = FALSE,
                 gate_designer = FALSE, report = FALSE)
  )

  expect_s3_class(run$value, "data.table")
  expect_true(any(grepl("Bead sample not found", run$messages, fixed = TRUE)))
  expect_true(any(grepl("No control samples found", run$messages, fixed = TRUE)))
})

test_that("non-interactive mode requires a pre-existing config", {
  test_dir <- setup_compat_experiment(include_config = FALSE)
  on.exit(unlink(test_dir, recursive = TRUE), add = TRUE)

  expect_error(
    run_pipeline(test_dir, view_config = FALSE, interactive = FALSE,
                 gate_designer = FALSE),
    "config.yml does not exist"
  )
})

test_that("legacy gating inspection option is rejected clearly", {
  test_dir <- setup_compat_experiment()
  on.exit(unlink(test_dir, recursive = TRUE), add = TRUE)

  expect_error(
    run_pipeline(test_dir, view_config = FALSE, interactive = FALSE,
                 gate_designer = FALSE,
                 gating_output = "inspect"),
    "no longer supported"
  )
})

test_that("report and save_gating_set write unified report and GatingSet", {
  test_dir <- setup_compat_experiment()
  on.exit(unlink(test_dir, recursive = TRUE), add = TRUE)

  result <- run_pipeline(
    test_dir,
    view_config = FALSE,
    interactive = FALSE,
    gate_designer = FALSE,
    report = TRUE,
    save_gating_set = TRUE
  )

  experiment_name <- basename(test_dir)
  expect_s3_class(result, "data.table")
  expect_true(file.exists(file.path(test_dir,
                                    paste0(experiment_name, "_report.html"))))
  expect_false(file.exists(file.path(test_dir,
                                     paste0(experiment_name, "_gating.html"))))
  expect_true(dir.exists(file.path(test_dir, "gs")))
})

test_that("legacy gating_output report alias maps to unified report", {
  test_dir <- setup_compat_experiment()
  on.exit(unlink(test_dir, recursive = TRUE), add = TRUE)

  result <- run_pipeline(
    test_dir,
    view_config = FALSE,
    interactive = FALSE,
    gate_designer = FALSE,
    gating_output = "report"
  )

  experiment_name <- basename(test_dir)
  expect_s3_class(result, "data.table")
  expect_true(file.exists(file.path(test_dir,
                                    paste0(experiment_name, "_report.html"))))
  expect_false(file.exists(file.path(test_dir,
                                     paste0(experiment_name, "_gating.html"))))

  html_text <- paste(readLines(file.path(test_dir,
                                         paste0(experiment_name, "_report.html")),
                               warn = FALSE),
                     collapse = "\n")
  expect_match(html_text, "Gate Plots", fixed = TRUE)
})

test_that("legacy positional gating_output report alias is preserved", {
  test_dir <- setup_compat_experiment()
  on.exit(unlink(test_dir, recursive = TRUE), add = TRUE)

  result <- run_pipeline(test_dir, FALSE, "report", interactive = FALSE,
                         gate_designer = FALSE)

  experiment_name <- basename(test_dir)
  expect_s3_class(result, "data.table")
  expect_true(file.exists(file.path(test_dir,
                                    paste0(experiment_name, "_report.html"))))
})

test_that("unsupported gating_output values are rejected clearly", {
  test_dir <- setup_compat_experiment()
  on.exit(unlink(test_dir, recursive = TRUE), add = TRUE)

  expect_error(
    run_pipeline(test_dir, view_config = FALSE, interactive = FALSE,
                 gate_designer = FALSE,
                 gating_output = TRUE),
    "gating_output must be one of"
  )

  expect_error(
    run_pipeline(test_dir, view_config = FALSE, interactive = FALSE,
                 gate_designer = FALSE,
                 gating_output = "save"),
    "gating_output must be one of"
  )
})

test_that("legacy gating_output set option is rejected clearly", {
  test_dir <- setup_compat_experiment()
  on.exit(unlink(test_dir, recursive = TRUE), add = TRUE)

  expect_error(
    run_pipeline(test_dir, view_config = FALSE, interactive = FALSE,
                 gate_designer = FALSE,
                 gating_output = "set"),
    "no longer supported"
  )
})

test_that("non-interactive mode rejects pipeline-launched gate designer", {
  test_dir <- setup_compat_experiment()
  on.exit(unlink(test_dir, recursive = TRUE), add = TRUE)

  expect_error(
    run_pipeline(test_dir, view_config = FALSE, interactive = FALSE,
                 gate_designer = TRUE),
    "gate_designer"
  )
})

test_that("non-interactive mode rejects manual config and cutoff prompts", {
  test_dir <- setup_compat_experiment(
    config_overrides = "  adjust_gating: TRUE"
  )
  on.exit(unlink(test_dir, recursive = TRUE), add = TRUE)

  expect_error(
    run_pipeline(test_dir, view_config = FALSE, interactive = FALSE,
                 gate_designer = FALSE),
    "adjust_gating"
  )
})

test_that("removed run_pipeline options are rejected clearly", {
  test_dir <- setup_compat_experiment()
  on.exit(unlink(test_dir, recursive = TRUE), add = TRUE)

  removed <- c("qc", "qc_output", "gating_report", "channel_report",
               "provenance", "plots", "use_gui")
  for (arg in removed) {
    args <- list(test_dir, view_config = FALSE, interactive = FALSE,
                 gate_designer = FALSE)
    args[[arg]] <- if (arg == "qc_output") "csv" else TRUE
    expect_error(
      do.call(run_pipeline, args),
      "no longer accepts"
    )
  }
})

test_that("summary mode writes additive population and marker tables", {
  test_dir <- setup_compat_experiment()
  on.exit(unlink(test_dir, recursive = TRUE), add = TRUE)

  result <- run_pipeline(test_dir, view_config = FALSE, interactive = FALSE,
                         gate_designer = FALSE, report = FALSE,
                         summaries = TRUE)

  experiment_name <- basename(test_dir)
  populations_path <- file.path(test_dir, paste0(experiment_name, "_populations.csv"))
  markers_path <- file.path(test_dir, paste0(experiment_name, "_markers.csv"))

  expect_s3_class(result, "data.table")
  expect_true(file.exists(populations_path))
  expect_true(file.exists(markers_path))

  populations <- data.table::fread(populations_path)
  markers <- data.table::fread(markers_path)

  expect_named(
    populations,
    c("experiment", "sample_id", "population", "parent_population",
      "n_events", "n_parent_events", "frequency_parent", "frequency_total",
      "gate_method", "gate_confidence", "review_required")
  )
  expect_true("singlets" %in% populations$population)
  expect_equal(sum(populations[population == "singlets"]$n_events), nrow(result))

  expect_named(
    markers,
    c("experiment", "sample_id", "population", "channel", "marker",
      "median", "mean", "sd", "iqr", "mad", "min", "max",
      "percent_positive", "scale")
  )
  expect_true(nrow(markers) > 0)
  expect_true(all(markers$population == "singlets"))
})

test_that("run_pipeline can apply multiple explicit gating templates", {
  test_dir <- setup_compat_experiment()
  on.exit(unlink(test_dir, recursive = TRUE), add = TRUE)
  file.copy(file.path(test_dir, "gt_samples.csv"),
            file.path(test_dir, "gt_HEK.csv"))
  file.copy(file.path(test_dir, "gt_samples.csv"),
            file.path(test_dir, "gt_Tcells.csv"))
  tcells_template <- data.table::fread(file.path(test_dir, "gt_Tcells.csv"))
  tcells_template <- tcells_template[alias != "singlets"]
  data.table::fwrite(tcells_template, file.path(test_dir, "gt_Tcells.csv"))

  result <- run_pipeline(
    test_dir,
    view_config = FALSE,
    interactive = FALSE,
    gate_designer = FALSE,
    gating_templates = c(
      HEK = "gt_HEK.csv",
      Tcells = "gt_Tcells.csv"
    ),
    report = TRUE,
    summaries = TRUE
  )

  expect_s3_class(result, "data.table")
  expect_true(all(c("gating_template", "gated_population") %in% names(result)))
  expect_equal(sort(unique(result$gating_template)),
               c("gt_HEK", "gt_Tcells"))
  expect_equal(sort(unique(result$gated_population)),
               c("HEK", "Tcells"))

  experiment_name <- basename(test_dir)
  expect_true(file.exists(file.path(test_dir,
                                    paste0(experiment_name, "_report.html"))))
  expect_false(file.exists(file.path(test_dir,
                                     paste0(experiment_name,
                                            "_gt_HEK_gating.html"))))
  expect_false(file.exists(file.path(test_dir,
                                     paste0(experiment_name,
                                            "_gt_Tcells_gating.html"))))
  expect_false(file.exists(file.path(test_dir,
                                     paste0(experiment_name, "_qc.csv"))))
  populations <- data.table::fread(
    file.path(test_dir, paste0(experiment_name, "_populations.csv"))
  )
  markers <- data.table::fread(
    file.path(test_dir, paste0(experiment_name, "_markers.csv"))
  )
  expect_equal(sort(unique(populations$gating_template)),
               c("gt_HEK", "gt_Tcells"))
  expect_equal(sort(unique(populations$gated_population)),
               c("HEK", "Tcells"))
  expect_equal(sort(unique(markers$gating_template)),
               c("gt_HEK", "gt_Tcells"))
  expect_equal(sort(unique(markers$gated_population)),
               c("HEK", "Tcells"))
  expect_equal(sort(unique(markers$population)),
               c("HEK", "Tcells"))
  expect_equal(
    unique(markers[gated_population == "Tcells"]$population),
    "Tcells"
  )
})
