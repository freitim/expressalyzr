test_that("compute_population_table returns expected schema from gated data", {
  test_path <- system.file("extdata", "example_fcs_files",
                           package = "expressalyzr", mustWork = TRUE)
  cs <- load_fcs(test_path)
  gs <- flowWorkspace::GatingSet(cs)

  gt_file <- system.file("tools", "gt_samples.csv",
                         package = "expressalyzr", mustWork = TRUE)
  openCyto::register_plugins(fun = density_gate, "density_gate", dep = NA, "gating")
  openCyto::register_plugins(fun = mixture_gate, "mixture_gate", dep = NA, "gating")
  gt <- openCyto::gatingTemplate(gt_file)
  openCyto::gt_gating(gt, gs)

  result <- compute_population_table(gs, experiment = "example")

  expect_s3_class(result, "data.table")
  expect_named(
    result,
    c("experiment", "sample_id", "population", "parent_population",
      "n_events", "n_parent_events", "frequency_parent", "frequency_total",
      "gate_method", "gate_confidence", "review_required")
  )
  expect_true("root" %in% result$population)
  expect_true("singlets" %in% result$population)
  expect_true(all(result$n_events >= 0))
  expect_true(all(result[population == "root"]$frequency_total == 100))
})

test_that("compute_marker_table returns expected schema and summaries", {
  dt <- data.table::data.table(
    File = rep(c("sample1.fcs", "sample2.fcs"), each = 4),
    FL1 = c(1, 2, 3, 4, 10, 20, 30, 40),
    FL1_pos = c(FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE),
    FL2 = c(5, 6, 7, 8, 50, 60, 70, 80)
  )

  result <- compute_marker_table(dt, c("FL1", "FL2"),
                                 experiment = "example",
                                 population = "singlets",
                                 scale = "raw")

  expect_s3_class(result, "data.table")
  expect_named(
    result,
    c("experiment", "sample_id", "population", "channel", "marker",
      "median", "mean", "sd", "iqr", "mad", "min", "max",
      "percent_positive", "scale")
  )
  expect_equal(nrow(result), 4)
  expect_equal(result[channel == "FL1" & sample_id == "sample1.fcs"]$median, 2.5)
  expect_equal(result[channel == "FL1" & sample_id == "sample1.fcs"]$percent_positive, 50)
  expect_true(is.na(result[channel == "FL2" & sample_id == "sample1.fcs"]$percent_positive))
})
