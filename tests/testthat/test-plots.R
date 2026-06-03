test_that("summary plot helpers return ggplot objects", {
  population_table <- data.table::data.table(
    experiment = "example",
    sample_id = rep(c("sample1", "sample2"), each = 2),
    population = rep(c("cells", "singlets"), 2),
    parent_population = rep(c("nondebris", "cells"), 2),
    n_events = c(100, 80, 120, 90),
    n_parent_events = c(150, 100, 160, 120),
    frequency_parent = c(66.7, 80, 75, 75),
    frequency_total = c(50, 40, 60, 45),
    gate_method = "opencyto",
    gate_confidence = NA_real_,
    review_required = FALSE
  )

  marker_table <- data.table::data.table(
    experiment = "example",
    sample_id = rep(c("sample1", "sample2"), each = 2),
    population = "singlets",
    channel = rep(c("FL1-A", "FL3-A"), 2),
    marker = rep(c("FL1-A", "FL3-A"), 2),
    median = c(10, 20, 30, 40),
    mean = c(12, 22, 32, 42),
    sd = c(1, 2, 3, 4),
    iqr = c(2, 3, 4, 5),
    mad = c(1.5, 2.5, 3.5, 4.5),
    min = c(1, 2, 3, 4),
    max = c(20, 30, 40, 50),
    percent_positive = c(NA_real_, NA_real_, NA_real_, NA_real_),
    scale = "raw"
  )

  event_counts <- data.table::data.table(
    sample = rep(c("sample1", "sample2"), each = 3),
    population = rep(c("root", "cells", "singlets"), 2),
    count = c(200, 100, 80, 200, 120, 90),
    parent_pct = c(100, 50, 80, 100, 60, 75)
  )

  expect_s3_class(plot_population_summary(population_table), "ggplot")
  expect_s3_class(plot_marker_summary(marker_table), "ggplot")
  expect_s3_class(plot_qc_overview(list(event_counts = event_counts)), "ggplot")
  expect_s3_class(theme_expressalyzr(), "theme")
  expect_true(all(c("blue", "teal", "ink") %in% names(expressalyzr_palette())))
})

test_that("plot_gate_hierarchy returns a ggplot object", {
  test_dir <- setup_compat_experiment()
  on.exit(unlink(test_dir, recursive = TRUE), add = TRUE)

  cs <- load_fcs(test_dir)
  gs <- flowWorkspace::GatingSet(cs)
  openCyto::register_plugins(fun = density_gate, "density_gate", dep = NA, "gating")
  openCyto::register_plugins(fun = mixture_gate, "mixture_gate", dep = NA, "gating")
  gt <- openCyto::gatingTemplate(file.path(test_dir, "gt_samples.csv"))
  openCyto::gt_gating(gt, gs)

  expect_s3_class(plot_gate_hierarchy(gs), "ggplot")
})

test_that("template gate fallback plot renders polygon gates", {
  test_dir <- setup_compat_experiment()
  on.exit(unlink(test_dir, recursive = TRUE), add = TRUE)

  template_file <- file.path(test_dir, "gt_polygon.csv")
  data.table::fwrite(
    data.table::data.table(
      alias = "nondebris",
      pop = "+",
      parent = "root",
      dims = "FSC-A,SSC-A",
      gating_method = "polygon_gate",
      gating_args = paste0(
        "x = c(0, 2e7, 2e7, 0), ",
        "y = c(0, 0, 2e7, 2e7)"
      ),
      collapseDataForGating = "",
      groupBy = "",
      preprocessing_method = "",
      preprocessing_args = ""
    ),
    template_file
  )

  cs <- load_fcs(test_dir)
  gate_designer_register_plugins()
  result <- run_gating_template(cs, template_file,
                                extraction_population = "nondebris")
  p <- plot_template_gate_page(
    result$gs,
    template_file = template_file,
    pop = "nondebris",
    samples = flowWorkspace::sampleNames(result$gs)[1:2],
    parent = "root",
    max_events = 100L
  )

  expect_s3_class(p, "ggplot")
})

test_that("template gate fallback renders automated gates without overlay", {
  test_dir <- setup_compat_experiment()
  on.exit(unlink(test_dir, recursive = TRUE), add = TRUE)

  template_file <- file.path(test_dir, "gt_mixture.csv")
  data.table::fwrite(
    data.table::data.table(
      alias = "cells",
      pop = "+",
      parent = "root",
      dims = "FSC-A,SSC-A",
      gating_method = "mixture_gate",
      gating_args = "n_clusters = 4:5, n_samples = 2000",
      collapseDataForGating = "",
      groupBy = "",
      preprocessing_method = "",
      preprocessing_args = ""
    ),
    template_file
  )

  cs <- load_fcs(test_dir)
  gs <- flowWorkspace::GatingSet(cs)
  p <- plot_template_gate_page(
    gs,
    template_file = template_file,
    pop = "cells",
    samples = flowWorkspace::sampleNames(gs)[1:2],
    parent = "root",
    max_events = 100L
  )

  expect_s3_class(p, "ggplot")
})

test_that("run_pipeline no longer accepts the old plots shortcut", {
  test_dir <- setup_compat_experiment()
  on.exit(unlink(test_dir, recursive = TRUE), add = TRUE)

  expect_error(
    run_pipeline(test_dir, view_config = FALSE, interactive = FALSE,
                 gate_designer = FALSE,
                 plots = TRUE),
    "no longer accepts"
  )
})
