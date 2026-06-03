#### test compute_channel_stats ####

test_that("compute_channel_stats returns correct structure", {
  dt <- data.table::data.table(
    File = rep(c("sample1.fcs", "sample2.fcs"), each = 100),
    FL1 = c(rnorm(100, 500, 50), rnorm(100, 800, 80)),
    FL2 = c(rnorm(100, 300, 30), rnorm(100, 600, 60))
  )

  result <- compute_channel_stats(dt, c("FL1", "FL2"))

  expect_s3_class(result, "data.table")
  expect_equal(nrow(result), 4)
  expect_named(result, c("sample", "channel", "mean", "median", "sd", "cv",
                          "pct_positive", "dynamic_range"))
  expect_true(all(result$pct_positive >= 0 & result$pct_positive <= 100))
})

test_that("compute_channel_stats handles missing channels gracefully", {
  dt <- data.table::data.table(File = "s1.fcs", FL1 = 1:10)
  result <- compute_channel_stats(dt, c("FL1", "NONEXISTENT"))

  expect_s3_class(result, "data.table")
  expect_equal(nrow(result), 1)
  expect_equal(result$channel, "FL1")
})

test_that("compute_channel_stats returns empty table for no matching channels", {
  dt <- data.table::data.table(File = "s1.fcs", X = 1:10)
  result <- compute_channel_stats(dt, c("FL1", "FL2"))

  expect_s3_class(result, "data.table")
  expect_equal(nrow(result), 0)
})

test_that("compute_channel_stats handles all-negative values", {
  dt <- data.table::data.table(File = "s1.fcs", FL1 = -10:-1)
  result <- compute_channel_stats(dt, "FL1")

  expect_equal(result$pct_positive, 0)
})

#### test compute_scatter_stats ####

test_that("compute_scatter_stats returns correct structure from example data", {
  test_path <- system.file("extdata", "example_fcs_files",
                           package = "expressalyzr", mustWork = TRUE)
  cs <- load_fcs(test_path)

  result <- compute_scatter_stats(cs)

  expect_s3_class(result, "data.table")
  expect_true(nrow(result) > 0)
  expect_named(result, c("sample", "channel", "mean", "median", "sd", "cv"))
  expect_true(all(result$mean > 0))
  expect_true(all(result$sd > 0))
})

test_that("compute_scatter_stats works with custom channels", {
  test_path <- system.file("extdata", "example_fcs_files",
                           package = "expressalyzr", mustWork = TRUE)
  cs <- load_fcs(test_path)
  avail <- flowCore::colnames(cs[[1]])
  ch <- avail[1]

  result <- compute_scatter_stats(cs, channels = ch)

  expect_equal(unique(result$channel), ch)
})

#### test compute_time_check ####

test_that("compute_time_check returns NULL when no Time channel", {
  test_path <- system.file("extdata", "example_fcs_files",
                           package = "expressalyzr", mustWork = TRUE)
  cs <- load_fcs(test_path)
  avail <- flowCore::colnames(cs[[1]])

  if (!"Time" %in% avail) {
    expect_message(result <- compute_time_check(cs), "No 'Time' channel")
    expect_null(result)
  } else {
    result <- compute_time_check(cs)
    expect_s3_class(result, "data.table")
    expect_named(result, c("sample", "time_bin", "event_count", "rate_anomaly"))
  }
})

#### test compute_event_counts ####

test_that("compute_event_counts returns correct structure", {
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

  result <- compute_event_counts(gs)

  expect_s3_class(result, "data.table")
  expect_named(result, c("sample", "population", "count", "parent_pct"))
  expect_true(all(result$count >= 0))
  expect_true("root" %in% result$population)
  expect_true(all(result[population == "root"]$parent_pct == 100))
})

#### test collect_qc ####

test_that("collect_qc returns a named list with all components", {
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

  data_cs <- flowWorkspace::gs_pop_get_data(gs, y = "singlets")
  dt <- cs_to_dt(data_cs)
  fl_chs <- colnames(dt)[grepl("FL|488|445|640", colnames(dt))]

  result <- collect_qc(gs, data_cs, dt, fl_chs)

  expect_type(result, "list")
  expect_named(result, c("event_counts", "scatter_stats", "channel_stats", "time_check"))
  expect_s3_class(result$event_counts, "data.table")
  expect_s3_class(result$scatter_stats, "data.table")
  expect_s3_class(result$channel_stats, "data.table")
})
