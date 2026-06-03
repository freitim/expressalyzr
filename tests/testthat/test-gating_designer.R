test_that("boundary gate arguments are formatted and validated", {
  args <- build_boundary_gating_args(min = c(1, 2), max = c(10, 20))

  expect_identical(args, "min = c(1, 2), max = c(10, 20)")
  expect_error(build_boundary_gating_args(c(1, 2), c(1, 3)),
               "smaller")
  expect_error(build_boundary_gating_args(c(1, NA), c(2, 3)),
               "finite")
})

test_that("polygon gate arguments are formatted and parsed", {
  args <- expressalyzr:::build_polygon_gating_args(
    x = c(1, 10, 10, 1),
    y = c(1, 1, 10, 10)
  )
  parsed <- expressalyzr:::parse_polygon_gating_args(args)

  expect_equal(args, "x = c(1, 10, 10, 1), y = c(1, 1, 10, 10)")
  expect_equal(unname(parsed[, 1]), c(1, 10, 10, 1))
  expect_equal(unname(parsed[, 2]), c(1, 1, 10, 10))
  expect_error(expressalyzr:::build_polygon_gating_args(c(1, 2), c(1, 2)),
               "at least three")
  expect_null(expressalyzr:::parse_polygon_gating_args("not valid"))
})

test_that("boundary gate template rows can be added and replaced", {
  template <- data.table::fread(system.file("tools", "gt_samples.csv",
                                           package = "expressalyzr",
                                           mustWork = TRUE))

  updated <- update_boundary_gate_template(
    template,
    alias = "nondebris",
    parent = "root",
    dims = c("FSC-A", "SSC-A"),
    min = c(100, 200),
    max = c(1000, 2000)
  )

  expect_s3_class(updated, "data.table")
  expect_equal(sum(updated$alias == "nondebris"), 1)
  expect_equal(updated[alias == "nondebris", gating_method], "boundary")
  expect_equal(updated[alias == "nondebris", dims], "FSC-A,SSC-A")
  expect_equal(updated[alias == "nondebris", gating_args],
               "min = c(100, 200), max = c(1000, 2000)")
  expect_true(all(c("cells", "singlets") %in% updated$alias))
})

test_that("gate template rows can be edited without reordering", {
  template <- data.table::fread(system.file("tools", "gt_samples.csv",
                                           package = "expressalyzr",
                                           mustWork = TRUE))
  row <- expressalyzr:::build_gate_template_row(
    alias = "cells",
    pop = "+",
    parent = "nondebris",
    dims = c("FSC-A", "SSC-A"),
    gating_method = "boundary",
    gating_args = build_boundary_gating_args(c(1, 2), c(3, 4)),
    collapseDataForGating = "FALSE",
    groupBy = "sample"
  )

  updated <- expressalyzr:::update_gate_template_row(template, 2L, row)

  expect_equal(updated$alias, template$alias)
  expect_equal(updated[2L, gating_method], "boundary")
  expect_equal(updated[2L, collapseDataForGating], "FALSE")
  expect_equal(updated[2L, groupBy], "sample")
})

test_that("gate template row choices and defaults are stable", {
  template <- data.table::data.table(
    alias = c("nondebris", "cells"),
    gating_method = c("boundary", "mixture_gate")
  )
  choices <- expressalyzr:::gate_template_row_choices(template)

  expect_equal(unname(choices), c("1", "2"))
  expect_equal(expressalyzr:::gate_template_selected_index(template, "cells"), 2)
  expect_equal(expressalyzr:::split_gate_dims("FSC-A, SSC-A"),
               c("FSC-A", "SSC-A"))
  expect_equal(expressalyzr:::gate_designer_default_dims(
    c("Time", "FSC-A", "SSC-A"), "missing"
  ), c("FSC-A", "SSC-A"))
})

test_that("axis layout supports linear and log views", {
  linear <- expressalyzr:::gate_designer_axis_layout("FSC-A", c(1, 10), "linear")
  log <- expressalyzr:::gate_designer_axis_layout("FL1-A", c(1, 10, 100), "log10")

  expect_equal(linear$type, "linear")
  expect_equal(linear$tickformat, ".2e")
  expect_equal(log$type, "log")
  expect_equal(log$exponentformat, "e")
  expect_equal(log$tickmode, "array")
  expect_equal(log$ticktext, c("1.0e+00", "1.0e+01", "1.0e+02"))
  expect_true(all(is.finite(log$range)))
})

test_that("log axes use sparse scientific ticks", {
  ticks <- expressalyzr:::gate_designer_log_ticks(c(9e4, 2e5, 2e7))

  expect_equal(ticks$tickvals, c(1e5, 1e6, 1e7))
  expect_equal(ticks$ticktext, c("1.0e+05", "1.0e+06", "1.0e+07"))
})

test_that("scatter point density highlights locally dense events", {
  density <- expressalyzr:::gate_designer_point_density(
    x = c(rep(1, 10), 100),
    y = c(rep(1, 10), 100),
    bins = 10L
  )

  expect_length(density, 11)
  expect_true(all(is.finite(density)))
  expect_true(density[[1L]] > density[[11L]])
})

test_that("boundary gate arguments can be parsed from template strings", {
  parsed <- expressalyzr:::parse_boundary_gating_args(
    "min = c(9e5, 1e5), max = c(16777215, 16777215)"
  )

  expect_equal(parsed$min, c(9e5, 1e5))
  expect_equal(parsed$max, c(16777215, 16777215))
  expect_null(expressalyzr:::parse_boundary_gating_args(""))
  expect_null(expressalyzr:::parse_boundary_gating_args("not valid"))
})

test_that("plotly polygon paths are converted to vertices", {
  vertices <- expressalyzr:::parse_plotly_polygon_path(
    "M 1,2 L 10,2 L 10,20 L 1,2 Z"
  )

  expect_equal(unname(vertices[, 1]), c(1, 10, 10))
  expect_equal(unname(vertices[, 2]), c(2, 2, 20))

  relayout <- list("shapes[0].path" = "M 1,1 L 5,1 L 5,5 Z")
  extracted <- expressalyzr:::extract_plotly_polygon_vertices(relayout)
  expect_equal(unname(extracted[, 1]), c(1, 5, 5))
  expect_equal(unname(extracted[, 2]), c(1, 1, 5))
})

test_that("pixel polygon vertices can be converted to data coordinates", {
  pixels <- cbind(c(0, 100, 100, 0), c(0, 0, 100, 100))
  data_vertices <- expressalyzr:::gate_designer_pixels_to_data(
    pixels,
    width = 100,
    height = 100,
    xaxis = list(type = "linear", range = c(0, 10)),
    yaxis = list(type = "linear", range = c(0, 20))
  )

  expect_equal(unname(data_vertices[, 1]), c(0, 10, 10, 0))
  expect_equal(unname(data_vertices[, 2]), c(20, 20, 0, 0))

  log_vertices <- expressalyzr:::gate_designer_pixels_to_data(
    cbind(c(0, 100, 100), c(0, 100, 0)),
    width = 100,
    height = 100,
    xaxis = list(type = "log", range = c(1, 3)),
    yaxis = list(type = "log", range = c(2, 4))
  )
  expect_equal(unname(round(log_vertices[, 1])), c(10, 1000, 1000))
  expect_equal(unname(round(log_vertices[, 2])), c(10000, 100, 10000))
})

test_that("gate bounds helpers validate and compare bounds", {
  valid <- list(min = c(1, 2), max = c(10, 20))
  same <- list(min = c(1, 2), max = c(10, 20))
  invalid <- list(min = c(10, 2), max = c(1, 20))

  expect_true(expressalyzr:::gate_designer_bounds_are_valid(valid))
  expect_false(expressalyzr:::gate_designer_bounds_are_valid(invalid))
  expect_true(expressalyzr:::gate_designer_bounds_equal(valid, same))
  expect_false(expressalyzr:::gate_designer_bounds_equal(valid, invalid))
})

test_that("polygon helpers validate, compare, and summarize gates", {
  dt <- data.table::data.table(
    x = c(1, 2, 8, 12),
    y = c(1, 2, 8, 12)
  )
  vertices <- cbind(c(0, 10, 10, 0), c(0, 0, 10, 10))

  expect_true(expressalyzr:::gate_designer_polygon_is_valid(vertices))
  expect_true(expressalyzr:::gate_designer_polygons_equal(vertices, vertices))
  summary <- expressalyzr:::gate_designer_polygon_summary(dt, vertices)
  expect_equal(summary$total, 4)
  expect_equal(summary$inside, 3)
  expect_equal(summary$percent, 75)
})

test_that("lasso-selected points are converted to a polygon hull", {
  selected <- data.frame(
    x = c(0, 10, 10, 0, 5),
    y = c(0, 0, 10, 10, 5)
  )

  vertices <- expressalyzr:::gate_designer_lasso_vertices(selected)

  expect_true(expressalyzr:::gate_designer_polygon_is_valid(vertices))
  expect_equal(nrow(vertices), 4)
  expect_true(all(vertices[, 1] %in% c(0, 10)))
  expect_true(all(vertices[, 2] %in% c(0, 10)))
})

test_that("control samples are inferred from config controls_pattern", {
  test_dir <- tempfile("gate_controls_")
  dir.create(test_dir)
  on.exit(unlink(test_dir, recursive = TRUE), add = TRUE)
  writeLines(c(
    "default:",
    "  controls_pattern: \"Controls\""
  ), file.path(test_dir, "config.yml"))

  controls <- expressalyzr:::gate_designer_control_sample_names(
    test_dir,
    c("01-Controls-A1.fcs", "01-Sample-B1.fcs", "02-Controls-A2.fcs")
  )

  expect_equal(controls, c("01-Controls-A1.fcs", "02-Controls-A2.fcs"))
  expect_equal(
    expressalyzr:::gate_designer_control_sample_names(test_dir, "sample.fcs"),
    character()
  )
})

test_that("gate designer template specs support multiple output files", {
  specs <- expressalyzr:::gate_designer_template_specs(
    output_file = "gt_samples.csv",
    template_file = "gt_samples.csv",
    template_files = c("gt_HEK.csv", "gt_Tcells.csv"),
    template_labels = c("HEK", "Tcells")
  )

  expect_equal(specs$template_label, c("HEK", "Tcells"))
  expect_equal(basename(specs$template_file), c("gt_HEK.csv", "gt_Tcells.csv"))
  expect_equal(basename(specs$output_file), c("gt_HEK.csv", "gt_Tcells.csv"))
  expect_error(
    expressalyzr:::gate_designer_template_specs(
      output_file = "gt_samples.csv",
      template_file = "gt_samples.csv",
      template_files = c("gt_HEK.csv", "gt_Tcells.csv"),
      template_labels = "HEK"
    ),
    "one non-empty label"
  )
})

test_that("polygon_gate builds a flowCore polygon gate", {
  gate <- expressalyzr:::polygon_gate(
    fr = NULL,
    pp_res = NULL,
    channels = c("FSC-A", "SSC-A"),
    x = c(0, 10, 10, 0),
    y = c(0, 0, 10, 10),
    filterId = "freeform"
  )

  expect_s4_class(gate, "polygonGate")
})

test_that("gate summary reports current gate coverage", {
  dt <- data.table::data.table(
    x = c(1, 2, 3, 4),
    y = c(1, 2, 3, 4)
  )
  bounds <- list(min = c(1.5, 1.5), max = c(3.5, 3.5))

  summary <- expressalyzr:::gate_designer_gate_summary(dt, bounds)
  expect_equal(summary$total, 4)
  expect_equal(summary$inside, 2)
  expect_equal(summary$percent, 50)
})

test_that("plotly shape relayout events are converted to gate bounds", {
  bounds <- expressalyzr:::extract_plotly_rect_bounds(list(
    "shapes[0].x0" = 10,
    "shapes[0].x1" = 2,
    "shapes[0].y0" = 20,
    "shapes[0].y1" = 5
  ))

  expect_equal(bounds$min, c(2, 5))
  expect_equal(bounds$max, c(10, 20))

  full_shape <- expressalyzr:::extract_plotly_rect_bounds(list(
    shapes = list(list(x0 = 100, x1 = 20, y0 = 200, y1 = 50))
  ))
  expect_equal(full_shape$min, c(20, 50))
  expect_equal(full_shape$max, c(100, 200))
  active_shape <- expressalyzr:::extract_plotly_rect_bounds(list(
    shapes = list(
      list(x0 = 1, x1 = 2, y0 = 1, y1 = 2, editable = FALSE),
      list(x0 = 10, x1 = 20, y0 = 30, y1 = 40, editable = TRUE)
    )
  ))
  expect_equal(active_shape$min, c(10, 30))
  expect_equal(active_shape$max, c(20, 40))
  expect_null(expressalyzr:::extract_plotly_rect_bounds(NULL))
})

test_that("gate designer overlay colors are valid rgba inputs", {
  colors <- expressalyzr:::gate_designer_overlay_colors(4)
  rgba <- expressalyzr:::gate_designer_rgba(colors[[1L]], 0.2)

  expect_length(colors, 4)
  expect_true(all(grepl("^#", colors)))
  expect_match(rgba, "^rgba\\([0-9]+, [0-9]+, [0-9]+, 0.2\\)$")
})

test_that("gate designer can be constructed without launching", {
  testthat::skip_if_not_installed("shiny")
  testthat::skip_if_not_installed("plotly")
  testthat::skip_if_not_installed("DT")

  test_dir <- setup_compat_experiment()
  on.exit(unlink(test_dir, recursive = TRUE), add = TRUE)

  app <- design_gating_template(test_dir, launch = FALSE, max_events = 500)
  expect_s3_class(app, "shiny.appobj")
})

test_that("combined multi-template gate designer can be constructed", {
  testthat::skip_if_not_installed("shiny")
  testthat::skip_if_not_installed("plotly")
  testthat::skip_if_not_installed("DT")

  test_dir <- setup_compat_experiment()
  on.exit(unlink(test_dir, recursive = TRUE), add = TRUE)
  file.copy(file.path(test_dir, "gt_samples.csv"),
            file.path(test_dir, "gt_HEK.csv"))
  file.copy(file.path(test_dir, "gt_samples.csv"),
            file.path(test_dir, "gt_Tcells.csv"))

  app <- design_gating_template(
    test_dir,
    template_files = file.path(test_dir, c("gt_HEK.csv", "gt_Tcells.csv")),
    output_files = file.path(test_dir, c("gt_HEK.csv", "gt_Tcells.csv")),
    template_labels = c("HEK", "Tcells"),
    launch = FALSE,
    max_events = 500
  )
  expect_s3_class(app, "shiny.appobj")
})
