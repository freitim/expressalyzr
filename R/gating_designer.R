gating_template_columns <- function() {
  c(
    "alias", "pop", "parent", "dims", "gating_method", "gating_args",
    "collapseDataForGating", "groupBy", "preprocessing_method",
    "preprocessing_args"
  )
}

read_gating_template <- function(template_file) {
  if (!file.exists(template_file)) {
    cols <- gating_template_columns()
    out <- data.table::as.data.table(stats::setNames(rep(list(character()), length(cols)), cols))
    return(out)
  }

  template <- data.table::fread(template_file)
  missing_cols <- setdiff(gating_template_columns(), names(template))
  for (col in missing_cols) {
    template[, (col) := ""]
  }
  data.table::setcolorder(template, gating_template_columns())
  for (col in gating_template_columns()) {
    template[, (col) := as.character(get(col))]
    template[is.na(get(col)), (col) := ""]
  }
  template
}

format_gate_number <- function(x) {
  format(signif(x, 10), scientific = FALSE, trim = TRUE)
}

#' Build OpenCyto boundary gate arguments.
#'
#' @param min Numeric vector of lower bounds for two channels.
#' @param max Numeric vector of upper bounds for two channels.
#'
#' @return A character string suitable for the \code{gating_args} column of an
#'   OpenCyto gating template using \code{gating_method = "boundary"}.
#' @export
build_boundary_gating_args <- function(min, max) {
  if (!is.numeric(min) || !is.numeric(max) ||
      length(min) != 2L || length(max) != 2L ||
      any(!is.finite(min)) || any(!is.finite(max))) {
    stop("min and max must be finite numeric vectors of length 2.",
         call. = FALSE)
  }
  if (any(min >= max)) {
    stop("Each min value must be smaller than its matching max value.",
         call. = FALSE)
  }

  paste0(
    "min = c(", paste(format_gate_number(min), collapse = ", "),
    "), max = c(", paste(format_gate_number(max), collapse = ", "), ")"
  )
}

parse_boundary_gating_args <- function(gating_args) {
  if (is.na(gating_args) || !nzchar(gating_args)) {
    return(NULL)
  }

  parsed <- tryCatch(
    eval(parse(text = paste0("list(", gating_args, ")")), envir = baseenv()),
    error = function(e) NULL
  )
  if (is.null(parsed) || is.null(parsed$min) || is.null(parsed$max)) {
    return(NULL)
  }
  min <- as.numeric(parsed$min)
  max <- as.numeric(parsed$max)
  if (length(min) != 2L || length(max) != 2L ||
      any(!is.finite(min)) || any(!is.finite(max))) {
    return(NULL)
  }
  list(min = min, max = max)
}

build_polygon_gating_args <- function(x, y) {
  if (!is.numeric(x) || !is.numeric(y) ||
      length(x) != length(y) || length(x) < 3L ||
      any(!is.finite(c(x, y)))) {
    stop("x and y must be finite numeric vectors with at least three points.",
         call. = FALSE)
  }

  paste0(
    "x = c(", paste(format_gate_number(x), collapse = ", "),
    "), y = c(", paste(format_gate_number(y), collapse = ", "), ")"
  )
}

parse_polygon_gating_args <- function(gating_args) {
  if (is.na(gating_args) || !nzchar(gating_args)) {
    return(NULL)
  }

  parsed <- tryCatch(
    eval(parse(text = paste0("list(", gating_args, ")")), envir = baseenv()),
    error = function(e) NULL
  )
  if (is.null(parsed)) {
    return(NULL)
  }

  if (!is.null(parsed$vertices)) {
    vertices <- as.matrix(parsed$vertices)
    if (ncol(vertices) != 2L) {
      return(NULL)
    }
  } else if (!is.null(parsed$x) && !is.null(parsed$y)) {
    x <- as.numeric(parsed$x)
    y <- as.numeric(parsed$y)
    if (length(x) != length(y)) {
      return(NULL)
    }
    vertices <- cbind(x, y)
  } else {
    return(NULL)
  }

  storage.mode(vertices) <- "numeric"
  if (nrow(vertices) < 3L || any(!is.finite(vertices))) {
    return(NULL)
  }
  colnames(vertices) <- c("x", "y")
  vertices
}

#' Update a gating template with a drawn boundary gate.
#'
#' @param template A data frame/data.table with OpenCyto gating-template
#'   columns.
#' @param alias Population alias to create or replace.
#' @param parent Parent population for the boundary gate.
#' @param dims Character vector of two channel names.
#' @param min Numeric vector of lower bounds for \code{dims}.
#' @param max Numeric vector of upper bounds for \code{dims}.
#' @param replace Whether to replace an existing row with the same alias.
#'
#' @return A \code{data.table} gating template.
#' @export
update_boundary_gate_template <- function(template, alias = "nondebris",
                                          parent = "root",
                                          dims = c("FSC-A", "SSC-A"),
                                          min, max, replace = TRUE) {
  if (!is.character(alias) || length(alias) != 1L || !nzchar(alias)) {
    stop("alias must be a non-empty character scalar.", call. = FALSE)
  }
  if (!is.character(parent) || length(parent) != 1L || !nzchar(parent)) {
    stop("parent must be a non-empty character scalar.", call. = FALSE)
  }
  if (!is.character(dims) || length(dims) != 2L || any(!nzchar(dims))) {
    stop("dims must be a character vector of two channel names.",
         call. = FALSE)
  }

  template <- data.table::as.data.table(data.table::copy(template))
  missing_cols <- setdiff(gating_template_columns(), names(template))
  for (col in missing_cols) {
    template[, (col) := ""]
  }
  data.table::setcolorder(template, gating_template_columns())
  for (col in gating_template_columns()) {
    template[, (col) := as.character(get(col))]
    template[is.na(get(col)), (col) := ""]
  }

  row <- data.table::data.table(
    alias = alias,
    pop = "+",
    parent = parent,
    dims = paste(dims, collapse = ","),
    gating_method = "boundary",
    gating_args = build_boundary_gating_args(min, max),
    collapseDataForGating = "",
    groupBy = "",
    preprocessing_method = "",
    preprocessing_args = ""
  )

  if (isTRUE(replace) && alias %in% template$alias) {
    template <- template[alias != row$alias]
    template <- data.table::rbindlist(list(row, template), use.names = TRUE,
                                      fill = TRUE)
  } else {
    template <- data.table::rbindlist(list(row, template), use.names = TRUE,
                                      fill = TRUE)
  }

  data.table::setcolorder(template, gating_template_columns())
  template
}

gate_template_row_label <- function(template) {
  if (nrow(template) == 0L) {
    return(character())
  }
  paste0(seq_len(nrow(template)), ". ", template$alias, " (",
         template$gating_method, ")")
}

gate_template_row_choices <- function(template) {
  labels <- gate_template_row_label(template)
  stats::setNames(as.character(seq_along(labels)), labels)
}

gate_template_selected_index <- function(template, gate_alias = NULL) {
  if (nrow(template) == 0L) {
    return(NA_integer_)
  }
  if (!is.null(gate_alias) && gate_alias %in% template$alias) {
    return(match(gate_alias, template$alias))
  }
  1L
}

split_gate_dims <- function(dims) {
  dims <- trimws(unlist(strsplit(as.character(dims), ",")))
  dims[nzchar(dims)]
}

gate_designer_row_bounds <- function(row) {
  if (nrow(row) != 1L || !identical(row$gating_method, "boundary")) {
    return(NULL)
  }
  parse_boundary_gating_args(row$gating_args)
}

gate_designer_row_polygon <- function(row) {
  if (nrow(row) != 1L || !identical(row$gating_method, "polygon_gate")) {
    return(NULL)
  }
  parse_polygon_gating_args(row$gating_args)
}

build_gate_template_row <- function(alias, pop, parent, dims, gating_method,
                                    gating_args = "",
                                    collapseDataForGating = "",
                                    groupBy = "",
                                    preprocessing_method = "",
                                    preprocessing_args = "") {
  data.table::data.table(
    alias = alias,
    pop = pop,
    parent = parent,
    dims = paste(split_gate_dims(dims), collapse = ","),
    gating_method = gating_method,
    gating_args = gating_args,
    collapseDataForGating = collapseDataForGating,
    groupBy = groupBy,
    preprocessing_method = preprocessing_method,
    preprocessing_args = preprocessing_args
  )
}

update_gate_template_row <- function(template, row_index = NULL, row,
                                     append = FALSE) {
  template <- data.table::as.data.table(data.table::copy(template))
  missing_cols <- setdiff(gating_template_columns(), names(template))
  for (col in missing_cols) {
    template[, (col) := ""]
  }
  data.table::setcolorder(template, gating_template_columns())
  for (col in gating_template_columns()) {
    template[, (col) := as.character(get(col))]
    template[is.na(get(col)), (col) := ""]
  }

  row <- data.table::as.data.table(row)
  missing_row_cols <- setdiff(gating_template_columns(), names(row))
  for (col in missing_row_cols) {
    row[, (col) := ""]
  }
  data.table::setcolorder(row, gating_template_columns())
  for (col in gating_template_columns()) {
    row[, (col) := as.character(get(col))]
    row[is.na(get(col)), (col) := ""]
  }

  if (isTRUE(append) || is.null(row_index) || is.na(row_index) ||
      row_index < 1L || row_index > nrow(template)) {
    out <- data.table::rbindlist(list(template, row), use.names = TRUE,
                                 fill = TRUE)
  } else {
    template[row_index] <- row
    out <- template
  }

  data.table::setcolorder(out, gating_template_columns())
  out
}

gate_designer_default_dims <- function(channels, dims = NULL) {
  dims <- split_gate_dims(dims)
  if (length(dims) >= 2L && all(dims[1:2] %in% channels)) {
    return(dims[1:2])
  }
  preferred <- c("FSC-A", "SSC-A")
  if (all(preferred %in% channels)) {
    return(preferred)
  }
  channels[seq_len(min(2L, length(channels)))]
}

gate_designer_method_choices <- function() {
  c("boundary", "polygon_gate", "singletGate", "mixture_gate", "density_gate")
}

gate_designer_axis_type <- function(scale) {
  if (identical(scale, "log10")) {
    "log"
  } else {
    "linear"
  }
}

gate_designer_scientific_labels <- function(values) {
  formatC(values, format = "e", digits = 1)
}

gate_designer_log_ticks <- function(values) {
  positive_range <- range(values[values > 0], finite = TRUE)
  if (!all(is.finite(positive_range)) || positive_range[[1L]] <= 0) {
    return(list(tickvals = NULL, ticktext = NULL))
  }

  powers <- seq(floor(log10(positive_range[[1L]])),
                ceiling(log10(positive_range[[2L]])))
  tickvals <- 10^powers
  tickvals <- tickvals[tickvals >= positive_range[[1L]] &
                         tickvals <= positive_range[[2L]]]
  if (length(tickvals) == 0L) {
    tickvals <- 10^powers
  }

  list(
    tickvals = tickvals,
    ticktext = gate_designer_scientific_labels(tickvals)
  )
}

gate_designer_axis_layout <- function(title, values, scale) {
  axis_range <- range(values, finite = TRUE)
  axis_pad <- diff(axis_range) * 0.05
  if (!is.finite(axis_pad) || axis_pad == 0) axis_pad <- 1

  if (identical(scale, "log10")) {
    positive_range <- range(values[values > 0], finite = TRUE)
    if (all(is.finite(positive_range)) && positive_range[[1L]] > 0) {
      positive_pad <- diff(log10(positive_range)) * 0.05
      if (!is.finite(positive_pad) || positive_pad == 0) positive_pad <- 0.1
      ticks <- gate_designer_log_ticks(values)
      return(list(
        title = title,
        type = "log",
        range = c(log10(positive_range[[1L]]) - positive_pad,
                  log10(positive_range[[2L]]) + positive_pad),
        exponentformat = "e",
        tickmode = "array",
        tickvals = ticks$tickvals,
        ticktext = ticks$ticktext,
        zeroline = FALSE
      ))
    }
  }

  list(
    title = title,
    type = "linear",
    range = c(axis_range[[1L]] - axis_pad, axis_range[[2L]] + axis_pad),
    exponentformat = "e",
    tickformat = ".2e",
    zeroline = FALSE
  )
}

gate_designer_point_density <- function(x, y, bins = 120L) {
  if (length(x) != length(y)) {
    stop("x and y must have the same length.", call. = FALSE)
  }
  density <- rep(NA_real_, length(x))
  finite <- is.finite(x) & is.finite(y)
  if (!any(finite)) {
    return(density)
  }

  x_finite <- x[finite]
  y_finite <- y[finite]
  x_breaks <- seq(min(x_finite), max(x_finite), length.out = bins + 1L)
  y_breaks <- seq(min(y_finite), max(y_finite), length.out = bins + 1L)
  if (length(unique(x_breaks)) < 2L || length(unique(y_breaks)) < 2L) {
    density[finite] <- 1
    return(density)
  }

  bin_dt <- data.table::data.table(
    row_id = which(finite),
    x_bin = findInterval(x_finite, x_breaks, all.inside = TRUE),
    y_bin = findInterval(y_finite, y_breaks, all.inside = TRUE)
  )
  bin_dt[, density := .N, by = .(x_bin, y_bin)]
  density[bin_dt$row_id] <- log1p(bin_dt$density)
  density
}

check_gate_designer_dependencies <- function() {
  missing <- c("shiny", "plotly", "DT")[
    !vapply(c("shiny", "plotly", "DT"), requireNamespace, logical(1),
            quietly = TRUE)
  ]
  if (length(missing) > 0L) {
    stop("design_gating_template() requires suggested packages: ",
         paste(missing, collapse = ", "), ".", call. = FALSE)
  }
}

gate_designer_sample_data <- function(cs, sample, dims, max_events) {
  sample_names <- flowCore::sampleNames(cs)
  sample_index <- match(sample, sample_names)
  if (is.na(sample_index)) {
    stop("Sample '", sample, "' is not present in the cytoset.", call. = FALSE)
  }

  fr <- cs[[sample_index]]
  mat <- flowCore::exprs(fr)
  missing_dims <- setdiff(dims, colnames(mat))
  if (length(missing_dims) > 0L) {
    stop("Missing gating channel(s): ", paste(missing_dims, collapse = ", "),
         call. = FALSE)
  }

  dt <- data.table::as.data.table(mat[, dims, drop = FALSE])
  data.table::setnames(dt, dims, c("x", "y"))
  dt <- dt[is.finite(x) & is.finite(y)]
  if (is.finite(max_events) && nrow(dt) > max_events) {
    dt <- dt[sample.int(.N, max_events)]
  }
  dt
}

gate_designer_register_plugins <- function() {
  openCyto::register_plugins(fun = density_gate, "density_gate", dep = NA,
                             "gating")
  openCyto::register_plugins(fun = mixture_gate, "mixture_gate", dep = NA,
                             "gating")
  openCyto::register_plugins(fun = polygon_gate, "polygon_gate", dep = NA,
                             "gating")
}

gate_designer_parent_sample_data <- function(cs, template, row_index, sample,
                                             parent, dims, max_events,
                                             use_parent = TRUE) {
  if (!isTRUE(use_parent) || is.na(row_index) || row_index <= 1L ||
      is.null(parent) || is.na(parent) || identical(parent, "root") ||
      !nzchar(parent)) {
    return(gate_designer_sample_data(cs, sample, dims, max_events))
  }

  prior_template <- template[seq_len(row_index - 1L)]
  if (!parent %in% prior_template$alias) {
    return(gate_designer_sample_data(cs, sample, dims, max_events))
  }

  template_file <- tempfile(fileext = ".csv")
  on.exit(unlink(template_file), add = TRUE)
  data.table::fwrite(prior_template, template_file)

  sample_index <- match(sample, flowCore::sampleNames(cs))
  if (is.na(sample_index)) {
    return(gate_designer_sample_data(cs, sample, dims, max_events))
  }

  gated_data <- tryCatch({
    gate_designer_register_plugins()
    gt <- openCyto::gatingTemplate(template_file)
    gs <- flowWorkspace::GatingSet(cs[sample_index])
    openCyto::gt_gating(gt, gs)
    parent_cs <- flowWorkspace::gs_pop_get_data(gs, y = parent)
    gate_designer_sample_data(parent_cs, sample, dims, max_events)
  }, error = function(e) NULL)

  if (is.null(gated_data)) {
    return(gate_designer_sample_data(cs, sample, dims, max_events))
  }
  gated_data
}

extract_plotly_rect_bounds <- function(relayout) {
  if (is.null(relayout) || length(relayout) == 0L) {
    return(NULL)
  }

  relayout_shapes <- relayout[["shapes"]]
  if (!is.null(relayout_shapes) && length(relayout_shapes) > 0L) {
    shape_index <- which(vapply(relayout_shapes, function(shape) {
      isTRUE(shape$editable) &&
        (is.null(shape$type) || identical(shape$type, "rect"))
    }, logical(1L)))
    if (length(shape_index) == 0L) {
      shape_index <- length(relayout_shapes)
    } else {
      shape_index <- tail(shape_index, 1L)
    }
    shape <- relayout_shapes[[shape_index]]
    x <- as.numeric(c(shape$x0, shape$x1))
    y <- as.numeric(c(shape$y0, shape$y1))
  } else {
    names_relayout <- names(relayout)
    shape_names <- grep("^shapes\\[[0-9]+\\]\\.(x0|x1|y0|y1)$",
                        names_relayout, value = TRUE)
    if (length(shape_names) < 4L) {
      return(NULL)
    }

    shape_ids <- sub("^shapes\\[([0-9]+)\\].*$", "\\1", shape_names)
    shape_id <- tail(unique(shape_ids), 1L)
    get_value <- function(axis) {
      as.numeric(relayout[[paste0("shapes[", shape_id, "].", axis)]])
    }

    x <- c(get_value("x0"), get_value("x1"))
    y <- c(get_value("y0"), get_value("y1"))
  }

  if (any(!is.finite(c(x, y)))) {
    return(NULL)
  }

  list(min = c(min(x), min(y)), max = c(max(x), max(y)))
}

gate_designer_polygon_path <- function(vertices) {
  vertices <- as.matrix(vertices)
  if (nrow(vertices) < 3L || ncol(vertices) != 2L) {
    return(NULL)
  }
  paste0(
    "M ", paste(vertices[1L, ], collapse = ","),
    paste0(" L ", apply(vertices[-1L, , drop = FALSE], 1L,
                        paste, collapse = ","), collapse = ""),
    " Z"
  )
}

gate_designer_rgba <- function(hex, alpha = 1) {
  rgb <- grDevices::col2rgb(hex)
  paste0(
    "rgba(",
    paste(as.integer(rgb[, 1L]), collapse = ", "),
    ", ", alpha, ")"
  )
}

gate_designer_overlay_colors <- function(n) {
  pal <- expressalyzr_palette()
  colors <- c(pal[["blue"]], pal[["green"]], pal[["amber"]],
              pal[["plum"]], pal[["ink"]], pal[["teal"]])
  rep(colors, length.out = n)
}

parse_svg_polygon_path <- function(path) {
  parse_plotly_polygon_path(path)
}

parse_plotly_polygon_path <- function(path) {
  if (is.null(path) || length(path) != 1L || is.na(path) || !nzchar(path)) {
    return(NULL)
  }

  cleaned <- gsub("[MLZmlz]", " ", path)
  nums <- suppressWarnings(as.numeric(
    unlist(regmatches(cleaned, gregexpr(
      "[-+]?(?:[0-9]*\\.?[0-9]+)(?:[eE][-+]?[0-9]+)?",
      cleaned
    )))
  ))
  if (length(nums) < 6L || length(nums) %% 2L != 0L ||
      any(!is.finite(nums))) {
    return(NULL)
  }

  vertices <- matrix(nums, ncol = 2L, byrow = TRUE)
  if (nrow(vertices) > 3L &&
      isTRUE(all.equal(vertices[1L, ], vertices[nrow(vertices), ],
                       tolerance = 1e-8, check.attributes = FALSE))) {
    vertices <- vertices[-nrow(vertices), , drop = FALSE]
  }
  if (nrow(vertices) < 3L) {
    return(NULL)
  }
  colnames(vertices) <- c("x", "y")
  vertices
}

gate_designer_pixels_to_data <- function(vertices, width, height, xaxis, yaxis) {
  if (!gate_designer_polygon_is_valid(vertices) ||
      !is.finite(width) || !is.finite(height) ||
      width <= 0 || height <= 0 ||
      length(xaxis$range) != 2L || length(yaxis$range) != 2L) {
    return(NULL)
  }

  x_axis <- xaxis$range[[1L]] +
    vertices[, 1L] / width * diff(xaxis$range)
  y_axis <- yaxis$range[[2L]] -
    vertices[, 2L] / height * diff(yaxis$range)

  x <- if (identical(xaxis$type, "log")) 10^x_axis else x_axis
  y <- if (identical(yaxis$type, "log")) 10^y_axis else y_axis
  out <- cbind(x, y)
  colnames(out) <- c("x", "y")
  out
}

extract_plotly_polygon_vertices <- function(relayout) {
  if (is.null(relayout) || length(relayout) == 0L) {
    return(NULL)
  }

  relayout_shapes <- relayout[["shapes"]]
  if (!is.null(relayout_shapes) && length(relayout_shapes) > 0L) {
    editable_shapes <- relayout_shapes[vapply(relayout_shapes, function(shape) {
      isTRUE(shape$editable) && !is.null(shape$path)
    }, logical(1L))]
    if (length(editable_shapes) > 0L) {
      relayout_shapes <- editable_shapes
    }
    paths <- vapply(relayout_shapes, function(shape) {
      if (!is.null(shape$path)) shape$path else NA_character_
    }, character(1L))
    paths <- paths[nzchar(paths) & !is.na(paths)]
    if (length(paths) == 0L) {
      return(NULL)
    }
    return(parse_plotly_polygon_path(tail(paths, 1L)))
  }

  path_names <- grep("^shapes\\[[0-9]+\\]\\.path$", names(relayout),
                     value = TRUE)
  if (length(path_names) == 0L) {
    return(NULL)
  }
  parse_plotly_polygon_path(as.character(relayout[[tail(path_names, 1L)]]))
}

gate_designer_data_bounds <- function(dt) {
  list(
    min = c(
      as.numeric(stats::quantile(dt$x, 0.01, na.rm = TRUE)),
      as.numeric(stats::quantile(dt$y, 0.01, na.rm = TRUE))
    ),
    max = c(
      as.numeric(stats::quantile(dt$x, 0.99, na.rm = TRUE)),
      as.numeric(stats::quantile(dt$y, 0.99, na.rm = TRUE))
    )
  )
}

gate_designer_bounds_are_valid <- function(bounds) {
  !is.null(bounds) &&
    is.list(bounds) &&
    length(bounds$min) == 2L &&
    length(bounds$max) == 2L &&
    all(is.finite(c(bounds$min, bounds$max))) &&
    all(bounds$min < bounds$max)
}

gate_designer_bounds_equal <- function(x, y) {
  if (!gate_designer_bounds_are_valid(x) ||
      !gate_designer_bounds_are_valid(y)) {
    return(FALSE)
  }
  isTRUE(all.equal(c(x$min, x$max), c(y$min, y$max),
                   tolerance = 1e-8, check.attributes = FALSE))
}

gate_designer_polygon_is_valid <- function(vertices) {
  !is.null(vertices) &&
    is.matrix(vertices) &&
    ncol(vertices) == 2L &&
    nrow(vertices) >= 3L &&
    all(is.finite(vertices))
}

gate_designer_polygons_equal <- function(x, y) {
  if (is.null(x) && is.null(y)) {
    return(TRUE)
  }
  if (!gate_designer_polygon_is_valid(x) ||
      !gate_designer_polygon_is_valid(y) ||
      !identical(dim(x), dim(y))) {
    return(FALSE)
  }
  isTRUE(all.equal(x, y, tolerance = 1e-8, check.attributes = FALSE))
}

gate_designer_points_in_polygon <- function(x, y, vertices) {
  if (!gate_designer_polygon_is_valid(vertices)) {
    return(rep(FALSE, length(x)))
  }

  px <- vertices[, 1L]
  py <- vertices[, 2L]
  inside <- rep(FALSE, length(x))
  j <- length(px)
  for (i in seq_along(px)) {
    intersects <- ((py[[i]] > y) != (py[[j]] > y)) &
      (x < (px[[j]] - px[[i]]) * (y - py[[i]]) /
         (py[[j]] - py[[i]] + .Machine$double.eps) + px[[i]])
    inside <- xor(inside, intersects)
    j <- i
  }
  inside
}

gate_designer_gate_summary <- function(dt, bounds) {
  if (!gate_designer_bounds_are_valid(bounds)) {
    return(list(total = nrow(dt), inside = NA_integer_, percent = NA_real_))
  }

  inside <- sum(
    dt$x >= bounds$min[[1L]] & dt$x <= bounds$max[[1L]] &
      dt$y >= bounds$min[[2L]] & dt$y <= bounds$max[[2L]],
    na.rm = TRUE
  )
  total <- nrow(dt)
  percent <- if (total > 0L) 100 * inside / total else NA_real_

  list(total = total, inside = inside, percent = percent)
}

gate_designer_polygon_summary <- function(dt, vertices) {
  if (!gate_designer_polygon_is_valid(vertices)) {
    return(list(total = nrow(dt), inside = NA_integer_, percent = NA_real_))
  }

  inside <- sum(gate_designer_points_in_polygon(dt$x, dt$y, vertices),
                na.rm = TRUE)
  total <- nrow(dt)
  percent <- if (total > 0L) 100 * inside / total else NA_real_
  list(total = total, inside = inside, percent = percent)
}

gate_designer_lasso_vertices <- function(selected) {
  if (is.null(selected) || nrow(selected) < 3L ||
      !all(c("x", "y") %in% names(selected))) {
    return(NULL)
  }
  vertices <- as.matrix(selected[, c("x", "y"), drop = FALSE])
  vertices <- vertices[stats::complete.cases(vertices), , drop = FALSE]
  if (nrow(vertices) < 3L) {
    return(NULL)
  }
  hull <- grDevices::chull(vertices[, 1L], vertices[, 2L])
  hull_vertices <- vertices[hull, , drop = FALSE]
  colnames(hull_vertices) <- c("x", "y")
  hull_vertices
}

gate_designer_control_sample_names <- function(data_path, sample_names) {
  config_file <- file.path(data_path, "config.yml")
  if (!file.exists(config_file)) {
    return(character())
  }

  config <- tryCatch(config::get(file = config_file), error = function(e) NULL)
  controls_pattern <- config$controls_pattern
  if (is.null(controls_pattern) || is.na(controls_pattern) ||
      !nzchar(controls_pattern)) {
    return(character())
  }

  sample_names[grepl(controls_pattern, sample_names)]
}

gate_designer_update_bound_inputs <- function(session, bounds) {
  if (!gate_designer_bounds_are_valid(bounds)) {
    return(invisible(NULL))
  }

  shiny::updateNumericInput(session, "x_min", value = bounds$min[[1L]])
  shiny::updateNumericInput(session, "x_max", value = bounds$max[[1L]])
  shiny::updateNumericInput(session, "y_min", value = bounds$min[[2L]])
  shiny::updateNumericInput(session, "y_max", value = bounds$max[[2L]])
  invisible(NULL)
}

gate_designer_template_specs <- function(output_file, template_file,
                                         template_files = NULL,
                                         output_files = NULL,
                                         template_labels = NULL) {
  if (is.null(template_files)) {
    template_files <- template_file
    output_files <- output_file
  } else if (is.null(output_files)) {
    output_files <- template_files
  }

  if (!is.character(template_files) || length(template_files) == 0L ||
      any(!nzchar(template_files))) {
    stop("template_files must be a non-empty character vector.",
         call. = FALSE)
  }
  if (!is.character(output_files) ||
      length(output_files) != length(template_files) ||
      any(!nzchar(output_files))) {
    stop("output_files must have one non-empty path per template file.",
         call. = FALSE)
  }

  if (is.null(template_labels)) {
    template_labels <- tools::file_path_sans_ext(basename(output_files))
  }
  if (!is.character(template_labels) ||
      length(template_labels) != length(template_files) ||
      any(!nzchar(template_labels))) {
    stop("template_labels must have one non-empty label per template file.",
         call. = FALSE)
  }

  data.table::data.table(
    template_id = seq_along(template_files),
    template_label = make.unique(template_labels, sep = "_"),
    template_file = normalizePath(template_files, mustWork = FALSE),
    output_file = normalizePath(output_files, mustWork = FALSE)
  )
}

#' Launch a visual gating-template designer.
#'
#' Opens a Shiny/plotly app for editing \code{gt_samples.csv}. Boundary gates
#' can be adjusted visually on any selected template row; automated methods such
#' as \code{mixture_gate} and \code{singletGate} can be edited as template rows.
#' The plot can switch between two-channel scatter and one-channel histogram
#' counts views, and saved gates remain overlaid while switching samples so
#' users can check whether the gate generalizes across wells.
#'
#' @param data_path Experiment directory containing a \code{data/} folder of
#'   FCS files.
#' @param output_file Gating template path to write. Defaults to
#'   \code{gt_samples.csv} in \code{data_path}.
#' @param template_file Existing template path to seed the editor. Defaults to
#'   \code{output_file} when it exists, otherwise the package template.
#' @param template_files Optional character vector of template files to edit in
#'   one combined app. When provided, users can switch between templates in the
#'   sidebar and each template is saved to its matching \code{output_files}
#'   path.
#' @param output_files Optional output paths matching \code{template_files}.
#'   Defaults to \code{template_files}.
#' @param template_labels Optional user-facing labels for
#'   \code{template_files}, such as population names.
#' @param gate_alias Initial gate row alias to select.
#' @param parent Parent population used when creating a new boundary gate.
#' @param dims Initial two channels to plot and gate.
#' @param sample Optional sample name to show initially.
#' @param max_events Maximum events to plot per sample.
#' @param launch Whether to run the app. Set \code{FALSE} to return the Shiny
#'   app object for testing or embedding.
#'
#' @return Runs the app when \code{launch = TRUE}; otherwise returns a
#'   \code{shiny.appobj}.
#' @export
design_gating_template <- function(data_path,
                                   output_file = file.path(data_path, "gt_samples.csv"),
                                   template_file = if (file.exists(output_file)) {
                                     output_file
                                   } else {
                                     system.file("tools", "gt_samples.csv",
                                                 package = "expressalyzr",
                                                 mustWork = TRUE)
                                   },
                                   template_files = NULL,
                                   output_files = NULL,
                                   template_labels = NULL,
                                   gate_alias = "nondebris",
                                   parent = "root",
                                   dims = c("FSC-A", "SSC-A"),
                                   sample = NULL,
                                   max_events = 10000L,
                                   launch = interactive()) {
  check_gate_designer_dependencies()

  if (!dir.exists(data_path)) {
    stop("data_path does not exist.", call. = FALSE)
  }
  if (!is.character(dims) || length(dims) != 2L) {
    stop("dims must be a character vector of two channel names.",
         call. = FALSE)
  }
  designer_templates <- gate_designer_template_specs(
    output_file = output_file,
    template_file = template_file,
    template_files = template_files,
    output_files = output_files,
    template_labels = template_labels
  )
  template_list <- lapply(designer_templates$template_file, read_gating_template)

  config_file <- file.path(data_path, "config.yml")
  channel_map <- "cytoflex"
  if (file.exists(config_file)) {
    config <- config::get(file = config_file)
    channel_map <- config_channel_map(config)
  }
  cs <- normalize_channels(load_fcs(data_path), channel_map)
  sample_names <- flowCore::sampleNames(cs)
  control_sample_names <- gate_designer_control_sample_names(data_path,
                                                             sample_names)
  channels <- flowWorkspace::colnames(cs)
  if (is.null(sample)) {
    sample <- sample_names[1L]
  }
  template <- template_list[[1L]]
  selected_index <- gate_template_selected_index(template, gate_alias)
  selected_row <- if (!is.na(selected_index)) template[selected_index] else NULL
  selected_dims <- if (!is.null(selected_row)) {
    gate_designer_default_dims(channels, selected_row$dims)
  } else {
    gate_designer_default_dims(channels, dims)
  }
  template_bounds <- if (!is.null(selected_row)) {
    gate_designer_row_bounds(selected_row)
  } else {
    NULL
  }
  template_polygon <- if (!is.null(selected_row)) {
    gate_designer_row_polygon(selected_row)
  } else {
    NULL
  }
  initial_dt <- gate_designer_sample_data(cs, sample, selected_dims, max_events)
  initial_bounds <- if (gate_designer_bounds_are_valid(template_bounds)) {
    template_bounds
  } else {
    gate_designer_data_bounds(initial_dt)
  }
  initial_polygon <- if (gate_designer_polygon_is_valid(template_polygon)) {
    template_polygon
  } else {
    NULL
  }
  initial_alias <- if (!is.null(selected_row)) selected_row$alias else gate_alias
  initial_parent <- if (!is.null(selected_row)) selected_row$parent else parent
  initial_method <- if (!is.null(selected_row)) selected_row$gating_method else "boundary"
  if (is.na(initial_method) || !nzchar(initial_method)) initial_method <- "boundary"
  method_choices <- unique(c(gate_designer_method_choices(), initial_method))
  initial_args <- if (!is.null(selected_row)) selected_row$gating_args else ""
  initial_pop <- if (!is.null(selected_row)) selected_row$pop else "+"
  initial_collapse <- if (!is.null(selected_row)) selected_row$collapseDataForGating else ""
  initial_group <- if (!is.null(selected_row)) selected_row$groupBy else ""
  initial_preprocess <- if (!is.null(selected_row)) selected_row$preprocessing_method else ""
  initial_preprocess_args <- if (!is.null(selected_row)) selected_row$preprocessing_args else ""

  ui <- shiny::fluidPage(
    shiny::tags$head(
      shiny::tags$style(shiny::HTML("
        .gate-help {
          background: #F7F9FB;
          border-left: 4px solid #2A9D8F;
          padding: 10px 12px;
          margin-bottom: 12px;
        }
        .gate-help ol {
          padding-left: 18px;
          margin-bottom: 0;
        }
        .gate-bounds {
          background: #243447;
          color: #FFFFFF;
          border-radius: 4px;
          padding: 8px 10px;
          margin-top: 10px;
        }
        .gate-save {
          width: 100%;
          font-weight: 700;
        }
        .gate-actions {
          display: flex;
          gap: 6px;
          flex-wrap: wrap;
        }
        .gate-actions .btn {
          flex: 1 1 auto;
        }
        .gate-status {
          display: flex;
          gap: 8px;
          flex-wrap: wrap;
          margin-bottom: 12px;
        }
        .gate-pill {
          border-radius: 999px;
          padding: 4px 10px;
          background: #E8EDF2;
          color: #243447;
          font-size: 12px;
          font-weight: 700;
        }
        .gate-pill-unsaved {
          background: #FBE7E2;
          color: #8A3B31;
        }
        .gate-pill-saved {
          background: #E6F3EF;
          color: #1F6E63;
        }
        .gate-metric {
          background: #FFFFFF;
          border: 1px solid #E8EDF2;
          border-radius: 4px;
          padding: 8px 10px;
          margin-bottom: 10px;
        }
        .gate-metric strong {
          display: block;
          font-size: 18px;
          color: #243447;
        }
        .gate-file {
          font-size: 12px;
          word-break: break-all;
          color: #6C757D;
          margin-bottom: 10px;
        }
        .gate-warning {
          background: #FFF8E8;
          border-left: 4px solid #E9A23B;
          padding: 8px 10px;
          margin-bottom: 10px;
        }
        .gate-plot-wrap {
          position: relative;
        }
        .gate-channel-row {
          display: flex;
          gap: 8px;
          align-items: flex-start;
        }
        .gate-channel-field {
          flex: 1 1 auto;
          min-width: 0;
        }
        .gate-scale-field {
          flex: 0 0 92px;
        }
        .gate-advanced {
          border: 1px solid #E8EDF2;
          border-radius: 4px;
          padding: 8px 10px;
          margin-bottom: 12px;
        }
        .gate-advanced summary {
          color: #44515E;
          cursor: pointer;
          font-weight: 700;
        }
        .js-plotly-plot .shapelayer path {
          cursor: move;
        }
        .js-plotly-plot .select-outline,
        .js-plotly-plot .selectionlayer path {
          stroke: #D76A5A !important;
          stroke-width: 2.5px !important;
          stroke-dasharray: none !important;
          stroke-opacity: 1 !important;
        }
        .js-plotly-plot .draglayer .nsewdrag {
          cursor: move;
        }
        .js-plotly-plot .draglayer .ewdrag {
          cursor: ew-resize;
        }
        .js-plotly-plot .draglayer .nsdrag {
          cursor: ns-resize;
        }
        .js-plotly-plot .draglayer .nwdrag,
        .js-plotly-plot .draglayer .sedrag {
          cursor: nwse-resize;
        }
        .js-plotly-plot .draglayer .nedrag,
        .js-plotly-plot .draglayer .swdrag {
          cursor: nesw-resize;
        }
      "))
    ),
    shiny::titlePanel("expressalyzr Gate Designer"),
    shiny::sidebarLayout(
      shiny::sidebarPanel(
        shiny::tags$div(
          class = "gate-help",
          shiny::tags$strong("Edit the gating template"),
          shiny::tags$ol(
            shiny::tags$li("Select a gate row or add a new one."),
            shiny::tags$li("Draw the gate and edit gate parameters."),
            shiny::tags$li("Switch samples to check whether the gate generalizes."),
            shiny::tags$li("Click Update selected gate to write the selected row.")
          )
        ),
        if (nrow(designer_templates) > 1L) {
          shiny::selectInput(
            "template_set",
            "Population template",
            choices = stats::setNames(
              as.character(designer_templates$template_id),
              designer_templates$template_label
            ),
            selected = "1"
          )
        },
        shiny::uiOutput("save_status"),
        shiny::uiOutput("template_context"),
        shiny::uiOutput("gate_summary"),
        shiny::selectInput("gate_row", "Gate row",
                           choices = gate_template_row_choices(template),
                           selected = as.character(selected_index)),
        shiny::tags$div(
          class = "gate-actions",
          shiny::actionButton("add_gate", "Add new row"),
          shiny::actionButton("delete_gate", "Delete selected")
        ),
        shiny::textInput("gate_alias", "Alias", initial_alias),
        shiny::textInput("pop", "Population sign", initial_pop),
        shiny::textInput("parent", "Parent population", initial_parent),
        shiny::selectizeInput("gating_method", "Gating method",
                              choices = method_choices,
                              selected = initial_method,
                              options = list(create = TRUE)),
        shiny::conditionalPanel(
          condition = "input.gating_method != 'boundary' && input.gating_method != 'polygon_gate'",
          shiny::textAreaInput("gating_args", "Gating args",
                               value = initial_args, rows = 3)
        ),
        shiny::tags$details(
          class = "gate-advanced",
          shiny::tags$summary("Advanced OpenCyto fields"),
          shiny::textInput("collapseDataForGating", "collapseDataForGating",
                           initial_collapse),
          shiny::textInput("groupBy", "groupBy", initial_group),
          shiny::textInput("preprocessing_method", "preprocessing method",
                           initial_preprocess),
          shiny::textAreaInput("preprocessing_args", "preprocessing args",
                               value = initial_preprocess_args, rows = 2)
        ),
        shiny::selectInput("sample", "Sample", sample_names, selected = sample),
        shiny::tags$div(
          class = "gate-actions",
          shiny::actionButton("show_all_samples", "All samples"),
          shiny::actionButton("show_control_samples", "Controls only")
        ),
        shiny::radioButtons("plot_mode", "Plot",
                            choices = c("Scatter" = "scatter",
                                        "Histogram counts" = "histogram"),
                            selected = "scatter", inline = TRUE),
        shiny::conditionalPanel(
          condition = "input.plot_mode == 'scatter'",
          shiny::checkboxInput("density_overlay", "Color by density",
                               value = TRUE)
        ),
        shiny::tags$div(
          class = "gate-channel-row",
          shiny::tags$div(
            class = "gate-channel-field",
            shiny::selectInput("x_channel", "X channel", channels,
                               selected = selected_dims[[1L]])
          ),
          shiny::tags$div(
            class = "gate-scale-field",
            shiny::selectInput("x_scale", "Scale",
                               choices = c("linear", "log10"),
                               selected = "linear")
          )
        ),
        shiny::conditionalPanel(
          condition = "input.plot_mode == 'scatter'",
          shiny::tags$div(
            class = "gate-channel-row",
            shiny::tags$div(
              class = "gate-channel-field",
              shiny::selectInput("y_channel", "Y channel", channels,
                                 selected = selected_dims[[2L]])
            ),
            shiny::tags$div(
              class = "gate-scale-field",
              shiny::selectInput("y_scale", "Scale",
                                 choices = c("linear", "log10"),
                                 selected = "linear")
            )
          )
        ),
        shiny::checkboxInput("preview_parent", "Preview parent population",
                             value = TRUE),
        shiny::numericInput("max_events", "Events to plot", max_events,
                            min = 100L, step = 500L),
        shiny::tags$hr(),
        shiny::conditionalPanel(
          condition = "input.gating_method == 'boundary'",
          shiny::tags$h4("Boundary gate"),
          shiny::uiOutput("histogram_notice"),
          shiny::numericInput("x_min", "X min", initial_bounds$min[[1L]]),
          shiny::numericInput("x_max", "X max", initial_bounds$max[[1L]]),
          shiny::conditionalPanel(
            condition = "input.plot_mode == 'scatter'",
            shiny::numericInput("y_min", "Y min", initial_bounds$min[[2L]]),
            shiny::numericInput("y_max", "Y max", initial_bounds$max[[2L]])
          ),
          shiny::tags$div(
            class = "gate-actions",
            shiny::actionButton("fit_gate", "Fit to visible data"),
            shiny::actionButton("reset_saved_gate", "Reset saved gate"),
            shiny::actionButton("undo_gate", "Undo")
          )
        ),
        shiny::conditionalPanel(
          condition = "input.gating_method == 'polygon_gate'",
          shiny::tags$h4("Freeform gate"),
          shiny::uiOutput("polygon_notice"),
          shiny::tags$div(
            class = "gate-actions",
            shiny::actionButton("clear_polygon", "Clear drawing"),
            shiny::actionButton("reset_saved_polygon", "Reset saved polygon")
          )
        ),
        shiny::conditionalPanel(
          condition = "input.gating_method != 'boundary' && input.gating_method != 'polygon_gate'",
          shiny::uiOutput("boundary_notice")
        ),
        shiny::tags$hr(),
        shiny::actionButton("save_gate", "Update selected gate",
                            class = "btn-primary gate-save"),
        shiny::actionButton("finish_design", "Done and continue pipeline",
                            class = "btn-default gate-save"),
        shiny::tags$div(class = "gate-bounds",
                        shiny::verbatimTextOutput("bounds_text"))
      ),
      shiny::mainPanel(
        shiny::tags$div(
          class = "gate-plot-wrap",
          plotly::plotlyOutput("gate_plot", height = "650px")
        ),
        shiny::tags$h4("Current template"),
        DT::DTOutput("template_table")
      )
    )
  )

  server <- function(input, output, session) {
    current_bounds <- shiny::reactiveVal(initial_bounds)
    saved_bounds <- shiny::reactiveVal(initial_bounds)
    current_polygon <- shiny::reactiveVal(initial_polygon)
    saved_polygon <- shiny::reactiveVal(initial_polygon)
    has_unsaved_changes <- shiny::reactiveVal(FALSE)
    bounds_history <- shiny::reactiveVal(list())
    template_states <- shiny::reactiveVal(template_list)
    selected_row_indices <- shiny::reactiveVal(
      rep(selected_index, nrow(designer_templates))
    )
    syncing_fields <- shiny::reactiveVal(FALSE)

    current_template_index <- function() {
      if (nrow(designer_templates) == 1L || is.null(input$template_set)) {
        return(1L)
      }
      index <- suppressWarnings(as.integer(input$template_set))
      if (is.na(index) || index < 1L || index > nrow(designer_templates)) {
        return(1L)
      }
      index
    }

    current_template_spec <- function() {
      designer_templates[current_template_index()]
    }

    template_state <- function(value) {
      index <- current_template_index()
      states <- template_states()
      if (missing(value)) {
        return(states[[index]])
      }
      states[[index]] <- value
      template_states(states)
      invisible(value)
    }

    selected_row_index <- function(value) {
      index <- current_template_index()
      indices <- selected_row_indices()
      if (missing(value)) {
        return(indices[[index]])
      }
      indices[[index]] <- value
      selected_row_indices(indices)
      invisible(value)
    }

    update_sample_choices <- function(choices) {
      if (length(choices) == 0L) {
        return(invisible(NULL))
      }
      selected <- if (!is.null(input$sample) && input$sample %in% choices) {
        input$sample
      } else {
        choices[[1L]]
      }
      shiny::updateSelectInput(session, "sample", choices = choices,
                               selected = selected)
      invisible(NULL)
    }

    refresh_gate_row_choices <- function(selected = selected_row_index()) {
      choices <- gate_template_row_choices(template_state())
      shiny::updateSelectInput(session, "gate_row", choices = choices,
                               selected = as.character(selected))
    }

    output$template_context <- shiny::renderUI({
      spec <- current_template_spec()
      shiny::tags$div(
        class = "gate-file",
        shiny::tags$strong(spec$template_label),
        shiny::tags$br(),
        spec$output_file
      )
    })

    set_bounds <- function(bounds, record = TRUE, update_inputs = TRUE,
                           mark_unsaved = TRUE) {
      if (!gate_designer_bounds_are_valid(bounds)) {
        return(FALSE)
      }
      old_bounds <- current_bounds()
      if (gate_designer_bounds_equal(old_bounds, bounds)) {
        return(TRUE)
      }
      if (isTRUE(record) && gate_designer_bounds_are_valid(old_bounds)) {
        bounds_history(c(bounds_history(), list(old_bounds)))
      }
      current_bounds(bounds)
      if (isTRUE(update_inputs)) {
        gate_designer_update_bound_inputs(session, bounds)
      }
      if (isTRUE(mark_unsaved)) {
        has_unsaved_changes(!gate_designer_bounds_equal(bounds, saved_bounds()))
      }
      TRUE
    }

    set_polygon <- function(vertices, mark_unsaved = TRUE) {
      if (!is.null(vertices) && !gate_designer_polygon_is_valid(vertices)) {
        return(FALSE)
      }
      current_polygon(vertices)
      if (isTRUE(mark_unsaved)) {
        has_unsaved_changes(
          !gate_designer_polygons_equal(vertices, saved_polygon())
        )
      }
      TRUE
    }

    mark_gate_saved <- function() {
      has_unsaved_changes(FALSE)
      session$onFlushed(function() {
        has_unsaved_changes(FALSE)
      }, once = TRUE)
      invisible(NULL)
    }

    load_row_into_inputs <- function(row_index) {
      template_now <- template_state()
      if (is.na(row_index) || row_index < 1L || row_index > nrow(template_now)) {
        return(invisible(NULL))
      }

      row <- template_now[row_index]
      row_dims <- gate_designer_default_dims(channels, row$dims)
      row_method <- row$gating_method
      if (is.na(row_method) || !nzchar(row_method)) row_method <- "boundary"

      syncing_fields(TRUE)
      shiny::updateTextInput(session, "gate_alias", value = row$alias)
      shiny::updateTextInput(session, "pop", value = row$pop)
      shiny::updateTextInput(session, "parent", value = row$parent)
      shiny::updateSelectizeInput(
        session,
        "gating_method",
        choices = unique(c(gate_designer_method_choices(), row_method)),
        selected = row_method,
        server = FALSE
      )
      shiny::updateTextAreaInput(session, "gating_args", value = row$gating_args)
      shiny::updateTextInput(session, "collapseDataForGating",
                             value = row$collapseDataForGating)
      shiny::updateTextInput(session, "groupBy", value = row$groupBy)
      shiny::updateTextInput(session, "preprocessing_method",
                             value = row$preprocessing_method)
      shiny::updateTextAreaInput(session, "preprocessing_args",
                                 value = row$preprocessing_args)
      shiny::updateSelectInput(session, "x_channel", selected = row_dims[[1L]])
      shiny::updateSelectInput(session, "y_channel", selected = row_dims[[2L]])
      syncing_fields(FALSE)

      row_bounds <- gate_designer_row_bounds(row)
      if (gate_designer_bounds_are_valid(row_bounds)) {
        saved_bounds(row_bounds)
        set_bounds(row_bounds, record = FALSE, mark_unsaved = FALSE)
      } else {
        bounds <- gate_designer_data_bounds(
          gate_designer_parent_sample_data(
            cs,
            template_state(),
            row_index,
            input$sample,
            row$parent,
            row_dims,
            input$max_events,
            use_parent = TRUE
          )
        )
        saved_bounds(bounds)
        set_bounds(bounds, record = FALSE, mark_unsaved = FALSE)
      }
      row_polygon <- gate_designer_row_polygon(row)
      if (gate_designer_polygon_is_valid(row_polygon)) {
        saved_polygon(row_polygon)
        set_polygon(row_polygon, mark_unsaved = FALSE)
      } else {
        saved_polygon(NULL)
        set_polygon(NULL, mark_unsaved = FALSE)
      }
      has_unsaved_changes(FALSE)
      invisible(NULL)
    }

    plot_data <- shiny::reactive({
      gate_designer_parent_sample_data(
        cs,
        template_state(),
        selected_row_index(),
        input$sample,
        input$parent,
        c(input$x_channel, input$y_channel),
        input$max_events,
        use_parent = input$preview_parent
      )
    })

    shiny::observeEvent(input$template_set, {
      row_index <- selected_row_index()
      template_now <- template_state()
      if (is.na(row_index) || row_index < 1L || row_index > nrow(template_now)) {
        row_index <- gate_template_selected_index(template_now, gate_alias)
        selected_row_index(row_index)
      }
      refresh_gate_row_choices(row_index)
      load_row_into_inputs(row_index)
      bounds_history(list())
      has_unsaved_changes(FALSE)
    }, ignoreInit = TRUE)

    shiny::observeEvent(input$gate_row, {
      row_index <- suppressWarnings(as.integer(input$gate_row))
      if (!is.na(row_index)) {
        selected_row_index(row_index)
        load_row_into_inputs(row_index)
      }
    }, ignoreInit = TRUE)

    shiny::observeEvent(input$show_all_samples, {
      update_sample_choices(sample_names)
      shiny::showNotification("Showing all samples.", type = "message")
    })

    shiny::observeEvent(input$show_control_samples, {
      if (length(control_sample_names) == 0L) {
        shiny::showNotification(
          "No control samples matched controls_pattern in config.yml.",
          type = "warning"
        )
        return(invisible(NULL))
      }
      update_sample_choices(control_sample_names)
      shiny::showNotification("Showing control samples only.", type = "message")
    })

    shiny::observeEvent(input$add_gate, {
      template_now <- template_state()
      new_alias <- paste0("gate_", nrow(template_now) + 1L)
      new_row <- build_gate_template_row(
        alias = new_alias,
        pop = "+",
        parent = if (nrow(template_now) > 0L) tail(template_now$alias, 1L) else "root",
        dims = c(input$x_channel, input$y_channel),
        gating_method = "boundary",
        gating_args = build_boundary_gating_args(current_bounds()$min,
                                                 current_bounds()$max)
      )
      updated <- update_gate_template_row(template_now, row = new_row,
                                          append = TRUE)
      template_state(updated)
      selected_row_index(nrow(updated))
      refresh_gate_row_choices(nrow(updated))
      load_row_into_inputs(nrow(updated))
      has_unsaved_changes(TRUE)
    })

    shiny::observeEvent(input$delete_gate, {
      template_now <- template_state()
      row_index <- selected_row_index()
      if (is.na(row_index) || nrow(template_now) <= 1L) {
        shiny::showNotification("Keep at least one gate row in the template.",
                                type = "warning")
        return(invisible(NULL))
      }
      updated <- template_now[-row_index]
      template_state(updated)
      new_index <- min(row_index, nrow(updated))
      selected_row_index(new_index)
      refresh_gate_row_choices(new_index)
      load_row_into_inputs(new_index)
      has_unsaved_changes(TRUE)
    })

    shiny::observeEvent(
      list(input$gate_alias, input$pop, input$parent, input$gating_method,
           input$gating_args, input$collapseDataForGating, input$groupBy,
           input$preprocessing_method, input$preprocessing_args,
           input$x_channel, input$y_channel),
      {
        if (isTRUE(syncing_fields())) {
          return(invisible(NULL))
        }
        has_unsaved_changes(TRUE)
      },
      ignoreInit = TRUE
    )

    output$save_status <- shiny::renderUI({
      if (isTRUE(has_unsaved_changes())) {
        shiny::tags$div(
          class = "gate-status",
          shiny::tags$span(class = "gate-pill gate-pill-unsaved",
                           "Unsaved gate changes")
        )
      } else {
        shiny::tags$div(
          class = "gate-status",
          shiny::tags$span(class = "gate-pill gate-pill-saved",
                           "Gate matches gt_samples.csv")
        )
      }
    })

    output$gate_summary <- shiny::renderUI({
      if (identical(input$gating_method, "polygon_gate")) {
        summary <- gate_designer_polygon_summary(plot_data(), current_polygon())
        percent <- if (is.finite(summary$percent)) {
          paste0(format(round(summary$percent, 1), nsmall = 1), "%")
        } else {
          "Draw gate"
        }
        return(shiny::tags$div(
          class = "gate-metric",
          shiny::tags$span("Events inside freeform gate"),
          shiny::tags$strong(percent),
          shiny::tags$span(
            paste0(
              if (is.na(summary$inside)) {
                "Draw at least three points"
              } else {
                paste0(format(summary$inside, big.mark = ","), " of ",
                       format(summary$total, big.mark = ","),
                       " plotted events")
              },
              if (isTRUE(input$preview_parent) &&
                  !identical(input$parent, "root")) {
                paste0(" from parent '", input$parent, "'")
              } else {
                ""
              }
            )
          )
        ))
      }

      if (!identical(input$gating_method, "boundary")) {
        return(shiny::tags$div(
          class = "gate-metric",
          shiny::tags$span("Selected template row"),
          shiny::tags$strong(input$gating_method),
          shiny::tags$span("Saved as OpenCyto template fields")
        ))
      }

      summary <- gate_designer_gate_summary(plot_data(), current_bounds())
      percent <- if (is.finite(summary$percent)) {
        paste0(format(round(summary$percent, 1), nsmall = 1), "%")
      } else {
        "NA"
      }

      shiny::tags$div(
        class = "gate-metric",
        shiny::tags$span("Events inside current gate"),
        shiny::tags$strong(percent),
        shiny::tags$span(
          paste0(format(summary$inside, big.mark = ","), " of ",
                 format(summary$total, big.mark = ","), " plotted events",
                 if (isTRUE(input$preview_parent) &&
                     !identical(input$parent, "root")) {
                   paste0(" from parent '", input$parent, "'")
                 } else {
                   ""
                 })
        )
      )
    })

    output$boundary_notice <- shiny::renderUI({
      if (identical(input$gating_method, "boundary")) {
        return(NULL)
      }
      shiny::tags$div(
        class = "gate-warning",
        "Visual drag editing is currently available for boundary and polygon gates. ",
        "Other methods are saved as template rows and can be previewed by running the pipeline/report."
      )
    })

    output$histogram_notice <- shiny::renderUI({
      if (!identical(input$plot_mode, "histogram")) {
        return(NULL)
      }
      shiny::tags$div(
        class = "gate-warning",
        "Histogram mode shows one-dimensional counts for the X channel. ",
        "Saving still writes the selected two-channel template row."
      )
    })

    output$polygon_notice <- shiny::renderUI({
      shiny::tags$div(
        class = "gate-warning",
        "Lasso events directly on the scatter plot. ",
        "The selected events are converted to a polygon_gate hull; click Update selected gate to save it."
      )
    })

    output$gate_plot <- plotly::renderPlotly({
      dt <- plot_data()
      bounds <- current_bounds()
      pal <- expressalyzr_palette()
      shapes <- list()
      annotations <- list()
      is_boundary <- identical(input$gating_method, "boundary")
      is_polygon <- identical(input$gating_method, "polygon_gate")
      is_histogram <- identical(input$plot_mode, "histogram")

      active_index <- current_template_index()
      all_templates <- template_states()
      overlay_colors <- gate_designer_overlay_colors(length(all_templates))
      active_alias <- input$gate_alias
      active_dims <- c(input$x_channel, input$y_channel)
      for (template_index in seq_along(all_templates)) {
        if (template_index == active_index) next
        template_now <- all_templates[[template_index]]
        if (nrow(template_now) == 0L) next

        match_index <- match(active_alias, template_now$alias)
        if (is.na(match_index)) {
          selected_indices <- selected_row_indices()
          candidate <- selected_indices[[template_index]]
          if (!is.na(candidate) && candidate >= 1L &&
              candidate <= nrow(template_now)) {
            match_index <- candidate
          }
        }
        if (is.na(match_index)) next

        overlay_row <- template_now[match_index]
        overlay_dims <- gate_designer_default_dims(channels, overlay_row$dims)
        if (!identical(overlay_dims, active_dims)) next

        overlay_color <- overlay_colors[[template_index]]
        overlay_label <- designer_templates$template_label[[template_index]]
        if (identical(overlay_row$gating_method, "boundary")) {
          overlay_bounds <- gate_designer_row_bounds(overlay_row)
          if (!gate_designer_bounds_are_valid(overlay_bounds)) next
          if (is_histogram) {
            shapes <- c(shapes, list(list(
              type = "rect",
              xref = "x",
              yref = "paper",
              x0 = overlay_bounds$min[[1L]],
              x1 = overlay_bounds$max[[1L]],
              y0 = 0,
              y1 = 1,
              editable = FALSE,
              line = list(color = overlay_color, width = 2),
              fillcolor = gate_designer_rgba(overlay_color, 0.06)
            )))
            annotations <- c(annotations, list(list(
              x = mean(c(overlay_bounds$min[[1L]], overlay_bounds$max[[1L]])),
              y = 1,
              xref = "x",
              yref = "paper",
              text = overlay_label,
              showarrow = FALSE,
              yanchor = "bottom",
              font = list(color = overlay_color, size = 11)
            )))
          } else {
            shapes <- c(shapes, list(list(
              type = "rect",
              xref = "x",
              yref = "y",
              x0 = overlay_bounds$min[[1L]],
              x1 = overlay_bounds$max[[1L]],
              y0 = overlay_bounds$min[[2L]],
              y1 = overlay_bounds$max[[2L]],
              editable = FALSE,
              line = list(color = overlay_color, width = 2),
              fillcolor = gate_designer_rgba(overlay_color, 0.06)
            )))
            annotations <- c(annotations, list(list(
              x = overlay_bounds$max[[1L]],
              y = overlay_bounds$max[[2L]],
              xref = "x",
              yref = "y",
              text = overlay_label,
              showarrow = FALSE,
              xanchor = "left",
              yanchor = "bottom",
              font = list(color = overlay_color, size = 11)
            )))
          }
        } else if (!is_histogram &&
                   identical(overlay_row$gating_method, "polygon_gate")) {
          overlay_polygon <- gate_designer_row_polygon(overlay_row)
          if (!gate_designer_polygon_is_valid(overlay_polygon)) next
          shapes <- c(shapes, list(list(
            type = "path",
            xref = "x",
            yref = "y",
            path = gate_designer_polygon_path(overlay_polygon),
            editable = FALSE,
            line = list(color = overlay_color, width = 2),
            fillcolor = gate_designer_rgba(overlay_color, 0.06)
          )))
          annotations <- c(annotations, list(list(
            x = max(overlay_polygon[, 1L], na.rm = TRUE),
            y = max(overlay_polygon[, 2L], na.rm = TRUE),
            xref = "x",
            yref = "y",
            text = overlay_label,
            showarrow = FALSE,
            xanchor = "left",
            yanchor = "bottom",
            font = list(color = overlay_color, size = 11)
          )))
        }
      }

      if (is_boundary && gate_designer_bounds_are_valid(bounds) &&
          is_histogram) {
        shapes <- c(shapes, list(list(
          type = "rect",
          xref = "x",
          yref = "paper",
          x0 = bounds$min[[1L]],
          x1 = bounds$max[[1L]],
          y0 = 0,
          y1 = 1,
          editable = TRUE,
          line = list(color = pal[["coral"]], width = 2),
          fillcolor = "rgba(215, 106, 90, 0.12)"
        )))
      } else if (is_boundary && gate_designer_bounds_are_valid(bounds)) {
        shapes <- c(shapes, list(list(
          type = "rect",
          xref = "x",
          yref = "y",
          x0 = bounds$min[[1L]],
          x1 = bounds$max[[1L]],
          y0 = bounds$min[[2L]],
          y1 = bounds$max[[2L]],
          editable = TRUE,
          line = list(color = pal[["coral"]], width = 2),
          fillcolor = "rgba(215, 106, 90, 0.10)"
        )))
      } else if (is_polygon &&
                 gate_designer_polygon_is_valid(current_polygon()) &&
                 !is_histogram) {
        shapes <- c(shapes, list(list(
          type = "path",
          xref = "x",
          yref = "y",
          path = gate_designer_polygon_path(current_polygon()),
          editable = TRUE,
          line = list(color = pal[["coral"]], width = 2),
          fillcolor = "rgba(215, 106, 90, 0.12)"
        )))
      }

      x_values <- c(dt$x, bounds$min[[1L]], bounds$max[[1L]])
      y_values <- c(dt$y, bounds$min[[2L]], bounds$max[[2L]])

      if (is_histogram) {
        p <- plotly::plot_ly(
          dt,
          x = ~x,
          type = "histogram",
          source = "gate_designer",
          nbinsx = 80,
          marker = list(
            color = pal[["teal"]],
            opacity = 0.55,
            line = list(color = "rgba(255, 255, 255, 0.35)", width = 1)
          )
        )
      } else {
        x_density <- if (identical(input$x_scale, "log10")) {
          ifelse(dt$x > 0, log10(dt$x), NA_real_)
        } else {
          dt$x
        }
        y_density <- if (identical(input$y_scale, "log10")) {
          ifelse(dt$y > 0, log10(dt$y), NA_real_)
        } else {
          dt$y
        }
        density_values <- if (isTRUE(input$density_overlay)) {
          gate_designer_point_density(x_density, y_density)
        } else {
          rep(1, nrow(dt))
        }
        marker <- list(
          color = if (isTRUE(input$density_overlay)) {
            density_values
          } else {
            pal[["teal"]]
          },
          size = 4,
          opacity = 0.42
        )
        if (isTRUE(input$density_overlay)) {
          marker$colorscale <- list(
            list(0, "#DCEAE6"),
            list(0.45, pal[["teal"]]),
            list(1, "#243447")
          )
          marker$showscale <- FALSE
        }

        p <- plotly::plot_ly(source = "gate_designer")
        p <- plotly::add_trace(
          p,
          data = dt,
          x = ~x,
          y = ~y,
          type = "scatter",
          mode = "markers",
          marker = marker,
          selected = list(
            marker = list(
              opacity = 0.42
            )
          ),
          unselected = list(
            marker = list(
              opacity = 0.42
            )
          ),
          showlegend = FALSE,
          inherit = FALSE
        )
      }
      p <- plotly::layout(
        p,
        title = list(text = if (is_boundary && is_histogram) {
          "Drag the red interval to adjust X bounds"
        } else if (is_boundary) {
          "Drag the red gate or its edges to adjust bounds"
      } else if (is_polygon && is_histogram) {
          "Switch to Scatter to draw a freeform gate"
        } else if (is_polygon) {
          "Lasso events to define a freeform gate"
        } else {
          paste("Viewing", input$gating_method, "row")
        }),
        xaxis = gate_designer_axis_layout(input$x_channel, x_values,
                                          input$x_scale),
        yaxis = if (is_histogram) {
          list(title = "Counts", type = "linear", exponentformat = "e",
               tickformat = ".2e", zeroline = FALSE)
        } else {
          gate_designer_axis_layout(input$y_channel, y_values, input$y_scale)
        },
        dragmode = if (is_polygon && !is_histogram) "lasso" else "pan",
        newshape = list(
          line = list(color = pal[["coral"]], width = 2),
          fillcolor = "rgba(215, 106, 90, 0.12)"
        ),
        shapes = shapes,
        annotations = annotations,
        showlegend = FALSE
      )
      plotly::config(
        p,
        displaylogo = FALSE,
        editable = TRUE,
        edits = list(
          shapePosition = TRUE,
          annotationPosition = FALSE,
          annotationTail = FALSE,
          annotationText = FALSE,
          axisTitleText = FALSE,
          colorbarPosition = FALSE,
          legendPosition = FALSE,
          titleText = FALSE
        ),
        modeBarButtonsToRemove = if (is_polygon && !is_histogram) {
          list("autoScale2d", "toggleSpikelines")
        } else {
          list("select2d", "lasso2d", "autoScale2d", "toggleSpikelines")
        },
        toImageButtonOptions = list(filename = "expressalyzr_gate_designer")
      )
    })

    handle_gate_relayout <- function(relayout) {
      if (identical(input$gating_method, "boundary")) {
        bounds <- extract_plotly_rect_bounds(relayout)
        if (!is.null(bounds)) {
          if (identical(input$plot_mode, "histogram")) {
            old_bounds <- current_bounds()
            if (gate_designer_bounds_are_valid(old_bounds)) {
              bounds$min[[2L]] <- old_bounds$min[[2L]]
              bounds$max[[2L]] <- old_bounds$max[[2L]]
            }
          }
          set_bounds(bounds)
        }
      } else if (identical(input$gating_method, "polygon_gate")) {
        vertices <- extract_plotly_polygon_vertices(relayout)
        if (!is.null(vertices)) {
          set_polygon(vertices)
        }
      }
      invisible(NULL)
    }

    shiny::observeEvent(plotly::event_data("plotly_relayout",
                                           source = "gate_designer"), {
      handle_gate_relayout(
        plotly::event_data("plotly_relayout", source = "gate_designer")
      )
    })

    shiny::observeEvent(input$gate_plot_relayout_js, {
      handle_gate_relayout(input$gate_plot_relayout_js)
    })

    shiny::observeEvent(input$gate_plot_polygon_pixels, {
      payload <- input$gate_plot_polygon_pixels
      if (!identical(input$gating_method, "polygon_gate")) {
        return(invisible(NULL))
      }
      if (!is.null(payload$error)) {
        shiny::showNotification(payload$error, type = "warning")
        return(invisible(NULL))
      }

      pixel_vertices <- parse_svg_polygon_path(payload$path)
      if (!gate_designer_polygon_is_valid(pixel_vertices)) {
        shiny::showNotification("Draw a closed polygon before using the drawing.",
                                type = "warning")
        return(invisible(NULL))
      }

      dt <- plot_data()
      bounds <- current_bounds()
      x_values <- c(dt$x, bounds$min[[1L]], bounds$max[[1L]])
      y_values <- c(dt$y, bounds$min[[2L]], bounds$max[[2L]])
      vertices <- gate_designer_pixels_to_data(
        pixel_vertices,
        width = as.numeric(payload$width),
        height = as.numeric(payload$height),
        xaxis = gate_designer_axis_layout(input$x_channel, x_values,
                                          input$x_scale),
        yaxis = gate_designer_axis_layout(input$y_channel, y_values,
                                          input$y_scale)
      )
      if (!gate_designer_polygon_is_valid(vertices)) {
        shiny::showNotification("Could not convert the drawing into gate coordinates.",
                                type = "warning")
        return(invisible(NULL))
      }
      set_polygon(vertices)
      shiny::showNotification("Freeform drawing captured.", type = "message")
    })

    shiny::observeEvent(plotly::event_data("plotly_selected",
                                           source = "gate_designer"), {
      if (!identical(input$gating_method, "polygon_gate")) {
        return(invisible(NULL))
      }
      selected <- plotly::event_data("plotly_selected",
                                     source = "gate_designer")
      if (is.null(selected) || nrow(selected) < 3L ||
          !all(c("x", "y") %in% names(selected))) {
        return(invisible(NULL))
      }
      hull_vertices <- gate_designer_lasso_vertices(selected)
      if (!gate_designer_polygon_is_valid(hull_vertices)) {
        return(invisible(NULL))
      }
      set_polygon(hull_vertices)
      shiny::showNotification("Freeform lasso captured.", type = "message")
    })

    shiny::observeEvent(
      list(input$x_min, input$x_max, input$y_min, input$y_max),
      {
        bounds <- list(
          min = c(input$x_min, input$y_min),
          max = c(input$x_max, input$y_max)
        )
        set_bounds(bounds, update_inputs = FALSE)
      },
      ignoreInit = TRUE
    )

    shiny::observeEvent(
      list(input$x_channel, input$y_channel, input$sample, input$max_events,
           input$preview_parent, input$parent),
      {
        if (!identical(input$gating_method, "boundary")) {
          return(invisible(NULL))
        }
        bounds <- current_bounds()
        if (!gate_designer_bounds_are_valid(bounds)) {
          bounds <- gate_designer_data_bounds(plot_data())
          set_bounds(bounds, record = FALSE)
        }
      },
      ignoreInit = TRUE
    )

    shiny::observeEvent(input$undo_gate, {
      history <- bounds_history()
      if (length(history) == 0L) {
        shiny::showNotification("No previous gate bounds to undo.",
                                type = "message")
        return(invisible(NULL))
      }
      previous <- history[[length(history)]]
      bounds_history(history[-length(history)])
      set_bounds(previous, record = FALSE)
    })

    shiny::observeEvent(input$fit_gate, {
      bounds <- gate_designer_data_bounds(plot_data())
      set_bounds(bounds)
    })

    shiny::observeEvent(input$reset_saved_gate, {
      set_bounds(saved_bounds(), record = TRUE)
    })

    shiny::observeEvent(input$clear_polygon, {
      set_polygon(NULL)
    })

    shiny::observeEvent(input$reset_saved_polygon, {
      set_polygon(saved_polygon())
    })

    shiny::observeEvent(input$keyboard_undo, {
      if (isTRUE(input$keyboard_undo > 0L)) {
        shiny::freezeReactiveValue(input, "keyboard_undo")
        history <- bounds_history()
        if (length(history) > 0L) {
          previous <- history[[length(history)]]
          bounds_history(history[-length(history)])
          set_bounds(previous, record = FALSE)
        }
      }
    })

    session$onFlushed(function() {
      shiny::insertUI(
        selector = "body",
        where = "beforeEnd",
        ui = shiny::tags$script(shiny::HTML("
          document.addEventListener('keydown', function(event) {
            var isUndo = (event.metaKey || event.ctrlKey) &&
              event.key && event.key.toLowerCase() === 'z';
            if (!isUndo || event.shiftKey) return;
            event.preventDefault();
            var shiny = window.Shiny || globalThis.Shiny;
            if (shiny && shiny.setInputValue) {
              shiny.setInputValue('keyboard_undo', Date.now(),
                {priority: 'event'});
            }
          });
          var gateRelayoutTimer = setInterval(function() {
            var plot = document.getElementById('gate_plot');
            if (!plot || typeof plot.on !== 'function' ||
                plot.dataset.expressalyzrRelayoutBound === 'true') {
              return;
            }
            plot.dataset.expressalyzrRelayoutBound = 'true';
            plot.on('plotly_relayout', function(eventData) {
              var shiny = window.Shiny || globalThis.Shiny;
              if (shiny && shiny.setInputValue) {
                shiny.setInputValue('gate_plot_relayout_js', eventData,
                  {priority: 'event'});
              }
            });
            clearInterval(gateRelayoutTimer);
          }, 500);
          function sendExpressalyzrPolygonDrawing() {
            var shiny = window.Shiny || globalThis.Shiny;
            if (!shiny || !shiny.setInputValue) return;
            var plot = document.getElementById('gate_plot');
            if (!plot) return;
            var drag = plot.querySelector('.nsewdrag');
            var paths = Array.prototype.slice.call(
              plot.querySelectorAll('.shapelayer path')
            );
            var path = null;
            for (var i = paths.length - 1; i >= 0; i--) {
              var d = paths[i].getAttribute('d');
              if (d && d.charAt(0) === 'M') {
                path = d;
                break;
              }
            }
            if (!path || !drag) {
              shiny.setInputValue('gate_plot_polygon_pixels',
                {error: 'No drawn polygon found.', nonce: Date.now()},
                {priority: 'event'});
              return;
            }
            shiny.setInputValue('gate_plot_polygon_pixels', {
              path: path,
              width: parseFloat(drag.getAttribute('width')),
              height: parseFloat(drag.getAttribute('height')),
              nonce: Date.now()
            }, {priority: 'event'});
          }
          var gatePolygonButtonTimer = setInterval(function() {
            var button = document.getElementById('use_polygon_drawing');
            if (!button || button.dataset.expressalyzrPolygonBound === 'true') {
              return;
            }
            button.dataset.expressalyzrPolygonBound = 'true';
            button.addEventListener('click', sendExpressalyzrPolygonDrawing);
            clearInterval(gatePolygonButtonTimer);
          }, 500);
        "))
      )
    }, once = TRUE)

    output$bounds_text <- shiny::renderText({
      bounds <- current_bounds()
      if (is.null(bounds)) {
        return("Draw a rectangle to define the boundary gate.")
      }
      paste0(
        "min = [", paste(format_gate_number(bounds$min), collapse = ", "),
        "]\nmax = [", paste(format_gate_number(bounds$max), collapse = ", "),
        "]"
      )
    })

    output$template_table <- DT::renderDT({
      DT::datatable(template_state(), rownames = FALSE,
                    options = list(pageLength = 6, scrollX = TRUE))
    })

    shiny::observeEvent(input$save_gate, {
      row_index <- selected_row_index()
      template_now <- template_state()
      dims_now <- c(input$x_channel, input$y_channel)
      gating_args <- input$gating_args

      if (!nzchar(input$gate_alias) || !nzchar(input$parent) ||
          !nzchar(input$gating_method)) {
        shiny::showNotification("Alias, parent, and gating method are required.",
                                type = "warning")
        return(invisible(NULL))
      }

      bounds <- current_bounds()
      if (identical(input$gating_method, "boundary")) {
        if (!gate_designer_bounds_are_valid(bounds)) {
          shiny::showNotification("Draw a rectangle or enter valid bounds before saving.",
                                  type = "warning")
          return(invisible(NULL))
        }
        gating_args <- build_boundary_gating_args(bounds$min, bounds$max)
      } else if (identical(input$gating_method, "polygon_gate")) {
        vertices <- current_polygon()
        if (!gate_designer_polygon_is_valid(vertices)) {
          shiny::showNotification("Draw a closed freeform gate before saving.",
                                  type = "warning")
          return(invisible(NULL))
        }
        gating_args <- build_polygon_gating_args(vertices[, 1L], vertices[, 2L])
      } else if (is.na(gating_args)) {
        gating_args <- ""
      }

      updated_row <- build_gate_template_row(
        alias = input$gate_alias,
        pop = input$pop,
        parent = input$parent,
        dims = dims_now,
        gating_method = input$gating_method,
        gating_args = gating_args,
        collapseDataForGating = input$collapseDataForGating,
        groupBy = input$groupBy,
        preprocessing_method = input$preprocessing_method,
        preprocessing_args = input$preprocessing_args
      )

      updated <- update_gate_template_row(template_now, row_index,
                                          updated_row)
      duplicate_alias <- duplicated(updated$alias) & nzchar(updated$alias)
      if (any(duplicate_alias)) {
        shiny::showNotification("Gate aliases must be unique before saving.",
                                type = "warning")
        return(invisible(NULL))
      }

      active_output_file <- current_template_spec()$output_file
      data.table::fwrite(updated, active_output_file)
      template_state(updated)
      if (identical(input$gating_method, "boundary")) {
        saved_bounds(bounds)
        saved_polygon(NULL)
        shiny::updateTextAreaInput(session, "gating_args", value = gating_args)
      } else if (identical(input$gating_method, "polygon_gate")) {
        saved_bounds(current_bounds())
        saved_polygon(current_polygon())
      } else {
        saved_bounds(current_bounds())
        saved_polygon(current_polygon())
      }
      refresh_gate_row_choices(row_index)
      mark_gate_saved()
      shiny::showNotification(paste("Saved", basename(active_output_file)),
                              type = "message")
    })

    shiny::observeEvent(input$finish_design, {
      if (isTRUE(has_unsaved_changes())) {
        shiny::showNotification(
          "Save the current gate first, then continue the pipeline.",
          type = "warning"
        )
        return(invisible(NULL))
      }
      shiny::stopApp(TRUE)
    })
  }

  app <- shiny::shinyApp(ui, server)
  if (isTRUE(launch)) {
    return(shiny::runApp(app))
  }
  app
}
