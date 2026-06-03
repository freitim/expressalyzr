#' expressalyzr plot palette.
#'
#' @return A named character vector of hex colours.
#' @export
expressalyzr_palette <- function() {
  c(
    blue = "#2F6F9F",
    teal = "#2A9D8F",
    green = "#59A14F",
    amber = "#E9A23B",
    coral = "#D76A5A",
    plum = "#8E6C8A",
    ink = "#243447",
    muted = "#6C757D",
    grid = "#E8EDF2",
    panel = "#F7F9FB"
  )
}

#' expressalyzr ggplot theme.
#'
#' @param base_size Base font size.
#'
#' @return A ggplot2 theme.
#' @export
theme_expressalyzr <- function(base_size = 11) {
  pal <- expressalyzr_palette()
  ggplot2::theme_minimal(base_size = base_size) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", colour = pal[["ink"]],
                                         hjust = 0, size = base_size * 1.15),
      plot.subtitle = ggplot2::element_text(colour = pal[["muted"]],
                                            size = base_size * 0.9),
      axis.title = ggplot2::element_text(colour = pal[["ink"]]),
      axis.text = ggplot2::element_text(colour = pal[["ink"]]),
      panel.grid.major = ggplot2::element_line(colour = pal[["grid"]],
                                               linewidth = 0.35),
      panel.grid.minor = ggplot2::element_blank(),
      strip.text = ggplot2::element_text(face = "bold", colour = pal[["ink"]]),
      strip.background = ggplot2::element_rect(fill = pal[["panel"]],
                                               colour = NA),
      legend.position = "bottom"
    )
}

#' Plot population summaries.
#'
#' @param population_table A table returned by \code{\link{compute_population_table}}.
#' @param populations Optional character vector of populations to include.
#'   By default, all non-root populations are shown.
#' @param value Summary value to plot: \code{"frequency_parent"},
#'   \code{"frequency_total"}, or \code{"n_events"}.
#'
#' @return A \code{ggplot} object.
#' @export
plot_population_summary <- function(population_table, populations = NULL,
                                    value = c("frequency_parent",
                                              "frequency_total",
                                              "n_events")) {
  value <- match.arg(value)
  dt <- data.table::as.data.table(data.table::copy(population_table))

  if (is.null(populations)) {
    dt <- dt[population != "root"]
  } else {
    dt <- dt[population %in% populations]
  }

  if (nrow(dt) == 0L) {
    stop("No population rows available to plot.", call. = FALSE)
  }

  dt[, plot_value := get(value)]
  y_label <- switch(value,
                    frequency_parent = "% of Parent",
                    frequency_total = "% of Total",
                    n_events = "Events")
  pal <- expressalyzr_palette()

  ggplot2::ggplot(dt, ggplot2::aes(x = sample_id, y = plot_value,
                                   fill = population)) +
    ggplot2::geom_col(width = 0.75) +
    ggplot2::scale_fill_manual(values = rep(
      c(pal[["blue"]], pal[["teal"]], pal[["green"]], pal[["amber"]],
        pal[["coral"]], pal[["plum"]]),
      length.out = length(unique(dt$population))
    ), guide = "none") +
    ggplot2::facet_wrap(ggplot2::vars(population), scales = "free_y") +
    ggplot2::labs(x = NULL, y = y_label, title = "Population Summary") +
    theme_expressalyzr(base_size = 11) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)
    )
}

#' Plot marker/channel summaries.
#'
#' @param marker_table A table returned by \code{\link{compute_marker_table}}.
#' @param channels Optional character vector of channels to include.
#' @param value Summary value to plot: \code{"median"}, \code{"mean"}, or
#'   \code{"percent_positive"}.
#'
#' @return A \code{ggplot} object.
#' @export
plot_marker_summary <- function(marker_table, channels = NULL,
                                value = c("median", "mean",
                                          "percent_positive")) {
  value <- match.arg(value)
  dt <- data.table::as.data.table(data.table::copy(marker_table))

  if (!is.null(channels)) {
    dt <- dt[channel %in% channels]
  }

  if (nrow(dt) == 0L) {
    stop("No marker rows available to plot.", call. = FALSE)
  }

  dt[, plot_value := get(value)]
  dt <- dt[is.finite(plot_value)]

  if (nrow(dt) == 0L) {
    stop("No finite marker values available to plot.", call. = FALSE)
  }

  y_label <- switch(value,
                    median = "Median Intensity",
                    mean = "Mean Intensity",
                    percent_positive = "% Positive")
  pal <- expressalyzr_palette()

  ggplot2::ggplot(dt, ggplot2::aes(x = sample_id, y = plot_value,
                                   group = channel)) +
    ggplot2::geom_line(color = pal[["muted"]], linewidth = 0.4) +
    ggplot2::geom_point(color = pal[["teal"]], size = 2.2) +
    ggplot2::facet_wrap(ggplot2::vars(channel), scales = "free_y") +
    ggplot2::labs(x = NULL, y = y_label, title = "Marker Summary") +
    theme_expressalyzr(base_size = 11) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)
    )
}

#' Plot a QC overview.
#'
#' @param qc_results A QC list returned by \code{\link{collect_qc}}, or an
#'   event-count table returned by \code{\link{compute_event_counts}}.
#'
#' @return A \code{ggplot} object.
#' @export
plot_qc_overview <- function(qc_results) {
  event_counts <- if (is.list(qc_results) && !is.null(qc_results$event_counts)) {
    qc_results$event_counts
  } else {
    qc_results
  }

  dt <- data.table::as.data.table(data.table::copy(event_counts))
  dt <- dt[population != "root" & !is.na(parent_pct)]

  if (nrow(dt) == 0L) {
    stop("No event-count QC rows available to plot.", call. = FALSE)
  }
  pal <- expressalyzr_palette()

  ggplot2::ggplot(dt, ggplot2::aes(x = population, y = parent_pct)) +
    ggplot2::geom_boxplot(fill = pal[["panel"]], color = pal[["muted"]],
                          outlier.shape = NA) +
    ggplot2::geom_jitter(color = pal[["blue"]], width = 0.15,
                         height = 0, size = 1.8, alpha = 0.8) +
    ggplot2::labs(x = NULL, y = "% of Parent",
                  title = "QC Overview: Gate Retention") +
    theme_expressalyzr(base_size = 11) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)
    )
}

#' Plot the gating hierarchy.
#'
#' @param gs A \code{GatingSet} object after gating.
#'
#' @return A \code{ggplot} object.
#' @export
plot_gate_hierarchy <- function(gs) {
  pops <- flowWorkspace::gs_get_pop_paths(gs, path = "auto")
  if (length(pops) == 0L) {
    stop("No populations available to plot.", call. = FALSE)
  }

  rows <- lapply(pops, function(pop) {
    parent <- if (identical(pop, "root")) {
      NA_character_
    } else {
      flowWorkspace::gs_pop_get_parent(gs, pop, path = "auto")
    }

    data.table::data.table(
      population = pop,
      parent_population = parent,
      depth = if (identical(pop, "root")) 0L else NA_integer_
    )
  })

  nodes <- data.table::rbindlist(rows)
  if (nrow(nodes) > 1L) {
    for (i in seq_len(nrow(nodes))) {
      for (j in which(is.na(nodes$depth))) {
        parent <- nodes$parent_population[j]
        parent_depth <- nodes[population == parent, depth]
        if (length(parent_depth) == 1L && !is.na(parent_depth)) {
          nodes$depth[j] <- parent_depth + 1L
        }
      }
    }
  }

  nodes[is.na(depth), depth := 0L]
  nodes[, y := rev(seq_len(.N))]
  nodes[, x := depth]

  edges <- nodes[!is.na(parent_population),
                 .(population, x, y, parent_population)]
  edges <- merge(
    edges,
    nodes[, .(parent_population = population, parent_x = x, parent_y = y)],
    by = "parent_population",
    all.x = TRUE,
    sort = FALSE
  )

  pal <- expressalyzr_palette()
  ggplot2::ggplot(nodes, ggplot2::aes(x = x, y = y)) +
    ggplot2::geom_segment(
      data = edges,
      ggplot2::aes(x = parent_x, y = parent_y, xend = x, yend = y),
      inherit.aes = FALSE,
      linewidth = 0.6,
      colour = pal[["muted"]]
    ) +
    ggplot2::geom_point(
      ggplot2::aes(fill = population == "root"),
      shape = 21,
      size = 5,
      stroke = 0.8,
      colour = pal[["ink"]]
    ) +
    ggplot2::geom_text(
      ggplot2::aes(label = population),
      hjust = 0,
      nudge_x = 0.08,
      size = 3.6,
      colour = pal[["ink"]]
    ) +
    ggplot2::scale_fill_manual(values = c("FALSE" = pal[["teal"]],
                                          "TRUE" = pal[["blue"]]),
                               guide = "none") +
    ggplot2::scale_x_continuous(
      breaks = seq(0, max(nodes$depth, na.rm = TRUE)),
      expand = ggplot2::expansion(mult = c(0.05, 0.4))
    ) +
    ggplot2::scale_y_continuous(NULL, breaks = NULL) +
    ggplot2::labs(x = NULL, title = "Gating Hierarchy") +
    theme_expressalyzr(base_size = 11) +
    ggplot2::theme(
      panel.grid.major.y = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_blank(),
      axis.ticks.x = ggplot2::element_blank()
    )
}

gate_plot_parent_data <- function(gs, parent) {
  if (is.null(parent) || is.na(parent) || !nzchar(parent) ||
      identical(parent, "root")) {
    return(tryCatch(
      flowWorkspace::gs_pop_get_data(gs, y = "root"),
      error = function(e) flowWorkspace::gs_pop_get_data(gs)
    ))
  }

  flowWorkspace::gs_pop_get_data(gs, y = parent)
}

gate_plot_template_row <- function(template_file, pop) {
  template <- read_gating_template(template_file)
  target_pop <- pop
  row <- template[which(template[["alias"]] == target_pop)]
  if (nrow(row) != 1L) {
    stop("Could not find one template row for population '", pop, "'.",
         call. = FALSE)
  }
  row
}

template_gate_spec <- function(row) {
  dims <- split_gate_dims(row$dims)
  if (length(dims) != 2L) {
    stop("Template row '", row$alias,
         "' must have exactly two dimensions for fallback plotting.",
         call. = FALSE)
  }

  if (identical(row$gating_method, "boundary")) {
    bounds <- parse_boundary_gating_args(row$gating_args)
    if (is.null(bounds)) {
      stop("Could not parse boundary gate arguments for '", row$alias, "'.",
           call. = FALSE)
    }
    vertices <- data.table::data.table(
      .x_plot = c(bounds$min[[1L]], bounds$max[[1L]], bounds$max[[1L]],
                  bounds$min[[1L]], bounds$min[[1L]]),
      .y_plot = c(bounds$min[[2L]], bounds$min[[2L]], bounds$max[[2L]],
                  bounds$max[[2L]], bounds$min[[2L]])
    )
  } else if (identical(row$gating_method, "polygon_gate")) {
    polygon <- parse_polygon_gating_args(row$gating_args)
    if (is.null(polygon)) {
      stop("Could not parse polygon gate arguments for '", row$alias, "'.",
           call. = FALSE)
    }
    vertices <- data.table::data.table(
      .x_plot = polygon[, 1L],
      .y_plot = polygon[, 2L]
    )
    vertices <- data.table::rbindlist(
      list(vertices, vertices[1L]),
      use.names = TRUE
    )
  } else {
    vertices <- NULL
  }

  list(
    dims = dims,
    vertices = vertices,
    has_overlay = !is.null(vertices)
  )
}

plot_template_gate_page <- function(gs, template_file, pop, samples,
                                    parent = "root", title = NULL,
                                    max_events = 5000L, ncol = 4L) {
  row <- gate_plot_template_row(template_file, pop)
  gate <- template_gate_spec(row)
  parent_cs <- gate_plot_parent_data(gs, parent)
  available_samples <- flowWorkspace::sampleNames(parent_cs)
  samples <- intersect(samples, available_samples)
  if (length(samples) == 0L) {
    stop("No requested samples are available for fallback gate plotting.",
         call. = FALSE)
  }

  parent_cs <- parent_cs[samples]
  dt <- cs_to_dt(parent_cs)
  dims <- gate$dims
  missing_dims <- setdiff(dims, names(dt))
  if (length(missing_dims) > 0L) {
    stop("Gate plot channels are missing from the cytometry data: ",
         paste(missing_dims, collapse = ", "), call. = FALSE)
  }

  dt <- dt[, .(
    name = File,
    .x_plot = as.numeric(get(dims[[1L]])),
    .y_plot = as.numeric(get(dims[[2L]]))
  )]
  dt <- dt[is.finite(.x_plot) & is.finite(.y_plot)]
  if (nrow(dt) == 0L) {
    stop("No finite events are available for fallback gate plotting.",
         call. = FALSE)
  }

  if (is.finite(max_events) && max_events > 0L) {
    dt <- dt[, .SD[seq_len(min(.N, as.integer(max_events)))], by = name]
  }
  dt[, name := factor(name, levels = samples)]

  pal <- expressalyzr_palette()
  if (is.null(title) || is.na(title) || !nzchar(title)) {
    title <- paste(parent, pop, sep = " -> ")
  }

  p <- ggplot2::ggplot(dt, ggplot2::aes(x = .x_plot, y = .y_plot)) +
    ggplot2::geom_bin2d(bins = 35, alpha = 0.9) +
    ggplot2::facet_wrap(ggplot2::vars(name), ncol = ncol, scales = "free") +
    ggplot2::scale_x_continuous(labels = compact_scientific_labels,
                                breaks = scales::breaks_extended(n = 3)) +
    ggplot2::scale_y_continuous(labels = compact_scientific_labels,
                                breaks = scales::breaks_extended(n = 3)) +
    ggplot2::scale_fill_gradientn(
      colours = c(pal[["panel"]], pal[["teal"]], pal[["ink"]]),
      name = "count"
    ) +
    ggplot2::labs(x = dims[[1L]], y = dims[[2L]], title = title) +
    theme_expressalyzr(base_size = 10) +
    ggplot2::theme(
      aspect.ratio = 1,
      strip.text = ggplot2::element_text(size = 8),
      axis.text = ggplot2::element_text(size = 7),
      axis.title = ggplot2::element_text(size = 9),
      panel.spacing = grid::unit(1.1, "lines")
    )

  if (isTRUE(gate$has_overlay)) {
    p <- p +
      ggplot2::geom_path(
        data = gate$vertices,
        ggplot2::aes(x = .x_plot, y = .y_plot),
        inherit.aes = FALSE,
        colour = "#FF3B30",
        linewidth = 0.55
      )
  }

  p
}

save_run_plots <- function(data_path, experiment_name, population_table = NULL,
                           marker_table = NULL, qc_results = NULL,
                           gs = NULL) {
  paths <- list()

  if (!is.null(gs)) {
    path <- file.path(data_path, paste0(experiment_name, "_gate_hierarchy.png"))
    ggplot2::ggsave(path, plot_gate_hierarchy(gs),
                    width = 8, height = 5, dpi = 150)
    message("Gate hierarchy plot written to ", path)
    paths$gate_hierarchy <- path
  }

  if (!is.null(population_table) && nrow(population_table) > 0L) {
    path <- file.path(data_path, paste0(experiment_name, "_population_summary.png"))
    ggplot2::ggsave(path, plot_population_summary(population_table),
                    width = 10, height = 6, dpi = 150)
    message("Population plot written to ", path)
    paths$population <- path
  }

  if (!is.null(marker_table) && nrow(marker_table) > 0L) {
    path <- file.path(data_path, paste0(experiment_name, "_marker_summary.png"))
    ggplot2::ggsave(path, plot_marker_summary(marker_table),
                    width = 10, height = 6, dpi = 150)
    message("Marker plot written to ", path)
    paths$marker <- path
  }

  if (!is.null(qc_results)) {
    path <- file.path(data_path, paste0(experiment_name, "_qc_overview.png"))
    ggplot2::ggsave(path, plot_qc_overview(qc_results),
                    width = 8, height = 5, dpi = 150)
    message("QC overview plot written to ", path)
    paths$qc <- path
  }

  invisible(paths)
}
