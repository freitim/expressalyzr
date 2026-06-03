#' Generate an HTML analysis report
#'
#' Renders one self-contained HTML report with QC, gating, channel,
#' compensation, and summary sections.
#'
#' @param data_path Path to the experiment data directory.
#' @param gs A \code{GatingSet} object (after gating).
#' @param qc_results Named list returned by \code{\link{collect_qc}}.
#' @param so_mat Optional spillover matrix (or \code{NULL}).
#' @param gating_results Optional list of per-template gating results. Used to
#'   render multi-template gating sections in the unified report.
#' @param data_dt Optional event table used for channel distribution detail.
#' @param population_table Optional population summary table.
#' @param marker_table Optional marker/channel summary table.
#' @param save_gating_set Whether GatingSet outputs were written.
#' @param summaries_written Whether summary CSV outputs were written.
#' @param include_gating_plots If \code{TRUE}, include per-sample gate scatter
#'   grids in the unified HTML report. These plots are useful for gate review
#'   but can make reports large for plate-scale experiments.
#' @param html_theme Bootstrap theme for the HTML report.
#'   One of \code{"flatly"} (default), \code{"cosmo"}, \code{"journal"},
#'   \code{"lumen"}, \code{"sandstone"}, \code{"united"}, or any valid
#'   \code{bslib} theme name.
#' @param gg_theme Base ggplot2 theme: \code{"minimal"} (default),
#'   \code{"classic"}, \code{"bw"}, \code{"light"}, or \code{"gray"}.
#' @param max_events Maximum events per sample to use in scatter and
#'   distribution plots. Lower values render faster; set to \code{Inf}
#'   to disable downsampling.
#'
#' @return Path to the generated HTML file (invisibly).
#' @export
generate_qc_report <- function(data_path, gs, qc_results, so_mat = NULL,
                                gating_results = NULL, data_dt = NULL,
                                population_table = NULL, marker_table = NULL,
                                save_gating_set = FALSE,
                                summaries_written = FALSE,
                                include_gating_plots = FALSE,
                                html_theme = "flatly", gg_theme = "minimal",
                                max_events = 5000L) {
  check_report_dependencies()

  experiment_name <- basename(data_path)

  template <- system.file("tools", "qc_report.Rmd",
                           package = "expressalyzr", mustWork = TRUE)

  output_file <- file.path(data_path, paste0(experiment_name, "_report.html"))

  rmarkdown::render(
    input = template,
    output_file = output_file,
    output_options = list(theme = html_theme),
    params = list(
      experiment_name = experiment_name,
      data_path       = data_path,
      gs              = gs,
      qc_results      = qc_results,
      so_mat          = so_mat,
      gating_results  = gating_results,
      data_dt         = data_dt,
      population_table = population_table,
      marker_table    = marker_table,
      save_gating_set = save_gating_set,
      summaries_written = summaries_written,
      include_gating_plots = include_gating_plots,
      gg_theme        = gg_theme,
      max_events      = max_events
    ),
    envir = new.env(parent = globalenv()),
    quiet = TRUE
  )

  message("Report written to ", output_file)
  invisible(output_file)
}

check_report_dependencies <- function() {
  required <- c("ggcyto", "hexbin", "knitr", "scales")
  missing <- required[!vapply(required, requireNamespace, logical(1),
                              quietly = TRUE)]

  if (length(missing) > 0) {
    stop(
      "Report generation requires missing package(s): ",
      paste(missing, collapse = ", "),
      call. = FALSE
    )
  }

  invisible(TRUE)
}
