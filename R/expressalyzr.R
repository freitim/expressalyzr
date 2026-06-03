#' expressalyzr: Automated gating and compensation of flow cytometry data
#'
#' The expressalyzr package provides a reproducible pipeline for gating,
#' compensation, quality control, and event-level export of flow cytometry data.
#'
#' @importFrom graphics title
#' @importFrom grDevices chull contourLines
#' @importFrom stats IQR median na.omit optim quantile runif sd
#' @importFrom utils file.edit flush.console head tail
#' @import flowStats
"_PACKAGE"

utils::globalVariables(c(
  ".", "..density..", "alias", "Auto", "Channel", "Cluster",
  "count_per_bin", "density", "depth", "File", "Fit", "gating_method",
  "Intensity", "Intercept", "MEFL", "mefl_dt", "no_negative", "parent_pct",
  "parent_population", "parent_x", "parent_y", "Peak", "plot_value",
  "population", "positive", "sample_id", "Slope", "channel", "x", "y",
  "control_channel", "control_channel_normalized", "desc",
  "gated_population", "gating_template", "inferred_channel", "name",
  "sample_match_index", "x_bin", "y_bin", ".x_plot", ".y_plot"
))
