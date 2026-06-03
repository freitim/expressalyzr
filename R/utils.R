#' Load .fcs files
#'
#' Load .fcs files into a cytoset.
#'
#' @param file_path The path to the .fcs files that are to be loaded.
#'
#' @return A GatingSet object.
#' @export
#'
#' @examples
#' path <- system.file("extdata", "example_fcs_files", package = "expressalyzr", mustWork = TRUE)
#' gating_set <- load_fcs(path)
load_fcs <- function(file_path) {
  if (!dir.exists(file_path)) {
    stop("The specified file path does either not exist or does not contain any .fcs files.")
  }

  data_dir <- file.path(file_path, "data")
  fcs_files <- if (dir.exists(data_dir)) {
    list.files(path = data_dir, pattern = "^.*\\.fcs$", full.names = TRUE)
  } else {
    character(0)
  }
  if (identical(fcs_files, character(0))) {
    fcs_files <- list.files(path = file_path, pattern = "^.*\\.fcs$",
                            full.names = TRUE)
  }

  for_order <- gsub("^(.*-.)(\\d{1}\\.fcs)$", "\\10\\2", fcs_files)
  fcs_files <- fcs_files[order(for_order)]

  if (identical(fcs_files, character(0))) {
    stop("The specified file path does either not exist or does not contain any .fcs files.")
  }

  flowWorkspace::load_cytoset_from_fcs(files = fcs_files)
}

cytoflex_channel_map <- function() {
  c(
    "445 [B]-A" = "FL6-A",
    "445 [B]-H" = "FL6-H",
    "445 [B]-W" = "FL6-Width",
    "488 [C]-A" = "FL1-A",
    "488 [C]-H" = "FL1-H",
    "488 [C]-W" = "FL1-Width",
    "640 [C]-A" = "FL3-A",
    "640 [C]-H" = "FL3-H",
    "640 [C]-W" = "FL3-Width"
  )
}

normalize_cytoflex_channels <- function(cs) {
  normalize_channels(cs, "cytoflex")
}

config_channel_map <- function(config) {
  if ("channel_map" %in% names(config)) {
    return(config$channel_map)
  }
  "cytoflex"
}

normalize_channels <- function(cs, channel_map = "cytoflex") {
  if (is.null(channel_map) ||
      identical(channel_map, FALSE) ||
      (length(channel_map) == 1L && is.na(channel_map)) ||
      (is.character(channel_map) &&
       length(channel_map) == 1L &&
       tolower(channel_map) %in% c("none", "null"))) {
    return(cs)
  }

  if (is.character(channel_map) &&
      length(channel_map) == 1L &&
      identical(tolower(channel_map), "cytoflex")) {
    channel_map <- cytoflex_channel_map()
  } else {
    channel_map <- unlist(channel_map, use.names = TRUE)
  }

  if (is.null(names(channel_map)) || any(!nzchar(names(channel_map)))) {
    stop("channel_map must be 'cytoflex', NULL, or a named map from FCS channel names to output channel names.",
         call. = FALSE)
  }

  chs <- flowWorkspace::colnames(cs)
  map_index <- match(chs, names(channel_map))
  mapped <- !is.na(map_index)
  chs[mapped] <- unname(channel_map[map_index[mapped]])
  flowWorkspace::colnames(cs) <- chs
  cs
}

fluorescence_channel_names <- function(channels, channel_pattern = "-A") {
  channels <- as.character(channels)

  if (!is.null(channel_pattern) &&
      !identical(channel_pattern, FALSE) &&
      length(channel_pattern) == 1L &&
      !is.na(channel_pattern) &&
      nzchar(channel_pattern)) {
    channels <- channels[grepl(channel_pattern, channels)]
  }

  excluded <- grepl("^(File|Time)$", channels) |
    grepl("^(FSC|SSC)(-|$)", channels) |
    grepl("(-W$|-Width$|Width$)", channels)

  channels[!excluded]
}

#' Create data subdirectory.
#'
#' This is a utility function for creating a subdirectory for the data in the
#' data path and moving the data files into the new directory. This way the
#' analysis results and data stay separated.
#'
#' @param data_path The path to the original data directory.
#'
#' @return \code{TRUE} if a data directory was created and files were moved,
#'   \code{FALSE} if the move failed, or \code{NULL} if the data directory
#'   already existed.
#'
create_data_subdir <- function(data_path) {

  new_data_path <- file.path(data_path, "data")

  if (!dir.exists(new_data_path)) {
    # create new folder
    data_files <- list.files(data_path, pattern = "\\.fcs$",
                             ignore.case = TRUE)
    create_out <- dir.create(new_data_path)

    # move files
    data_files_from <- file.path(data_path, data_files)
    data_files_to <- file.path(new_data_path, data_files)

    rename_out <- file.rename(data_files_from, data_files_to)

    output <- create_out & all(rename_out)
  } else {
    output <- NULL
  }

  return(output)
}

config_template_path <- function(config_profile = c("default", "legacy", "advanced")) {
  config_profile <- match.arg(config_profile)
  template_file <- paste0("config_", config_profile, ".yml")

  system.file("tools", template_file,
              package = "expressalyzr",
              mustWork = TRUE)
}

#' Create and/or load configuration file.
#'
#' A configuration file is created from a template and opened for editing.
#' After editing is complete or if a configuration file already exists
#' the file is loaded.
#'
#' @param config_file_path Path to the experiment \code{config.yml}.
#' @param view_config Whether to open the config file for review before
#'   loading it.
#' @param interactive Whether opening editors and waiting for console input is
#'   allowed.
#' @param config_profile Starter config profile to copy when
#'   \code{config.yml} does not exist. Existing config files are loaded
#'   unchanged. One of \code{"default"}, \code{"legacy"}, or
#'   \code{"advanced"}.
#'
#' @export
load_config <- function(config_file_path, view_config, interactive = TRUE,
                        config_profile = c("default", "legacy", "advanced")) {
  config_profile <- match.arg(config_profile)

  if (!file.exists(config_file_path)) {
    if (!interactive) {
      stop("config.yml does not exist and cannot be created/reviewed in non-interactive mode. Create config.yml first or run with interactive = TRUE.",
           call. = FALSE)
    }

    template_path <- config_template_path(config_profile)
    file.copy(template_path, config_file_path)

    template_gating <- system.file("tools", "gt_samples.csv",
                                   package = "expressalyzr",
                                   mustWork = TRUE)

    gating_file_path <- file.path(dirname(config_file_path), "gt_samples.csv")
    file.copy(template_gating, gating_file_path)

    file_created <- TRUE
  } else {
    file_created <- FALSE
  }

  if (view_config || file_created) {
    if (!interactive) {
      stop("view_config = TRUE requires interactive = TRUE.",
           call. = FALSE)
    }

    utils::file.edit(config_file_path)

    cat("\n")
    readline(prompt = "Press [Enter] to continue.")
  }

  return(config::get(file = config_file_path))
}

#'
#'
write_config <- function(config_file_path, active_config, value) {

  config_file <- readLines(config_file_path)
  config_ind <- which(!grepl("^  .*$", config_file))
  start_config <- which(grepl(active_config, config_file[config_ind]))

  end_config <- first(config_ind > start_config)
  if (!end_config) {
    end_config <- length(config_file)
  }

  active_section <- config_file[start_config:end_config]

  value_name <- gsub("^(.*):.*$", "\\1", value)
  value_ind <- which(grepl(value_name, active_section))

  active_section[value_ind] <- paste0("  ", value)
  config_file[start_config:end_config] <- active_section

  writeLines(config_file, config_file_path)

  return(TRUE)
}

#' Compute the coefficient of variation
#'
#' @param x Numeric vector.
#' @param log_t Whether to compute the coefficient of variation after log
#'   transformation of positive values.
#'
comp_cv <- function(x, log_t = FALSE) {

  if (log_t) {
    x <- log(x[x > 0])
  }

  return(stats::sd(x) / mean(x))
}

#' Binning function.
#'
#' @param x Numeric vector to bin.
#' @param bins Number of bins.
#' @param lower Lower quantile or range bound.
#' @param upper Upper quantile or range bound.
#' @param s_fun Optional summary function to apply within each bin.
#' @param use_quantiles Whether to use quantile breaks instead of evenly spaced
#'   breaks.
#'
bin <- function(x, bins, lower = 0, upper = 1, s_fun = NULL, use_quantiles = FALSE) {

  qs <- seq(lower, upper, length.out = bins + 1)
  breaks <- quantile(x, qs)

  if (!use_quantiles) {
    breaks <- seq(min(breaks), max(breaks), length.out = bins + 1)
  }

  bins <- cut(x, breaks, labels = FALSE, include.lowest = TRUE)

  if (!is.null(s_fun)) {
    t_dt <- data.table::data.table(bin = bins, x)[, s_fun(x), by = .(bin)]
    bins <- merge(data.table::data.table(bin = bins), t_dt,
                  by = "bin", all.x = TRUE, sort = FALSE)$V1
  }

  return(bins)
}

#'
#'
cf_to_dt <- function(cf) {
  values <- flowCore::exprs(cf)
  return(data.table::as.data.table(values))
}

#'
#'
cs_to_dt <- function(cs) {
  data_dt <- flowWorkspace::lapply(cs, cf_to_dt)
  return(data.table::rbindlist(data_dt, idcol = "File"))
}

shell_single_quote <- function(x) {
  paste0("'", gsub("'", "'\"'\"'", x, fixed = TRUE), "'")
}

read_terminal_prompt <- function(prompt) {
  if (.Platform$OS.type == "windows" || !file.exists("/dev/tty")) {
    return(NULL)
  }

  script <- paste(
    paste("printf '%s'", shell_single_quote(prompt), "> /dev/tty"),
    "IFS= read -r value < /dev/tty",
    "printf '%s\\n' \"$value\"",
    sep = "; "
  )
  value <- tryCatch(
    system2("sh", c("-c", script), stdout = TRUE, stderr = FALSE),
    error = function(e) NULL,
    warning = function(w) NULL
  )
  if (is.null(value)) {
    return(NULL)
  }
  if (length(value) == 0L) {
    return("")
  }
  value[[1L]]
}

compact_log10_breaks <- function(x, n = 3L) {
  x <- x[is.finite(x) & x > 0]
  if (length(x) == 0L) {
    return(numeric())
  }
  10^pretty(log10(range(x)), n = n)
}

compact_scientific_labels <- function(x) {
  format(x, scientific = TRUE, digits = 2, trim = TRUE)
}

read_prompt_input <- function(prompt, input = NULL) {
  if (interactive()) {
    return(readline(prompt = prompt))
  }

  if (!is.null(input)) {
    cat(prompt)
    flush.console()
    value <- tryCatch(
      readLines(input, n = 1L, warn = FALSE),
      error = function(e) character()
    )
    if (length(value) == 0L) {
      return("")
    }
    return(value[[1L]])
  }

  terminal_value <- read_terminal_prompt(prompt)
  if (!is.null(terminal_value)) {
    return(terminal_value)
  }

  cat(prompt)
  flush.console()
  value <- tryCatch(
    readLines(stdin(), n = 1L, warn = FALSE),
    error = function(e) character()
  )
  if (length(value) == 0L) {
    return("")
  }
  value[[1L]]
}

show_manual_plot <- function(plot, prefix = "expressalyzr_manual_plot",
                             output_dir = tempdir()) {
  if (interactive()) {
    print(plot)
    return(invisible(NULL))
  }

  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  path <- tempfile(pattern = paste0(prefix, "_"),
                   tmpdir = output_dir,
                   fileext = ".png")
  ggplot2::ggsave(path, plot = plot, width = 8, height = 6, dpi = 120)
  message("Manual adjustment plot written to ", path)

  opened <- try(utils::browseURL(normalizePath(path, mustWork = FALSE)),
                silent = TRUE)
  if (inherits(opened, "try-error")) {
    message("Open the PNG above to inspect the current adjustment plot.")
  }
  invisible(path)
}

#'
#'
adjust_threshold <- function(cont_dt, th) {

  tall_dt <- data.table::melt(cont_dt, id.vars = "File",
                              variable.name = "Channel",
                              value.name = "Intensity")

  new_th <- th

  while (!identical(new_th, "")) {
    th <- as.numeric(as.character(new_th))

    bg_dt <- tall_dt[File == unique(File)[1],
                     .(Intensity = quantile(Intensity, th)),
                     by = .(Channel)]

    pl <- ggplot2::ggplot(data = tall_dt[Intensity > 0],
                          ggplot2::aes(x = Intensity, fill = Channel)) +
      ggplot2::geom_histogram(bins = 64) +
      ggplot2::geom_vline(data = bg_dt,
                          ggplot2::aes(xintercept = Intensity),
                          linetype = "dashed") +
      ggplot2::scale_x_continuous(
        trans = "log10",
        breaks = compact_log10_breaks,
        labels = compact_scientific_labels,
        guide = ggplot2::guide_axis(check.overlap = TRUE)
      ) +
      ggplot2::facet_grid(File ~ Channel, scales = "free_x") +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))

    show_manual_plot(pl, prefix = "expressalyzr_background_cutoff")

    new_th <- read_prompt_input(
      "Adjust background cutoff (numeric + Enter redraws; blank Enter accepts): "
    )
  }

  return(th)
}
