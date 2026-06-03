#' Load .fcs files.
#'
#' Load .fcs files into GatingSet object using read.ncdfFlowSet function.
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

  data_dir <- file.path(file_path, "data")
  fcs_files <- list.files(path = data_dir, pattern = "^.*\\.fcs$", full.names = TRUE)

  for_order <- gsub("^(.*-.)(\\d{1}\\.fcs)$", "\\10\\2", fcs_files)
  fcs_files <- fcs_files[order(for_order)]

  if (identical(fcs_files, character(0))) {

    stop("The specified file path does either not exist or does not contain any .fcs files.")

  } else {
    return(flowWorkspace::load_cytoset_from_fcs(files = fcs_files))
  }
}

# this function causes a segfault form C stack overflow on my computer
load_fcs_ncdf <- function(file_path) {
  fcs_files <- list.files(path = file_path, pattern = "^.*\\.fcs$", full.names = TRUE)
  if (identical(fcs_files, character(0))) {
    stop("The specified file path does either not exist or does not contain any .fcs files.")
  } else {
    gating_set <- ncdfFlow::read.ncdfFlowSet(fcs_files)
  }
  return(gating_set)
}

#' Create data subdirectory.
#'
#' This is a utility function for creating a subdirectory for the data in the
#' data path and moving the data files into the new directory. This way the
#' analysis results and data stay separated.
#'
#' @param data_path The path to the original data directory.
#'
#' @return
#'
create_data_subdir <- function(data_path) {

  new_data_path <- file.path(data_path, "data")

  if (!dir.exists(new_data_path)) {
    # create new folder
    data_files <- list.files(data_path)
    data_files <- data_files[!grepl("\\.R$|\\.yml$|\\.csv$", data_files)]
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

#' Create and/or load configuration file.
#'
#' A configuration file is created from a template and opened for editing.
#' After editing is complete or if a configuration file already exists
#' the file is loaded.
#'
#' @param data_path The path to the original data directory.
#'
#' @export
load_config <- function(config_file_path, view_config) {

  if (!file.exists(config_file_path)) {

    template_path <- system.file("tools", "config_template.yml",
                                 package = "expressalyzr",
                                 mustWork = TRUE)
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
    file.edit(config_file_path)

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
comp_cv <- function(x, log_t = FALSE) {

  if (log_t) {
    x <- log(x[x > 0])
  }

  return(stats::sd(x) / mean(x))
}

#' Binning function.
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
      ggplot2::scale_x_continuous(trans = "log") +
      ggplot2::facet_grid(File ~ Channel, scales = "free_x")

    print(pl)

    new_th <- readline("Adjust background cutoff: ")
  }

  return(th)
}
