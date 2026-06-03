#' Compute the spillover matrix based on single color transfection controls.
#'
spillover_matrix <- function(data, cont_ind, comp_pattern, threshold, manual_comp) {

  if (length(data) != length(cont_ind)) {
    stop("Number of control samples and index of control samples are of different length.")
  }

  # data <- data[order(flowCore::sampleNames(data))]
  data <- data[cont_ind]

  comp_pattern <- paste0("FL\\d{1,2}", comp_pattern)

  chs_all <- flowCore::colnames(data)
  chs <- chs_all[grepl(comp_pattern, chs_all)]

  n_bins <- 32

  data <- filter_density(data, chs, n_bins, threshold)

  so_mat <- flowCore::spillover(data,
                                unstained = 1,
                                patt = comp_pattern,
                                fsc = "FSC-A",
                                ssc = "SSC-A",
                                method = "mean",
                                stain_match = "ordered")

  if (manual_comp) {
    so_mat <- manual_compensation(data[-cont_ind[1], chs], so_mat, chs)
  }

  return(so_mat)
}

#'
#'
manual_compensation <- function(data, so_mat, chs) {

  a_df <- flowCore::parameters(data[1, chs[1]][[1]])

  ch_r <- c(a_df[["minRange"]], a_df[["maxRange"]])

  for (i in chs) {
    for (j in chs[!chs %in% i]) {

      new_v <- so_mat[i, j]

      while (!identical(new_v, "")) {
        comp_data <- flowWorkspace::realize_view(data)
        comp_data <- flowWorkspace::compensate(comp_data, so_mat)

        local_cf <- comp_data[chs %in% i]

        pl <- ggplot2::ggplot(data = local_cf,
                              ggplot2::aes(x = get(j), y = get(i))) +
          ggplot2::geom_hex(ggplot2::aes(fill = ..density..),
                            bins = 64) +
          ggplot2::scale_x_continuous(trans = "log10", limits = c(NA, ch_r[2])) +
          ggplot2::scale_y_continuous(trans = "log10", limits = c(NA, ch_r[2])) +
          ggplot2::xlab(j) +
          ggplot2::ylab(i) +
          ggplot2::scale_fill_viridis_c() +
          ggplot2::ggtitle(paste("Current value:", so_mat[i, j]))

        print(pl)

        new_v <- readline(prompt = "Adjust spillover: ")

        if (suppressWarnings(!is.na(as.numeric(new_v)))) {
          so_mat[i, j] <- as.numeric(as.character(new_v))
        }
      }
    }
  }
  return(so_mat)
}

#' Filter data according to a density threshold.
#'
filter_density <- function(data, channels, bins, th) {

  l_data <- length(data)

  data_dt <- cs_to_dt(data)

  chs <- paste0(channels, "_bin")

  data_dt[, (channels) := lapply(.SD, log),
          .SDcols = channels,
          by = .(File)]

  data_dt <- na.omit(data_dt)

  for (ch in channels) {
    data_dt <- data_dt[!is.infinite(get(ch))]
  }

  data_dt[, (chs) := lapply(.SD, bin,
                            bins = bins, s_fun = mean),
          .SDcols = channels,
          by = .(File)]

  data_dt[, count_per_bin := .N, by = c(chs, "File")]
  data_dt[, density := count_per_bin / .N, by = .(File)]

  data_dt <- data_dt[density > th]

  data_dt[, (channels) := lapply(.SD, exp),
          .SDcols = channels,
          by = .(File)]

  col_remove <- c("File", chs, "count_per_bin", "density")

  to_data.frame <- function(file) {

    dat <- data_dt[File == file, !col_remove, with = FALSE]
    dat <- as.matrix(dat)

    return(flowCore::flowFrame(dat))
  }

  out_data <- sapply(unique(data_dt$File), to_data.frame,
                     simplify = FALSE, USE.NAMES = TRUE)

  out_data <- flowWorkspace::flowSet_to_cytoset(flowCore::flowSet(out_data))

  if (l_data != length(out_data)) {
    stop("Choose a lower threshold, it appears that one or more samples got filtered out.")
  }

  return(out_data)
}
