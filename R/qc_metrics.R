#' Compute event counts at each gating stage.
#'
#' Extracts event counts for every population in a GatingSet, per sample,
#' and computes the percentage relative to the parent population.
#'
#' @param gs A \code{GatingSet} object with applied gates.
#'
#' @return A \code{data.table} with columns: \code{sample}, \code{population},
#'   \code{count}, \code{parent_pct}.
#' @export
compute_event_counts <- function(gs) {
  pops <- flowWorkspace::gs_get_pop_paths(gs, path = "auto")
  samples <- flowWorkspace::sampleNames(gs)

  rows <- vector("list", length(samples) * length(pops))
  idx <- 1L

  for (samp in samples) {
    for (pop in pops) {
      count <- flowWorkspace::gh_pop_get_count(gs[[samp]], y = pop)

      if (pop == "root") {
        parent_count <- count
      } else {
        parent <- flowWorkspace::gs_pop_get_parent(gs, pop, path = "auto")
        parent_count <- flowWorkspace::gh_pop_get_count(gs[[samp]], y = parent)
      }

      pct <- if (parent_count > 0) count / parent_count * 100 else NA_real_
      rows[[idx]] <- list(
        sample = samp,
        population = pop,
        count = count,
        parent_pct = round(pct, 2)
      )
      idx <- idx + 1L
    }
  }

  data.table::rbindlist(rows)
}

#' Compute scatter statistics per sample.
#'
#' Calculates mean, median, standard deviation, and coefficient of variation
#' for the specified scatter channels in a cytoset or GatingSet.
#'
#' @param cs A \code{cytoset} or \code{GatingSet} object.
#' @param channels Character vector of channel names (default: FSC-A, SSC-A).
#'
#' @return A \code{data.table} with columns: \code{sample}, \code{channel},
#'   \code{mean}, \code{median}, \code{sd}, \code{cv}.
#' @export
compute_scatter_stats <- function(cs, channels = c("FSC-A", "SSC-A")) {
  samples <- flowCore::sampleNames(cs)

  rows <- vector("list", length(samples) * length(channels))
  idx <- 1L

  for (samp in samples) {
    mat <- flowCore::exprs(cs[[samp]])
    available <- intersect(channels, colnames(mat))

    for (ch in available) {
      vals <- mat[, ch]
      m <- mean(vals, na.rm = TRUE)
      rows[[idx]] <- list(
        sample = samp,
        channel = ch,
        mean = round(m, 2),
        median = round(stats::median(vals, na.rm = TRUE), 2),
        sd = round(stats::sd(vals, na.rm = TRUE), 2),
        cv = round(if (m != 0) stats::sd(vals, na.rm = TRUE) / m * 100 else NA_real_, 2)
      )
      idx <- idx + 1L
    }
  }

  data.table::rbindlist(rows)
}

#' Compute fluorescence channel statistics per sample.
#'
#' Calculates summary statistics for fluorescence channels from a data.table
#' (typically produced by \code{cs_to_dt}).
#'
#' @param dt A \code{data.table} with a \code{File} column and fluorescence channels.
#' @param channels Character vector of fluorescence channel names to summarise.
#'
#' @return A \code{data.table} with columns: \code{sample}, \code{channel},
#'   \code{mean}, \code{median}, \code{sd}, \code{cv}, \code{pct_positive},
#'   \code{dynamic_range}.
#' @export
compute_channel_stats <- function(dt, channels) {
  available <- intersect(channels, colnames(dt))
  if (length(available) == 0L) {
    return(data.table::data.table(
      sample = character(), channel = character(),
      mean = numeric(), median = numeric(), sd = numeric(),
      cv = numeric(), pct_positive = numeric(), dynamic_range = numeric()
    ))
  }

  rows <- vector("list", data.table::uniqueN(dt$File) * length(available))
  idx <- 1L
  files <- unique(dt$File)

  for (f in files) {
    sub <- dt[File == f]
    for (ch in available) {
      vals <- sub[[ch]]
      m <- mean(vals, na.rm = TRUE)
      s <- stats::sd(vals, na.rm = TRUE)
      pct_pos <- sum(vals > 0, na.rm = TRUE) / length(vals) * 100
      rng <- range(vals, na.rm = TRUE)
      dyn_range <- if (rng[1] > 0) log10(rng[2] / rng[1]) else rng[2] - rng[1]

      rows[[idx]] <- list(
        sample = f,
        channel = ch,
        mean = round(m, 2),
        median = round(stats::median(vals, na.rm = TRUE), 2),
        sd = round(s, 2),
        cv = round(if (m != 0) s / m * 100 else NA_real_, 2),
        pct_positive = round(pct_pos, 2),
        dynamic_range = round(dyn_range, 2)
      )
      idx <- idx + 1L
    }
  }

  data.table::rbindlist(rows)
}

#' Check for fluidics anomalies using event rate over time.
#'
#' Divides acquisition time into bins and flags bins where the event rate
#' deviates more than 2 standard deviations from the sample mean rate.
#' Requires a \code{Time} parameter in the FCS data.
#'
#' @param cs A \code{cytoset} object.
#' @param n_bins Number of time bins (default: 20).
#'
#' @return A \code{data.table} with columns: \code{sample}, \code{time_bin},
#'   \code{event_count}, \code{rate_anomaly}. Returns \code{NULL} if no
#'   Time channel is present.
#' @export
compute_time_check <- function(cs, n_bins = 20L) {
  samples <- flowCore::sampleNames(cs)

  first_mat <- flowCore::exprs(cs[[samples[1]]])
  if (!"Time" %in% colnames(first_mat)) {
    message("No 'Time' channel found in data; skipping time check.")
    return(NULL)
  }

  rows <- vector("list", length(samples) * n_bins)
  idx <- 1L

  for (samp in samples) {
    mat <- flowCore::exprs(cs[[samp]])
    time_vals <- mat[, "Time"]
    breaks <- seq(min(time_vals), max(time_vals), length.out = n_bins + 1L)
    bin_ids <- cut(time_vals, breaks = breaks, labels = FALSE, include.lowest = TRUE)
    counts <- tabulate(bin_ids, nbins = n_bins)
    mean_rate <- mean(counts)
    sd_rate <- stats::sd(counts)

    for (b in seq_len(n_bins)) {
      rows[[idx]] <- list(
        sample = samp,
        time_bin = b,
        event_count = counts[b],
        rate_anomaly = abs(counts[b] - mean_rate) > 2 * sd_rate
      )
      idx <- idx + 1L
    }
  }

  data.table::rbindlist(rows)
}

#' Flag outlier values using median absolute deviation.
#'
#' Values more than \code{k} MADs from the median are flagged as outliers.
#'
#' @param x Numeric vector.
#' @param k Number of MADs to use as threshold (default: 2).
#'
#' @return Logical vector the same length as \code{x}; \code{TRUE} for outliers.
#' @export
flag_outliers <- function(x, k = 2) {
  med <- stats::median(x, na.rm = TRUE)
  mad_val <- stats::mad(x, na.rm = TRUE)
  if (mad_val == 0) mad_val <- stats::sd(x, na.rm = TRUE)
  if (is.na(mad_val) || mad_val == 0) return(rep(FALSE, length(x)))
  abs(x - med) > k * mad_val
}

#' Collect all QC metrics.
#'
#' Convenience wrapper that runs all QC computations and returns a named list.
#'
#' @param gs A \code{GatingSet} object (after gating).
#' @param cs A \code{cytoset} object (post-gating singlets).
#' @param dt A \code{data.table} from \code{cs_to_dt} (post-compensation).
#' @param channels Character vector of fluorescence channel names.
#'
#' @return A named list with elements: \code{event_counts}, \code{scatter_stats},
#'   \code{channel_stats}, \code{time_check}.
#' @export
collect_qc <- function(gs, cs, dt, channels) {
  list(
    event_counts  = compute_event_counts(gs),
    scatter_stats = compute_scatter_stats(cs),
    channel_stats = compute_channel_stats(dt, channels),
    time_check    = compute_time_check(cs)
  )
}
