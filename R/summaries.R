#' Compute a population summary table.
#'
#' Builds one row per sample and gated population from a \code{GatingSet}.
#'
#' @param gs A \code{GatingSet} object after gating.
#' @param experiment Experiment name to include in the output table.
#'
#' @return A \code{data.table} with population counts and frequencies.
#' @export
compute_population_table <- function(gs, experiment = NA_character_) {
  pops <- flowWorkspace::gs_get_pop_paths(gs, path = "auto")
  samples <- flowWorkspace::sampleNames(gs)

  rows <- vector("list", length(samples) * length(pops))
  idx <- 1L

  for (samp in samples) {
    total_count <- flowWorkspace::gh_pop_get_count(gs[[samp]], y = "root")

    for (pop in pops) {
      n_events <- flowWorkspace::gh_pop_get_count(gs[[samp]], y = pop)

      if (identical(pop, "root")) {
        parent <- NA_character_
        parent_count <- NA_real_
        freq_parent <- NA_real_
      } else {
        parent <- flowWorkspace::gs_pop_get_parent(gs, pop, path = "auto")
        parent_count <- flowWorkspace::gh_pop_get_count(gs[[samp]], y = parent)
        freq_parent <- if (parent_count > 0) n_events / parent_count * 100 else NA_real_
      }

      freq_total <- if (total_count > 0) n_events / total_count * 100 else NA_real_

      rows[[idx]] <- list(
        experiment = experiment,
        sample_id = samp,
        population = pop,
        parent_population = parent,
        n_events = n_events,
        n_parent_events = parent_count,
        frequency_parent = round(freq_parent, 4),
        frequency_total = round(freq_total, 4),
        gate_method = if (identical(pop, "root")) NA_character_ else "opencyto",
        gate_confidence = NA_real_,
        review_required = FALSE
      )
      idx <- idx + 1L
    }
  }

  data.table::rbindlist(rows)
}

#' Compute a marker/channel summary table.
#'
#' Summarises fluorescence channel values in an event-level table, typically
#' for the extracted population exported by \code{\link{run_pipeline}}.
#'
#' @param dt Event-level \code{data.table} with a \code{File} column.
#' @param channels Character vector of fluorescence channels to summarise.
#' @param experiment Experiment name to include in the output table.
#' @param population Population represented by \code{dt}.
#' @param scale Label describing the value scale in \code{dt}.
#'
#' @return A \code{data.table} with marker/channel summaries.
#' @export
compute_marker_table <- function(dt, channels, experiment = NA_character_,
                                 population = "singlets",
                                 scale = "processed") {
  available <- intersect(channels, colnames(dt))
  if (length(available) == 0L || nrow(dt) == 0L) {
    return(data.table::data.table(
      experiment = character(), sample_id = character(),
      population = character(), channel = character(), marker = character(),
      median = numeric(), mean = numeric(), sd = numeric(), iqr = numeric(),
      mad = numeric(), min = numeric(), max = numeric(),
      percent_positive = numeric(), scale = character()
    ))
  }

  files <- unique(dt$File)
  rows <- vector("list", length(files) * length(available))
  idx <- 1L

  for (file in files) {
    sub <- dt[File == file]

    for (channel in available) {
      vals <- sub[[channel]]
      positive_col <- paste0(channel, "_pos")
      percent_positive <- if (positive_col %in% colnames(sub)) {
        round(mean(sub[[positive_col]], na.rm = TRUE) * 100, 4)
      } else {
        NA_real_
      }

      rows[[idx]] <- list(
        experiment = experiment,
        sample_id = file,
        population = population,
        channel = channel,
        marker = channel,
        median = round(stats::median(vals, na.rm = TRUE), 4),
        mean = round(mean(vals, na.rm = TRUE), 4),
        sd = round(stats::sd(vals, na.rm = TRUE), 4),
        iqr = round(stats::IQR(vals, na.rm = TRUE), 4),
        mad = round(stats::mad(vals, na.rm = TRUE), 4),
        min = round(min(vals, na.rm = TRUE), 4),
        max = round(max(vals, na.rm = TRUE), 4),
        percent_positive = percent_positive,
        scale = scale
      )
      idx <- idx + 1L
    }
  }

  data.table::rbindlist(rows)
}
