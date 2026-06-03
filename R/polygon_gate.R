polygon_gate <- function(fr, pp_res, channels = NA, filterId = "",
                         x = NULL, y = NULL, vertices = NULL) {
  if (!is.null(vertices)) {
    gate <- as.matrix(vertices)
    if (ncol(gate) != 2L) {
      stop("vertices must have two columns.", call. = FALSE)
    }
  } else if (!is.null(x) && !is.null(y)) {
    x <- as.numeric(x)
    y <- as.numeric(y)
    if (length(x) != length(y)) {
      stop("x and y must have the same length.", call. = FALSE)
    }
    gate <- cbind(x, y)
  } else {
    stop("polygon_gate requires x/y vectors or a two-column vertices matrix.",
         call. = FALSE)
  }

  storage.mode(gate) <- "numeric"
  if (nrow(gate) < 3L || any(!is.finite(gate))) {
    stop("polygon_gate requires at least three finite vertices.",
         call. = FALSE)
  }
  colnames(gate) <- channels

  flowCore::polygonGate(.gate = gate, filterId = filterId)
}
