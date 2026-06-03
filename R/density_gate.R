#'
density_gate <- function(fr, pp_res, channels = NA, filterId = "",
                         res = 100, th = 1e-14, conv = FALSE, ellipse = FALSE) {
  max_r <- range(fr)[2, ]
  dat <- as.data.frame(flowCore::exprs(fr))
  dat <- dat[dat$`FSC-A` <= max_r[[1]] & dat$`SSC-A` <= max_r[[2]], ]
  kde <- MASS::kde2d(x = dat$`FSC-A`, y = dat$`SSC-A`, n = res)
  cl <- contourLines(kde, levels = th)

  a <- lapply(cl, function(li) abs(pracma::polyarea(li$x, li$y)))
  max_area <- which.max(unlist(a))

  gate <- cbind(cl[[max_area]]$x, cl[[max_area]]$y)

  if (conv) {
    gate <- gate[chull(gate),]
  } else if (ellipse) {
    og <- gate

    fit <- conicfit::EllipseDirectFit(og / 1000)
    fit_g <- conicfit::AtoG(fit)$ParG
    gate <- conicfit::calculateEllipse(fit_g[1], fit_g[2], fit_g[3], fit_g[4], 180 / pi * fit_g[5])
    gate <- gate * 1000
#
#     plot(og)
#     lines(gate, col = "red")
#
#     readline()
  }

  colnames(gate) <- channels
  return(flowCore::polygonGate(.gate = gate, filterId = filterId))
}
