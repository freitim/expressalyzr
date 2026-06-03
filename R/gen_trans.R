#'
generate_transformation <- function(data) {

  gt_file <- system.file("tools", "gt_beads.csv", package = "expressalyzr",
                         mustWork = TRUE)
  gt <- openCyto::gatingTemplate(gt_file)

  gs <- flowWorkspace::GatingSet(data)

  openCyto::register_plugins(fun = mixture_gate, methodName = "mixture_gate",
                             dep = c("Rmixmod"), "gating")

  openCyto::gt_gating(gt, gs)

  data_dt <- flowWorkspace::gs_pop_get_data(gs, y = "singlets")
  data_dt <- flowCore::exprs(data_dt[[1]])
  data_dt <- data.table::as.data.table(data_dt)

  chs <- colnames(data_dt)
  chs <- chs[grepl("FL", chs)]

  data_dt <- data_dt[, chs, with = FALSE]
  data_dt <- data_dt[rowSums(data_dt <= 0) == 0]

  mods <- Rmixmod::mixmodGaussianModel(listModels = "Gaussian_pk_Lk_Ck")
  strat <- Rmixmod::mixmodStrategy(algo = "CEM",
                                   nbTry = 50,
                                   initMethod = "random")
  fit <- Rmixmod::mixmodCluster(log(data_dt),
                                nbCluster = 6:9,
                                models = mods,
                                strategy = strat)

  data_dt[, Cluster := fit["partition"]]
  data.table::setkey(data_dt, Cluster)

  met_dt <- data_dt[, lapply(.SD, function(x) IQR(x)/median(x)), by = .(Cluster)]

  select_i <- rowSums(met_dt[, !"Cluster"] > 1.275) == 0
  select_cl <- met_dt[select_i]$Cluster

  # pl <- ggplot2::ggplot(data = data_dt, ggplot2::aes(x = `FL1-A`, y = `FL3-A`,
  #                                                    color = Cluster %in% select_cl)) +
  #   ggplot2::geom_point(alpha = 0.1) +
  #   ggplot2::scale_x_continuous(trans = "log") +
  #   ggplot2::scale_y_continuous(trans = "log") +
  #   ggplot2::facet_wrap(Cluster~.)
  #
  # print(pl)

  data_dt <- data_dt[Cluster %in% select_cl]

  data_dt <- data_dt[, lapply(.SD, mean), by = .(Cluster)]

  setkeyv(data_dt, tail(colnames(data_dt), 1))

  peaks <- nrow(mefl_dt)
  data_dt[, Peak := (1 + peaks - length(Cluster)) : peaks]
  data_dt <- merge(data_dt, mefl_dt, by = "Peak")

  data_dt <- data_dt[!is.na(MEFL)]

  tall_dt <- melt(data_dt,
                  id.vars = c("Peak", "Cluster", "MEFL"),
                  variable.name = "Channel",
                  value.name = "Intensity")

  tall_dt[, c("Slope", "Intercept", "Auto") := as.list(run_fit(MEFL, Intensity, sos)),
          by = .(Channel)]

  tall_dt[, Fit := bead_model(MEFL, Slope, Intercept, Auto), by = .(Channel)]

  pl <- ggplot2::ggplot(tall_dt,
                        ggplot2::aes(x = Intensity,
                                     y = MEFL)) +
    ggplot2::geom_line(ggplot2::aes(x = Fit)) +
    ggplot2::geom_point(size = 3,
                        ggplot2::aes(color = as.factor(Peak))) +
    ggplot2::scale_x_continuous(trans = "log") +
    ggplot2::scale_y_continuous(trans = "log") +
    ggplot2::labs(color = "Peak") +
    ggplot2::facet_wrap(~Channel)

  print(pl)

  trans_par <- unique(tall_dt[, .(Channel, Slope, Intercept)])
  trans_fun <- function(value, channel) {

    m <- trans_par[Channel == channel]$Slope
    b <- trans_par[Channel == channel]$Intercept

    return(sign(value) * exp(b) * (abs(value) ^ m))
  }

  return(trans_fun)
}

#' Sum of squares cost function.
#'
sos <- function(f) sum(f ^ 2)

#' Fit bead model to mean peak data.
#'
run_fit <- function(x, y, f) {

  fit_fun <- function(p) f(log(x + p[3]) - (p[1] * log(y) + p[2]))

  # estimate start values
  m <- diff(log(tail(x, 2))) / diff(log(tail(y, 2)))
  b <- log(tail(x, 1)) - m * log(tail(y, 1))
  a <- exp(m * log(y[1]) + b) - x[1]

  fit <- optim(c(m, b, a),
               fit_fun,
               method = "L-BFGS-B",
               lower = c(-Inf, -Inf, 0))

  return(fit$par)
}

#'
#'
bead_model <- function(x, m, b, a) exp((log(x + a) - b) / m)

#' Transform arbitrary values.
#'
apply_transform <- function(data, t_fun) {

  chs <- flowCore::colnames(data)
  chs <- chs[grepl("FL", chs)]

  set_fun <- function(channel) {
    force(channel)
    out_f <- function(values) t_fun(values, channel)
    return(out_f)
  }

  t_list <- lapply(chs, set_fun)
  t_list <- flowCore::transformList(chs, t_list)

  return(flowWorkspace::transform(data, t_list))
}
