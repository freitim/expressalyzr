#' ---
#' title: "Gating report"
#' output: html_document
#' ---

#+ echo=F, warnings=F, message=T, fig.width = 15, fig.height = 15
cat(experiment_name)

pops <- flowWorkspace::gs_get_pop_paths(gs, path = "auto")

for (pop in pops[2:length(pops)]) {
  p <- ggcyto::autoplot(gs, pop, bins = 30) +
    ggcyto::ggcyto_par_set(limits = "data") +
    ggplot2::scale_fill_distiller(palette = "Spectral") +
    ggplot2::theme_minimal()
  print(p)
}
