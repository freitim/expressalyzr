#### helper: set up gated data for report test ####
setup_gated_data <- function() {
  test_path <- system.file("extdata", "example_fcs_files",
                           package = "expressalyzr", mustWork = TRUE)

  tmp_dir <- file.path(tempdir(), paste0("report_test_", Sys.getpid()))
  data_dir <- file.path(tmp_dir, "data")
  dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)
  fcs_files <- list.files(file.path(test_path, "data"), full.names = TRUE)
  file.copy(fcs_files, data_dir)

  cs <- load_fcs(tmp_dir)
  gs <- flowWorkspace::GatingSet(cs)

  gt_file <- system.file("tools", "gt_samples.csv",
                         package = "expressalyzr", mustWork = TRUE)
  openCyto::register_plugins(fun = density_gate, "density_gate", dep = NA, "gating")
  openCyto::register_plugins(fun = mixture_gate, "mixture_gate", dep = NA, "gating")
  gt <- openCyto::gatingTemplate(gt_file)
  openCyto::gt_gating(gt, gs)

  data_cs <- flowWorkspace::gs_pop_get_data(gs, y = "singlets")
  dt <- cs_to_dt(data_cs)
  fl_chs <- colnames(dt)[grepl("FL|488|445|640", colnames(dt))]

  list(tmp_dir = tmp_dir, gs = gs, data_cs = data_cs, dt = dt, fl_chs = fl_chs)
}

#### test generate_qc_report ####

test_that("generate_qc_report produces an HTML file", {
  env <- setup_gated_data()
  on.exit(unlink(env$tmp_dir, recursive = TRUE), add = TRUE)

  qc_results <- collect_qc(env$gs, env$data_cs, env$dt, env$fl_chs)

  out_file <- generate_qc_report(env$tmp_dir, env$gs, qc_results)

  expect_true(file.exists(out_file))
  expect_match(out_file, "\\.html$")

  html_content <- readLines(out_file, warn = FALSE)
  html_text <- paste(html_content, collapse = "\n")
  expect_match(html_text, "Event Counts", fixed = TRUE)
  expect_match(html_text, "Scatter Quality", fixed = TRUE)
  expect_match(html_text, "Per-well gate scatter plots are omitted",
               fixed = TRUE)
  expect_false(grepl("Gate Plots", html_text, fixed = TRUE))
  expect_false(grepl("Could not render gate plot", html_text,
                     fixed = TRUE))
})

test_that("generate_qc_report renders separate gated population sections", {
  env <- setup_gated_data()
  on.exit(unlink(env$tmp_dir, recursive = TRUE), add = TRUE)

  base_qc <- collect_qc(env$gs, env$data_cs, env$dt, env$fl_chs)
  tag_table <- function(tbl, population, template) {
    if (is.null(tbl)) {
      return(NULL)
    }
    population_label <- population
    template_label <- template
    tbl <- data.table::copy(tbl)
    tbl[, gated_population := population_label]
    tbl[, gating_template := template_label]
    tbl
  }
  tag_qc <- function(population, template) {
    lapply(base_qc, tag_table, population = population, template = template)
  }
  hek_qc <- tag_qc("HEK", "gt_HEK")
  tcells_qc <- tag_qc("Tcells", "gt_Tcells")
  combine_metric <- function(metric) {
    tables <- Filter(Negate(is.null), list(hek_qc[[metric]], tcells_qc[[metric]]))
    if (length(tables) == 0L) {
      return(NULL)
    }
    data.table::rbindlist(tables, use.names = TRUE, fill = TRUE)
  }
  qc_results <- list(
    event_counts = combine_metric("event_counts"),
    scatter_stats = combine_metric("scatter_stats"),
    channel_stats = combine_metric("channel_stats"),
    time_check = combine_metric("time_check")
  )

  tagged_dt <- data.table::rbindlist(lapply(c("HEK", "Tcells"), function(pop) {
    dt <- data.table::copy(env$dt)
    dt[, gated_population := pop]
    dt[, gating_template := paste0("gt_", pop)]
    dt
  }), use.names = TRUE, fill = TRUE)

  population_table <- data.table::rbindlist(lapply(c("HEK", "Tcells"), function(pop) {
    table <- compute_population_table(env$gs, experiment = basename(env$tmp_dir))
    table[, gated_population := pop]
    table[, gating_template := paste0("gt_", pop)]
    table
  }), use.names = TRUE, fill = TRUE)

  marker_table <- data.table::rbindlist(lapply(c("HEK", "Tcells"), function(pop) {
    table <- compute_marker_table(env$dt, env$fl_chs,
                                  experiment = basename(env$tmp_dir),
                                  population = pop)
    table[, gated_population := pop]
    table[, gating_template := paste0("gt_", pop)]
    table
  }), use.names = TRUE, fill = TRUE)

  gating_results <- list(
    list(label = "gt_HEK", population_label = "HEK", gs = env$gs),
    list(label = "gt_Tcells", population_label = "Tcells", gs = env$gs)
  )

  out_file <- generate_qc_report(
    env$tmp_dir,
    env$gs,
    qc_results,
    gating_results = gating_results,
    data_dt = tagged_dt,
    population_table = population_table,
    marker_table = marker_table,
    summaries_written = TRUE,
    include_gating_plots = TRUE
  )

  html_text <- paste(readLines(out_file, warn = FALSE), collapse = "\n")
  expect_match(html_text, "HEK Gate Plots", fixed = TRUE)
  expect_match(html_text, "Tcells Gate Plots", fixed = TRUE)
  expect_match(html_text, "Gating Templates", fixed = TRUE)
  expect_false(grepl("Population and Marker Summary", html_text,
                     fixed = TRUE))
  expect_match(html_text, "root -&gt; nondebris", fixed = TRUE)
  expect_false(grepl("Could not render gate plot", html_text,
                     fixed = TRUE))
  expect_false(grepl("Per-well gate scatter plots are omitted", html_text,
                     fixed = TRUE))
})
