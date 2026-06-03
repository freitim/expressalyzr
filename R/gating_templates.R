ensure_default_gating_template <- function(template_file, interactive = TRUE) {
  if (file.exists(template_file)) {
    return(invisible(template_file))
  }

  if (!interactive) {
    stop("gt_samples.csv does not exist and cannot be created in non-interactive mode. Create gt_samples.csv first or run with interactive = TRUE.",
         call. = FALSE)
  }

  template_gating <- system.file("tools", "gt_samples.csv",
                                 package = "expressalyzr",
                                 mustWork = TRUE)
  ok <- file.copy(template_gating, template_file, overwrite = FALSE)
  if (!ok) {
    stop("Could not create gt_samples.csv from the package template.",
         call. = FALSE)
  }

  message("Created starter gating template: ", template_file)
  invisible(template_file)
}

starter_gating_template_path <- function() {
  system.file("tools", "gt_samples.csv",
              package = "expressalyzr",
              mustWork = TRUE)
}

copy_starter_gating_template <- function(template_file, source_template = NULL,
                                         interactive = TRUE,
                                         label = "gating template") {
  if (file.exists(template_file)) {
    return(invisible(template_file))
  }

  if (!interactive) {
    stop(label, " does not exist and cannot be created in non-interactive mode: ",
         template_file, ". Create the template first or run with interactive = TRUE.",
         call. = FALSE)
  }

  if (is.null(source_template) || !file.exists(source_template)) {
    source_template <- starter_gating_template_path()
  }

  dir.create(dirname(template_file), recursive = TRUE, showWarnings = FALSE)
  ok <- file.copy(source_template, template_file, overwrite = FALSE)
  if (!ok) {
    stop("Could not create ", label, " from starter template: ",
         template_file, call. = FALSE)
  }

  message("Created starter ", label, ": ", template_file)
  invisible(template_file)
}

resolve_pipeline_gating_templates <- function(data_path, default_template,
                                              gating_templates = NULL,
                                              gated_populations = NULL,
                                              interactive = TRUE) {
  if (is.null(gating_templates)) {
    return(data.table::data.table(
      gating_template = NA_character_,
      gated_population = NA_character_,
      extraction_population = infer_template_extraction_population(default_template),
      template_file = default_template,
      explicit = FALSE
    ))
  }

  if (is.list(gating_templates) && !is.data.frame(gating_templates)) {
    gating_templates <- unlist(gating_templates, use.names = TRUE)
  }
  if (!is.character(gating_templates) || length(gating_templates) == 0L) {
    stop("gating_templates must be NULL or a character vector/list of template files.",
         call. = FALSE)
  }
  if (any(!nzchar(gating_templates))) {
    stop("gating_templates cannot contain empty paths.", call. = FALSE)
  }

  population_labels <- names(gating_templates)
  missing_labels <- is.null(population_labels) || any(!nzchar(population_labels))
  if (missing_labels) {
    population_labels <- tools::file_path_sans_ext(basename(gating_templates))
  }

  if (!is.null(gated_populations)) {
    if (!is.character(gated_populations) ||
        length(gated_populations) != length(gating_templates) ||
        any(!nzchar(gated_populations))) {
      stop("gated_populations must be NULL or a non-empty character vector with one label per gating template.",
           call. = FALSE)
    }
    population_labels <- gated_populations
  }
  population_labels <- make.unique(population_labels, sep = "_")

  template_files <- vapply(gating_templates, function(path) {
    if (grepl("^~|^/", path)) {
      normalizePath(path, mustWork = FALSE)
    } else {
      normalizePath(file.path(data_path, path), mustWork = FALSE)
    }
  }, character(1L), USE.NAMES = FALSE)

  missing_files <- template_files[!file.exists(template_files)]
  if (length(missing_files) > 0L) {
    for (missing_file in missing_files) {
      copy_starter_gating_template(
        missing_file,
        source_template = default_template,
        interactive = interactive,
        label = "population gating template"
      )
    }
  }

  data.table::data.table(
    gating_template = tools::file_path_sans_ext(basename(template_files)),
    gated_population = population_labels,
    extraction_population = vapply(template_files,
                                   infer_template_extraction_population,
                                   character(1L)),
    template_file = template_files,
    explicit = TRUE
  )
}

infer_template_extraction_population <- function(template_file) {
  template <- read_gating_template(template_file)
  aliases <- template$alias[nzchar(template$alias)]
  if (length(aliases) == 0L) {
    stop("Gating template has no populations: ", template_file, call. = FALSE)
  }
  aliases[[length(aliases)]]
}

run_gating_template <- function(cs, template_file, extraction_population = NULL) {
  gt <- openCyto::gatingTemplate(template_file)
  gs <- flowWorkspace::GatingSet(cs)
  openCyto::gt_gating(gt, gs)
  if (is.null(extraction_population)) {
    extraction_population <- infer_template_extraction_population(template_file)
  }
  list(
    gt = gt,
    gs = gs,
    extraction_population = extraction_population,
    data_cs = flowWorkspace::gs_pop_get_data(gs, y = extraction_population)
  )
}
