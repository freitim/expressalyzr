test_that("pipeline gating template resolver preserves legacy default", {
  test_dir <- tempfile("template_resolver_")
  dir.create(test_dir)
  on.exit(unlink(test_dir, recursive = TRUE), add = TRUE)
  default_template <- file.path(test_dir, "gt_samples.csv")
  writeLines(c(
    "alias,pop,parent,dims,gating_method,gating_args",
    "nondebris,+,root,\"FSC-A,SSC-A\",boundary,",
    "cells,+,nondebris,\"FSC-A,SSC-A\",density_gate,"
  ), default_template)

  specs <- expressalyzr:::resolve_pipeline_gating_templates(
    test_dir,
    default_template = default_template
  )

  expect_equal(nrow(specs), 1)
  expect_false(specs$explicit)
  expect_true(is.na(specs$gating_template))
  expect_true(is.na(specs$gated_population))
  expect_equal(specs$extraction_population, "cells")
  expect_equal(specs$template_file, default_template)
})

test_that("missing default gating template is created only interactively", {
  test_dir <- tempfile("template_resolver_")
  dir.create(test_dir)
  on.exit(unlink(test_dir, recursive = TRUE), add = TRUE)
  default_template <- file.path(test_dir, "gt_samples.csv")

  expect_error(
    expressalyzr:::ensure_default_gating_template(
      default_template,
      interactive = FALSE
    ),
    "gt_samples.csv does not exist"
  )

  expect_message(
    expressalyzr:::ensure_default_gating_template(
      default_template,
      interactive = TRUE
    ),
    "Created starter gating template"
  )
  expect_true(file.exists(default_template))
  expect_silent(
    expressalyzr:::ensure_default_gating_template(
      default_template,
      interactive = FALSE
    )
  )
})

test_that("pipeline gating template resolver names explicit templates", {
  test_dir <- tempfile("template_resolver_")
  dir.create(test_dir)
  on.exit(unlink(test_dir, recursive = TRUE), add = TRUE)
  writeLines(c(
    "alias,pop,parent,dims,gating_method,gating_args",
    "nondebris,+,root,\"FSC-A,SSC-A\",boundary,",
    "singlets,+,nondebris,\"FSC-A,FSC-H\",singletGate,"
  ), file.path(test_dir, "gt_population_a.csv"))
  writeLines(c(
    "alias,pop,parent,dims,gating_method,gating_args",
    "nondebris,+,root,\"FSC-A,SSC-A\",boundary,",
    "cells,+,nondebris,\"FSC-A,FSC-H\",singletGate,"
  ), file.path(test_dir, "gt_population_b.csv"))

  specs <- expressalyzr:::resolve_pipeline_gating_templates(
    test_dir,
    default_template = file.path(test_dir, "gt_samples.csv"),
    gating_templates = c(
      HEK = "gt_population_a.csv",
      Tcells = "gt_population_b.csv"
    )
  )

  expect_equal(specs$gating_template,
               c("gt_population_a", "gt_population_b"))
  expect_equal(specs$gated_population, c("HEK", "Tcells"))
  expect_equal(specs$extraction_population, c("singlets", "cells"))
  expect_true(all(specs$explicit))
  expect_true(all(file.exists(specs$template_file)))
})

test_that("pipeline gating template resolver creates missing explicit templates interactively", {
  test_dir <- tempfile("template_resolver_")
  dir.create(test_dir)
  on.exit(unlink(test_dir, recursive = TRUE), add = TRUE)

  default_template <- file.path(test_dir, "gt_samples.csv")
  writeLines(c(
    "alias,pop,parent,dims,gating_method,gating_args",
    "nondebris,+,root,\"FSC-A,SSC-A\",boundary,",
    "singlets,+,nondebris,\"FSC-A,FSC-H\",singletGate,"
  ), default_template)

  expect_message(
    specs <- expressalyzr:::resolve_pipeline_gating_templates(
      test_dir,
      default_template = default_template,
      gating_templates = c(
        HEK = "gt_HEK.csv",
        Tcells = "gt_Tcells.csv"
      ),
      interactive = TRUE
    ),
    "Created starter population gating template"
  )

  expect_equal(specs$gating_template, c("gt_HEK", "gt_Tcells"))
  expect_equal(specs$gated_population, c("HEK", "Tcells"))
  expect_equal(specs$extraction_population, c("singlets", "singlets"))
  expect_true(all(file.exists(specs$template_file)))
  expect_equal(readLines(file.path(test_dir, "gt_HEK.csv")),
               readLines(default_template))
  expect_equal(readLines(file.path(test_dir, "gt_Tcells.csv")),
               readLines(default_template))
})

test_that("pipeline gating template resolver does not create explicit templates non-interactively", {
  test_dir <- tempfile("template_resolver_")
  dir.create(test_dir)
  on.exit(unlink(test_dir, recursive = TRUE), add = TRUE)

  expect_error(
    expressalyzr:::resolve_pipeline_gating_templates(
      test_dir,
      default_template = file.path(test_dir, "gt_samples.csv"),
      gating_templates = c(HEK = "gt_HEK.csv"),
      interactive = FALSE
    ),
    "cannot be created in non-interactive mode"
  )
  expect_false(file.exists(file.path(test_dir, "gt_HEK.csv")))
})

test_that("pipeline gating template resolver accepts explicit population labels", {
  test_dir <- tempfile("template_resolver_")
  dir.create(test_dir)
  on.exit(unlink(test_dir, recursive = TRUE), add = TRUE)
  writeLines(c(
    "alias,pop,parent,dims,gating_method,gating_args",
    "nondebris,+,root,\"FSC-A,SSC-A\",boundary,"
  ), file.path(test_dir, "gt_a.csv"))
  writeLines(c(
    "alias,pop,parent,dims,gating_method,gating_args",
    "cells,+,root,\"FSC-A,SSC-A\",boundary,"
  ), file.path(test_dir, "gt_b.csv"))

  specs <- expressalyzr:::resolve_pipeline_gating_templates(
    test_dir,
    default_template = file.path(test_dir, "gt_samples.csv"),
    gating_templates = c("gt_a.csv", "gt_b.csv"),
    gated_populations = c("HEK", "Tcells")
  )

  expect_equal(specs$gating_template, c("gt_a", "gt_b"))
  expect_equal(specs$gated_population, c("HEK", "Tcells"))
  expect_error(
    expressalyzr:::resolve_pipeline_gating_templates(
      test_dir,
      default_template = file.path(test_dir, "gt_samples.csv"),
      gating_templates = c("gt_a.csv", "gt_b.csv"),
      gated_populations = "HEK"
    ),
    "one label per gating template"
  )
})
