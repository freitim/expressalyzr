#### test load_fcs function ####

## errors if non-existant file path or no .fcs files in the directory

test_that("non-existing file paths throw an error", {
  test_dir_path <- file.path(".", "does_not_exist")
  # create dummy directory just to remove again
  expect_true(dir.create(file.path(test_dir_path, "data"), recursive = TRUE))
  expect_equal(unlink(test_dir_path, recursive = TRUE), 0)
  expect_false(dir.exists(test_dir_path))
  expect_error(load_fcs(test_dir_path), "^.* not exist .*$")
})

test_that("existing file paths with no .fcs files throw an error", {
  test_dir_path <- file.path(".", "no_fcs_files")
  # create dummy directory for testing
  expect_true(dir.create(file.path(test_dir_path, "data"), recursive = TRUE))
  expect_true(dir.exists(test_dir_path))
  expect_error(load_fcs(test_dir_path), "^.* not contain any .*$")
  # clean up
  expect_equal(unlink(test_dir_path, recursive = TRUE), 0)
  expect_false(dir.exists(test_dir_path))
})

## return value

test_that("a cytoset is returned", {
  test_file_path <- system.file("extdata", "example_fcs_files",
                                package = "expressalyzr", mustWork = TRUE)
  expect_equal(is(load_fcs(test_file_path)), c("cytoset", "flowSet"))
})

test_that("load_fcs falls back to experiment-root .fcs files", {
  test_dir_path <- file.path(tempdir(), "root_level_fcs")
  dir.create(test_dir_path, showWarnings = FALSE)
  on.exit(unlink(test_dir_path, recursive = TRUE), add = TRUE)
  fcs_source <- list.files(
    system.file("extdata", "example_fcs_files", "data",
                package = "expressalyzr", mustWork = TRUE),
    pattern = "\\.fcs$",
    full.names = TRUE
  )[1L]
  file.copy(fcs_source, file.path(test_dir_path, basename(fcs_source)))
  expect_equal(is(load_fcs(test_dir_path)), c("cytoset", "flowSet"))
})


#### test create_data_subdir function ####

test_that("a new directory with the name 'data' is created", {
  # setup
  test_dir_path <- file.path(".", "dir_without_data_subdir")
  dummy_file <- file.path(test_dir_path, "dummy.fcs")
  dummy_subdir <- file.path(test_dir_path, "dummy_subdir")
  expect_true(dir.create(test_dir_path))
  expect_true(file.create(dummy_file))
  expect_true(dir.create(dummy_subdir))
  # test
  expect_true(create_data_subdir(test_dir_path))
  expect_true(dir.exists(file.path(test_dir_path, "data")))
  # clean up
  expect_equal(unlink(test_dir_path, recursive = TRUE), 0)
  expect_false(dir.exists(test_dir_path))
})

test_that("only FCS files are moved from the old to the new data directory", {
  # setup
  test_dir_path <- file.path(".", "dir_for_moving_data")
  dummy_file <- file.path(test_dir_path, "dummy.fcs")
  notes_file <- file.path(test_dir_path, "notes.txt")
  dummy_subdir <- file.path(test_dir_path, "dummy_subdir")
  expect_true(dir.create(test_dir_path))
  expect_true(file.create(dummy_file))
  expect_true(file.create(notes_file))
  expect_true(dir.create(dummy_subdir))
  # test
  expect_true(create_data_subdir(test_dir_path))
  expect_true(file.exists(file.path(test_dir_path, "data", "dummy.fcs")))
  expect_true(file.exists(notes_file))
  expect_true(dir.exists(dummy_subdir))
  expect_false(file.exists(file.path(test_dir_path, "data", "notes.txt")))
  expect_false(dir.exists(file.path(test_dir_path, "data", "dummy_subdir")))
  # clean up
  expect_equal(unlink(test_dir_path, recursive = TRUE), 0)
  expect_false(dir.exists(test_dir_path))
})

test_that("nothing happens if data subdirectory already exists", {
  # setup
  test_dir_path <- file.path(".", "dir_with_data_subdir")
  expect_true(dir.create(test_dir_path))
  expect_true(dir.create(file.path(test_dir_path, "data")))
  # test
  expect_null(create_data_subdir(test_dir_path))
  # clean up
  expect_equal(unlink(test_dir_path, recursive = TRUE), 0)
  expect_false(dir.exists(test_dir_path))
})


#### test config profiles ####

test_that("config profile templates exist and parse", {
  for (profile in c("default", "legacy", "advanced")) {
    template_path <- config_template_path(profile)
    expect_true(file.exists(template_path))

    parsed <- config::get(file = template_path)
    expect_type(parsed, "list")
    expect_false(is.null(parsed$controls_pattern))
    expect_equal(parsed$channel_map, "cytoflex")
  }
})

test_that("default and legacy config profiles encode different UX defaults", {
  default_config <- config::get(file = config_template_path("default"))
  legacy_config <- config::get(file = config_template_path("legacy"))
  advanced_config <- config::get(file = config_template_path("advanced"))

  expect_true(default_config$manual_comp)
  expect_true(default_config$manual_cutoff)
  expect_null(default_config$bg_channels)

  expect_true(legacy_config$manual_comp)
  expect_true(legacy_config$manual_cutoff)
  expect_equal(legacy_config$bg_channels, c("FL1-H", "FL3-H", "FL11-H"))

  expect_true(advanced_config$manual_comp)
  expect_true(advanced_config$manual_cutoff)
  expect_null(advanced_config$bg_channels)
})

test_that("load_config does not overwrite existing configs with profile defaults", {
  test_dir <- tempfile("existing_config_")
  dir.create(test_dir)
  on.exit(unlink(test_dir, recursive = TRUE), add = TRUE)

  config_file <- file.path(test_dir, "config.yml")
  writeLines(c(
    "default:",
    "  manual_comp: FALSE",
    "  manual_cutoff: FALSE",
    "  bg_channels: NULL",
    "  controls_pattern: \"NOMATCH\""
  ), config_file)

  parsed <- load_config(config_file, view_config = FALSE,
                        interactive = FALSE,
                        config_profile = "legacy")

  expect_false(parsed$manual_comp)
  expect_false(parsed$manual_cutoff)
  expect_null(parsed$bg_channels)
  expect_equal(parsed$controls_pattern, "NOMATCH")
})

test_that("missing config channel_map preserves CytoFLEX compatibility", {
  expect_equal(config_channel_map(list()), "cytoflex")
  expect_null(config_channel_map(list(channel_map = NULL)))
})

test_that("normalize_channels leaves unmapped FCS channels unchanged", {
  test_file_path <- system.file("extdata", "example_fcs_files",
                                package = "expressalyzr", mustWork = TRUE)
  cs <- load_fcs(test_file_path)
  before <- flowWorkspace::colnames(cs)

  expect_no_error(normalized <- normalize_channels(cs, "cytoflex"))
  after <- flowWorkspace::colnames(normalized)

  unmapped <- setdiff(before, names(cytoflex_channel_map()))
  expect_true(all(unmapped %in% after))
})

test_that("fluorescence_channel_names avoids CytoFLEX-only FL matching", {
  channels <- c("File", "Time", "FSC-A", "SSC-A", "FITC-A", "APC-A",
                "FL1-A", "FITC-W", "SSC-H")

  expect_equal(
    fluorescence_channel_names(channels, "-A"),
    c("FITC-A", "APC-A", "FL1-A")
  )
})

test_that("read_prompt_input reads from supplied input outside interactive R", {
  con <- textConnection("0.95")
  on.exit(close(con), add = TRUE)

  expect_equal(expressalyzr:::read_prompt_input("Adjust: ", input = con),
               "0.95")
})

test_that("read_prompt_input returns blank when input is exhausted", {
  con <- textConnection(character())
  on.exit(close(con), add = TRUE)

  expect_equal(expressalyzr:::read_prompt_input("Adjust: ", input = con), "")
})

test_that("shell_single_quote safely quotes terminal prompts", {
  quoted <- expressalyzr:::shell_single_quote("Adjust Bob's value: ")
  expect_equal(quoted, "'Adjust Bob'\"'\"'s value: '")
})

test_that("show_manual_plot writes a PNG outside interactive R", {
  old_options <- options(browser = function(url) invisible(TRUE))
  on.exit(options(old_options), add = TRUE)
  plot_dir <- tempfile("manual_plot_")
  on.exit(unlink(plot_dir, recursive = TRUE), add = TRUE)

  plot <- ggplot2::ggplot(data.frame(x = 1:3, y = 1:3),
                          ggplot2::aes(x = x, y = y)) +
    ggplot2::geom_point()
  path <- expressalyzr:::show_manual_plot(
    plot,
    prefix = "manual_test",
    output_dir = plot_dir
  )

  expect_true(file.exists(path))
  expect_match(basename(path), "^manual_test_.*\\.png$")
})
