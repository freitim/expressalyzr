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

test_that("files and folders are moved from the old to the new data directory", {
  # setup
  test_dir_path <- file.path(".", "dir_for_moving_data")
  dummy_file <- file.path(test_dir_path, "dummy.fcs")
  dummy_subdir <- file.path(test_dir_path, "dummy_subdir")
  expect_true(dir.create(test_dir_path))
  expect_true(file.create(dummy_file))
  expect_true(dir.create(dummy_subdir))
  # test
  expect_true(create_data_subdir(test_dir_path))
  expect_true(file.exists(file.path(test_dir_path, "data", "dummy.fcs")))
  expect_true(dir.exists(file.path(test_dir_path, "data", "dummy_subdir")))
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
