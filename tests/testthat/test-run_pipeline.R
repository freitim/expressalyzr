test_that("remove later", {
  test_file_path <- system.file("extdata", "example_fcs_files",
                                package = "expressalyzr", mustWork = TRUE)
  fcs_files <- list.files(file.path(test_file_path, "data"), full.names = TRUE)
  test_dir <- file.path(test_file_path, "test")
  test_data_dir <- file.path(test_dir, "data")

  expect_true(dir.create(test_data_dir, recursive = TRUE))
  expect_true(all(file.copy(fcs_files, test_data_dir)))

  # run_pipeline(test_dir)

  expect_equal(unlink(test_dir, recursive = TRUE), 0)
  expect_false(dir.exists(test_dir))
})
