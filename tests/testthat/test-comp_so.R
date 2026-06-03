test_that("spillover_matrix computes compensation from cytoset controls", {
  make_control <- function(fl1, fl3) {
    n <- length(fl1)
    flowCore::flowFrame(cbind(
      `FSC-A` = seq(1e5, 2e5, length.out = n),
      `SSC-A` = seq(8e4, 1.6e5, length.out = n),
      `FL1-A` = fl1,
      `FL3-A` = fl3
    ))
  }

  n <- 80L
  jitter <- seq(-5, 5, length.out = n)
  controls <- flowCore::flowSet(list(
    unstained = make_control(100 + jitter, 100 - jitter),
    fl1 = make_control(1000 + jitter, 200 - jitter),
    fl3 = make_control(150 + jitter, 900 - jitter)
  ))
  controls <- flowWorkspace::flowSet_to_cytoset(controls)

  result <- spillover_matrix(
    controls,
    cont_ind = c(1, 2, 3),
    comp_pattern = "-A",
    threshold = 0,
    manual_comp = FALSE,
    interactive = FALSE
  )

  expect_true(is.matrix(result))
  expect_equal(rownames(result), c("FL1-A", "FL3-A"))
  expect_equal(colnames(result), c("FL1-A", "FL3-A"))
  expect_true(all(is.finite(result)))
})
