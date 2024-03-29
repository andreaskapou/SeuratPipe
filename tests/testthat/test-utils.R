context("Checking utils functions work")

test_that("plot dimensions works", {
  feat_len <- 2
  expect_identical(.plot_dims(feat_len), list(width = 15, height = 6, ncols = 2))
  feat_len <- 4
  expect_identical(.plot_dims(feat_len), list(width = 30, height = 6, ncols = 4))
})
