## Target-constructor error guards and the with_samples sampling branches.

test_that("target constructors reject bad input", {
  expect_error(gmm_target_from_samples(1:5), "numeric matrix")
  expect_error(epanechnikov_target(half_width = -1), "positive")
})

test_that("fit_proxymix rejects a non-target", {
  expect_error(fit_proxymix(42), "gmm_target")
})

test_that("banana_target attaches exact change-of-variables samples", {
  b <- banana_target(with_samples = TRUE, n = 300L, seed = 1L)
  expect_equal(dim(b@samples), c(300L, 2L))
  expect_true(all(is.finite(b@samples)))
})

test_that("donut_target attaches polar rejection samples", {
  d <- donut_target(with_samples = TRUE, n = 300L, seed = 1L)
  expect_equal(dim(d@samples), c(300L, 2L))
  expect_true(all(is.finite(d@samples)))
  ## Samples lie near the annulus radius.
  r <- sqrt(rowSums(d@samples^2))
  expect_equal(mean(r), 2.5, tolerance = 0.3)
})
