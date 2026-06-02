## Tests for autoplot.gmm_fit() and its internal helpers.

.helper_sample_fit <- function(p = 2L, N = 2L, seed = 1L) {
  samples <- .helper_iid_normal(n = 300, p = p, seed = seed)
  tgt <- gmm_target_from_samples(samples)
  fit_proxymix(tgt, N = N, regime = "sample", max_iter = 20L)
}

test_that("autoplot() returns a ggplot for a two-dimensional fit", {
  skip_if_not_installed("ggplot2")
  fit <- .helper_sample_fit(p = 2L)
  plt <- ggplot2::autoplot(fit)
  expect_s3_class(plt, "ggplot")
  expect_gt(length(plt$layers), 0L)
})

test_that("autoplot() renders a one-dimensional marginal on request", {
  skip_if_not_installed("ggplot2")
  fit <- .helper_sample_fit(p = 2L)
  plt <- ggplot2::autoplot(fit, dims = 1L)
  expect_s3_class(plt, "ggplot")
})

test_that("autoplot() marginalises a higher-dimensional fit to the asked dims", {
  skip_if_not_installed("ggplot2")
  fit <- .helper_sample_fit(p = 3L)
  expect_s3_class(ggplot2::autoplot(fit, dims = c(1L, 3L)), "ggplot")
  expect_s3_class(ggplot2::autoplot(fit, dims = 2L), "ggplot")
})

test_that("autoplot() honours the overlay switches", {
  skip_if_not_installed("ggplot2")
  fit <- .helper_sample_fit(p = 2L)
  plt <- ggplot2::autoplot(
    fit,
    show_components = FALSE,
    show_data = FALSE
  )
  expect_s3_class(plt, "ggplot")
})

test_that("autoplot() rejects malformed `dims`", {
  skip_if_not_installed("ggplot2")
  fit <- .helper_sample_fit(p = 2L)
  expect_error(ggplot2::autoplot(fit, dims = c(1L, 1L)), "unique")
  expect_error(ggplot2::autoplot(fit, dims = 5L), "indices")
  expect_error(ggplot2::autoplot(fit, dims = c(1L, 2L, 3L)), "one or two")
})

test_that(".gmm_ellipse() traces a closed bivariate contour", {
  ell <- proxymix:::.gmm_ellipse(c(1, -1), diag(c(2, 0.5)), level = 0.9)
  expect_named(ell, c("x", "y"))
  expect_equal(
    as.numeric(ell[1L, ]),
    as.numeric(ell[nrow(ell), ]),
    tolerance = 1e-8
  )
})
