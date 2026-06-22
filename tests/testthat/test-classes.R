test_that("gmm constructs and validates", {
  g <- gmm(weights = c(0.4, 0.6),
           means = list(c(-1, 0), c(1, 0)),
           covariances = list(diag(2), diag(2)))
  expect_s7_class(g, gmm)
  expect_equal(gmm_dim(g), 2L)
  expect_equal(gmm_n_components(g), 2L)
})

test_that("gmm rejects non-summing-to-one weights", {
  expect_error(
    gmm(weights = c(0.4, 0.4),
        means = list(c(0, 0), c(1, 1)),
        covariances = list(diag(2), diag(2))),
    "sum to 1"
  )
})

test_that("gmm rejects mismatched component shapes", {
  expect_error(
    gmm(weights = c(0.5, 0.5),
        means = list(c(0, 0), c(1, 1, 1)),
        covariances = list(diag(2), diag(2))),
    "length-2 numeric vector"
  )
})

test_that("gmm rejects a covariance that is not positive-definite", {
  expect_error(
    gmm(weights = c(0.5, 0.5),
        means = list(0, 1),
        covariances = list(matrix(-1), matrix(1))),
    "positive-definite"
  )
  ## an indefinite 2x2 (one negative eigenvalue) is rejected too
  expect_error(
    gmm(weights = 1, means = list(c(0, 0)),
        covariances = list(matrix(c(1, 0, 0, -2), 2))),
    "positive-definite"
  )
  ## a valid (near-singular but PSD) covariance still constructs
  expect_s7_class(
    gmm(weights = 1, means = list(c(0, 0)),
        covariances = list(matrix(c(1, 0, 0, 1e-10), 2))),
    gmm
  )
})

test_that("gmm_fit inherits from gmm", {
  x <- matrix(stats::rnorm(200), ncol = 2)
  tgt <- gmm_target_from_samples(x)
  fit <- fit_em_samples(tgt, N = 2L, max_iter = 25L, n_starts = 2L)
  expect_s7_class(fit, gmm_fit)
  expect_s7_class(fit, gmm)
  expect_equal(fit@regime, "sample")
})

test_that("gmm_target requires log_density or samples", {
  expect_error(
    gmm_target(n_dim = 2L),
    "log_density.*samples"
  )
})

test_that("gmm_target_from_samples carries the matrix verbatim", {
  x <- matrix(stats::rnorm(200), ncol = 2)
  tgt <- gmm_target_from_samples(x)
  expect_identical(tgt@samples, x)
  expect_equal(tgt@n_dim, 2L)
})

test_that("is_proposal builds correctly", {
  q <- is_mvn(n_dim = 2L)
  expect_s7_class(q, is_proposal)
  expect_equal(q@n_dim, 2L)
  expect_true(is.function(q@sample))
  expect_true(is.function(q@log_density))
})

test_that("print methods do not error", {
  g <- gmm(weights = 1, means = list(c(0, 0)),
           covariances = list(diag(2)))
  expect_invisible(print(g))
  tgt <- banana_target()
  expect_invisible(print(tgt))
  q <- is_mvn(n_dim = 2L)
  expect_invisible(print(q))
})
