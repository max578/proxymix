## Class validator, accessor, and print-branch contracts.
## These exercise the defensive error paths and the display branches that the
## happy-path tests do not reach.

test_that("gmm validator rejects malformed mixtures", {
  expect_error(gmm(weights = numeric(0), means = list(), covariances = list()),
               "at least one")
  expect_error(gmm(weights = c(0.5, 0.5), means = list(c(0, 0)),
                   covariances = list(diag(2), diag(2))),
               "same length")
  expect_error(gmm(weights = c(0.5, 0.5), means = list(c(0, 0), c(1, 1)),
                   covariances = list(diag(2))),
               "same length")
  expect_error(gmm(weights = c(-0.5, 1.5), means = list(c(0, 0), c(1, 1)),
                   covariances = list(diag(2), diag(2))),
               "non-negative")
  expect_error(gmm(weights = c(0.3, 0.3), means = list(c(0, 0), c(1, 1)),
                   covariances = list(diag(2), diag(2))),
               "sum to 1")
  expect_error(gmm(weights = c(0.5, 0.5), means = list(c(0, 0), c(1, 1, 1)),
                   covariances = list(diag(2), diag(2))),
               "length-2 numeric")
  expect_error(gmm(weights = c(0.5, 0.5), means = list(c(0, 0), c(1, 1)),
                   covariances = list(diag(2), diag(3))),
               "2x2 matrix")
})

test_that("gmm_dim and gmm_n_components require a gmm", {
  expect_error(gmm_dim(42), "must be a")
  expect_error(gmm_n_components("x"), "must be a")
})

test_that("gmm_target validator rejects malformed targets", {
  expect_error(gmm_target(n_dim = 0L, samples = matrix(1, 1, 1)),
               "positive integer")
  expect_error(gmm_target(n_dim = 1L), "at least one of")
  expect_error(gmm_target(n_dim = 1L, log_density = 42), "must be a function")
  expect_error(gmm_target(n_dim = 2L, samples = matrix(1:6, ncol = 3L)),
               "must have 2 columns")
  expect_error(gmm_target(n_dim = 1L, samples = matrix(1:4, ncol = 1L),
                          normalised = c(TRUE, FALSE)),
               "length-1 logical")
})

test_that("gmm_fit validator checks the regime tag", {
  expect_error(
    gmm_fit(weights = 1, means = list(0), covariances = list(matrix(1, 1, 1)),
            target = NULL, regime = c("a", "b")),
    "single string")
  expect_error(
    gmm_fit(weights = 1, means = list(0), covariances = list(matrix(1, 1, 1)),
            target = NULL, regime = "bogus"),
    "must be one of")
})

test_that("is_proposal validator rejects malformed proposals", {
  expect_error(
    is_proposal(n_dim = 0L, sample = function(n) 1, log_density = function(x) 1),
    "positive integer")
  expect_error(
    is_proposal(n_dim = 1L, sample = 42, log_density = function(x) 1),
    "function of one argument")
  expect_error(
    is_proposal(n_dim = 1L, sample = function(n) 1, log_density = 42),
    "function of one argument")
})

test_that("print methods cover fit, target-with-support, and overflow branches", {
  fit <- fit_proxymix(gmm_target_from_samples(
    withr::with_seed(69, matrix(stats::rnorm(200), ncol = 2))),
                      N = 1L, regime = "moment")
  expect_output(print(fit), "gmm_fit")

  expect_output(print(epanechnikov_target()), "support")

  K <- 7L
  g <- gmm(weights = rep(1 / K, K),
           means = lapply(seq_len(K), function(k) c(k, 0)),
           covariances = rep(list(diag(2)), K))
  expect_output(print(g), "more components")
})
