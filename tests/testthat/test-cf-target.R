## Expected values are exact: the Cauchy distribution function is closed form.

gamma_noise_cf <- function(t) (1 - 1i * t[, 1L])^(-3) * exp(-0.125 * t[, 1L]^2)

test_that("cf_target() checks the characteristic function", {
  tgt <- cf_target(gamma_noise_cf, n_dim = 1L)
  expect_true(S7::S7_inherits(tgt, gmm_target))
  expect_identical(tgt@cf, gamma_noise_cf)
  expect_null(tgt@log_density)
  expect_null(tgt@samples)
  expect_error(cf_target(gamma_noise_cf, n_dim = 0L), "positive integer")
  expect_error(cf_target("not a function", n_dim = 1L), "must be a function")
  expect_error(
    cf_target(function(t) 2 * gamma_noise_cf(t), n_dim = 1L),
    class = "proxymix_invalid_cf"
  )
  expect_error(gmm_target(n_dim = 1L, cf = 1), "`cf` must be a function")
})

test_that("cf_target() leaves the random number generator unchanged", {
  set.seed(5L)
  before <- .Random.seed
  cf_target(gamma_noise_cf, n_dim = 1L)
  expect_identical(.Random.seed, before)
})

test_that("existing targets print as before and a cf target says so", {
  expect_false(any(grepl("cf ", capture.output(print(banana_target())))))
  expect_true(any(grepl(
    "cf +: supplied", capture.output(print(cf_target(gamma_noise_cf, 1L)))
  )))
})

test_that("fit_proxymix() picks regime cf for a cf-only target", {
  tgt <- cf_target(gamma_noise_cf, n_dim = 1L)
  expect_message(
    fit <- fit_proxymix(tgt, N = 2L, n_freq = 1024L, seed = 1L),
    "cf"
  )
  direct <- from_cf(gamma_noise_cf, 1L, 2L, n_freq = 1024L, seed = 1L)
  expect_identical(fit@regime, "cf")
  expect_identical(fit@target, tgt)
  expect_equal(fit@weights, direct@weights)
  expect_equal(fit@means, direct@means)
  expect_equal(fit@covariances, direct@covariances)
  explicit <- fit_proxymix(tgt, N = 2L, regime = "cf", n_freq = 1024L,
                           seed = 1L)
  expect_equal(explicit@weights, direct@weights)
})

test_that("regime cf accepts and ignores canonicalise", {
  tgt <- cf_target(gamma_noise_cf, n_dim = 1L)
  fit <- fit_proxymix(tgt, N = 2L, regime = "cf", canonicalise = FALSE,
                      n_freq = 1024L, seed = 1L)
  direct <- from_cf(gamma_noise_cf, 1L, 2L, n_freq = 1024L, seed = 1L)
  expect_equal(fit@means, direct@means)
})

test_that("regime cf needs a characteristic function", {
  expect_error(
    fit_proxymix(banana_target(), N = 2L, regime = "cf"),
    "characteristic function"
  )
})

test_that("auto keeps the existing regimes when a target has more than a cf", {
  x <- withr::with_seed(2L, matrix(stats::rnorm(400L), ncol = 2L))
  tgt <- gmm_target(n_dim = 2L, samples = x,
                    cf = function(t) exp(-0.5 * rowSums(t^2)))
  expect_message(fit <- fit_proxymix(tgt, N = 2L, max_iter = 10L), "sample")
  expect_identical(fit@regime, "sample")
})

test_that("the bracket is on by default, changes nothing else, and can be skipped", {
  truth <- gmm(weights = c(0.4, 0.6), means = list(c(-2, 0), c(2, 1)),
               covariances = list(diag(2), matrix(c(1, 0.5, 0.5, 2), 2L)))
  cf <- function(t) gmm_cf(truth, t)
  with_b <- from_cf(cf, 2L, 2L, n_freq = 1024L, seed = 1L)
  without <- from_cf(cf, 2L, 2L, n_freq = 1024L, seed = 1L, bracket = FALSE)
  expect_identical(with_b@weights, without@weights)
  expect_identical(with_b@means, without@means)
  tab <- with_b@diagnostics$bracket
  expect_identical(nrow(tab), 8L)
  expect_equal(with_b@diagnostics$bracket_directions[1:2, ], diag(2))
  expect_equal(rowSums(with_b@diagnostics$bracket_directions^2), rep(1, 8L))
  expect_true(all(tab$distance <= tab$upper))
  q <- gmm_fit_quality(with_b)
  expect_identical(q$cf_distance, max(tab$distance))
  expect_identical(q$cf_upper, max(tab$upper))
  expect_identical(
    with_b@diagnostics$bracket_directions,
    from_cf(cf, 2L, 2L, n_freq = 1024L, seed = 1L)@diagnostics$bracket_directions
  )
  expect_null(without@diagnostics$bracket)
  expect_true(is.na(gmm_fit_quality(without)$cf_upper))
  expect_error(from_cf(cf, 2L, 2L, bracket = NA), "TRUE")
})

test_that("the bracket leaves the random number generator unchanged", {
  cf <- function(t) exp(-0.5 * rowSums(t^2))
  set.seed(3L)
  from_cf(cf, 2L, 1L, n_freq = 512L, bracket = FALSE)
  r_without <- stats::runif(1L)
  set.seed(3L)
  from_cf(cf, 2L, 1L, n_freq = 512L)
  expect_identical(stats::runif(1L), r_without)
})

test_that("the bracket contains the exact Kolmogorov distance for a Cauchy target", {
  fit <- suppressWarnings(
    from_cf(function(t) exp(-abs(t[, 1L])), 1L, 4L, n_freq = 2048L, seed = 1L),
    classes = "proxymix_cf_heavy_tail"
  )
  q <- gmm_fit_quality(fit)
  x <- seq(-60, 60, length.out = 120001L)
  exact <- max(abs(stats::pcauchy(x) - pgmm(x, fit)))
  ## the lower bound is a grid maximum inside the 0.005 and 0.995
  ## quantiles; the upper bound covers the whole line
  expect_lte(q$cf_distance, exact + 1e-6)
  expect_gte(q$cf_upper, exact)
  ## and the bracket is narrow, so a trivial one would fail
  expect_gt(q$cf_distance, 0.5 * exact)
  expect_lt(q$cf_upper - q$cf_distance, 0.01)
})

test_that("glance() and print() report the bracket", {
  skip_if_not_installed("generics")
  fit <- from_cf(gamma_noise_cf, 1L, 2L, n_freq = 1024L, seed = 1L)
  gl <- generics::glance(fit)
  expect_identical(gl$regime, "cf")
  expect_identical(gl$cf_upper, gmm_fit_quality(fit)$cf_upper)
  expect_identical(gl$cf_heldout_rms, fit@diagnostics$heldout_rms)
  out <- capture.output(print(fit))
  expect_true(any(grepl("Kolmogorov distance +: at least", out)))
  out_off <- capture.output(print(
    from_cf(gamma_noise_cf, 1L, 2L, n_freq = 1024L, seed = 1L, bracket = FALSE)
  ))
  expect_true(any(grepl("not computed", out_off)))
  expect_true(is.na(generics::glance(
    fit_proxymix(banana_target(), N = 2L, regime = "kld", is_size = 500L,
                 max_iter = 5L, seed = 1L)
  )$cf_upper))
})

test_that("the structured fits carry the bracket", {
  S <- matrix(c(2, 0.5, 0.5, 1), 2L)
  cf <- function(t) exp(1i * t[, 1] - sqrt(rowSums((t %*% S) * t)))
  fit <- suppressWarnings(from_cf_elliptical(cf, n_dim = 2L, seed = 1L))
  expect_identical(fit@regime, "cf")
  expect_identical(nrow(fit@diagnostics$bracket), 8L)
  expect_true(is.finite(gmm_fit_quality(fit)$cf_upper))
  A <- matrix(c(0.7, 0.3, 0, -0.7), 2L)
  cf2 <- function(t) {
    u <- t %*% A
    (1 - 1i * u[, 1])^(-3) * (1 - 1i * u[, 2])^(-3) * exp(-0.125 * rowSums(t^2))
  }
  src <- from_cf_sources(cf2, n_dim = 2L, N = 2L, n_freq = 1024L, seed = 1L)
  expect_identical(nrow(src$sources[[1L]]@diagnostics$bracket), 1L)
  src_off <- from_cf_sources(cf2, n_dim = 2L, N = 2L, n_freq = 1024L,
                             seed = 1L, bracket = FALSE)
  expect_null(src_off$sources[[1L]]@diagnostics$bracket)
  expect_equal(src_off$mixing, src$mixing)
})

test_that("operator results do not carry the bounds of the fit", {
  fit <- from_cf(gamma_noise_cf, 1L, 2L, n_freq = 1024L, seed = 1L)
  fit2 <- from_cf(gamma_noise_cf, 1L, 2L, n_freq = 1024L, seed = 2L)
  q_fit <- gmm_fit_quality(fit)
  expect_true(is.finite(q_fit$cf_upper))
  q <- gmm_fit_quality(gmm_marginalise(fit, keep = 1L))
  expect_identical(q$regime, "cf")
  expect_identical(q$cf_upper, NA_real_)
  expect_identical(q$cf_distance, NA_real_)
  expect_identical(q$cf_heldout_rms, NA_real_)
  qp <- gmm_fit_quality(gmm_product(fit, fit2))
  expect_identical(qp$cf_upper, NA_real_)
  expect_identical(qp$quality_sources[[1L]]$quality$cf_upper, q_fit$cf_upper)
  expect_identical(qp$quality_sources[[2L]]$quality$cf_upper,
                   gmm_fit_quality(fit2)$cf_upper)
})
