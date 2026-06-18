## E4: maximum-entropy targets.
##
## Independent oracles: the analytic Gaussian density (stats::dnorm and the
## hand-coded multivariate normal), numerical quadrature for normalisation and
## differential entropy, and the maximum-entropy ordering (the uniform is the
## maximum-entropy density on a compact support, so any other density on the
## same box has strictly smaller entropy).

# ---- full support: the Gaussian ------------------------------------------

test_that("full-support second-moment maxent equals the Gaussian to 1e-8", {
  g <- maxent_target(moments = list(mean = 0, cov = matrix(2)))
  xs <- matrix(seq(-3, 3, length.out = 21L), ncol = 1L)
  expect_equal(g@log_density(xs),
               stats::dnorm(xs[, 1L], 0, sqrt(2), log = TRUE),
               tolerance = 1e-10)
  expect_true(isTRUE(g@normalised))
  expect_equal(g@metadata$family, "maxent_gaussian")
  ## 2-D against a hand-coded multivariate normal.
  mu <- c(1, -1)
  S <- matrix(c(2, 0.5, 0.5, 1), 2, 2)
  g2 <- maxent_target(moments = list(mean = mu, cov = S))
  X <- withr::with_seed(1L, matrix(stats::rnorm(40), ncol = 2L))
  oracle <- apply(X, 1L, function(r) {
    as.numeric(-0.5 * t(r - mu) %*% solve(S) %*% (r - mu) -
                 0.5 * log((2 * pi)^2 * det(S)))
  })
  expect_equal(g2@log_density(X), oracle, tolerance = 1e-10)
})

test_that("regime (i) moment matching recovers the Gaussian maxent target", {
  mu <- c(1, -1)
  S <- matrix(c(2, 0.5, 0.5, 1), 2, 2)
  g <- maxent_target(moments = list(mean = mu, cov = S))
  fit <- fit_moment_match(g, N = 1L, ridge_eps = 1e-10)
  expect_equal(fit@means[[1L]], mu, tolerance = 1e-8)
  expect_equal(fit@covariances[[1L]], S, tolerance = 1e-7)
})

# ---- support alone: the uniform ------------------------------------------

test_that("support-only maxent is the uniform with entropy log(volume)", {
  u <- maxent_target(support = list(lower = c(0, 0), upper = c(2, 3)))
  expect_equal(exp(u@log_density(matrix(c(1, 1), nrow = 1L))), 1 / 6)
  expect_equal(u@log_density(matrix(c(3, 1), nrow = 1L)), -Inf)
  expect_true(isTRUE(u@normalised))
  expect_equal(u@metadata$log_volume, log(6))
  ## 1-D: integrates to one and has differential entropy log(volume).
  u1 <- maxent_target(support = list(lower = -1, upper = 2))
  dens <- function(z) exp(u1@log_density(matrix(z, ncol = 1L)))
  expect_equal(stats::integrate(dens, -1, 2)$value, 1, tolerance = 1e-8)
  h <- -stats::integrate(function(z) {
    p <- dens(z)
    ifelse(p > 0, p * log(p), 0)
  }, -1, 2)$value
  expect_equal(h, log(3), tolerance = 1e-8)
})

test_that("the uniform maxent requires a finite box", {
  expect_error(
    maxent_target(support = list(lower = 0, upper = Inf)),
    "finite box"
  )
})

# ---- bounded support with moments: the truncated Gaussian ----------------

test_that("box second-moment maxent is a normalised truncated Gaussian (diagonal)", {
  tg <- maxent_target(moments = list(mean = 0, cov = matrix(4)),
                      support = list(lower = -3, upper = 3))
  expect_true(isTRUE(tg@normalised))
  expect_equal(tg@metadata$family, "maxent_truncated_gaussian")
  dens <- function(z) exp(tg@log_density(matrix(z, ncol = 1L)))
  ## integrates to one (closed-form box normaliser).
  expect_equal(stats::integrate(dens, -3, 3)$value, 1, tolerance = 1e-7)
  ## maximum-entropy ordering: the uniform on the same box dominates.
  h_tg <- -stats::integrate(function(z) {
    p <- dens(z)
    ifelse(p > 0, p * log(p), 0)
  }, -3, 3)$value
  expect_lt(h_tg, log(6))
  ## outside the support the log-density is -Inf.
  expect_equal(tg@log_density(matrix(5, ncol = 1L)), -Inf)
})

test_that("a non-diagonal box truncated Gaussian is declared unnormalised", {
  S <- matrix(c(1, 0.5, 0.5, 1), 2, 2)
  tg <- maxent_target(moments = list(mean = c(0, 0), cov = S),
                      support = list(lower = c(-2, -2), upper = c(2, 2)))
  expect_false(isTRUE(tg@normalised))
  expect_true(is.na(tg@log_normalizer))
  expect_false(is.null(tg@support))
})

test_that("a bounded maxent target fits via regime (iii) under the auto proposal", {
  tg <- maxent_target(moments = list(mean = 0, cov = matrix(4)),
                      support = list(lower = -3, upper = 3))
  fit <- fit_kld_em(tg, N = 3L, is_size = 4000L, seed = 1L, max_iter = 40L)
  ## the declared support drives automatic support-matched proposal selection.
  expect_match(fit@diagnostics$proposal_name, "uniform")
  expect_equal(fit@diagnostics$support_fraction, 1, tolerance = 1e-9)
  expect_true(is.finite(fit@diagnostics$kld_final))
  expect_false(anyNA(fit@weights))
})

# ---- input validation -----------------------------------------------------

test_that("maxent_target rejects ill-posed constraint sets", {
  expect_error(maxent_target(), "at least one constraint")
  expect_error(maxent_target(moments = list(mean = 0)),
               "both `mean` and `cov`")
  expect_error(
    maxent_target(moments = list(mean = 0, cov = matrix(-1))),
    "positive-definite"
  )
  expect_error(
    maxent_target(moments = list(mean = c(0, 0), cov = matrix(1))),
    "2-by-2"
  )
})
