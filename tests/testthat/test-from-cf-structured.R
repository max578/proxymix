## from_cf_elliptical() fits concentric normals to an elliptical
## characteristic function; from_cf_sources() fits a linear mix of
## independent sources. Expected values are exact quantiles and
## distribution functions where the target has them, and the target's own
## characteristic function otherwise.

shape_2d <- function() matrix(c(2, 0.5, 0.5, 1), 2L)

gamma_mix_cf <- function(A) {
  function(t) {
    u <- t %*% A
    (1 - 1i * u[, 1L])^(-3) * (1 - 1i * u[, 2L])^(-3) *
      exp(-0.125 * rowSums(t^2))
  }
}

# ---------------------------------------------------------------------------
# non-negative least squares
# ---------------------------------------------------------------------------

test_that(".cf_nnls returns non-negative x with zero gradient on the positive entries", {
  withr::local_seed(3L)
  A <- matrix(abs(rnorm(200L * 12L)), 200L, 12L)
  b <- as.vector(A %*% c(1, 0, 0, 2, 0, 0.5, rep(0, 6L))) + 0.1 * rnorm(200L)
  x <- .cf_nnls(A, b)
  slope <- as.vector(crossprod(A, b - A %*% x))
  expect_true(all(x >= 0))
  expect_lt(max(abs(slope[x > 0])), 1e-8)
  expect_lt(max(slope[x == 0]), 1e-8)
})

# ---------------------------------------------------------------------------
# from_cf_elliptical
# ---------------------------------------------------------------------------

test_that("a bivariate Cauchy is matched to its 0.999 quantile", {
  S <- shape_2d()
  cf <- function(t) exp(1i * t[, 1L] - sqrt(rowSums((t %*% S) * t)))
  fit <- from_cf_elliptical(cf, n_dim = 2L, seed = 1L)
  expect_equal(fit@diagnostics$shape_ratio, c(1, 1), tolerance = 1e-4)
  expect_equal(fit@diagnostics$tail_index, 1, tolerance = 1e-3)
  expect_equal(do.call(rbind, fit@means),
               matrix(c(1, 0), length(fit@weights), 2L, byrow = TRUE),
               tolerance = 1e-6)
  first <- gmm_marginalise(fit, 1L)
  probs <- c(0.9, 0.99, 0.999)
  expect_equal(
    vapply(probs, function(q) qgmm(q, first), numeric(1L)),
    1 + sqrt(2) * qcauchy(probs),
    tolerance = 0.02
  )
  expect_lt(gmm_cf_distance(fit, cf)$distance, 1e-3)
})

test_that("a bivariate Laplace is matched in distribution and covariance", {
  S <- shape_2d()
  cf <- function(t) 1 / (1 + 0.5 * rowSums((t %*% S) * t))
  fit <- from_cf_elliptical(cf, n_dim = 2L, seed = 1L)
  expect_equal(gmm_cov(fit), S, tolerance = 0.01)
  ## the first coordinate is Laplace with variance 2, so scale 1
  x <- seq(-8, 8, length.out = 401L)
  exact <- ifelse(x < 0, 0.5 * exp(x), 1 - 0.5 * exp(-x))
  expect_lt(max(abs(pgmm(x, gmm_marginalise(fit, 1L)) - exact)), 1e-3)
  expect_identical(gmm_fit_quality(fit)$regime, "cf")
})

test_that("a target that is not elliptical raises a warning", {
  two <- gmm(
    weights = c(0.4, 0.6),
    means = list(c(-2, 0), c(2, 1)),
    covariances = list(diag(2), matrix(c(1, 0.5, 0.5, 2), 2L))
  )
  expect_warning(
    from_cf_elliptical(function(t) gmm_cf(two, t), n_dim = 2L, seed = 1L),
    class = "proxymix_cf_not_elliptical"
  )
})

test_that("from_cf_elliptical is reproducible and validates its arguments", {
  cf <- function(t) exp(-sqrt(rowSums(t^2)))
  withr::local_seed(11L)
  before <- .Random.seed
  a <- from_cf_elliptical(cf, n_dim = 2L, N = 12L, seed = 5L)
  expect_identical(.Random.seed, before)
  b <- from_cf_elliptical(cf, n_dim = 2L, N = 12L, seed = 5L)
  expect_identical(a@weights, b@weights)
  expect_lte(length(a@weights), 12L)
  expect_error(from_cf_elliptical("nope", 2L), "must be a function")
  expect_error(from_cf_elliptical(cf, 2L, N = 1L), "at least 2")
  expect_error(from_cf_elliptical(cf, 2L, scale_range = c(2, 1)),
               "scale_range")
})

# ---------------------------------------------------------------------------
# from_cf_sources
# ---------------------------------------------------------------------------

test_that("two mixed Gamma sources are separated and fitted", {
  skip_on_cran()
  A <- matrix(c(0.7, 0.3, 0, -0.7), 2L)
  cf <- gamma_mix_cf(A)
  fit <- from_cf_sources(cf, n_dim = 2L, N = 3L, n_freq = 1024L, seed = 1L)
  expect_s3_class(fit, "proxymix_cf_sources")
  expect_lt(fit$diagnostics$off_diagonal, 0.01)
  expect_lt(fit$diagnostics$dependence, 1e-3)
  unit <- function(M) M / rep(sqrt(colSums(M^2)), each = 2L)
  cosine <- abs(crossprod(unit(fit$mixing), unit(A)))
  expect_gt(min(apply(cosine, 1L, max)), 1 - 1e-4)
  expect_lt(fit$diagnostics$heldout_rms, 0.01)
  ## the fit's characteristic function against the target's
  withr::local_seed(2L)
  Tm <- matrix(rnorm(40L), 20L, 2L)
  expect_lt(max(Mod(fit$cf(Tm) - cf(Tm))), 0.02)
  ## Kolmogorov distance of the first coordinate
  first <- fit$project(c(1, 0))
  expect_equal(gmm_dim(first), 1L)
  d <- gmm_cf_distance(first, function(s) cf(cbind(s[, 1L], 0)))
  expect_lt(d$distance, 0.01)
  ## a source with a coefficient below 1e-3 of the largest is replaced by
  ## one normal, which keeps the mean of the projection
  a <- as.vector(c(1, 5e-4) %*% solve(fit$mixing))
  small <- fit$project(a)
  expect_length(small@weights, 3L)
  expect_length(fit$project(as.vector(c(1, 2e-3) %*% solve(fit$mixing)))@weights, 9L)
  expect_equal(sum(small@weights * unlist(small@means)),
               sum(a * as.vector(A %*% c(3, 3))), tolerance = 0.01)
  ## draws have the target's mean, A %*% c(3, 3)
  X <- withr::with_seed(3L, fit$draw(20000L))
  expect_equal(dim(X), c(20000L, 2L))
  expect_equal(colMeans(X), as.vector(A %*% c(3, 3)), tolerance = 0.03)
})

test_that("the pairwise dependence is 0 at the true directions and grows with a rotation", {
  A <- matrix(c(0.7, 0.3, 0, -0.7), 2L)
  phi <- .cf_checked(gamma_mix_cf(A), 2L)
  covariance <- 3 * A %*% t(A) + 0.25 * diag(2L)
  dependence <- function(angle) {
    R <- matrix(c(cos(angle), sin(angle), -sin(angle), cos(angle)), 2L)
    B <- R %*% solve(A)
    B <- B / sqrt(diag(B %*% covariance %*% t(B)))
    noise <- B %*% covariance %*% t(B)
    diag(noise) <- 0
    .cf_pair_dependence(phi, B, noise)
  }
  expect_lt(dependence(0), 1e-10)
  expect_gt(dependence(0.1), 0.01)
  expect_gt(dependence(0.3), dependence(0.1))
})

test_that("a target that is not a mix of independent sources is flagged by its dependence", {
  skip_on_cran()
  cf <- function(t) {
    0.5 * exp(1i * (t %*% c(2, 1))[, 1L] - 0.5 * rowSums(t^2)) +
      0.5 * exp(-1i * (t %*% c(2, -1))[, 1L] -
        0.5 * (t[, 1L]^2 + 4 * t[, 2L]^2 + 1.6 * t[, 1L] * t[, 2L]))
  }
  expect_warning(
    fit <- from_cf_sources(cf, n_dim = 2L, N = 3L, n_freq = 1024L, seed = 1L),
    class = "proxymix_cf_not_sources"
  )
  expect_gt(fit$diagnostics$dependence, 0.01)
})

test_that("from_cf_sources errors on a Gaussian target and on a target without variance", {
  ## Gaussian: no direction can be read
  S <- shape_2d()
  expect_error(
    from_cf_sources(function(t) exp(-0.5 * rowSums((t %*% S) * t)),
                    n_dim = 2L, N = 2L, n_freq = 1024L, seed = 1L),
    class = "proxymix_cf_not_sources"
  )
  ## a Cauchy has no variance
  expect_error(
    suppressWarnings(from_cf_sources(
      function(t) exp(-sqrt(rowSums(t^2))),
      n_dim = 2L, N = 2L, n_freq = 1024L, seed = 1L
    )),
    "positive definite covariance"
  )
})

test_that("from_cf_sources validates its arguments", {
  cf <- gamma_mix_cf(diag(2))
  expect_error(from_cf_sources("nope", 2L, 2L), "must be a function")
  expect_error(from_cf_sources(cf, 1L, 2L), "at least 2")
  expect_error(from_cf_sources(cf, 2L, 0L), "positive integer")
})
