## from_cf() fits a mixture to a characteristic function; gmm_cf() evaluates
## a mixture's characteristic function; gmm_cf_distance() is the Kolmogorov
## distance along a direction. Every expected value below is exact or comes
## from a density integrated numerically.

two_cluster <- function() {
  gmm(
    weights = c(0.4, 0.6),
    means = list(c(-2, 0), c(2, 1)),
    covariances = list(diag(2), matrix(c(1, 0.5, 0.5, 2), 2L))
  )
}

# ---------------------------------------------------------------------------
# gmm_cf
# ---------------------------------------------------------------------------

test_that("gmm_cf matches the normal characteristic function", {
  mu <- c(1, -2)
  S <- matrix(c(2, 0.6, 0.6, 1), 2L)
  g <- gmm(weights = 1, means = list(mu), covariances = list(S))
  t <- rbind(c(0, 0), c(0.3, -0.7), c(-1.2, 0.4))
  expected <- exp(1i * as.vector(t %*% mu) - 0.5 * rowSums((t %*% S) * t))
  expect_equal(gmm_cf(g, t), expected)
  expect_equal(gmm_cf(g, -t), Conj(expected))
})

test_that("gmm_cf agrees with the sample mean of exp(i t'X)", {
  g <- two_cluster()
  x <- withr::with_seed(11L, rgmm(2e5L, g))
  t <- c(0.4, -0.3)
  empirical <- mean(exp(1i * as.vector(x %*% t)))
  expect_lt(Mod(gmm_cf(g, t) - empirical), 0.01)
})

test_that("gmm_cf reads a vector by the mixture's dimension", {
  g1 <- gmm(weights = 1, means = list(0), covariances = list(matrix(1)))
  expect_equal(gmm_cf(g1, c(0, 1, 2)), as.complex(exp(-0.5 * c(0, 1, 4))))
  expect_length(gmm_cf(two_cluster(), c(0.1, 0.2)), 1L)
  expect_error(gmm_cf(two_cluster(), matrix(0, 2L, 3L)), "2 columns")
  expect_error(gmm_cf(list(), 0), class = "rlang_error")
})

# ---------------------------------------------------------------------------
# loss derivatives
# ---------------------------------------------------------------------------

test_that("the loss gradient matches central differences", {
  K <- 3L
  p <- 2L
  withr::local_seed(3L)
  W <- matrix(rnorm(400L), 200L, p)
  ph <- exp(-0.5 * rowSums(W^2) + 1i * W[, 1L])
  par <- rnorm(K - 1L + K * p + K * 3L, sd = 0.3)
  value_at <- function(x) .cf_loss_grad(x, W, ph, K, p)$value
  numeric_grad <- vapply(seq_along(par), function(j) {
    e <- replace(numeric(length(par)), j, 1e-6)
    (value_at(par + e) - value_at(par - e)) / 2e-6
  }, numeric(1L))
  analytic <- .cf_loss_grad(par, W, ph, K, p)
  expect_equal(analytic$gradient, numeric_grad, tolerance = 1e-7)

  ## the same loss written as a sum of squares
  rj <- .cf_resid_jac(par, W, ph, K, p)
  expect_equal(sum(rj$r^2), analytic$value)
  expect_equal(as.vector(2 * crossprod(rj$J, rj$r)), analytic$gradient)
})

# ---------------------------------------------------------------------------
# from_cf: recovery
# ---------------------------------------------------------------------------

test_that("one component returns the target's mean and covariance", {
  mu <- c(1, -2)
  S <- matrix(c(2, 0.6, 0.6, 1), 2L)
  cf <- function(t) exp(1i * as.vector(t %*% mu) - 0.5 * rowSums((t %*% S) * t))
  fit <- from_cf(cf, n_dim = 2L, N = 1L, n_freq = 1024L, seed = 1L)
  expect_s7_class(fit, gmm_fit)
  expect_equal(fit@means[[1L]], mu, tolerance = 1e-6)
  expect_equal(fit@covariances[[1L]], S, tolerance = 1e-6)
  expect_equal(fit@diagnostics$cf_mean, mu, tolerance = 1e-6)
  expect_equal(fit@diagnostics$tail_index, c(2, 2), tolerance = 1e-3)
})

test_that("a two-component mixture is recovered", {
  truth <- two_cluster()
  cf <- function(t) gmm_cf(truth, t)
  fit <- from_cf(cf, n_dim = 2L, N = 2L, n_freq = 2048L, seed = 1L)
  expect_true(fit@converged)
  expect_lt(fit@diagnostics$loss, 1e-20)
  expect_lt(fit@diagnostics$heldout_max, 1e-10)
  ## components come back in canonical order, so match by weight
  ord <- order(fit@weights)
  expect_equal(fit@weights[ord], c(0.4, 0.6), tolerance = 1e-8)
  expect_equal(fit@means[ord], truth@means, tolerance = 1e-8)
  expect_equal(fit@covariances[ord], truth@covariances, tolerance = 1e-8)
  expect_lt(gmm_cf_distance(fit, cf, c(1, 1))$distance, 1e-10)
})

test_that("a skewed one-dimensional target is matched in distribution", {
  skip_on_cran()
  ## normal inverse Gaussian, alpha = 2, beta = 1, delta = 1, mu = 0
  al <- 2
  be <- 1
  ga <- sqrt(al^2 - be^2)
  cf <- function(t) exp(ga - sqrt(al^2 - (be + 1i * t[, 1L])^2))
  dens <- function(x) {
    r <- sqrt(1 + x^2)
    al * besselK(al * r, 1, expon.scaled = TRUE) / (pi * r) *
      exp(ga + be * x - al * r)
  }
  fit <- from_cf(cf, n_dim = 1L, N = 5L, n_freq = 2048L, seed = 1L)
  cert <- gmm_cf_distance(fit, cf)
  expect_lt(cert$distance, 1e-3)

  ## the reported difference equals the one computed from the density
  at <- cert$x[c(20L, 60L, 100L, 140L, 180L)]
  cdf <- vapply(at, function(x) {
    stats::integrate(dens, -Inf, x, rel.tol = 1e-10)$value
  }, numeric(1L))
  reported <- cert$difference[c(20L, 60L, 100L, 140L, 180L)]
  expect_lt(max(abs(reported - (cdf - pgmm(at, fit)))), 1e-7)
})

# ---------------------------------------------------------------------------
# from_cf: returned object
# ---------------------------------------------------------------------------

test_that("the fit carries its diagnostics and works with the operators", {
  truth <- two_cluster()
  cf <- function(t) gmm_cf(truth, t)
  fit <- from_cf(cf, n_dim = 2L, N = 2L, n_freq = 1024L, seed = 1L)
  expect_identical(fit@regime, "cf")
  expect_named(
    fit@diagnostics,
    c("loss", "loss_by_start", "loss_path", "heldout_rms", "heldout_max",
      "relocations", "tail_index", "cf_mean", "cf_cov", "n_freq",
      "bandwidths", "bracket", "bracket_directions")
  )
  expect_length(fit@diagnostics$loss_path, 2L)
  q <- gmm_fit_quality(fit)
  expect_identical(q$regime, "cf")
  expect_identical(q$cf_heldout_rms, fit@diagnostics$heldout_rms)
  total <- gmm_affine(fit, matrix(c(1, 1), 1L))
  expect_equal(gmm_dim(total), 1L)
  expect_identical(gmm_fit_quality(total)$regime, "cf")
})

test_that("extra starts are reported and never worsen the fit", {
  cf <- function(t) (1 - 1i * t[, 1L])^(-3) * exp(-0.125 * t[, 1L]^2)
  one <- from_cf(cf, n_dim = 1L, N = 2L, n_freq = 1024L, seed = 1L)
  two <- from_cf(cf, n_dim = 1L, N = 2L, n_freq = 1024L, n_start = 2L,
                 seed = 1L)
  expect_length(two@diagnostics$loss_by_start, 2L)
  expect_equal(two@diagnostics$loss, min(two@diagnostics$loss_by_start))
  expect_lte(two@diagnostics$loss, one@diagnostics$loss)
})

test_that("a fit can be started from a given mixture", {
  cf <- function(t) (1 - 1i * t[, 1L])^(-3) * exp(-0.125 * t[, 1L]^2)
  first <- from_cf(cf, n_dim = 1L, N = 3L, n_freq = 1024L, seed = 1L)
  again <- from_cf(cf, n_dim = 1L, N = 3L, n_freq = 1024L, seed = 1L,
                   init = first)
  expect_null(again@diagnostics$loss_path)
  expect_lte(again@diagnostics$loss, first@diagnostics$loss * (1 + 1e-6))
  ## a two-cluster target from a rough two-component start
  truth <- two_cluster()
  start <- gmm(weights = c(0.5, 0.5), means = list(c(-1, 0), c(1, 0)),
               covariances = list(diag(2), diag(2)))
  fit <- from_cf(function(t) gmm_cf(truth, t), n_dim = 2L, N = 2L,
                 n_freq = 1024L, seed = 1L, init = start)
  expect_lt(fit@diagnostics$loss, 1e-20)
  expect_error(from_cf(cf, 1L, 2L, init = first), "must have")
  expect_error(from_cf(cf, 1L, 3L, init = list()), "init")
})

test_that("from_cf is reproducible and leaves the RNG state untouched", {
  cf <- function(t) (1 - 1i * t[, 1L])^(-3) * exp(-0.125 * t[, 1L]^2)
  a <- from_cf(cf, n_dim = 1L, N = 2L, n_freq = 1024L, seed = 7L)
  b <- from_cf(cf, n_dim = 1L, N = 2L, n_freq = 1024L, seed = 7L)
  expect_identical(a@means, b@means)
  expect_identical(a@diagnostics$heldout_rms, b@diagnostics$heldout_rms)

  set.seed(9L)
  expected <- runif(1L)
  set.seed(9L)
  invisible(from_cf(cf, n_dim = 1L, N = 2L, n_freq = 1024L, seed = 7L))
  expect_identical(runif(1L), expected)
})

# ---------------------------------------------------------------------------
# from_cf: heavy tails and invalid input
# ---------------------------------------------------------------------------

test_that("a target far from the origin is centred on its mean", {
  ## the principal phase at the differencing frequency wraps past pi near a
  ## mean of 130 for a unit normal
  for (mu in c(140, -3000)) {
    fit <- from_cf(function(t) exp(1i * mu * t[, 1L] - 0.5 * t[, 1L]^2),
                   1L, 1L, n_freq = 1024L, seed = 1L)
    expect_equal(fit@means[[1L]], mu, tolerance = 1e-8)
    expect_equal(fit@covariances[[1L]], matrix(1), tolerance = 1e-6)
    expect_lt(gmm_fit_quality(fit)$cf_distance, 1e-6)
  }
  m <- c(300, -250)
  S <- matrix(c(2, 0.6, 0.6, 1), 2L)
  fit <- from_cf(function(t) exp(1i * drop(t %*% m) - 0.5 * rowSums((t %*% S) * t)),
                 2L, 1L, n_freq = 1024L, seed = 1L)
  expect_equal(fit@means[[1L]], m, tolerance = 1e-8)
  expect_equal(fit@covariances[[1L]], S, tolerance = 1e-6)
})

test_that("a heavy-tailed target far from the origin is centred on its location", {
  fit <- suppressWarnings(
    from_cf(function(t) exp(1i * 500 * t[, 1L] - abs(t[, 1L])), 1L, 3L,
            n_freq = 2048L, seed = 1L),
    classes = "proxymix_cf_heavy_tail"
  )
  med <- stats::uniroot(function(x) pgmm(x, fit) - 0.5, c(400, 600))$root
  expect_equal(med, 500, tolerance = 1e-3)
  expect_lt(gmm_fit_quality(fit)$cf_distance, 0.02)
})

test_that("a Cauchy target raises the heavy-tail warning", {
  cf <- function(t) exp(-abs(t[, 1L]))
  expect_warning(
    fit <- from_cf(cf, n_dim = 1L, N = 3L, n_freq = 1024L, seed = 1L),
    class = "proxymix_cf_heavy_tail"
  )
  expect_equal(fit@diagnostics$tail_index, 1, tolerance = 1e-3)
  expect_null(fit@diagnostics$cf_mean)
  ## the centre is matched: the Cauchy quartiles are -1 and 1
  expect_equal(qgmm(c(0.25, 0.75), fit), c(-1, 1), tolerance = 0.05)
})

test_that("functions that are not characteristic functions are rejected", {
  expect_error(
    from_cf(function(t) 2 * exp(-0.5 * t[, 1L]^2), 1L, 1L),
    "equal 1 at zero", class = "proxymix_invalid_cf"
  )
  expect_error(
    from_cf(function(t) exp(-0.5 * t[, 1L]^2 + t[, 1L]^3 / 6), 1L, 1L),
    class = "proxymix_invalid_cf"
  )
  expect_error(
    from_cf(function(t) exp(-0.5 * t[, 1L]^2) * (1 + 0.5i * tanh(t[, 1L]^2)),
            1L, 1L),
    "conjugate", class = "proxymix_invalid_cf"
  )
  expect_error(
    from_cf(function(t) 1, 1L, 1L),
    "one complex", class = "proxymix_invalid_cf"
  )
  expect_error(
    from_cf(function(t) rep(NaN, nrow(t)), 1L, 1L),
    "not finite", class = "proxymix_invalid_cf"
  )
  ## a point mass at 1: the modulus never decays
  expect_error(
    from_cf(function(t) exp(1i * t[, 1L]), 1L, 1L),
    "does not fall", class = "proxymix_invalid_cf"
  )
})

test_that("from_cf validates its arguments", {
  cf <- function(t) exp(-0.5 * t[, 1L]^2)
  expect_error(from_cf("nope", 1L, 1L), "must be a function")
  expect_error(from_cf(cf, 0L, 1L), "n_dim")
  expect_error(from_cf(cf, 1L, 1.5), "N")
  expect_error(from_cf(cf, 1L, 1L, n_start = 0L), "n_start")
  expect_error(from_cf(cf, 1L, 1L, max_iter = NA), "max_iter")
  expect_error(from_cf(cf, 1L, 1L, bandwidths = c(2, 1)), "increasing")
  expect_error(from_cf(cf, 1L, 1L, n_freq = 100L), "at least 64")
})

# ---------------------------------------------------------------------------
# gmm_cf_distance
# ---------------------------------------------------------------------------

test_that("the distance between two shifted normals is exact", {
  cf <- function(t) exp(-0.5 * t[, 1L]^2)
  g <- gmm(weights = 1, means = list(0.1), covariances = list(matrix(1)))
  res <- gmm_cf_distance(g, cf)
  expect_equal(res$distance, 2 * pnorm(0.05) - 1, tolerance = 1e-4)
  expect_equal(res$at, 0.05, tolerance = 0.05)
  expect_lt(res$cutoff_gap, 1e-7)
  expect_length(res$difference, 201L)
})

test_that("the exact distance lies between the grid maximum and the bound", {
  cf <- function(t) exp(-0.5 * t[, 1L]^2)
  g <- gmm(weights = 1, means = list(0.1), covariances = list(matrix(1)))
  exact <- 2 * pnorm(0.05) - 1
  res <- gmm_cf_distance(g, cf)
  ## the inversion integral is accurate to about 1e-6
  expect_lte(res$distance, exact + 2e-6)
  expect_gte(res$upper_inside, exact)
  expect_gte(res$upper, res$upper_inside)
  expect_lt(res$upper - res$distance, 0.005)
  expect_lt(res$truncation, 1e-8)
  ## the grid is at equal steps of the mixture's distribution function
  expect_equal(pnorm(res$x, 0.1), seq(0.005, 0.995, length.out = 201L),
               tolerance = 1e-8)
  fine <- gmm_cf_distance(g, cf, n_grid = 2001L, probs = c(1e-4, 1 - 1e-4))
  expect_lt(fine$upper - exact, 5e-4)
})

test_that("the distance does not depend on the length of the direction", {
  truth <- two_cluster()
  cf <- function(t) gmm_cf(truth, t)
  g <- gmm(weights = 1, means = list(c(0.4, 0.6)),
           covariances = list(matrix(c(5, 1, 1, 1.5), 2L)))
  a <- gmm_cf_distance(g, cf, c(1, 2))$distance
  b <- gmm_cf_distance(g, cf, c(3, 6))$distance
  expect_gt(a, 0.01)
  expect_equal(a, b, tolerance = 1e-5)
})

test_that("gmm_cf_distance validates its arguments", {
  cf <- function(t) exp(-0.5 * t[, 1L]^2)
  g <- gmm(weights = 1, means = list(0), covariances = list(matrix(1)))
  expect_error(gmm_cf_distance(list(), cf), class = "rlang_error")
  expect_error(gmm_cf_distance(g, "nope"), "must be a function")
  expect_error(gmm_cf_distance(g, cf, direction = 0), "not all zero")
  expect_error(gmm_cf_distance(g, cf, direction = c(1, 1)), "not all zero")
  expect_error(gmm_cf_distance(g, cf, probs = c(0.9, 0.1)), "probs")
  expect_error(gmm_cf_distance(g, cf, n_grid = 1L), "n_grid")
})
