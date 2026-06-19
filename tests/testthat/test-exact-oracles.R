## Cheap exact-oracle regression guards (Independent Oracle Principle).
##
## These promote the cheapest, fastest checks from the workspace end-to-end
## validation suite (tiers V1 and C1.1-C1.3) into the package's own `testthat`
## suite, so the closed-form exactness identities and the unifying-primitive
## equalities are re-checked on every `R CMD check` against oracles the package
## did not produce (analytic formulas, grid quadrature, `stats::lm` /
## `stats::prcomp` / a hand kernel smoother).
##
## The remaining V1 / C1 checks are already in the suite and are not duplicated
## here: V1.1 affine moments in test-operator-calculus.R (A1); V1.5-V1.7
## Renyi-2 / Cauchy-Schwarz / mutual-information in test-entropy.R; C1.1 the
## K = 1 == OLS reduction in test-gmr-k1-lm.R. The heavier cross-implementation
## oracles (grf, DoubleML, particle filter) stay in the workspace suite.

## Mixture moments by direct summation (independent of any package accessor).
.mix_moments <- function(g) {
  w <- g@weights
  mu <- Reduce(`+`, Map(function(a, b) a * b, w, g@means))
  m2 <- Reduce(`+`, Map(function(a, b, cc) a * (cc + tcrossprod(b)),
                        w, g@means, g@covariances))
  list(mean = as.numeric(mu), cov = m2 - tcrossprod(mu))
}

## Grid (quadrature) Bayesian posterior for a 1-D mixture prior under a scalar
## linear-Gaussian observation y ~ N(C x, R).
.grid_posterior_1d <- function(prior_g, C, y, R, lo, hi, n = 8000L) {
  xs <- seq(lo, hi, length.out = n)
  dx <- xs[2L] - xs[1L]
  pr <- dgmm(matrix(xs, ncol = 1L), prior_g)
  lik <- stats::dnorm(y, mean = C * xs, sd = sqrt(R))
  post <- pr * lik
  z <- sum(post) * dx
  post <- post / z
  list(mean = sum(xs * post) * dx, logevidence = log(z))
}

## Total mass of a 1-D mixture by grid integration of its density.
.gmm_mass_1d <- function(g, lo, hi, n = 8000L) {
  xs <- seq(lo, hi, length.out = n)
  sum(dgmm(matrix(xs, ncol = 1L), g)) * (xs[2L] - xs[1L])
}

# --- V1.2: conditioning equals the analytic Gaussian Schur complement -------

test_that("single-Gaussian conditioning equals the analytic Schur complement", {
  mu <- c(0.3, -0.4)
  s <- matrix(c(1, 0.5, 0.5, 1.2), 2L)
  g <- gmm(weights = 1, means = list(mu), covariances = list(s))
  yobs <- 0.7
  cond <- gmm_conditionalise(g, given = c(NA, yobs))

  cond_mean <- mu[1L] + s[1L, 2L] / s[2L, 2L] * (yobs - mu[2L])
  cond_var <- s[1L, 1L] - s[1L, 2L]^2 / s[2L, 2L]

  expect_equal(cond@means[[1L]][1L], cond_mean, tolerance = 1e-8)
  expect_equal(cond@covariances[[1L]][1L, 1L], cond_var, tolerance = 1e-8)
})

# --- V1.3: the K = 2 Bayesian update equals a grid quadrature posterior -----

test_that("gmm_observe at K = 2 equals an independent grid posterior", {
  prior <- gmm(weights = c(0.5, 0.5), means = list(-1.5, 1.5),
               covariances = list(matrix(0.5), matrix(0.5)))
  y <- 0.8
  post <- gmm_observe(prior, A = matrix(1, 1L, 1L), y = y,
                      noise_cov = matrix(0.25, 1L, 1L), ridge_eps = 0)
  pm <- .mix_moments(post)
  grid <- .grid_posterior_1d(prior, C = 1, y = y, R = 0.25, lo = -8, hi = 8)

  expect_equal(pm$mean[1L], grid$mean, tolerance = 1e-3)
  expect_equal(post@metadata$log_marginal_evidence, grid$logevidence,
               tolerance = 1e-3)
})

# --- V1.4: the mixture density integrates to one ----------------------------

test_that("dgmm integrates to one over its effective support", {
  g <- gmm(weights = c(0.3, 0.7), means = list(-2, 2),
           covariances = list(matrix(0.6), matrix(0.9)))
  mass <- .gmm_mass_1d(g, -12, 12)
  expect_equal(mass, 1, tolerance = 1e-4)
})

# --- C1.2: the K = 1 covariance eigenvector equals the principal component ---

test_that("K = 1 covariance eigenvector equals prcomp's first PC", {
  x <- withr::with_seed(2L,
    mvnfast::rmvn(600L, c(0, 0), matrix(c(2, 1.2, 1.2, 1.5), 2L)))
  fit <- fit_proxymix(gmm_target_from_samples(x), N = 1L, regime = "moment")

  ev <- eigen(fit@covariances[[1L]])$vectors[, 1L]
  pc <- stats::prcomp(x)$rotation[, 1L]
  ## eigenvectors are defined up to sign; compare the absolute cosine.
  expect_equal(abs(sum(ev * pc)), 1, tolerance = 1e-6)
})

# --- C1.3: the K = n conditional mean equals Nadaraya-Watson ----------------

test_that("K = n conditional mean equals the Nadaraya-Watson estimator", {
  set.seed(3L)
  n <- 60L
  xx <- sort(stats::rnorm(n))
  yy <- sin(xx) + stats::rnorm(n, 0, 0.2)
  h <- 0.4
  s2 <- 0.04
  ## one Gaussian component per datum: a kernel density estimate as a mixture.
  gn <- gmm(weights = rep(1 / n, n),
            means = lapply(seq_len(n), function(i) c(yy[i], xx[i])),
            covariances = rep(list(diag(c(s2, h^2))), n))
  x0 <- 0.3

  gmr <- .mix_moments(gmm_conditionalise(gn, c(NA, x0)))$mean[1L]
  nw <- sum(stats::dnorm(x0 - xx, 0, h) * yy) /
    sum(stats::dnorm(x0 - xx, 0, h))

  expect_equal(gmr, nw, tolerance = 1e-10)
})
