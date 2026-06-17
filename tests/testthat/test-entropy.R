## Independent oracles for the entropy layer: analytic Gaussian forms (to 1e-10)
## and numerical quadrature via stats::integrate on 1-D mixtures (to 1e-6).

.dens1 <- function(x, g) dgmm(matrix(x, ncol = 1L), g)

# --- SE1: closed-form Renyi-2 entropy ---------------------------------------
test_that("renyi2 entropy equals the analytic Gaussian form at K = 1 (1e-10)", {
  scales <- c(0.5, 1, 2)
  for (p in c(1L, 2L, 3L)) {
    S <- diag(scales[seq_len(p)], nrow = p)
    g <- gmm(weights = 1, means = list(rep(0, p)), covariances = list(S))
    ## H_2(N(mu, S)) = (p/2) log(4 pi) + 1/2 log|S|.
    analytic <- 0.5 * (p * log(4 * pi) +
                         as.numeric(determinant(S, logarithm = TRUE)$modulus))
    expect_equal(gmm_entropy(g, order = "renyi2"), analytic, tolerance = 1e-10)
  }
})

test_that("renyi2 entropy matches numerical quadrature on a 1-D mixture (1e-6)", {
  g <- gmm(weights = c(0.3, 0.7),
           means = list(-1.5, 2),
           covariances = list(matrix(0.6), matrix(1.3)))
  v <- stats::integrate(function(x) .dens1(x, g)^2, -Inf, Inf,
                        rel.tol = 1e-10)$value
  expect_equal(gmm_entropy(g, order = "renyi2"), -log(v), tolerance = 1e-6)
})

# --- SE2: Shannon Monte Carlo + analytic upper bound ------------------------
test_that("shannon entropy matches the Gaussian form at K = 1 (MC error)", {
  S <- matrix(c(1, 0.3, 0.3, 1), 2, 2)
  g <- gmm(weights = 1, means = list(c(0, 0)), covariances = list(S))
  analytic <- 0.5 * log((2 * pi * exp(1))^2 * det(S))
  out <- gmm_entropy(g, order = "shannon", n_mc = 20000L, seed = 1L)
  expect_lt(abs(out$mc - analytic), 0.05)
})

test_that("the analytic upper bound brackets the shannon MC estimate", {
  g <- gmm(weights = c(0.5, 0.5),
           means = list(c(-2, 0), c(2, 0)),
           covariances = list(diag(2), diag(2)))
  out <- gmm_entropy(g, order = "shannon", n_mc = 20000L, seed = 1L)
  expect_lte(out$mc, out$upper_bound + 3 * out$mc_se)
})

# --- SE3: Cauchy-Schwarz divergence -----------------------------------------
test_that("cs divergence is symmetric, non-negative, and zero at p = q (1e-10)", {
  p <- gmm(weights = c(0.4, 0.6), means = list(c(-1, 0), c(1, 1)),
           covariances = list(diag(2), diag(2) * 1.5))
  q <- gmm(weights = 1, means = list(c(0, 0)), covariances = list(diag(2) * 2))
  expect_equal(gmm_divergence(p, q), gmm_divergence(q, p), tolerance = 1e-10)
  expect_gte(gmm_divergence(p, q), 0)
  expect_equal(gmm_divergence(p, p), 0, tolerance = 1e-10)
})

test_that("cs divergence matches numerical quadrature on 1-D mixtures (1e-6)", {
  p <- gmm(weights = c(0.5, 0.5), means = list(-1, 1.5),
           covariances = list(matrix(0.5), matrix(0.8)))
  q <- gmm(weights = 1, means = list(0.2), covariances = list(matrix(1.4)))
  vpq <- stats::integrate(function(x) .dens1(x, p) * .dens1(x, q), -Inf, Inf,
                          rel.tol = 1e-10)$value
  vpp <- stats::integrate(function(x) .dens1(x, p)^2, -Inf, Inf,
                          rel.tol = 1e-10)$value
  vqq <- stats::integrate(function(x) .dens1(x, q)^2, -Inf, Inf,
                          rel.tol = 1e-10)$value
  oracle <- 0.5 * log(vpp) + 0.5 * log(vqq) - log(vpq)
  expect_equal(gmm_divergence(p, q, type = "cs"), oracle, tolerance = 1e-6)
})

test_that("gmm_divergence type = 'kl' delegates to gmm_kld", {
  p <- gmm(weights = 1, means = list(c(0, 0)), covariances = list(diag(2)))
  q <- gmm(weights = 1, means = list(c(0.5, 0)), covariances = list(diag(2)))
  a <- gmm_divergence(p, q, type = "kl", n_mc = 2000L)
  expect_true(is.list(a) && all(c("mc", "mc_se") %in% names(a)))
})

# --- input validation -------------------------------------------------------
test_that("entropy and divergence reject bad input", {
  g1 <- gmm(weights = 1, means = list(0), covariances = list(matrix(1)))
  g2 <- gmm(weights = 1, means = list(c(0, 0)), covariances = list(diag(2)))
  expect_error(gmm_entropy(42), "gmm")
  expect_error(gmm_divergence(g1, 42), "gmm")
  expect_error(gmm_divergence(g1, g2), "dimension")
})
