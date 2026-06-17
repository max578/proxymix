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

# --- SE4: Cauchy-Schwarz mutual information ---------------------------------
## Analytic CS mutual information of a standard correlated bivariate Gaussian
## (hand-derived from the Gaussian-product integral, independent of the code):
##   I_CS(rho) = -log 2 + 1/2 log(4 - rho^2) - 1/4 log(1 - rho^2).
.ics_bivariate <- function(rho) -log(2) + 0.5 * log(4 - rho^2) - 0.25 * log(1 - rho^2)

test_that("CS mutual information equals the analytic bivariate-Gaussian form (1e-10)", {
  for (rho in c(0, 0.3, 0.6, 0.9)) {
    s <- matrix(c(1, rho, rho, 1), 2, 2)
    g <- gmm(weights = 1, means = list(c(0, 0)), covariances = list(s))
    expect_equal(gmm_mutual_information(g, 1L, 2L), .ics_bivariate(rho),
                 tolerance = 1e-10)
  }
})

test_that("CS mutual information is non-negative, zero at rho = 0, monotone in |rho|", {
  ics <- function(rho) {
    s <- matrix(c(1, rho, rho, 1), 2, 2)
    gmm_mutual_information(
      gmm(weights = 1, means = list(c(0, 0)), covariances = list(s)), 1L, 2L)
  }
  expect_equal(ics(0), 0, tolerance = 1e-10)
  expect_gte(ics(0.5), 0)
  expect_lt(ics(0.3), ics(0.6))
  expect_lt(ics(0.6), ics(0.9))
})

test_that("CS mutual information is ~zero under planted independence", {
  ## coord 2 is N(0, 1) in BOTH components, so the joint factors as
  ## [bimodal in coord 1] x [N(0,1) in coord 2]: the blocks are independent.
  g <- gmm(weights = c(0.5, 0.5),
           means = list(c(-2, 0), c(2, 0)),
           covariances = list(diag(2), diag(2)))
  expect_equal(gmm_mutual_information(g, 1L, 2L), 0, tolerance = 1e-9)
})

test_that("gmm_mutual_information rejects bad input", {
  g <- gmm(weights = 1, means = list(c(0, 0)), covariances = list(diag(2)))
  expect_error(gmm_mutual_information(42, 1L, 2L), "gmm")
  expect_error(gmm_mutual_information(g, 1L, 1L), "disjoint")
  expect_error(gmm_mutual_information(g, 1L, 3L), "1.*2")
})

# --- SE5: conditional predictive entropy ------------------------------------
test_that("conditional entropy at K = 1 equals the analytic Gaussian conditional Renyi-2 (1e-10)", {
  ## joint N(0, s) over (Y = coord 1, X = coord 2); conditional var of Y is the
  ## scalar Schur complement and does not depend on the conditioning value.
  s <- matrix(c(2, 0.8, 0.8, 1), 2, 2)
  g <- gmm(weights = 1, means = list(c(0, 0)), covariances = list(s))
  cond_var <- s[1, 1] - s[1, 2]^2 / s[2, 2]
  analytic <- 0.5 * (log(4 * pi) + log(cond_var))
  ce <- gmm_conditional_entropy(g, given = c(NA, 0.7))
  expect_equal(ce, analytic, tolerance = 1e-10)
  expect_equal(gmm_conditional_entropy(g, given = c(NA, -3)), ce, tolerance = 1e-10)
})

test_that("a multimodal conditional has higher Renyi-2 entropy than its single component", {
  ## conditioning Y on X = 0 yields a well-separated equal-weight bimodal in Y.
  g <- gmm(weights = c(0.5, 0.5),
           means = list(c(-6, 0), c(6, 0)),
           covariances = list(diag(2), diag(2)))
  ce <- gmm_conditional_entropy(g, given = c(NA, 0))
  h_comp <- 0.5 * (log(4 * pi) + log(1))            # Renyi-2 of one N(., 1) component
  expect_gt(ce, h_comp)                              # multimodal => more uncertain
  expect_equal(ce, h_comp + log(2), tolerance = 0.05) # ~ +log K, well-separated equal weights
})

test_that("gmm_conditional_entropy vectorises over rows and rejects bad input", {
  s <- matrix(c(2, 0.8, 0.8, 1), 2, 2)
  g <- gmm(weights = 1, means = list(c(0, 0)), covariances = list(s))
  ce <- gmm_conditional_entropy(g, given = rbind(c(NA, 0), c(NA, 1), c(NA, 2)))
  expect_length(ce, 3L)
  expect_true(all(abs(diff(ce)) < 1e-10))           # constant across x for one Gaussian
  expect_error(gmm_conditional_entropy(42, c(NA, 0)), "gmm")
  expect_error(gmm_conditional_entropy(g, c(NA, 0, 0)), "one entry per coordinate")
})
