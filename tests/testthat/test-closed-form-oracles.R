## Independent closed-form / quadrature oracles for quantities that were
## previously graded only on self-consistency (positivity, finiteness).
## The reference formulas below are hand-written in this file so they do
## not share source with the package implementations.

## Hand-written Gaussian-to-Gaussian KL (test-local reference).
.ref_kl_gauss <- function(m0, S0, m1, S1) {
  p <- length(m0)
  S1i <- solve(S1)
  0.5 * (sum(diag(S1i %*% S0)) +
           drop(t(m1 - m0) %*% S1i %*% (m1 - m0)) - p +
           determinant(S1, logarithm = TRUE)$modulus -
           determinant(S0, logarithm = TRUE)$modulus)
}

## Hand-written squared Hellinger distance between two Gaussians.
.ref_hellinger2_gauss <- function(m0, S0, m1, S1) {
  Sm <- (S0 + S1) / 2
  bc <- exp(0.25 * determinant(S0, logarithm = TRUE)$modulus +
              0.25 * determinant(S1, logarithm = TRUE)$modulus -
              0.5 * determinant(Sm, logarithm = TRUE)$modulus -
              0.125 * drop(t(m0 - m1) %*% solve(Sm) %*% (m0 - m1)))
  1 - as.numeric(bc)
}

test_that("gmm_kld at K = 1 matches the closed-form Gaussian KL", {
  m0 <- c(0.3, -0.7); S0 <- matrix(c(1.2, 0.3, 0.3, 0.8), 2)
  m1 <- c(-0.4, 0.5); S1 <- matrix(c(0.9, -0.2, -0.2, 1.5), 2)
  p <- gmm(weights = 1, means = list(m0), covariances = list(S0))
  q <- gmm(weights = 1, means = list(m1), covariances = list(S1))
  truth <- as.numeric(.ref_kl_gauss(m0, S0, m1, S1))
  est <- withr::with_seed(5, gmm_kld(p, q, n_mc = 40000L))
  ## Monte Carlo estimate within a generous multiple of its own SE.
  expect_lt(abs(est$mc - truth), 5 * est$mc_se + 1e-6)
  ## At K = 1 the Hershey-Olsen variational term is exact.
  expect_equal(est$variational, truth, tolerance = 1e-10)
})

test_that("gmm_kld at K = 1 (one dimension) matches numerical quadrature", {
  p <- gmm(weights = 1, means = list(0.5), covariances = list(matrix(1.4)))
  q <- gmm(weights = 1, means = list(-0.2), covariances = list(matrix(0.7)))
  truth <- stats::integrate(function(x) {
    fp <- dgmm(matrix(x, ncol = 1L), p)
    fq <- dgmm(matrix(x, ncol = 1L), q)
    fp * (log(fp) - log(fq))
  }, -12, 12, rel.tol = 1e-10)$value
  est <- withr::with_seed(6, gmm_kld(p, q, n_mc = 40000L))
  expect_lt(abs(est$mc - truth), 5 * est$mc_se + 1e-6)
  expect_equal(est$variational, truth, tolerance = 1e-8)
})

test_that("hellinger_mc matches the closed-form Gaussian Hellinger", {
  ## A normalised Gaussian target and a deliberately offset single-Gaussian
  ## proxy: the true squared Hellinger distance is closed form.
  m_t <- c(0, 0); S_t <- diag(2)
  tgt <- gmm_target(
    n_dim = 2L,
    log_density = function(x) {
      if (is.null(dim(x))) x <- matrix(x, ncol = 2L)
      -0.5 * rowSums(x^2) - log(2 * pi)
    },
    normalised = TRUE, name = "std_normal"
  )
  fit <- fit_kld_em(tgt, N = 1L, is_size = 4000L, max_iter = 60L, seed = 2L,
                    proposal = is_mvt(2L, mean = c(0, 0),
                                      sigma = 4 * diag(2), df = 5))
  truth <- .ref_hellinger2_gauss(fit@means[[1L]], fit@covariances[[1L]],
                                 m_t, S_t)
  est <- hellinger_mc(fit, seed = 3L)
  expect_lt(abs(est$h2 - truth), 5 * est$se + 5e-3)
})

test_that("Renyi-2 entropy is affine-equivariant: H(AX + b) = H(X) + log|det A|", {
  g <- gmm(weights = c(0.3, 0.7),
           means = list(c(-1, 0.5), c(1.2, -0.3)),
           covariances = list(matrix(c(1, 0.4, 0.4, 1.1), 2),
                              matrix(c(0.8, -0.1, -0.1, 0.6), 2)))
  A <- matrix(c(2, 0.5, -0.3, 1.4), 2)
  ga <- gmm_affine(g, A, b = c(3, -2), ridge_eps = 0)
  expect_equal(
    gmm_entropy(ga, order = "renyi2"),
    gmm_entropy(g, order = "renyi2") + log(abs(det(A))),
    tolerance = 1e-8
  )
})

test_that("marginalise is consistent with integrating the joint (tower check)", {
  g <- gmm(weights = c(0.4, 0.6),
           means = list(c(-1, 1), c(2, -0.5)),
           covariances = list(matrix(c(1, 0.5, 0.5, 1.2), 2),
                              matrix(c(0.7, -0.2, -0.2, 0.9), 2)))
  marg <- gmm_marginalise(g, keep = 1L)
  for (a in c(-1.5, 0, 1.8)) {
    num <- stats::integrate(function(x2) {
      dgmm(cbind(rep(a, length(x2)), x2), g)
    }, -15, 15, rel.tol = 1e-9)$value
    expect_equal(dgmm(matrix(a), marg), num, tolerance = 1e-7)
  }
})

test_that("observe with a near-noiseless selection row agrees with conditionalise", {
  g <- gmm(weights = c(0.4, 0.6),
           means = list(c(-1, 1), c(2, -0.5)),
           covariances = list(matrix(c(1, 0.5, 0.5, 1.2), 2),
                              matrix(c(0.7, -0.2, -0.2, 0.9), 2)))
  y2 <- 0.4
  o <- gmm_observe(g, A = matrix(c(0, 1), nrow = 1L), y = y2,
                   noise_cov = matrix(1e-10), ridge_eps = 0)
  om <- gmm_marginalise(o, keep = 1L)
  cd <- gmm_conditionalise(g, given = c(NA, y2))
  expect_equal(om@weights, cd@weights, tolerance = 1e-5)
  expect_equal(unlist(om@means), unlist(cd@means), tolerance = 1e-4)
  expect_equal(unlist(om@covariances), unlist(cd@covariances),
               tolerance = 1e-4)
})

test_that("operator results still integrate to one (mass preservation)", {
  g <- gmm(weights = c(0.5, 0.5),
           means = list(c(-1, 0), c(1, 0.5)),
           covariances = list(diag(2), 0.5 * diag(2)))
  cond <- gmm_conditionalise(g, given = c(NA, 0.3))
  mass_c <- stats::integrate(function(x) dgmm(matrix(x), cond),
                             -15, 15, rel.tol = 1e-9)$value
  expect_equal(mass_c, 1, tolerance = 1e-7)
  aff <- gmm_affine(gmm_marginalise(g, keep = 1L), matrix(2), b = 1,
                    ridge_eps = 0)
  mass_a <- stats::integrate(function(x) dgmm(matrix(x), aff),
                             -25, 30, rel.tol = 1e-9)$value
  expect_equal(mass_a, 1, tolerance = 1e-7)
})

test_that("rgmm draws are distributed as dgmm claims (KS cross-check)", {
  g <- gmm(weights = c(0.35, 0.65), means = list(-2, 1.5),
           covariances = list(matrix(0.8), matrix(1.3)))
  x <- withr::with_seed(17, rgmm(3000L, g))
  cdf <- function(q) {
    0.35 * stats::pnorm(q, -2, sqrt(0.8)) +
      0.65 * stats::pnorm(q, 1.5, sqrt(1.3))
  }
  ks <- suppressWarnings(stats::ks.test(as.numeric(x), cdf))
  expect_gt(ks$p.value, 0.01)
})

test_that("gmm_impute on complete data returns finite completions (idempotence guard)", {
  x <- withr::with_seed(23, matrix(stats::rnorm(120), ncol = 2,
                                   dimnames = list(NULL, c("a", "b"))))
  imp <- suppressWarnings(gmm_impute(x, m = 2L, seed = 1L))
  for (comp in imp@completions) {
    expect_true(all(is.finite(comp)))
  }
})
