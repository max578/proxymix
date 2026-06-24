# Conditional-independence (Gaussian graphical model) structure of a mixture.

test_that("a tridiagonal precision recovers the chain skeleton (exact)", {
  ## A single Gaussian whose precision is tridiagonal is a Markov chain
  ## x1 - x2 - x3 - x4: partial correlation is nonzero only between neighbours.
  omega <- diag(4)
  for (i in 1:3) {
    omega[i, i + 1L] <- omega[i + 1L, i] <- -0.5
  }
  covx <- solve(omega)
  g <- gmm(weights = 1, means = list(rep(0, 4)), covariances = list(covx))
  adj <- gmm_independence_graph(g)
  truth <- matrix(0L, 4, 4)
  for (i in 1:3) truth[i, i + 1L] <- truth[i + 1L, i] <- 1L
  bare <- adj
  attr(bare, "pcor") <- NULL
  dimnames(bare) <- NULL
  expect_equal(bare, truth)
  expect_true(is.matrix(attr(adj, "pcor")))
})

test_that("independent coordinates give the empty graph", {
  g <- gmm(weights = 1, means = list(c(0, 0, 0)), covariances = list(diag(3)))
  adj <- gmm_independence_graph(g)
  expect_equal(sum(adj), 0L)
})

test_that("a dense precision gives the complete graph", {
  s <- matrix(0.4, 3, 3)
  diag(s) <- 1
  g <- gmm(weights = 1, means = list(c(0, 0, 0)), covariances = list(s))
  adj <- gmm_independence_graph(g)
  expect_equal(sum(adj) / 2L, 3L)            # all 3 off-diagonal pairs present
})

test_that("threshold controls edge inclusion", {
  s <- matrix(c(1, 0.1, 0.1, 1), 2, 2)
  g <- gmm(weights = 1, means = list(c(0, 0)), covariances = list(s))
  expect_equal(sum(gmm_independence_graph(g, threshold = 0.05)), 2L)
  expect_equal(sum(gmm_independence_graph(g, threshold = 0.5)), 0L)
})

test_that("input validation", {
  g <- gmm(weights = 1, means = list(c(0, 0)), covariances = list(diag(2)))
  expect_error(gmm_independence_graph(42), "must be a")
  expect_error(gmm_independence_graph(g, threshold = -1), "non-negative")
  g1 <- gmm(weights = 1, means = list(0), covariances = list(matrix(1, 1, 1)))
  expect_error(gmm_independence_graph(g1), "at least")
})

test_that("regime (iii): structure recovered from an evaluable-only energy model", {
  skip_on_cran()
  energy <- function(X) {
    X <- matrix(X, ncol = 3)
    rowSums((X^2 - 1)^2) - 0.7 * (X[, 1] * X[, 2] + X[, 2] * X[, 3])
  }
  target <- gmm_target(n_dim = 3L, log_density = function(X) -energy(X))
  g <- fit_kld_em(target, N = 8L, proposal = is_uniform(3L, -3, 3),
                  is_size = 8000L, anneal = TRUE, seed = 1L, support_warn = FALSE)
  adj <- gmm_independence_graph(g)
  ## the chain x1 - x2 - x3: neighbours coupled, x1 - x3 conditionally independent
  expect_equal(adj["x1", "x2"], 1L)
  expect_equal(adj["x2", "x3"], 1L)
  expect_equal(adj["x1", "x3"], 0L)
})
