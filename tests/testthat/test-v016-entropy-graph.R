# gmm_independence_graph(): Fisher-z edge test when the sample size is known.

## Covariance of a 3-coordinate model whose precision has unit diagonal, so the
## partial correlations are exactly the negated off-diagonal precision entries.
.chain_cov <- function(r12, r23, r13) {
  omega <- diag(3)
  omega[1L, 2L] <- omega[2L, 1L] <- -r12
  omega[2L, 3L] <- omega[3L, 2L] <- -r23
  omega[1L, 3L] <- omega[3L, 1L] <- -r13
  solve(omega)
}

## A regime (ii) fit record carrying n data rows and a chosen overall covariance.
.fit_with_n <- function(covx, n, regime = "sample") {
  x <- withr::with_seed(1L, matrix(stats::rnorm(n * 3L), ncol = 3L))
  gmm_fit(
    weights = 1, means = list(c(0, 0, 0)), covariances = list(covx),
    target = gmm_target_from_samples(x), regime = regime,
    converged = TRUE, iterations = 1L
  )
}

test_that("a sample fit tests each partial correlation against its sampling error", {
  ## x1 - x3 partial correlation 0.069 at n = 500: z = 1.54, below 1.96.
  fit <- .fit_with_n(.chain_cov(0.4, 0.4, 0.069), n = 500L)
  adj <- gmm_independence_graph(fit)
  expect_equal(adj["x1", "x2"], 1L)
  expect_equal(adj["x2", "x3"], 1L)
  expect_equal(adj["x1", "x3"], 0L)
  expect_equal(attr(adj, "pcor")["x1", "x3"], 0.069, tolerance = 1e-6)
})

test_that("the same partial correlation is an edge once n is large enough", {
  ## At n = 2000, z = atanh(0.069) * sqrt(2000 - 4) = 3.09.
  fit <- .fit_with_n(.chain_cov(0.4, 0.4, 0.069), n = 2000L)
  expect_equal(gmm_independence_graph(fit)["x1", "x3"], 1L)
})

test_that("a moment-matched fit also takes n from its samples", {
  fit <- .fit_with_n(.chain_cov(0.4, 0.4, 0.069), n = 500L, regime = "moment")
  expect_equal(gmm_independence_graph(fit)["x1", "x3"], 0L)
})

test_that("the edge decision matches the Fisher-z critical partial correlation", {
  ## Two coordinates, n = 100: the critical |r| is tanh(qnorm(0.975) / sqrt(97)).
  r_crit <- tanh(stats::qnorm(0.975) / sqrt(97))
  one_pair <- function(r) {
    g <- gmm(weights = 1, means = list(c(0, 0)),
             covariances = list(matrix(c(1, r, r, 1), 2L, 2L)))
    gmm_independence_graph(g, n = 100L)[1L, 2L]
  }
  expect_equal(one_pair(r_crit - 0.005), 0L)
  expect_equal(one_pair(r_crit + 0.005), 1L)
  expect_equal(one_pair(-(r_crit + 0.005)), 1L)
})

test_that("alpha sets the per-edge test level", {
  g <- gmm(weights = 1, means = list(c(0, 0, 0)),
           covariances = list(.chain_cov(0.4, 0.4, 0.069)))
  ## z = 1.54 for x1 - x3 at n = 500: an edge at alpha = 0.2, none at 0.05.
  expect_equal(gmm_independence_graph(g, n = 500L, alpha = 0.2)[1L, 3L], 1L)
  expect_equal(gmm_independence_graph(g, n = 500L, alpha = 0.05)[1L, 3L], 0L)
})

test_that("an explicit threshold keeps the fixed-magnitude rule on a sample fit", {
  fit <- .fit_with_n(.chain_cov(0.4, 0.4, 0.069), n = 500L)
  adj <- gmm_independence_graph(fit, threshold = 0.05)
  expect_equal(sum(adj) / 2L, 3L)
})

test_that("without a sample size the rule is the 0.05 magnitude threshold", {
  g <- gmm(weights = 1, means = list(c(0, 0, 0)),
           covariances = list(.chain_cov(0.4, 0.4, 0.069)))
  expect_identical(gmm_independence_graph(g),
                   gmm_independence_graph(g, threshold = 0.05))
  expect_equal(sum(gmm_independence_graph(g)) / 2L, 3L)
})

test_that("the returned object keeps its shape under the Fisher-z rule", {
  fit <- .fit_with_n(.chain_cov(0.4, 0.4, 0.069), n = 500L)
  adj <- gmm_independence_graph(fit)
  expect_true(is.integer(adj))
  expect_equal(dim(adj), c(3L, 3L))
  expect_equal(dimnames(adj), list(c("x1", "x2", "x3"), c("x1", "x2", "x3")))
  expect_setequal(names(attributes(adj)), c("dim", "dimnames", "pcor"))
  expect_equal(adj, t(adj))
})

test_that("invalid n, alpha, or a conflicting threshold are rejected", {
  g <- gmm(weights = 1, means = list(c(0, 0, 0)), covariances = list(diag(3)))
  expect_error(gmm_independence_graph(g, n = 4L), "greater than")
  expect_error(gmm_independence_graph(g, n = c(100L, 200L)), "single")
  expect_error(gmm_independence_graph(g, n = 100L, alpha = 0), "alpha")
  expect_error(gmm_independence_graph(g, n = 100L, alpha = 1.5), "alpha")
  expect_error(gmm_independence_graph(g, n = 100L, threshold = 0.1), "not both")
})

test_that("a positional second argument is still the threshold", {
  s <- matrix(c(1, 0.3, 0.3, 1), 2L, 2L)
  g <- gmm(weights = 1, means = list(c(0, 0)), covariances = list(s))
  expect_equal(gmm_independence_graph(g, 0.5)[1L, 2L], 0L)
  expect_equal(gmm_independence_graph(g, 0.1)[1L, 2L], 1L)
})
