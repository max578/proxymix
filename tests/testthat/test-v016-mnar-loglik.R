## Observed-data log-likelihood and convergence in proxy_mnar_sensitivity().

mnar_toy <- function(n = 300L) {
  withr::with_seed(20260622L, {
    k <- sample.int(2L, n, replace = TRUE)
    ctr <- rbind(c(0, 0), c(1.5, 0.5))
    z <- matrix(stats::rnorm(2L * n), n, 2L) %*%
      chol(matrix(c(1, 0.6, 0.6, 1), 2L)) + ctr[k, ]
    colnames(z) <- c("x1", "y")
    z[stats::runif(n) < stats::plogis(-0.5 + 0.7 * z[, 2L]), "y"] <- NA
    z
  })
}

## Reference by stats::integrate, row by row, independent of the package's
## Gauss-Hermite rule and conditional-density helpers.
mnar_loglik_ref <- function(g, d, alpha, beta) {
  lik <- vapply(seq_len(nrow(d)), function(i) {
    x <- d[i, 1L]
    y <- d[i, 2L]
    sum(vapply(seq_along(g@weights), function(k) {
      m <- g@means[[k]]
      S <- g@covariances[[k]]
      fx <- stats::dnorm(x, m[1L], sqrt(S[1L, 1L]))
      cm <- m[2L] + S[2L, 1L] / S[1L, 1L] * (x - m[1L])
      cs <- sqrt(S[2L, 2L] - S[2L, 1L]^2 / S[1L, 1L])
      if (is.na(y)) {
        pm <- stats::integrate(function(t) stats::dnorm(t, cm, cs) *
                                 stats::plogis(alpha + beta * t),
                               -Inf, Inf, rel.tol = 1e-10)$value
        g@weights[k] * fx * pm
      } else {
        g@weights[k] * fx * stats::dnorm(y, cm, cs) *
          (1 - stats::plogis(alpha + beta * y))
      }
    }, numeric(1L)))
  }, numeric(1L))
  sum(log(lik))
}

test_that("loglik is the selection model's observed-data log-likelihood", {
  d <- mnar_toy()
  st <- proxy_mnar_sensitivity(d, "y", beta_grid = c(0, 0.7), N = 2L, m = 2L,
                               seed = 1L)
  expect_true(all(c("loglik", "converged") %in% names(st)))
  for (b in c(0, 0.7)) {
    imp <- gmm_impute(d, N = 2L, m = 2L, mechanism = mnar("y", beta = b),
                      seed = 1L, max_iter = 500L)
    ref <- mnar_loglik_ref(imp@point_fit, d, imp@diagnostics$mnar_alpha, b)
    expect_equal(st$loglik[st$beta == b], ref, tolerance = 1e-6)
  }
})

test_that("the data lean away from missing at random under the mixture", {
  d <- mnar_toy()
  st <- proxy_mnar_sensitivity(d, "y", beta_grid = c(0, 0.7), N = 2L, m = 2L,
                               seed = 1L)
  expect_true(all(st$converged))
  expect_gt(st$loglik[2L], st$loglik[1L])
})

test_that("an unconverged fit is flagged and warned about", {
  d <- mnar_toy()
  expect_warning(
    st <- proxy_mnar_sensitivity(d, "y", beta_grid = c(0, 0.7), N = 2L, m = 2L,
                                 seed = 1L, max_iter = 2L),
    "did not converge"
  )
  expect_false(any(st$converged))
})
