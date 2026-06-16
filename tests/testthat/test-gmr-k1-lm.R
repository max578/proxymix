## The OLS reduction (Independent Oracle Principle: graded against `stats::lm`).
##
## A single-component (K = 1) Gaussian-mixture proxy on the joint (Y, X) encodes
## ordinary least squares exactly: the Schur-complement conditional mean of Y
## given X reproduces `predict(lm)`, and the conditional variance reproduces the
## residual variance. `lm` is computed independently of any proxymix output, so
## this anchors the regression-read layer (regress spec R1) to an external
## authority rather than to a self-consistent expectation.

test_that("K = 1 conditional mean equals predict(lm) (exact OLS reduction)", {
  dat <- withr::with_seed(123L, {
    n <- 200L
    x1 <- stats::rnorm(n)
    x2 <- stats::rnorm(n)
    y <- 1.5 + 2 * x1 - 0.7 * x2 + stats::rnorm(n, sd = 0.3)
    data.frame(y = y, x1 = x1, x2 = x2)
  })

  ## proxymix: one Gaussian on the joint, then condition on X (Y is coord 1).
  joint <- gmm_target_from_samples(cbind(dat$y, dat$x1, dat$x2))
  fit <- fit_proxymix(joint, N = 1L, regime = "moment", ridge_eps = 0)

  ## Independent oracle.
  ols <- stats::lm(y ~ x1 + x2, data = dat)

  newx <- cbind(c(-1, 0, 1, 2), c(0.5, -0.5, 1, -1))
  gmr_mean <- vapply(seq_len(nrow(newx)), function(i) {
    cond <- gmm_conditionalise(fit, given = c(NA, newx[i, 1L], newx[i, 2L]))
    cond@means[[1L]]
  }, numeric(1L))
  ols_mean <- as.numeric(stats::predict(
    ols, newdata = data.frame(x1 = newx[, 1L], x2 = newx[, 2L])
  ))

  expect_equal(gmr_mean, ols_mean, tolerance = 1e-8)
})

test_that("K = 1 conditional variance equals the OLS residual variance", {
  dat <- withr::with_seed(456L, {
    n <- 300L
    x1 <- stats::rnorm(n)
    x2 <- stats::rnorm(n)
    y <- -0.4 + 1.1 * x1 + 0.9 * x2 + stats::rnorm(n, sd = 0.5)
    data.frame(y = y, x1 = x1, x2 = x2)
  })

  joint <- gmm_target_from_samples(cbind(dat$y, dat$x1, dat$x2))
  fit <- fit_proxymix(joint, N = 1L, regime = "moment", ridge_eps = 0)
  ols <- stats::lm(y ~ x1 + x2, data = dat)

  ## Conditional variance is constant in X for a single Gaussian.
  cond <- gmm_conditionalise(fit, given = c(NA, 0, 0))
  cond_var <- cond@covariances[[1L]][1L, 1L]

  expect_equal(cond_var, stats::var(stats::residuals(ols)), tolerance = 1e-8)
})
