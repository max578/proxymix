## Synthetic validation against known ground truth (spec section 9.1). Each
## data-generating process has an analytic CATE, computed independently of any
## proxymix output, so PEHE / bias are graded against the truth (Independent
## Oracle Principle). These are the absolute-accuracy checks; the grf-relative
## comparison (S1, S2) belongs to the gated benchmark pass.

.pehe <- function(est, truth) sqrt(mean((est - truth)^2))

test_that("crossing regimes: the effect changes sign with x", {
  skip_on_cran()
  dat <- withr::with_seed(701L, {
    n <- 3000L
    x <- stats::rnorm(n)
    t <- stats::rbinom(n, 1L, 0.5)
    tau_true <- 2 * x                       # crosses zero at x = 0
    y <- 0.5 + tau_true * t + stats::rnorm(n, sd = 0.5)
    data.frame(y = y, t = t, x = x)
  })
  m <- fit_uplift(dat, "y", "t", "x", N = 3L, regime = "sample",
                  max_iter = 300L, n_starts = 6L, seed = 11L)
  grid <- data.frame(x = seq(-1.5, 1.5, by = 0.5))
  ce <- proxy_cate(m, grid, se = FALSE)
  tau_true <- 2 * grid$x

  expect_lt(.pehe(ce$tau, tau_true), 0.5)
  expect_lt(ce$tau[grid$x == -1.5], 0)
  expect_gt(ce$tau[grid$x == 1.5], 0)
})

test_that("heteroscedastic outcome: a constant effect is recovered", {
  dat <- withr::with_seed(702L, {
    n <- 3000L
    x <- stats::rnorm(n)
    t <- stats::rbinom(n, 1L, 0.5)
    ## Effect is a constant 1; noise scale grows with |x|.
    y <- 0.3 * x + 1 * t + stats::rnorm(n, sd = 0.3 + 0.6 * abs(x))
    data.frame(y = y, t = t, x = x)
  })
  m <- fit_uplift(dat, "y", "t", "x", N = 3L, regime = "sample",
                  max_iter = 200L, n_starts = 6L, seed = 12L)
  grid <- data.frame(x = seq(-1, 1, by = 0.5))
  ce <- proxy_cate(m, grid, se = FALSE)
  expect_lt(.pehe(ce$tau, rep(1, nrow(grid))), 0.3)
})

test_that("non-linear gate: a smooth effect is tracked by the mixture", {
  dat <- withr::with_seed(703L, {
    n <- 4000L
    x <- stats::runif(n, -2, 2)
    t <- stats::rbinom(n, 1L, 0.5)
    tau_true <- 1 / (1 + exp(-2 * x))       # smooth logistic ramp in [0, 1]
    y <- 0.2 + tau_true * t + stats::rnorm(n, sd = 0.3)
    data.frame(y = y, t = t, x = x)
  })
  m <- fit_uplift(dat, "y", "t", "x", N = 4L, regime = "sample",
                  max_iter = 250L, n_starts = 6L, seed = 13L)
  grid <- data.frame(x = seq(-1.5, 1.5, by = 0.5))
  ce <- proxy_cate(m, grid, se = FALSE)
  tau_true <- 1 / (1 + exp(-2 * grid$x))
  expect_lt(.pehe(ce$tau, tau_true), 0.2)
  ## Monotone like the truth, with a small slack: strict positivity of every
  ## adjacent difference from a local-optimum EM is one BLAS rounding away
  ## from a flake on another platform.
  expect_true(all(diff(ce$tau) > -0.02))
  expect_gt(stats::cor(ce$tau, tau_true, method = "spearman"), 0.99)
})

test_that("null effect: the estimate is near zero with covering intervals", {
  dat <- withr::with_seed(704L, {
    n <- 2000L
    x <- stats::rnorm(n)
    t <- stats::rbinom(n, 1L, 0.5)
    y <- 1 + 0.5 * x + stats::rnorm(n, sd = 0.5)   # treatment does nothing
    data.frame(y = y, t = t, x = x)
  })
  m <- fit_uplift(dat, "y", "t", "x", N = 2L, regime = "sample",
                  max_iter = 150L, seed = 14L)
  ce <- proxy_cate(m, data.frame(x = c(-1, 0, 1)))
  expect_lt(max(abs(ce$tau)), 0.15)
  ## 95% intervals should cover zero for a true null.
  expect_true(all(ce$ci_lo < 0 & ce$ci_hi > 0))
})

test_that("binary outcome: response-scale CATE is a sensible probability shift", {
  dat <- withr::with_seed(705L, {
    n <- 3000L
    x <- stats::rnorm(n)
    t <- stats::rbinom(n, 1L, 0.5)
    ## Treatment lifts the success probability, more so at high x.
    p <- stats::plogis(-0.5 + 0.5 * x + (0.8 + 0.5 * x) * t)
    y <- stats::rbinom(n, 1L, p)
    data.frame(y = y, t = t, x = x)
  })
  m <- fit_uplift(dat, "y", "t", "x", N = 3L, regime = "sample",
                  outcome_type = "binary", max_iter = 200L, seed = 15L)
  grid <- data.frame(x = c(-1, 0, 1))
  ce <- proxy_cate(m, grid, scale = "response", se = FALSE)

  ## Effects are valid probability shifts and positive (treatment helps).
  expect_true(all(ce$tau >= -1 & ce$tau <= 1))
  expect_true(all(ce$tau > 0))
  ## Larger lift at higher x.
  expect_true(ce$tau[3L] > ce$tau[1L])
})
