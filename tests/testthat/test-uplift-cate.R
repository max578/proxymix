## The decision layer, graded against stats::lm (independent oracle) and against
## known synthetic ground truth.
##
## At K = 1 the joint Gaussian's conditional mean is ordinary least squares, so
## proxy_cate() must reproduce the lm treatment coefficient exactly and its
## delta standard error must match the lm coefficient SE asymptotically. lm is
## computed independently of any proxymix output (Independent Oracle Principle).

test_that("proxy_cate at K = 1 equals the lm treatment coefficient (OLS)", {
  dat <- withr::with_seed(11L, {
    n <- 1500L
    x1 <- stats::rnorm(n)
    x2 <- stats::rnorm(n)
    t <- stats::rbinom(n, 1L, 0.5)
    y <- 0.7 + 1.3 * t + 0.5 * x1 - 0.4 * x2 + stats::rnorm(n, sd = 0.6)
    data.frame(y = y, t = t, x1 = x1, x2 = x2)
  })

  m <- fit_uplift(dat, outcome = "y", treatment = "t",
                  covariates = c("x1", "x2"),
                  N = 1L, regime = "moment", ridge_eps = 0)
  ce <- proxy_cate(m, newdata = data.frame(x1 = c(-1, 0, 1), x2 = c(1, 0, -1)))

  ols <- stats::lm(y ~ t + x1 + x2, data = dat)
  b_t <- unname(stats::coef(ols)["t"])

  ## K = 1 effect is constant in x and equals the lm coefficient on t.
  expect_equal(ce$tau, rep(b_t, 3L), tolerance = 1e-6)
})

test_that("proxy_cate delta SE matches the lm coefficient SE asymptotically", {
  dat <- withr::with_seed(22L, {
    n <- 4000L
    x <- stats::rnorm(n)
    t <- stats::rbinom(n, 1L, 0.5)
    y <- 0.2 + 0.9 * t - 0.6 * x + stats::rnorm(n, sd = 0.8)
    data.frame(y = y, t = t, x = x)
  })
  m <- fit_uplift(dat, "y", "t", "x", N = 1L, regime = "moment", ridge_eps = 0)
  ce <- proxy_cate(m, newdata = data.frame(x = 0))

  ols <- stats::lm(y ~ t + x, data = dat)
  se_lm <- unname(summary(ols)$coefficients["t", "Std. Error"])

  ## Delta SE is the OLS coefficient SE up to df corrections (n vs n - p).
  expect_equal(ce$se, se_lm, tolerance = 0.02)
})

test_that("a two-regime model recovers heterogeneous effects (ground truth)", {
  skip_on_cran()
  ## Effect grows with x: tau(x) = 0.5 + 1.5 x. Two regimes split on x let the
  ## mixture capture the heterogeneity a single line cannot.
  dat <- withr::with_seed(33L, {
    n <- 3000L
    x <- stats::rnorm(n)
    t <- stats::rbinom(n, 1L, 0.5)
    tau_true <- 0.5 + 1.5 * x
    y <- 1 + 0.3 * x + tau_true * t + stats::rnorm(n, sd = 0.5)
    data.frame(y = y, t = t, x = x)
  })
  m <- fit_uplift(dat, "y", "t", "x", N = 3L, regime = "sample",
                  max_iter = 200L, n_starts = 4L, seed = 7L)

  grid <- data.frame(x = seq(-1.5, 1.5, by = 0.5))
  ce <- proxy_cate(m, newdata = grid, se = FALSE)
  tau_true <- 0.5 + 1.5 * grid$x

  ## Effect estimates increase with x (slack-tolerant: strict adjacent
  ## differences from a local-optimum EM flake across BLAS builds).
  expect_true(all(diff(ce$tau) > -0.02))
  expect_gt(stats::cor(ce$tau, tau_true, method = "spearman"), 0.99)
  pehe <- sqrt(mean((ce$tau - tau_true)^2))
  expect_lt(pehe, 0.4)
})

test_that("proxy_decide treats exactly where value * tau exceeds cost", {
  dat <- withr::with_seed(44L, {
    n <- 2500L
    x <- stats::rnorm(n)
    t <- stats::rbinom(n, 1L, 0.5)
    y <- 1 + (0.5 + x) * t + stats::rnorm(n, sd = 0.5)
    data.frame(y = y, t = t, x = x)
  })
  m <- fit_uplift(dat, "y", "t", "x", N = 2L, regime = "sample",
                  max_iter = 150L, seed = 3L)
  grid <- data.frame(x = seq(-2, 2, by = 0.5))
  dec <- proxy_decide(m, grid, value = 1, cost = 0.4)
  ce <- proxy_cate(m, grid, se = FALSE)

  ## Action matches the hand threshold tau > cost / value = 0.4.
  expect_equal(dec$action, as.integer(ce$tau > 0.4))
  ## Expected value is monotone in tau.
  expect_equal(order(dec$expected_value), order(ce$tau))
})

test_that("proxy_overlap flags out-of-support covariate configurations", {
  dat <- withr::with_seed(55L, {
    n <- 800L
    x <- stats::rnorm(n)
    t <- stats::rbinom(n, 1L, 0.5)
    y <- x + t + stats::rnorm(n, sd = 0.5)
    data.frame(y = y, t = t, x = x)
  })
  m <- fit_uplift(dat, "y", "t", "x", N = 2L, regime = "sample",
                  max_iter = 100L, seed = 1L)
  ov <- proxy_overlap(m, newdata = data.frame(x = c(0, 10)))

  expect_false(ov$overlap_flag[1L])   # x = 0 is well covered
  expect_true(ov$overlap_flag[2L])    # x = 10 is far outside the support
})

test_that("the mc standard error is finite and positive", {
  dat <- withr::with_seed(66L, {
    n <- 600L
    x <- stats::rnorm(n)
    t <- stats::rbinom(n, 1L, 0.5)
    y <- 0.5 * t + x + stats::rnorm(n, sd = 0.6)
    data.frame(y = y, t = t, x = x)
  })
  m <- fit_uplift(dat, "y", "t", "x", N = 1L, regime = "moment", ridge_eps = 0)
  ce <- proxy_cate(m, data.frame(x = c(-1, 0, 1)),
                   se_method = "mc", B = 40L)
  expect_true(all(is.finite(ce$se)))
  expect_true(all(ce$se > 0))
})
