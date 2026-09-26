## The audit layer: confounding gap, policy value, retrospective uplift,
## regime segments, identification report. Property tests plus a planted
## latent-confounder recovery (Independent Oracle: known data-generating truth).

test_that("the confounding gap is identically zero at K = 1", {
  dat <- withr::with_seed(101L, {
    n <- 800L
    x <- stats::rnorm(n)
    t <- stats::rbinom(n, 1L, 0.5)
    y <- 0.8 * t + x + stats::rnorm(n, sd = 0.5)
    data.frame(y = y, t = t, x = x)
  })
  m <- fit_uplift(dat, "y", "t", "x", N = 1L, regime = "moment", ridge_eps = 0)
  gap <- proxy_confounding_gap(m, data.frame(x = c(-1, 0, 1, 2)))

  ## One regime -> the gate cannot move with T -> observational == do.
  expect_equal(gap$gap, rep(0, 4L), tolerance = 1e-8)
  expect_equal(gap$tau_obs, gap$tau_do, tolerance = 1e-8)
})

test_that("the confounding gap recovers a planted latent confounder", {
  ## A latent class U drives BOTH treatment and outcome; the true causal effect
  ## of T is 0.5 per unit of T. T is continuous, so it varies within each
  ## regime. Across seeds 1-7 and 202 the do-mode slope fell in [0.475, 0.506]
  ## and the naive slope in [0.888, 0.961].
  dat <- withr::with_seed(202L, {
    n <- 2000L
    u <- stats::rbinom(n, 1L, 0.5)
    t <- 1.5 * u + stats::rnorm(n, sd = 0.5)
    x <- stats::rnorm(n)
    y <- 1 + 2 * u + 0.5 * t + 0.3 * x + stats::rnorm(n, sd = 0.4)
    data.frame(y = y, t = t, x = x)
  })
  expect_warning(
    m <- fit_uplift(dat, "y", "t", "x", N = 2L, regime = "sample",
                    max_iter = 300L, n_starts = 6L, seed = 5L),
    "not binary"
  )
  grid <- withr::with_seed(202L, data.frame(x = stats::rnorm(300)))
  expect_no_warning(gap <- proxy_confounding_gap(m, grid))
  dt <- diff(m@treatment_levels)

  ## Effects per unit of T: the naive estimate is confounded upward, the
  ## do-mode lands on the planted 0.5.
  expect_lt(abs(mean(gap$tau_do) / dt - 0.5), 0.05)
  expect_gt(mean(gap$tau_obs) / dt - 0.5, 0.3)
  expect_gt(mean(gap$gap) / dt, 0.3)
})

test_that("the gap warns when a binary confounder splits the arms across regimes", {
  dat <- withr::with_seed(202L, {
    n <- 4000L
    u <- stats::rbinom(n, 1L, 0.5)
    t <- stats::rbinom(n, 1L, 0.2 + 0.6 * u)
    x <- stats::rnorm(n)
    y <- 1 + 2 * u + 0.5 * t + 0.3 * x + stats::rnorm(n, sd = 0.4)
    data.frame(y = y, t = t, x = x)
  })
  m <- fit_uplift(dat, "y", "t", "x", N = 2L, regime = "sample",
                  max_iter = 300L, n_starts = 6L, seed = 5L)
  grid <- withr::with_seed(202L, data.frame(x = stats::rnorm(300)))

  expect_warning(
    gap <- proxy_confounding_gap(m, grid),
    class = "proxymix_constant_treatment"
  )
  expect_equal(gap$tau_do, rep(0, 300L), tolerance = 1e-8)
})

test_that("the optimal policy weakly dominates treat-all and treat-none", {
  dat <- withr::with_seed(303L, {
    n <- 1500L
    x <- stats::rnorm(n)
    t <- stats::rbinom(n, 1L, 0.5)
    y <- 1 + (0.3 + x) * t + stats::rnorm(n, sd = 0.5)
    data.frame(y = y, t = t, x = x)
  })
  m <- fit_uplift(dat, "y", "t", "x", N = 2L, regime = "sample",
                  max_iter = 120L, seed = 2L)
  nd <- data.frame(x = stats::rnorm(400))

  v_all <- proxy_policy_value(m, nd, "all", value = 1, cost = 0.3)$policy_value
  v_none <- proxy_policy_value(m, nd, "none", value = 1, cost = 0.3)$policy_value
  v_opt <- proxy_policy_value(m, nd, "optimal", value = 1, cost = 0.3)$policy_value

  expect_gte(v_opt, v_all - 1e-8)
  expect_gte(v_opt, v_none - 1e-8)

  ## treat-none treats nobody; treat-all treats everybody (no exclusions here).
  expect_equal(proxy_policy_value(m, nd, "none", 1, 0.3)$n_treated, 0L)
})

test_that("retrospective uplift equals proxy_cate at K = 1", {
  dat <- withr::with_seed(404L, {
    n <- 600L
    x <- stats::rnorm(n)
    t <- stats::rbinom(n, 1L, 0.5)
    y <- 0.9 * t + x + stats::rnorm(n, sd = 0.5)
    data.frame(y = y, t = t, x = x)
  })
  m <- fit_uplift(dat, "y", "t", "x", N = 1L, regime = "moment", ridge_eps = 0)

  ru <- proxy_retrospective_uplift(m, observed = dat[1:6, ])
  ce <- proxy_cate(m, newdata = dat[1:6, ], se = FALSE)
  ## One regime: the abduction gate is irrelevant -> both are the slope.
  expect_equal(ru$retro_uplift, ce$tau, tolerance = 1e-7)
})

test_that("regime segments report K rows that sum to a valid mixture", {
  dat <- withr::with_seed(505L, {
    n <- 800L
    x <- stats::rnorm(n)
    t <- stats::rbinom(n, 1L, 0.5)
    y <- 1 + (0.5 + x) * t + stats::rnorm(n, sd = 0.5)
    data.frame(y = y, t = t, x = x)
  })
  m <- fit_uplift(dat, "y", "t", "x", N = 3L, regime = "sample",
                  max_iter = 100L, seed = 1L)
  expect_warning(
    seg <- proxy_regime_segments(m),
    class = "proxymix_constant_treatment"
  )

  expect_equal(nrow(seg), 3L)
  expect_equal(sum(seg$weight), 1, tolerance = 1e-8)
  expect_true("x" %in% names(seg))
  expect_true(all(c("regime", "weight", "effect", "sigma") %in% names(seg)))
})

test_that("the identification report prints and carries an overlap rate", {
  dat <- withr::with_seed(606L, {
    n <- 600L
    x <- stats::rnorm(n)
    t <- stats::rbinom(n, 1L, 0.5)
    y <- 0.5 * t + x + stats::rnorm(n, sd = 0.5)
    data.frame(y = y, t = t, x = x)
  })
  m <- fit_uplift(dat, "y", "t", "x", N = 2L, regime = "sample",
                  max_iter = 80L, seed = 1L)
  expect_warning(
    rep <- proxy_identification_report(m, data.frame(x = stats::rnorm(100))),
    class = "proxymix_constant_treatment"
  )

  expect_s7_class(rep, uplift_identification)
  expect_true(rep@overlap_pct >= 0 && rep@overlap_pct <= 100)
  expect_equal(rep@K, 2L)
  expect_output(print(rep), "Identification report")
  expect_output(print(rep), "NOT identified")
})
