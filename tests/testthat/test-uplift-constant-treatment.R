## When each fitted regime holds a single treatment arm, the within-regime
## effect is zero by construction; the confounding gap and the segment table
## warn instead of reporting it silently.

.binary_arm_model <- function() {
  dat <- withr::with_seed(20260902L, {
    n <- 600L
    x <- stats::rnorm(n)
    t <- stats::rbinom(n, 1L, 0.5)
    y <- 1 + (0.5 + x) * t + stats::rnorm(n, sd = 0.5)
    data.frame(y = y, t = t, x = x)
  })
  fit_uplift(dat, "y", "t", "x", N = 2L, regime = "sample",
             max_iter = 80L, seed = 1L)
}

.continuous_treatment_model <- function() {
  dat <- withr::with_seed(2L, {
    n <- 600L
    x <- stats::rnorm(n)
    t <- stats::rnorm(n)
    y <- 1 + 0.5 * t + x + stats::rnorm(n, sd = 0.5)
    data.frame(y = y, t = t, x = x)
  })
  suppressWarnings(
    fit_uplift(dat, "y", "t", "x", N = 2L, regime = "sample",
               max_iter = 80L, seed = 1L)
  )
}

test_that("the gap and segments warn when each regime holds one treatment arm", {
  m <- .binary_arm_model()
  nd <- data.frame(x = c(-1, 0, 1))

  expect_warning(
    gap <- proxy_confounding_gap(m, nd),
    class = "proxymix_constant_treatment"
  )
  expect_equal(gap$tau_do, rep(0, 3L), tolerance = 1e-8)
  expect_warning(
    seg <- proxy_regime_segments(m),
    class = "proxymix_constant_treatment"
  )
  expect_equal(seg$effect, rep(0, 2L), tolerance = 1e-8)
})

test_that("the gap and segments do not warn when treatment varies within regimes", {
  m <- .continuous_treatment_model()
  nd <- data.frame(x = c(-1, 0, 1))

  expect_no_warning(proxy_confounding_gap(m, nd))
  expect_no_warning(proxy_regime_segments(m))
})
