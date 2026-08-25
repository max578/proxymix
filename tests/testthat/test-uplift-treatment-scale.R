## PX-01 -- the decision family must read the treatment contrast from
## `model@treatment_levels` (what fit_uplift() observed in the data), not from
## a silent `{0, 1}` default. A treatment coded on a different scale (e.g. a
## nitrogen rate `{0, 100}`) must not silently return an effect on the wrong
## scale.

test_that("proxy_cate defaults t1/t0 to the observed treatment levels, not {0, 1}", {
  dat <- withr::with_seed(33L, {
    n <- 1500L
    x <- stats::rnorm(n)
    t <- stats::rbinom(n, 1L, 0.5) * 100
    y <- 1 + 2.0 * (t / 100) + 0.5 * x + stats::rnorm(n, sd = 0.5)
    data.frame(y = y, t = t, x = x)
  })

  m <- fit_uplift(dat, outcome = "y", treatment = "t", covariates = "x",
                  N = 1L, regime = "moment", ridge_eps = 0)
  expect_equal(m@treatment_levels, c(0, 100))

  ## No t1/t0 supplied: must use the observed {0, 100} contrast, i.e. an
  ## effect near 2.0 -- not the {0, 1}-coded 0.02.
  ce <- proxy_cate(m, newdata = data.frame(x = 0))
  expect_equal(ce$tau, 2.0, tolerance = 0.15)

  ## Explicit t1/t0 within the observed range still works and is unchanged.
  ce_explicit <- proxy_cate(m, newdata = data.frame(x = 0), t1 = 100, t0 = 0)
  expect_equal(ce$tau, ce_explicit$tau, tolerance = 1e-10)
})

test_that("proxy_cate aborts with a classed condition when t1/t0 fall outside the observed range", {
  dat <- withr::with_seed(34L, {
    n <- 800L
    x <- stats::rnorm(n)
    t <- stats::rbinom(n, 1L, 0.5) * 100
    y <- 1 + 0.02 * t + 0.5 * x + stats::rnorm(n, sd = 0.5)
    data.frame(y = y, t = t, x = x)
  })
  m <- fit_uplift(dat, outcome = "y", treatment = "t", covariates = "x",
                  N = 1L, regime = "moment", ridge_eps = 0)

  expect_error(
    proxy_cate(m, newdata = data.frame(x = 0), t1 = 1, t0 = 0),
    class = "proxymix_treatment_scale_error"
  )
})

test_that("proxy_overlap, proxy_decide and proxy_confounding_gap also default to the observed arms", {
  dat <- withr::with_seed(35L, {
    n <- 1200L
    x <- stats::rnorm(n)
    t <- stats::rbinom(n, 1L, 0.5) * 100
    y <- 1 + 2.0 * (t / 100) + 0.5 * x + stats::rnorm(n, sd = 0.5)
    data.frame(y = y, t = t, x = x)
  })
  m <- fit_uplift(dat, outcome = "y", treatment = "t", covariates = "x",
                  N = 1L, regime = "moment", ridge_eps = 0)

  ov <- proxy_overlap(m, newdata = data.frame(x = 0))
  expect_true(all(!ov$overlap_flag))

  dec <- proxy_decide(m, newdata = data.frame(x = 0), value = 1, cost = 0.5)
  expect_equal(dec$tau, 2.0, tolerance = 0.15)

  gap <- proxy_confounding_gap(m, newdata = data.frame(x = 0))
  expect_equal(gap$tau_obs, 2.0, tolerance = 0.15)
})
