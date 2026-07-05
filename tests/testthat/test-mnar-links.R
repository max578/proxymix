## Link-symmetric oracles for the gated (MNAR) moment engine: the logit
## link is the shipped default but previously had no moment-level oracle
## (only probit had the closed form), and probit had no end-to-end
## recovery test.

test_that("logit-link gated moments match numerical quadrature", {
  gs <- proxymix:::.gated_smooth
  cases <- list(
    list(mu = 0.4, s2 = 1.3, a = -0.5, b = 1.2),
    list(mu = -1.0, s2 = 0.6, a = 0.3, b = -0.8),
    list(mu = 2.0, s2 = 2.5, a = -2.0, b = 0.5)
  )
  for (cs in cases) {
    dens <- function(y) {
      stats::dnorm(y, cs$mu, sqrt(cs$s2)) * stats::plogis(cs$a + cs$b * y)
    }
    I_num <- stats::integrate(dens, -30, 30, rel.tol = 1e-11)$value
    m1_num <- stats::integrate(function(y) y * dens(y), -30, 30,
                               rel.tol = 1e-11)$value / I_num
    m2_num <- stats::integrate(function(y) y^2 * dens(y), -30, 30,
                               rel.tol = 1e-11)$value / I_num
    got <- gs(cs$mu, cs$s2, cs$a, cs$b, "logit")
    expect_equal(got$I, I_num, tolerance = 1e-6)
    expect_equal(got$m1, m1_num, tolerance = 1e-6)
    expect_equal(got$v, m2_num - m1_num^2, tolerance = 1e-5)
  }
})

test_that("probit-link MNAR imputation recovers the mean end-to-end", {
  skip_on_cran()
  dat_full <- withr::with_seed(41, {
    n <- 800L
    x <- stats::rnorm(n)
    y <- 1 + 0.8 * x + stats::rnorm(n, sd = 0.7)
    miss <- stats::runif(n) < stats::pnorm(-1.5 + 1.0 * y)
    list(x = x, y = y, miss = miss)
  })
  truth <- mean(dat_full$y)
  dat <- cbind(x = dat_full$x,
               y = ifelse(dat_full$miss, NA_real_, dat_full$y))
  imp <- suppressMessages(
    gmm_impute(dat, mechanism = mnar("y", beta = 1.0, link = "probit"),
               m = 8L, seed = 2L)
  )
  est <- mean(vapply(imp@completions, function(cc) mean(cc[, "y"]),
                     numeric(1L)))
  ignorable <- mean(dat[, "y"], na.rm = TRUE)
  ## The mechanism-aware fit must beat the ignorable (complete-case) mean
  ## and land close to the truth.
  expect_lt(abs(est - truth), abs(ignorable - truth))
  expect_lt(abs(est - truth), 0.12)
})
