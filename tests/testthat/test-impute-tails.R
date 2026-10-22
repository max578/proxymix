## gmm_impute() under a detection limit. Expected truncated means are
## integrated numerically.

test_that("truncated draws far in a tail stay inside the interval", {
  cases <- list(c(29.2, 1, -0.98), c(0, 1, -40), c(10, 2, -5))
  for (cs in cases) {
    x <- withr::with_seed(2L, .draw_trunc(rep(cs[1], 2e4), cs[2], -Inf, cs[3]))
    expect_true(all(x <= cs[3]))
    dens <- function(y) {
      exp(stats::dnorm(y, cs[1], cs[2], log = TRUE) -
            stats::pnorm(cs[3], cs[1], cs[2], log.p = TRUE))
    }
    exact <- stats::integrate(function(y) y * dens(y), -Inf, cs[3])$value
    expect_lt(abs(mean(x) - exact), 0.01)
  }
  x <- withr::with_seed(3L, .draw_trunc(rep(0, 2e4), 1, 38, Inf))
  expect_true(all(x >= 38))
})

test_that("truncated draws away from the tails are unchanged", {
  mu <- seq(-3, 3, length.out = 500L)
  old <- withr::with_seed(1L, {
    pL <- stats::pnorm(-Inf)
    pU <- stats::pnorm((0.4 - mu) / 1.3)
    u <- pL + stats::runif(length(mu)) * (pU - pL)
    mu + 1.3 * stats::qnorm(pmin(pmax(u, 1e-12), 1 - 1e-12))
  })
  expect_identical(withr::with_seed(1L, .draw_trunc(mu, 1.3, -Inf, 0.4)), old)
})

test_that("a near-singular component does not stop imputation under a limit", {
  skip_on_cran()
  ## replicate 323 of the detection-limit pilot: t(3) data, p = 5, 20% censored
  S <- matrix(0.5, 5L, 5L)
  diag(S) <- 1
  Z <- withr::with_seed(1e6 * 4 + 1e5 + 1e4 * 5 + 323, {
    Z <- matrix(stats::rnorm(1000L), 200L, 5L) %*% chol(S)
    Z * sqrt(3 / stats::rchisq(200L, df = 3))
  })
  colnames(Z) <- c("y", paste0("x", 1:4))
  lod <- stats::qt(0.2, 3)
  D <- Z
  D[Z[, 1] < lod, 1] <- NA
  imp <- gmm_impute(D, m = 10L, mechanism = censored("y", upper = lod),
                    seed = 323L)
  filled <- unlist(lapply(gmm_complete(imp, "all"), function(d) {
    d[is.na(D[, 1]), 1]
  }))
  expect_true(all(filled <= lod))
})
