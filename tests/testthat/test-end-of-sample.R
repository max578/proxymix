## End-of-sample instability test (Independent Oracle Principle: graded against
## the known break presence/absence in a simulated series, and -- for the
## parametric method -- the chi-square reference; never against the function's
## own output).

.eos_model <- function() list(
  prior = gmm(weights = 1, means = list(0), covariances = list(matrix(10))),
  dynamics = list(A = matrix(1), Q = matrix(0.04)),
  measurement = list(C = matrix(1), R = matrix(1)))

.sim_ll <- function(seed, n = 120L, q = 0.04, r = 1, brk = 0, m = 1L, obs_df = Inf) {
  withr::with_seed(seed, {
    level <- numeric(n)
    for (t in 2:n) level[t] <- level[t - 1L] + stats::rnorm(1, 0, sqrt(q))
    noise <- if (is.finite(obs_df)) stats::rt(n, obs_df) * sqrt(r * (obs_df - 2) / obs_df)
             else stats::rnorm(n, 0, sqrt(r))
    y <- level + noise
    if (brk != 0) y[(n - m + 1L):n] <- y[(n - m + 1L):n] + brk
    y
  })
}

test_that("gmm_eos_test validates its inputs", {
  mdl <- .eos_model()
  y <- .sim_ll(1L)
  expect_error(gmm_eos_test(matrix(1), mdl$dynamics, mdl$measurement, y), "gmm")
  expect_error(gmm_eos_test(gmm(weights = c(.5, .5), means = list(0, 1),
                                covariances = list(matrix(1), matrix(1))),
                            mdl$dynamics, mdl$measurement, y), "single-component")
  expect_error(gmm_eos_test(mdl$prior, list(A = matrix(1)), mdl$measurement, y), "Q")
  expect_error(gmm_eos_test(mdl$prior, mdl$dynamics, mdl$measurement, y, m = 0L),
               "1 <= m")
  expect_error(gmm_eos_test(mdl$prior, mdl$dynamics, mdl$measurement, y, m = 500L),
               "1 <= m")
})

test_that("a stable series is not flagged; a large break is, by both methods", {
  mdl <- .eos_model()
  y_stable <- .sim_ll(2L, brk = 0)
  y_break <- .sim_ll(2L, brk = 8, m = 1L)
  for (meth in c("chisq", "andrews")) {
    s <- gmm_eos_test(mdl$prior, mdl$dynamics, mdl$measurement, y_stable,
                      m = 1L, method = meth)
    b <- gmm_eos_test(mdl$prior, mdl$dynamics, mdl$measurement, y_break,
                      m = 1L, method = meth)
    expect_false(s$reject)
    expect_true(b$reject)
    expect_gt(b$statistic, s$statistic)
    expect_true(s$p_value >= 0 && s$p_value <= 1)
  }
})

test_that("the statistic equals the chi-square reference by construction", {
  ## Independent oracle: re-derive the last-m squared standardised innovation via
  ## a hand-rolled textbook Kalman recursion (not the package's filter), then the
  ## chi-square tail. The function must match both.
  mdl <- .eos_model(); y <- .sim_ll(5L); m <- 2L
  A <- 1; Q <- 0.04; C <- 1; R <- 1; mt <- 0; Pt <- 10
  z2 <- numeric(length(y))
  for (t in seq_along(y)) {
    mp <- A * mt; Pp <- A * Pt * A + Q
    S <- C * Pp * C + R; e <- y[t] - C * mp
    z2[t] <- e^2 / S
    K <- Pp * C / S; mt <- mp + K * e; Pt <- Pp - K * C * Pp
  }
  ref_stat <- sum(z2[(length(y) - m + 1L):length(y)])
  ref_p <- stats::pchisq(ref_stat, df = m, lower.tail = FALSE)
  res <- gmm_eos_test(mdl$prior, mdl$dynamics, mdl$measurement, y, m = m,
                      method = "chisq")
  expect_equal(res$statistic, ref_stat, tolerance = 1e-8)
  expect_equal(res$p_value, ref_p, tolerance = 1e-8)
})

test_that("size is near nominal and power is high (vs known truth)", {
  skip_on_cran()
  mdl <- .eos_model(); R <- 150L; alpha <- 0.05
  for (m in 1:3) {
    rej_null <- rej_brk <- logical(R)
    for (rr in seq_len(R)) {
      yn <- .sim_ll(rr, brk = 0, m = m)
      yb <- .sim_ll(1000L + rr, brk = 3, m = m)
      rej_null[rr] <- gmm_eos_test(mdl$prior, mdl$dynamics, mdl$measurement, yn,
                                   m = m, method = "chisq")$reject
      rej_brk[rr] <- gmm_eos_test(mdl$prior, mdl$dynamics, mdl$measurement, yb,
                                  m = m, method = "chisq")$reject
    }
    size <- mean(rej_null); power <- mean(rej_brk)
    expect_lt(size, 0.12)          # size near nominal 0.05 (MC band)
    expect_gt(power, 0.45)         # detects a delta=3 break
  }
})

test_that("the Andrews subsampling method stays calibrated under heavy tails", {
  skip_on_cran()
  ## The Andrews (2003) distribution-free advantage: with t(3) observation noise
  ## the chi-square reference is mis-specified and over-rejects, while the
  ## subsampling P-test keeps size near nominal. Oracle: known no-break truth.
  mdl <- .eos_model(); R <- 200L; alpha <- 0.05
  rej_chisq <- rej_andrews <- logical(R)
  for (rr in seq_len(R)) {
    y <- .sim_ll(7000L + rr, brk = 0, m = 1L, obs_df = 3)
    rej_chisq[rr] <- gmm_eos_test(mdl$prior, mdl$dynamics, mdl$measurement, y,
                                  m = 1L, method = "chisq", alpha = alpha)$reject
    rej_andrews[rr] <- gmm_eos_test(mdl$prior, mdl$dynamics, mdl$measurement, y,
                                    m = 1L, method = "andrews", alpha = alpha)$reject
  }
  ## the subsampling test is better calibrated than the parametric one here
  expect_lt(mean(rej_andrews), mean(rej_chisq) + 1e-9)
  expect_lt(mean(rej_andrews), 0.12)
})
