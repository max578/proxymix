## Regression tests for the v0.11.6 correctness fixes. Each test pins a
## specific defect: the offset-dependent regime-(iii) stopping rule, the
## end-of-sample recursion dropping model offsets, the half-applied NaN
## guard in the log-sum-exp kernel, and the unfloored gate normaliser in
## the smooth-gated moments.

test_that("regime (iii) fitting is invariant to the target's normalising constant", {
  base <- banana_target()
  shifted <- gmm_target(
    n_dim = 2L,
    log_density = local({
      f <- base@log_density
      function(x) f(x) + 5000
    }),
    normalised = FALSE,
    name = "banana_shifted"
  )
  q <- is_mvt(n_dim = 2L, mean = c(0, 0), sigma = 4 * diag(2), df = 5)
  fit_a <- withr::with_seed(99,
    fit_kld_em(base, N = 3L, proposal = q, is_size = 1500L,
               max_iter = 40L, seed = 7L))
  fit_b <- withr::with_seed(99,
    fit_kld_em(shifted, N = 3L, proposal = q, is_size = 1500L,
               max_iter = 40L, seed = 7L))
  ## The stopping rule is judged on the offset-free EM objective, so a
  ## constant shift of the log-density must not change the trajectory.
  expect_identical(fit_b@iterations, fit_a@iterations)
  expect_identical(fit_b@converged, fit_a@converged)
  expect_equal(fit_b@weights, fit_a@weights, tolerance = 1e-12)
  expect_equal(fit_b@means, fit_a@means, tolerance = 1e-12)
  expect_equal(fit_b@covariances, fit_a@covariances, tolerance = 1e-12)
  ## The KLD trace itself shifts by exactly the constant (it is reported,
  ## not used for stopping).
  expect_equal(
    fit_b@diagnostics$kld_trace - fit_a@diagnostics$kld_trace,
    rep(5000, length(fit_a@diagnostics$kld_trace)),
    tolerance = 1e-8
  )
})

test_that("gmm_eos_test honours a measurement offset exactly", {
  prior <- gmm(weights = 1, means = list(0), covariances = list(matrix(10)))
  dyn <- list(A = matrix(1), Q = matrix(0.04))
  meas0 <- list(C = matrix(1), R = matrix(1))
  measd <- list(C = matrix(1), R = matrix(1), d = 3)
  y <- withr::with_seed(4, c(stats::rnorm(60), 5))
  t0 <- gmm_eos_test(prior, dyn, meas0, y, m = 1L, method = "chisq")
  td <- gmm_eos_test(prior, dyn, measd, y + 3, m = 1L, method = "chisq")
  ## Shifting the series and declaring the offset is the same problem.
  expect_equal(td$statistic, t0$statistic, tolerance = 1e-12)
  expect_equal(td$p_value, t0$p_value, tolerance = 1e-12)
})

test_that("gmm_eos_test honours a dynamics offset (previously silently dropped)", {
  prior <- gmm(weights = 1, means = list(0), covariances = list(matrix(10)))
  meas <- list(C = matrix(1), R = matrix(1))
  y <- withr::with_seed(4, c(stats::rnorm(60), 5))
  t0 <- gmm_eos_test(prior, list(A = matrix(1), Q = matrix(0.04)),
                     meas, y, m = 1L, method = "chisq")
  tb <- gmm_eos_test(prior, list(A = matrix(1), b = 0.5, Q = matrix(0.04)),
                     meas, y, m = 1L, method = "chisq")
  expect_false(isTRUE(all.equal(tb$statistic, t0$statistic)))
})

test_that("gmm_eos_test accepts a function-valued (time-varying) spec", {
  prior <- gmm(weights = 1, means = list(0), covariances = list(matrix(10)))
  y <- withr::with_seed(4, c(stats::rnorm(40), 5))
  t_const <- gmm_eos_test(prior, list(A = matrix(1), Q = matrix(0.04)),
                          list(C = matrix(1), R = matrix(1)),
                          y, m = 1L, method = "chisq")
  t_fun <- gmm_eos_test(prior,
                        function(t) list(A = matrix(1), Q = matrix(0.04)),
                        function(t) list(C = matrix(1), R = matrix(1)),
                        y, m = 1L, method = "chisq")
  expect_equal(t_fun$statistic, t_const$statistic, tolerance = 1e-12)
})

test_that("gmm_eos_test rejects Gaussian-sum noise", {
  prior <- gmm(weights = 1, means = list(0), covariances = list(matrix(10)))
  qmix <- gmm(weights = c(0.9, 0.1), means = list(0, 0),
              covariances = list(matrix(0.05), matrix(1)))
  y <- withr::with_seed(4, stats::rnorm(20))
  expect_error(
    gmm_eos_test(prior, list(A = matrix(1), Q = qmix),
                 list(C = matrix(1), R = matrix(1)), y, m = 1L),
    "Gaussian"
  )
})

test_that("logsumexp_rows treats NaN as excluded and +Inf as divergent", {
  lse <- proxymix:::logsumexp_rows
  m <- rbind(c(0, NaN), c(0, -Inf), c(-Inf, -Inf), c(Inf, 0))
  out <- lse(m)
  expect_equal(out[1L], 0)
  expect_equal(out[2L], 0)
  expect_identical(out[3L], -Inf)
  expect_identical(out[4L], Inf)
})

test_that(".gated_smooth survives a gate-normaliser underflow", {
  ## A conditional mean deep in the never-missing region underflows the
  ## gate normaliser I to zero; the floored division must stay finite
  ## rather than sending NaN into the M-step.
  gs <- proxymix:::.gated_smooth(mu = c(0, -60), s2 = 1, a = 0, b = 1,
                                 link = "probit")
  expect_true(all(is.finite(gs$m1)))
  expect_true(all(is.finite(gs$v)))
  expect_true(all(gs$I > 0))
})
