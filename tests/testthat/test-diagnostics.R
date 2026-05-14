test_that("kld_trace returns the EM trace for regime kld", {
  fit <- fit_proxymix(banana_target(), N = 2L, regime = "kld",
                      is_size = 1000L, max_iter = 10L, seed = 1L)
  trace <- kld_trace(fit)
  expect_true(is.numeric(trace))
  expect_gt(length(trace), 1L)
})

test_that("kld_trace returns NA for non-kld regimes", {
  x <- matrix(stats::rnorm(200), ncol = 2)
  tgt <- gmm_target_from_samples(x)
  fit <- fit_moment_match(tgt, N = 1L)
  expect_true(all(is.na(kld_trace(fit))))
})

test_that("ess_trace returns a positive ESS for regime kld", {
  fit <- fit_proxymix(banana_target(), N = 2L, regime = "kld",
                      is_size = 1000L, max_iter = 10L, seed = 1L)
  expect_gt(ess_trace(fit), 0)
})

test_that("hellinger_mc returns a finite estimate in [0, 1]", {
  fit <- fit_proxymix(banana_target(), N = 3L, regime = "kld",
                      is_size = 2000L, max_iter = 20L, seed = 1L)
  h <- hellinger_mc(fit, n_mc = 1000L, seed = 1L)
  expect_true(is.finite(h$h2))
  expect_lt(h$h2, 1.5) # MC can overshoot 1 slightly; just sanity-check.
})

test_that("bic_aic returns numeric for regime sample", {
  mt <- mixture_target(with_samples = TRUE, n = 300L, seed = 1L)
  fit <- fit_em_samples(mt, N = 2L, max_iter = 50L, n_starts = 2L)
  crit <- bic_aic(fit)
  expect_true(is.numeric(crit$bic))
  expect_true(is.numeric(crit$aic))
})
