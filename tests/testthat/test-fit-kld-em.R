test_that("regime (iii) gives a small KLD on a 3-mixture target", {
  withr::with_seed(2026, {
    mt <- mixture_target(with_samples = FALSE)
    fit <- fit_kld_em(mt, N = 3L, is_size = 3000L,
                      max_iter = 40L, seed = 5L)
    expect_equal(fit@regime, "kld")
    expect_lt(abs(fit@diagnostics$kld_final), 0.05)
  })
})

test_that("regime (iii) ESS is reported and positive", {
  fit <- fit_kld_em(banana_target(), N = 3L,
                    is_size = 1500L, max_iter = 15L, seed = 1L)
  expect_gt(fit@diagnostics$ess, 0)
  expect_lt(fit@diagnostics$ess, 1500L + 1)
})

test_that("regime (iii) refuses targets without log_density", {
  x <- matrix(stats::rnorm(200), ncol = 2)
  tgt <- gmm_target_from_samples(x)
  expect_error(fit_kld_em(tgt, N = 2L), "log_density")
})

test_that("regime (iii) accepts a user-supplied proposal", {
  q <- is_mvt(n_dim = 2L, mean = c(0, 0),
              sigma = 4 * diag(2), df = 5)
  fit <- fit_kld_em(banana_target(), N = 3L, proposal = q,
                    is_size = 1500L, max_iter = 15L, seed = 1L)
  expect_equal(fit@diagnostics$proposal_name, "is_mvt")
})

test_that("regime (iii) reduces KLD over iterations on average", {
  withr::with_seed(2026, {
    fit <- fit_kld_em(banana_target(), N = 3L,
                      is_size = 2000L, max_iter = 25L, seed = 13L)
    trace <- fit@diagnostics$kld_trace
    ## KLD trace need not be strictly monotone (Monte Carlo noise), but
    ## the final value should be no worse than the first by a wide
    ## margin once the chain settles.
    expect_lt(trace[length(trace)], trace[1L])
  })
})
