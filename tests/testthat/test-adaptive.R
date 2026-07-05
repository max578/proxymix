## PMC proposal adaptation in regime (iii): the refreshed proposal tracks
## the target, so the effective sample size recovers from a poor initial
## proposal, while `adapt = "none"` reproduces the fixed-proposal
## behaviour exactly.

test_that("adapt = 'none' is the default and matches the explicit call exactly", {
  q <- is_mvt(n_dim = 2L, mean = c(0, 0), sigma = 4 * diag(2), df = 5)
  f_default <- fit_kld_em(banana_target(), N = 2L, proposal = q,
                          is_size = 1200L, max_iter = 15L, seed = 1L,
                          validation_size = 0L)
  f_none <- fit_kld_em(banana_target(), N = 2L, proposal = q,
                       is_size = 1200L, max_iter = 15L, seed = 1L,
                       validation_size = 0L, adapt = "none")
  expect_equal(f_default@weights, f_none@weights, tolerance = 1e-14)
  expect_equal(f_default@means, f_none@means, tolerance = 1e-14)
  expect_equal(f_default@diagnostics$kld_final,
               f_none@diagnostics$kld_final, tolerance = 1e-14)
  expect_identical(f_default@diagnostics$n_refresh, 0L)
})

test_that("PMC recovers the effective sample size from a poor initial proposal", {
  ## A proposal four times too wide: the fixed-proposal ESS is mediocre,
  ## and every refresh should improve it as the proposal tracks the target.
  tgt <- banana_target()
  q_poor <- is_mvt(n_dim = 2L, mean = c(4, 4), sigma = 36 * diag(2), df = 5)
  f_fixed <- suppressWarnings(
    fit_kld_em(tgt, N = 3L, proposal = q_poor, is_size = 1500L,
               max_iter = 30L, seed = 2L, validation_size = 0L)
  )
  f_pmc <- suppressWarnings(
    fit_kld_em(tgt, N = 3L, proposal = q_poor, is_size = 1500L,
               max_iter = 30L, seed = 2L, validation_size = 0L,
               adapt = "pmc", refresh_every = 5L)
  )
  expect_gt(f_pmc@diagnostics$n_refresh, 0L)
  expect_length(f_pmc@diagnostics$ess_history,
                f_pmc@diagnostics$n_refresh + 1L)
  ## Final-batch ESS beats the fixed-proposal ESS decisively.
  expect_gt(f_pmc@diagnostics$ess, 5 * f_fixed@diagnostics$ess)
  ## The PMC fit's in-sample KLD is trustworthy AND small. (The fixed
  ## fit's kld_final is not comparable: at a collapsed ESS the in-sample
  ## estimate is biased optimistic -- it can even go negative on a
  ## normalised target -- which is exactly why the certificate flags it.)
  expect_lt(abs(f_pmc@diagnostics$kld_final), 0.15)
})

test_that("PMC fits are reproducible end-to-end under a seed", {
  tgt <- banana_target()
  q <- is_mvt(n_dim = 2L, mean = c(0, 0), sigma = 4 * diag(2), df = 5)
  f1 <- suppressWarnings(withr::with_seed(5,
    fit_kld_em(tgt, N = 2L, proposal = q, is_size = 1000L, max_iter = 12L,
               seed = 9L, validation_size = 0L,
               adapt = "pmc", refresh_every = 4L)))
  f2 <- suppressWarnings(withr::with_seed(777,
    fit_kld_em(tgt, N = 2L, proposal = q, is_size = 1000L, max_iter = 12L,
               seed = 9L, validation_size = 0L,
               adapt = "pmc", refresh_every = 4L)))
  expect_equal(f1@weights, f2@weights, tolerance = 1e-14)
  expect_equal(f1@means, f2@means, tolerance = 1e-14)
  expect_equal(f1@diagnostics$ess_history, f2@diagnostics$ess_history,
               tolerance = 1e-12)
})

test_that("PMC argument validation", {
  tgt <- banana_target()
  expect_error(fit_kld_em(tgt, N = 2L, adapt = "pmc", refresh_every = 0L),
               "refresh_every")
  expect_error(fit_kld_em(tgt, N = 2L, adapt = "pmc", defensive_gamma = 1),
               "defensive_gamma")
  expect_error(fit_kld_em(tgt, N = 2L, adapt = "pmc", inflate = 0.5),
               "inflate")
})
