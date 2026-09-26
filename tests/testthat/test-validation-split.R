test_that("validation_size = 0 leaves validation diagnostics NA", {
  fit <- fit_kld_em(banana_target(), N = 2L,
                    is_size = 1000L, max_iter = 10L, seed = 1L,
                    validation_size = 0L)
  expect_true(is.na(fit@diagnostics$validation_kld))
  expect_true(is.na(fit@diagnostics$validation_ess))
  expect_equal(fit@diagnostics$validation_size, 0L)
})

test_that("the default validation size draws until the standard error is small", {
  fit <- fit_kld_em(banana_target(), N = 3L,
                    proposal = is_mvt(n_dim = 2L, sigma = 4 * diag(2), df = 5),
                    is_size = 2000L, max_iter = 60L, seed = 1L)
  d <- fit@diagnostics
  expect_true(is.finite(d$validation_kld))
  expect_equal(d$validation_size %% 2000L, 0L)
  expect_lte(d$validation_size, 20000L)
  expect_lt(d$validation_mc_se, 0.15 * d$validation_kld)
})

test_that("an explicit validation_size draws exactly that many", {
  fit <- fit_kld_em(banana_target(), N = 2L,
                    is_size = 1000L, max_iter = 10L, seed = 1L,
                    validation_size = 250L)
  expect_equal(fit@diagnostics$validation_size, 250L)
  expect_equal(fit@diagnostics$n_target_evals, 1250L)
})

test_that("validation_mc_se matches the spread over repeated validation draws", {
  q <- is_mvt(n_dim = 2L, sigma = 4 * diag(2), df = 5)
  runs <- vapply(1:40, function(s1) {
    d <- fit_kld_em(banana_target(), N = 3L, proposal = q,
                    is_size = 1000L, max_iter = 30L, seed = 1L,
                    validation_size = 4000L,
                    validation_seed = s1)@diagnostics
    c(d$validation_kld, d$validation_mc_se)
  }, numeric(2L))
  ratio <- stats::sd(runs[1L, ]) / mean(runs[2L, ])
  expect_gt(ratio, 0.7)
  expect_lt(ratio, 1.4)
})

test_that("validation_size > 0 populates validation_kld / validation_ess", {
  fit <- fit_kld_em(banana_target(), N = 2L,
                    is_size = 1500L, max_iter = 15L, seed = 1L,
                    validation_size = 1500L)
  d <- fit@diagnostics
  expect_true(is.finite(d$validation_kld))
  expect_true(is.finite(d$validation_ess))
  expect_gt(d$validation_ess, 0)
  expect_lt(d$validation_ess, 1500L + 1)
})

test_that("validation_kld is close to in-sample kld_final at convergence", {
  withr::with_seed(2026, {
    fit <- fit_kld_em(banana_target(), N = 3L,
                      is_size = 3000L, max_iter = 40L, seed = 5L,
                      validation_size = 3000L, validation_seed = 17L)
    d <- fit@diagnostics
    ## Allow a few times the Monte Carlo SE.
    gap <- abs(d$validation_kld - d$kld_final)
    expect_lt(gap, 10 * d$mc_se_kld + 0.05)
  })
})

test_that("validation_proposal can differ from fitting proposal", {
  q_fit <- is_mvt(n_dim = 2L, sigma = 4 * diag(2), df = 5)
  q_val <- is_mvt(n_dim = 2L, sigma = 9 * diag(2), df = 5)
  fit <- fit_kld_em(banana_target(), N = 2L, proposal = q_fit,
                    is_size = 1500L, max_iter = 15L, seed = 1L,
                    validation_size = 1500L,
                    validation_proposal = q_val,
                    validation_seed = 3L)
  expect_equal(fit@diagnostics$validation_proposal_name, "is_mvt")
  expect_true(is.finite(fit@diagnostics$validation_kld))
})

test_that("ess_summary() exposes both fit and validation diagnostics", {
  fit <- fit_kld_em(banana_target(), N = 2L,
                    is_size = 1000L, max_iter = 10L, seed = 1L,
                    validation_size = 1000L)
  s <- ess_summary(fit)
  expect_true(is.finite(s$ess_relative))
  expect_true(is.finite(s$validation_ess_relative))
  expect_lte(s$max_weight, 1)
  expect_lte(s$validation_max_weight, 1)
  expect_gte(s$max_weight, 1 / 1000L)
})
