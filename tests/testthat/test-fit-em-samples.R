test_that("regime (ii) is monotone in log-likelihood across iterations", {
  withr::with_seed(2026, {
    mt <- mixture_target(with_samples = TRUE, n = 500L, seed = 7L)
    fit <- fit_em_samples(mt, N = 3L, max_iter = 50L, n_starts = 2L)
    ll <- fit@diagnostics$loglik_trace
    expect_true(length(ll) >= 2L)
    deltas <- diff(ll)
    expect_true(all(deltas >= -1e-6))
  })
})

test_that("regime (ii) approaches mclust's log-likelihood on a 3-mixture", {
  skip_if_not_installed("mclust")
  skip_on_cran()
  ## mclust uses internal cross-package lookups that resolve only when its
  ## namespace is attached, so we attach it here and detach on exit.
  withr::local_namespace("mclust")
  withr::with_seed(2026, {
    mt <- mixture_target(with_samples = TRUE, n = 800L, seed = 11L)
    fit <- fit_em_samples(mt, N = 3L, max_iter = 200L, n_starts = 8L)
    mm <- mclust::Mclust(mt@samples, G = 3L, verbose = FALSE)
    ## Allow a small gap due to EM local optima.
    expect_lt(abs(fit@diagnostics$loglik_final - mm$loglik) /
                abs(mm$loglik), 0.05)
  })
})

test_that("regime (ii) populates BIC and AIC", {
  withr::with_seed(2026, {
    mt <- mixture_target(with_samples = TRUE, n = 400L, seed = 3L)
    fit <- fit_em_samples(mt, N = 2L, max_iter = 50L, n_starts = 2L)
    crit <- bic_aic(fit)
    expect_false(is.na(crit$bic))
    expect_false(is.na(crit$aic))
    expect_true(crit$bic > crit$aic)
  })
})

test_that("regime (ii) refuses targets without samples", {
  b <- banana_target() # log_density only
  expect_error(fit_em_samples(b, N = 2L), "samples")
})
