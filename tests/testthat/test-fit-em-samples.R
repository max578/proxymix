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

test_that("a sample target with no rows is refused before initialisation", {
  tgt <- gmm_target_from_samples(matrix(numeric(0), ncol = 2))
  expect_error(fit_proxymix(tgt, N = 1L, regime = "sample"), "samples")
  expect_error(fit_proxymix(tgt, N = 1L, regime = "moment"), "samples")
})

test_that("regime (ii) keeps small-scale coordinates when scales differ widely", {
  sds <- c(0.026, 0.3, 5, 80, 5078)
  x <- withr::with_seed(1, sapply(sds, function(s) stats::rnorm(5000, 0, s)))
  fit <- fit_em_samples(gmm_target_from_samples(x), N = 1L, seed = 1L)
  ratio <- sqrt(diag(fit@covariances[[1L]])) / apply(x, 2L, stats::sd)
  expect_equal(ratio, rep(1, 5L), tolerance = 1e-3)
})
