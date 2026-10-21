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

test_that("a component on fewer than p + 1 samples is flagged as collapsed", {
  x <- withr::with_seed(3, rbind(
    matrix(stats::rnorm(400), ncol = 2),
    matrix(c(10, 10), nrow = 2, ncol = 2, byrow = TRUE) + 1e-4 * stats::rnorm(4)
  ))
  init <- gmm(weights = c(0.99, 0.01),
              means = list(c(0, 0), c(10, 10)),
              covariances = list(diag(2), diag(0.01, 2)))
  expect_warning(
    fit <- fit_em_samples(gmm_target_from_samples(x), N = 2L, init = init),
    class = "proxymix_component_collapse"
  )
  expect_identical(fit@diagnostics$collapsed_components, 2L)
  expect_true(gmm_fit_quality(fit)$degenerate)
})

test_that("a point mass of many identical rows is flagged as collapsed", {
  x <- withr::with_seed(1, rbind(
    matrix(stats::rnorm(400), ncol = 2),
    matrix(10, nrow = 30, ncol = 2)
  ))
  init <- gmm(weights = c(0.87, 0.13),
              means = list(c(0, 0), c(10, 10)),
              covariances = list(diag(2), diag(1e-8, 2)))
  expect_warning(
    fit <- fit_em_samples(gmm_target_from_samples(x), N = 2L, init = init,
                          max_iter = 40L),
    class = "proxymix_component_collapse"
  )
  expect_identical(fit@diagnostics$collapsed_components, 2L)
  expect_true(gmm_fit_quality(fit)$degenerate)
})

test_that("a well-posed fit records no collapse and no decrease", {
  withr::with_seed(2026, {
    mt <- mixture_target(with_samples = TRUE, n = 500L, seed = 7L)
    expect_no_warning(fit <- fit_em_samples(mt, N = 3L, max_iter = 50L, n_starts = 2L))
  })
  expect_identical(fit@diagnostics$collapsed_components, integer(0L))
  expect_false(fit@diagnostics$nonmonotone)
  expect_gte(fit@diagnostics$loglik_max_drop, 0)
  expect_false(gmm_fit_quality(fit)$degenerate)
})

test_that("select_N reports collapse per candidate and warns only for the chosen fit", {
  x <- withr::with_seed(5, rbind(
    matrix(stats::rnorm(300, mean = -4), ncol = 2),
    matrix(stats::rnorm(300, mean = 4), ncol = 2)
  ))
  expect_no_warning(sel <- select_N(gmm_target_from_samples(x), candidates = 1:6))
  expect_true("collapsed" %in% names(sel$table))
  expect_false(sel$table$collapsed[sel$table$chosen])
})
