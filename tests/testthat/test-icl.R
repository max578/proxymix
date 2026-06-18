## E4: integrated completed likelihood (ICL) in bic_aic().
##
## Independent oracles: ICL = BIC at K = 1 (the classification entropy is
## exactly zero when there is one component); the classification entropy
## recomputed by hand from the responsibilities; and the BCG (2000) ordering
## ICL >= BIC, strict when components overlap.

test_that("bic_aic keeps returning bic, aic, n_params (backward compatible)", {
  x <- withr::with_seed(1L, matrix(stats::rnorm(200), ncol = 2L))
  fit <- fit_em_samples(gmm_target_from_samples(x), N = 2L, max_iter = 50L)
  b <- bic_aic(fit)
  expect_true(all(c("bic", "aic", "icl", "classification_entropy",
                    "n_params") %in% names(b)))
  expect_true(is.finite(b$bic))
  expect_true(is.finite(b$aic))
})

test_that("ICL equals BIC at K = 1 (zero classification entropy)", {
  x <- withr::with_seed(2L, rbind(
    matrix(stats::rnorm(200, -4), ncol = 2L),
    matrix(stats::rnorm(200, 4), ncol = 2L)
  ))
  fit <- fit_em_samples(gmm_target_from_samples(x), N = 1L, max_iter = 50L)
  b <- bic_aic(fit)
  expect_equal(b$icl, b$bic, tolerance = 1e-10)
  expect_equal(b$classification_entropy, 0, tolerance = 1e-10)
})

test_that("classification entropy is non-negative and matches a hand recomputation", {
  x <- withr::with_seed(3L, rbind(
    matrix(stats::rnorm(200, -4), ncol = 2L),
    matrix(stats::rnorm(200, 4), ncol = 2L)
  ))
  tgt <- gmm_target_from_samples(x)
  fit <- fit_em_samples(tgt, N = 2L, max_iter = 100L)
  b <- bic_aic(fit)
  expect_gte(b$classification_entropy, 0)
  expect_gte(b$icl, b$bic)
  ## Independent recomputation of the responsibilities and their entropy.
  lu <- proxymix:::gmm_log_unnorm(x, fit@weights, fit@means, fit@covariances)
  lr <- lu - proxymix:::logsumexp_rows(lu)
  r <- exp(lr)
  en_manual <- -sum(ifelse(r > 0, r * log(r), 0))
  expect_equal(b$classification_entropy, en_manual, tolerance = 1e-10)
  expect_equal(b$icl, b$bic + 2 * en_manual, tolerance = 1e-10)
})

test_that("ICL strictly exceeds BIC when components overlap", {
  x <- withr::with_seed(4L, rbind(
    matrix(stats::rnorm(300, -0.5), ncol = 2L),
    matrix(stats::rnorm(300, 0.5), ncol = 2L)
  ))
  fit <- fit_em_samples(gmm_target_from_samples(x), N = 2L, max_iter = 100L)
  b <- bic_aic(fit)
  expect_gt(b$classification_entropy, 1)
  expect_gt(b$icl, b$bic + 1)
})

test_that("ICL is NA for regimes without an empirical likelihood", {
  ## moment regime
  x <- withr::with_seed(5L, matrix(stats::rnorm(200), ncol = 2L))
  fm <- fit_proxymix(gmm_target_from_samples(x), N = 1L, regime = "moment")
  expect_true(is.na(bic_aic(fm)$icl))
  expect_true(is.na(bic_aic(fm)$classification_entropy))
  ## kld regime
  fk <- fit_kld_em(banana_target(), N = 2L, is_size = 1500L, seed = 1L,
                   max_iter = 20L)
  expect_true(is.na(bic_aic(fk)$icl))
})
