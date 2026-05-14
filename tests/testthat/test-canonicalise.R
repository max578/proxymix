test_that("gmm_canonicalise sorts components by weight desc, ||mu|| desc", {
  g <- gmm(weights = c(0.1, 0.6, 0.3),
           means = list(c(0, 0), c(3, 0), c(-1, 1)),
           covariances = list(diag(2), diag(2), diag(2)))
  c1 <- gmm_canonicalise(g)
  expect_equal(c1@weights, c(0.6, 0.3, 0.1), tolerance = 1e-12)
  expect_equal(c1@means[[1L]], c(3, 0))
})

test_that("gmm_canonicalise tiebreaks by ||mu|| when weights are equal", {
  g <- gmm(weights = c(0.5, 0.5),
           means = list(c(0.5, 0), c(2, 0)),
           covariances = list(diag(2), diag(2)))
  c1 <- gmm_canonicalise(g)
  expect_equal(c1@means[[1L]], c(2, 0))
  expect_equal(c1@means[[2L]], c(0.5, 0))
})

test_that("gmm_canonicalise is a no-op on already-canonical input", {
  g <- gmm(weights = c(0.6, 0.4),
           means = list(c(2, 0), c(0, 0)),
           covariances = list(diag(2), diag(2)))
  c1 <- gmm_canonicalise(g)
  expect_identical(c1@weights, g@weights)
})

test_that("gmm_canonicalise preserves the mixture distribution", {
  g <- gmm(weights = c(0.1, 0.6, 0.3),
           means = list(c(0, 0), c(3, 0), c(-1, 1)),
           covariances = list(diag(2), diag(2), diag(2)))
  c1 <- gmm_canonicalise(g)
  x <- matrix(c(0, 0, 1, 0, -1, 0, 3, 0), ncol = 2, byrow = TRUE)
  expect_equal(dgmm(x, g), dgmm(x, c1), tolerance = 1e-12)
})

test_that("gmm_canonicalise preserves the gmm_fit subclass and slots", {
  mt <- mixture_target(with_samples = TRUE, n = 200L, seed = 1L)
  fit <- fit_em_samples(mt, N = 3L, max_iter = 20L, n_starts = 2L,
                        canonicalise = FALSE)
  c1 <- gmm_canonicalise(fit)
  expect_s7_class(c1, gmm_fit)
  expect_equal(c1@regime, fit@regime)
  expect_identical(c1@target, fit@target)
  expect_identical(c1@diagnostics, fit@diagnostics)
})

test_that("fit_*() canonicalises by default and respects canonicalise = FALSE", {
  set.seed(1)
  x <- matrix(stats::rnorm(400), ncol = 2)
  tgt <- gmm_target_from_samples(x)

  fit_default <- fit_em_samples(tgt, N = 2L, max_iter = 25L, n_starts = 2L)
  fit_raw <- fit_em_samples(tgt, N = 2L, max_iter = 25L, n_starts = 2L,
                            canonicalise = FALSE)
  ## After canonicalisation, weights are non-increasing.
  expect_true(all(diff(fit_default@weights) <= 1e-12))
  ## In the raw fit, weight order is whatever EM produced.
  expect_s7_class(fit_raw, gmm_fit)
})

test_that("regime kld also canonicalises by default", {
  fit <- fit_proxymix(banana_target(), N = 3L, regime = "kld",
                      is_size = 1000L, max_iter = 15L, seed = 1L)
  expect_true(all(diff(fit@weights) <= 1e-12))
})
