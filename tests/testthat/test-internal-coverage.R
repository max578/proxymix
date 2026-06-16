## Internal-helper and diagnostic error/edge branches.

test_that("logsumexp_rows guards orientation and handles the empty matrix", {
  expect_error(logsumexp_rows(1:3), "must be a matrix")
  expect_equal(logsumexp_rows(matrix(numeric(0), nrow = 0L, ncol = 2L)),
               numeric(0L))
})

test_that("kl_gauss falls back to determinant when S_a is not chol-able", {
  ## S_a indefinite -> chol() fails -> determinant() branch.
  val <- kl_gauss(c(0, 0), matrix(c(1, 2, 2, 1), 2, 2), c(0, 0), diag(2))
  expect_length(val, 1L)
  expect_true(is.finite(val))
})

test_that("is_pd distinguishes PD from indefinite matrices", {
  expect_true(is_pd(diag(2)))
  expect_false(is_pd(matrix(c(1, 2, 2, 1), 2, 2)))
})

test_that("as_input_matrix validates shape", {
  expect_error(as_input_matrix(matrix(1:4, ncol = 2L), p = 3L), "column")
  expect_error(as_input_matrix(1:3, p = 2L), "length 2")
  expect_error(as_input_matrix("x", p = 2L), "numeric matrix or vector")
})

test_that("diagnostic accessors reject a non-fit", {
  expect_error(kld_trace(42), "gmm_fit")
  expect_error(ess_trace(42), "gmm_fit")
  expect_error(hellinger_mc(42), "gmm_fit")
  expect_error(ess_summary(42), "gmm_fit")
  expect_error(bic_aic(42), "gmm_fit")
})

test_that("hellinger_mc errors without a target log-density", {
  ## A samples-only target has no log_density.
  fit <- fit_proxymix(gmm_target_from_samples(matrix(stats::rnorm(200), ncol = 2)),
                      N = 2L, regime = "sample", max_iter = 10L)
  expect_error(hellinger_mc(fit), "log_density")
})

test_that("hellinger_mc samples from the fit when the regime is not kld", {
  ## Target carries both samples and a normalised log-density: a sample-regime
  ## fit then exercises the sample-from-g_theta Hellinger branch.
  tgt <- banana_target(with_samples = TRUE, n = 600L, seed = 1L)
  fit <- fit_proxymix(tgt, N = 2L, regime = "sample", max_iter = 15L)
  h <- hellinger_mc(fit, n_mc = 800L, seed = 1L)
  expect_true(is.finite(h$h2))
  expect_true(is.finite(h$se))
  expect_gt(h$n_mc, 0)
})
