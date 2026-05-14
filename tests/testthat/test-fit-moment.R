test_that("regime (i) moment match recovers the global mean and covariance", {
  withr::with_seed(2026, {
    n <- 600L
    p <- 3L
    truth_mu <- c(1, -2, 0.5)
    truth_S <- matrix(c(2, 0.3, 0,
                        0.3, 1, -0.1,
                        0, -0.1, 1.5), 3, 3)
    samples <- mvnfast::rmvn(n, mu = truth_mu, sigma = truth_S)
    tgt <- gmm_target_from_samples(samples)
    fit <- fit_moment_match(tgt, N = 1L)
    expect_equal(fit@regime, "moment")
    expect_equal(fit@means[[1L]], truth_mu, tolerance = 0.15)
    expect_equal(fit@covariances[[1L]], truth_S, tolerance = 0.2)
  })
})

test_that("regime (i) N>1 returns a deterministic moment-seeded mixture", {
  set.seed(1)
  samples <- matrix(stats::rnorm(200), ncol = 2)
  tgt <- gmm_target_from_samples(samples)
  fit <- fit_moment_match(tgt, N = 3L)
  expect_equal(fit@regime, "moment")
  expect_equal(gmm_n_components(fit), 3L)
  expect_true(grepl("PC1", fit@diagnostics$note %||% ""))
})

test_that("regime (i) errors when neither samples nor moments are supplied", {
  tgt <- gmm_target(n_dim = 2L,
                    log_density = function(x) -0.5 * rowSums(x^2))
  expect_error(fit_moment_match(tgt, N = 1L), "moment-match")
})
