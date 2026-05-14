test_that("auto regime picks moment for N == 1 with samples", {
  x <- matrix(stats::rnorm(200), ncol = 2)
  tgt <- gmm_target_from_samples(x)
  fit <- suppressMessages(fit_proxymix(tgt, N = 1L))
  expect_equal(fit@regime, "moment")
})

test_that("auto regime picks sample for N >= 2 with samples", {
  x <- matrix(stats::rnorm(200), ncol = 2)
  tgt <- gmm_target_from_samples(x)
  fit <- suppressMessages(
    fit_proxymix(tgt, N = 2L, max_iter = 25L, n_starts = 2L)
  )
  expect_equal(fit@regime, "sample")
})

test_that("auto regime picks kld for log-density-only target", {
  fit <- suppressMessages(
    fit_proxymix(banana_target(), N = 3L,
                 is_size = 1500L, max_iter = 15L, seed = 1L)
  )
  expect_equal(fit@regime, "kld")
})

test_that("explicit regime overrides auto", {
  mt <- mixture_target(with_samples = TRUE, n = 200L, seed = 1L)
  fit <- fit_proxymix(mt, N = 3L, regime = "kld",
                      is_size = 1000L, max_iter = 10L, seed = 1L)
  expect_equal(fit@regime, "kld")
})

test_that("auto errors when target has neither samples nor log_density", {
  expect_error(
    fit_proxymix(suppressWarnings(
      gmm_target(n_dim = 2L, samples = matrix(0, 0, 2),
                 log_density = function(x) 0)
    ), N = 1L, regime = "auto"),
    NA
  )
})
