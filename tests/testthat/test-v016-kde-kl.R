## In-sample versus held-out divergence estimates for regime (iii) fits.

grid_hellinger2 <- function(log_f, g, lim1, lim2, n1 = 400L, n2 = 600L) {
  g1 <- seq(lim1[1L], lim1[2L], length.out = n1)
  g2 <- seq(lim2[1L], lim2[2L], length.out = n2)
  gm <- as.matrix(expand.grid(g1, g2))
  cell <- (g1[2L] - g1[1L]) * (g2[2L] - g2[1L])
  1 - sum(exp(0.5 * (log_f(gm) + dgmm(gm, g, log = TRUE)))) * cell
}

test_that("glance reports the held-out validation_kld", {
  skip_if_not_installed("generics")
  fit <- fit_kld_em(banana_target(), N = 2L, is_size = 1000L,
                    max_iter = 10L, seed = 1L, validation_size = 500L)
  gl <- generics::glance(fit)
  expect_true("validation_kld" %in% names(gl))
  expect_equal(gl$validation_kld, fit@diagnostics$validation_kld)
  fit0 <- fit_kld_em(banana_target(), N = 2L, is_size = 1000L,
                     max_iter = 10L, seed = 1L, validation_size = 0L)
  expect_true(is.na(generics::glance(fit0)$validation_kld))
})

test_that("hellinger_mc uses n_mc fresh draws and the seed for a regime kld fit", {
  fit <- fit_kld_em(banana_target(), N = 2L, is_size = 1000L,
                    max_iter = 15L, seed = 1L, validation_size = 0L)
  h1 <- hellinger_mc(fit, n_mc = 3000L, seed = 1L)
  h1b <- hellinger_mc(fit, n_mc = 3000L, seed = 1L)
  h2 <- hellinger_mc(fit, n_mc = 3000L, seed = 2L)
  expect_equal(h1$n_mc, 3000L)
  expect_identical(h1$h2, h1b$h2)
  expect_false(identical(h1$h2, h2$h2))
  expect_true(is.finite(h1$se) && h1$se > 0)
})

test_that("hellinger_mc on a regime kld fit agrees with grid quadrature", {
  tgt <- banana_target()
  fit <- fit_kld_em(tgt, N = 3L, is_size = 2000L, max_iter = 60L, seed = 1L,
                    validation_size = 0L,
                    proposal = is_mvt(2L, mean = c(0, 0),
                                      sigma = 4 * diag(2), df = 3))
  truth <- grid_hellinger2(tgt@log_density, fit, c(-7, 7), c(-5, 30))
  est <- hellinger_mc(fit, n_mc = 20000L, seed = 2L)
  expect_gt(est$h2, 0)
  expect_lt(abs(est$h2 - truth), 4 * est$se)
})

test_that("from_kde carries a held-out validation estimate by default", {
  x <- withr::with_seed(1L, matrix(stats::rnorm(160L), ncol = 2L))
  fit <- from_kde(x, N = 2L, is_size = 1000L, max_iter = 20L, seed = 1L)
  expect_equal(fit@diagnostics$validation_size, 250L)
  expect_true(is.finite(fit@diagnostics$validation_kld))
})

test_that("from_kde default is_size grows with N", {
  x <- withr::with_seed(2L, matrix(stats::rnorm(80L), ncol = 2L))
  for (N in c(1L, 4L)) {
    fit <- suppressWarnings(from_kde(x, N = N, max_iter = 2L, seed = 1L,
                                     validation_size = 0L))
    expect_equal(fit@diagnostics$is_size, 10000L + 2500L * N)
  }
})
