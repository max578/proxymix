## Tests for from_kde() — v0.2.0 graduation.
##
## C-tier coverage:
##   - end-to-end recovery on a known mixture (within MC tolerance);
##   - bandwidth sensitivity (component variance monotone in h);
##   - dimensional guard (warn at p > 5, error at p > 10);
##   - normalisation declaration matches contract;
##   - default proposal sanity (ESS > 50%);
##   - metadata pass-through.

test_that("from_kde recovers a bimodal Gaussian mixture", {
  skip_if_not_installed("mvnfast")
  set.seed(1L)
  x <- rbind(
    mvnfast::rmvn(200L, mu = c(-2, 0), sigma = diag(2)),
    mvnfast::rmvn(200L, mu = c( 2, 0), sigma = diag(2))
  )
  fit <- from_kde(
    x, N = 2L, bandwidth = "silverman",
    is_size = 2000L, max_iter = 60L, seed = 1L
  )

  expect_s7_class(fit, gmm_fit)
  expect_equal(gmm_n_components(fit), 2L)
  expect_equal(fit@regime, "kld")
  ## Means recovered to within ~0.5 of truth.
  mus <- do.call(rbind, fit@means)
  mu1 <- mus[which.min(mus[, 1L]), ]
  mu2 <- mus[which.max(mus[, 1L]), ]
  expect_lt(abs(mu1[1L] - (-2)), 0.5)
  expect_lt(abs(mu2[1L] - ( 2)), 0.5)
  expect_lt(abs(mu1[2L]),         0.5)
  expect_lt(abs(mu2[2L]),         0.5)
})

test_that("from_kde target declares normalised = TRUE", {
  set.seed(2L)
  x <- matrix(stats::rnorm(120L), ncol = 2L)
  fit <- from_kde(x, N = 1L, is_size = 600L, max_iter = 20L, seed = 2L)
  expect_true(isTRUE(fit@target@normalised))
  expect_equal(fit@target@log_normalizer, 0)
  expect_false(isTRUE(fit@diagnostics$kld_is_shifted))
})

test_that("from_kde produces a healthy default proposal (ESS > 50%)", {
  set.seed(3L)
  x <- matrix(stats::rnorm(200L), ncol = 2L)
  fit <- from_kde(x, N = 1L, is_size = 1500L, max_iter = 25L, seed = 3L)
  expect_gt(fit@diagnostics$ess_relative, 0.50)
  expect_lt(fit@diagnostics$max_weight, 0.05)
})

test_that("from_kde refuses p > 10 and warns 5 < p <= 10", {
  set.seed(4L)
  bad <- matrix(stats::rnorm(11L * 20L), ncol = 11L)
  expect_error(from_kde(bad, N = 2L, is_size = 200L),
               regexp = "supports.*p <= 10")
  warn <- matrix(stats::rnorm(60L * 6L), ncol = 6L)
  expect_warning(
    fit <- from_kde(warn, N = 1L, is_size = 400L,
                    max_iter = 8L, seed = 4L),
    regexp = "p <= 5"
  )
  expect_s7_class(fit, gmm_fit)
})

test_that("from_kde validates samples (matrix, no NA, n >= 5, positive p)", {
  expect_error(from_kde(as.data.frame(matrix(0, 10, 2)), N = 1L),
               regexp = "numeric matrix")
  bad_na <- matrix(stats::rnorm(20L), ncol = 2L); bad_na[1L, 1L] <- NA
  expect_error(from_kde(bad_na, N = 1L), regexp = "NA")
  expect_error(from_kde(matrix(0, 3L, 2L), N = 1L),
               regexp = "at least 5 rows")
})

test_that("from_kde bandwidth selection: scalar, vector, silverman, scott", {
  set.seed(5L)
  x <- matrix(stats::rnorm(80L), ncol = 2L)
  ## scalar
  fit_s <- from_kde(x, N = 1L, bandwidth = 0.5,
                    is_size = 400L, max_iter = 8L, seed = 5L)
  expect_equal(fit_s@target@metadata$bandwidth, c(0.5, 0.5))
  ## vector
  fit_v <- from_kde(x, N = 1L, bandwidth = c(0.4, 0.6),
                    is_size = 400L, max_iter = 8L, seed = 5L)
  expect_equal(fit_v@target@metadata$bandwidth, c(0.4, 0.6))
  ## silverman
  fit_si <- from_kde(x, N = 1L, bandwidth = "silverman",
                     is_size = 400L, max_iter = 8L, seed = 5L)
  expect_equal(fit_si@target@metadata$bandwidth_method, "silverman")
  ## scott
  fit_sc <- from_kde(x, N = 1L, bandwidth = "scott",
                     is_size = 400L, max_iter = 8L, seed = 5L)
  expect_equal(fit_sc@target@metadata$bandwidth_method, "scott")
  ## rejection
  expect_error(from_kde(x, N = 1L, bandwidth = c(0.4, 0.6, 0.8)),
               regexp = "length 1 or 2")
  expect_error(from_kde(x, N = 1L, bandwidth = -0.5),
               regexp = "positive")
})

test_that("from_kde bandwidth sensitivity: larger h -> smoother proxy", {
  set.seed(6L)
  x <- matrix(stats::rnorm(120L), ncol = 2L)
  fits <- lapply(c(0.3, 1.5), function(h) {
    from_kde(x, N = 1L, bandwidth = h,
             is_size = 1500L, max_iter = 30L, seed = 6L)
  })
  tr_small <- sum(diag(fits[[1L]]@covariances[[1L]]))
  tr_large <- sum(diag(fits[[2L]]@covariances[[1L]]))
  expect_gt(tr_large, tr_small)
})
