## The bootstrap ensemble (functional-space error budget) and the
## component-count selector.

test_that("the kld-regime ensemble covers a known functional", {
  ## Evaluable-only correlated Gaussian: the mean functional's truth is
  ## known exactly, and the 90% functional interval should contain it.
  Sig <- matrix(c(1.5, 0.6, 0.6, 1.0), 2)
  Si <- solve(Sig)
  logdet <- as.numeric(determinant(Sig, logarithm = TRUE)$modulus)
  tgt <- gmm_target(
    n_dim = 2L,
    log_density = function(x) {
      if (is.null(dim(x))) x <- matrix(x, ncol = 2L)
      xc <- sweep(x, 2L, c(1, -0.5))
      -0.5 * rowSums((xc %*% Si) * xc) - 0.5 * (2 * log(2 * pi) + logdet)
    },
    normalised = TRUE, name = "corr_gauss"
  )
  fit <- fit_kld_em(tgt, N = 1L, is_size = 2000L, max_iter = 40L, seed = 3L,
                    validation_size = 0L,
                    proposal = is_mvt(2L, sigma = 6 * diag(2), df = 5))
  ens <- gmm_fit_ensemble(fit, B = 120L, seed = 4L)
  ci <- proxy_functional_ci(ens, gmm_mean, level = 0.9)
  expect_equal(nrow(ci), 2L)
  expect_true(all(ci$conf.low < c(1, -0.5) & c(1, -0.5) < ci$conf.high))
  ## Intervals have positive width and bracket the point estimate.
  expect_true(all(ci$conf.low < ci$estimate & ci$estimate < ci$conf.high))
})

test_that("the sample-regime ensemble bootstraps the target's samples", {
  x <- withr::with_seed(6, matrix(stats::rnorm(600, mean = 2), ncol = 2))
  tgt <- gmm_target_from_samples(x)
  fit <- fit_em_samples(tgt, N = 1L, max_iter = 30L, seed = 1L)
  ens <- gmm_fit_ensemble(fit, B = 60L, seed = 2L)
  ci <- proxy_functional_ci(ens, function(g) gmm_mean(g)[1L], level = 0.9)
  expect_true(ci$conf.low < 2 && 2 < ci$conf.high)
  ## Width should be in the vicinity of the classical SE (sd/sqrt(n)).
  se <- stats::sd(x[, 1L]) / sqrt(nrow(x))
  width <- ci$conf.high - ci$conf.low
  expect_gt(width, 1.5 * se)
  expect_lt(width, 8 * se)
})

test_that("ensembles are reproducible under a seed and validate inputs", {
  fit <- fit_kld_em(banana_target(), N = 2L, is_size = 1000L,
                    max_iter = 12L, seed = 1L, validation_size = 0L)
  e1 <- gmm_fit_ensemble(fit, B = 20L, seed = 9L)
  e2 <- gmm_fit_ensemble(fit, B = 20L, seed = 9L)
  expect_equal(e1$members[[7L]]@means, e2$members[[7L]]@means,
               tolerance = 1e-14)
  expect_error(gmm_fit_ensemble(fit, B = 2L), "at least 10")
  expect_error(proxy_functional_ci(list(), gmm_mean), "gmm_ensemble")
  expect_error(proxy_functional_ci(e1, gmm_mean, level = 1.2), "inside")
})

test_that("select_N recovers the component count of an evaluable 3-mixture", {
  sel <- suppressWarnings(suppressMessages(
    select_N(mixture_target(with_samples = FALSE), candidates = 1:4,
             is_size = 3000L, max_iter = 40L, seed = 5L)
  ))
  expect_s3_class(sel, "proxymix_selection")
  expect_equal(sel$regime, "kld")
  expect_equal(sel$best_n, 3L)
  expect_equal(nrow(sel$table), 4L)
  expect_true(is.finite(sel$table$validation_mc_se[sel$table$chosen]))
})

test_that("select_N's sample arm agrees with mclust's BIC choice", {
  skip_if_not_installed("mclust")
  x <- withr::with_seed(11, rbind(
    mvnfast::rmvn(250L, mu = c(-3, 0), sigma = diag(2)),
    mvnfast::rmvn(250L, mu = c(3, 1), sigma = diag(2))
  ))
  sel <- select_N(gmm_target_from_samples(x), candidates = 1:4, seed = 2L,
                  max_iter = 60L)
  expect_equal(sel$regime, "sample")
  bicm <- mclust::mclustBIC(x, G = 1:4, modelNames = "VVV", verbose = FALSE)
  pick <- which(bicm == max(bicm, na.rm = TRUE), arr.ind = TRUE)
  expect_equal(sel$best_n, as.integer(rownames(bicm)[pick[1L, "row"]]))
})
