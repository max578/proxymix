test_that("one component equals the closed-form PPCA solution", {
  set.seed(1)
  p <- 6
  W <- matrix(rnorm(p * 2), p)
  x <- matrix(rnorm(400 * 2), 400) %*% t(W) + matrix(rnorm(400 * p, sd = 0.5), 400)
  fit <- fit_mppca_samples(gmm_target_from_samples(x), N = 1L, q = 2L, seed = 1L)

  ## Tipping and Bishop (1999): sigma^2 is the mean of the trailing eigenvalues
  ## of the sample covariance (divisor n); W W' adds lambda_j - sigma^2 along
  ## the leading eigenvectors.
  S <- stats::cov(x) * (nrow(x) - 1) / nrow(x)
  e <- eigen(S, symmetric = TRUE)
  s2 <- mean(e$values[3:p])
  U <- e$vectors[, 1:2]
  C <- U %*% diag(e$values[1:2] - s2) %*% t(U) + s2 * diag(p)

  expect_equal(fit@metadata$mppca$noise_var, s2, tolerance = 1e-10)
  expect_equal(gmm_covariances(fit)[[1]], C, tolerance = 1e-10)
  expect_equal(gmm_means(fit)[[1]], colMeans(x), tolerance = 1e-12)
})

test_that("q = p - 1 gives the full-covariance fit", {
  set.seed(11)
  S <- matrix(c(1.0, 0.8, 0.8, 1.0), 2)
  draw <- function(n, mu) sweep(matrix(rnorm(n * 2), n) %*% chol(S), 2, mu, "+")
  x <- rbind(draw(300, c(-2.5, 0)), draw(300, c(2.5, 0)))
  tgt <- gmm_target_from_samples(x)
  mp <- fit_mppca_samples(tgt, N = 2L, q = 1L, seed = 1L)
  full <- fit_em_samples(tgt, N = 2L, seed = 1L, max_iter = 500L, tol = 1e-10)
  for (k in 1:2) {
    j <- which.min(vapply(gmm_means(full), function(m) sum((m - gmm_means(mp)[[k]])^2), numeric(1)))
    expect_equal(gmm_covariances(mp)[[k]], gmm_covariances(full)[[j]], tolerance = 1e-4)
  }
  expect_equal(mp@diagnostics$loglik_final, full@diagnostics$loglik_final, tolerance = 1e-5)
})

test_that("clusters are recovered with fewer points per component than variables", {
  set.seed(22)
  p <- 30
  mus <- list(rnorm(p, 0, 0.6), rnorm(p, 0, 0.6))
  Ws <- lapply(1:2, function(k) qr.Q(qr(matrix(rnorm(p * 2), p))) %*% diag(c(3, 2)))
  z <- rep(1:2, each = 20)
  x <- t(vapply(z, function(k) {
    mus[[k]] + as.numeric(Ws[[k]] %*% rnorm(2)) + rnorm(p, sd = 0.5)
  }, numeric(p)))
  fit <- fit_mppca_samples(gmm_target_from_samples(x), N = 2L, q = 2L, seed = 1L)
  lp <- sapply(1:2, function(k) {
    d <- sweep(x, 2, gmm_means(fit)[[k]])
    R <- chol(gmm_covariances(fit)[[k]])
    log(gmm_weights(fit)[k]) - sum(log(diag(R))) - 0.5 * colSums(backsolve(R, t(d), transpose = TRUE)^2)
  })
  cl <- max.col(lp)
  expect_true(all(cl == z) || all(cl == 3 - z))
  expect_true(fit@converged)
  expect_true(all(diff(fit@diagnostics$loglik_trace) > -1e-8 * abs(fit@diagnostics$loglik_final)))
  expect_equal(fit@diagnostics$n_params, 1 + 2 * p + 2 * (p * 2 - 1 + 1))
})

test_that("the result works with the mixture operators", {
  set.seed(3)
  x <- matrix(rnorm(300 * 4), 300) %*% matrix(rnorm(16), 4)
  fit <- fit_mppca_samples(gmm_target_from_samples(x), N = 2L, q = 2L, seed = 1L)
  expect_true(S7::S7_inherits(fit, gmm_fit))
  expect_equal(gmm_dim(gmm_marginalise(fit, c(1L, 2L))), 2L)
  expect_equal(gmm_dim(gmm_missing(fit, observed = 1L, values = 0)), 3L)
  expect_length(fit@metadata$mppca$loadings, 2L)
  expect_equal(dim(fit@metadata$mppca$loadings[[1]]), c(4L, 2L))
})

test_that("invalid arguments are refused", {
  x <- matrix(rnorm(100), ncol = 4)
  tgt <- gmm_target_from_samples(x)
  expect_error(fit_mppca_samples(tgt, N = 2L, q = 4L), "q")
  expect_error(fit_mppca_samples(tgt, N = 2L, q = 0L), "q")
  expect_error(fit_mppca_samples(x, N = 2L, q = 1L), "gmm_target")
  expect_error(fit_mppca_samples(tgt, N = 0L, q = 1L), "N")
})

test_that("defaults fit two components with a one-dimensional subspace", {
  set.seed(5)
  x <- rbind(matrix(rnorm(100 * 3), 100), matrix(rnorm(100 * 3), 100) + 5)
  fit <- fit_mppca_samples(gmm_target_from_samples(x))
  expect_equal(gmm_n_components(fit), 2L)
  expect_equal(fit@metadata$mppca$q, 1L)
  expect_equal(sort(round(vapply(gmm_means(fit), mean, numeric(1)))), c(0, 5))
})
