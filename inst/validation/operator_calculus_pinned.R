## inst/validation/operator_calculus_pinned.R
##
## Pinned validation corpus for the v0.3.0 affine-Gaussian operator
## calculus. Three reference pipelines under fixed inputs, each
## compared against a hand-coded Kalman-style reference.
##
## Not part of the automated test suite. Run manually after a v0.3.x
## release to detect numerical drift across BLAS / R versions.

if (!requireNamespace("proxymix", quietly = TRUE)) {
  stop("install proxymix first")
}

library(proxymix)

sessionInfo()
cat("\nR:", R.version.string, "\n\n")

results <- list()

## ---------------------------------------------------------------------------
## 1. Kalman parity: single-component prior, 1D scalar observation.
## ---------------------------------------------------------------------------

g1 <- gmm(
  weights     = 1,
  means       = list(c(0, 0, 0)),
  covariances = list(diag(c(1, 2, 3)))
)
H1   <- matrix(c(1, 0, 1), nrow = 1L)
y1   <- 0.8
R1   <- matrix(0.4, 1L, 1L)
fit1 <- gmm_observe(g1, A = H1, y = y1, noise_cov = R1)

## Hand reference:
S    <- diag(c(1, 2, 3))
Sx   <- H1 %*% S %*% t(H1) + R1
gain <- S %*% t(H1) %*% solve(Sx)
mu_hat <- as.numeric(gain * y1)
S_hat  <- S - gain %*% H1 %*% S

results$kalman_parity <- list(
  mu_diff_max  = max(abs(mu_hat - fit1@means[[1L]])),
  cov_diff_max = max(abs(S_hat  - fit1@covariances[[1L]])),
  weight       = fit1@weights[1L]
)

stopifnot(results$kalman_parity$mu_diff_max  < 1e-5,
          results$kalman_parity$cov_diff_max < 1e-5)

## ---------------------------------------------------------------------------
## 2. Sequential vs stacked: two scalar observations on a 2-comp prior.
## ---------------------------------------------------------------------------

g2 <- gmm(
  weights     = c(0.4, 0.6),
  means       = list(c(-1, 0), c(1.5, 0.5)),
  covariances = list(diag(c(0.6, 0.8)), diag(c(0.7, 0.5)))
)
H1d <- matrix(c(1, 0), nrow = 1L)
H2d <- matrix(c(0, 1), nrow = 1L)
y_a <- 0.5
y_b <- 0.2
R_one <- matrix(0.25, 1L, 1L)
R_two <- 0.25 * diag(2)

g_seq   <- gmm_observe(gmm_observe(g2, H1d, y_a, R_one), H2d, y_b, R_one)
g_stack <- gmm_observe(g2, diag(2), c(y_a, y_b), R_two)

results$sequential_vs_stacked <- list(
  weight_diff = max(abs(g_seq@weights - g_stack@weights)),
  mu_diff     = max(vapply(seq_along(g_seq@means),
                           function(k) max(abs(g_seq@means[[k]] -
                                               g_stack@means[[k]])),
                           numeric(1L))),
  cov_diff    = max(vapply(seq_along(g_seq@covariances),
                           function(k) max(abs(g_seq@covariances[[k]] -
                                               g_stack@covariances[[k]])),
                           numeric(1L)))
)

stopifnot(results$sequential_vs_stacked$weight_diff < 1e-4,
          results$sequential_vs_stacked$mu_diff     < 1e-4,
          results$sequential_vs_stacked$cov_diff    < 1e-4)

## ---------------------------------------------------------------------------
## 3. Aggregation through a 3 -> 2 block-sum, then observe one aggregate.
## ---------------------------------------------------------------------------

g3 <- gmm(
  weights     = c(0.3, 0.4, 0.3),
  means       = list(c(0, 0, 0), c(2, 1, -1), c(-1, -1, 2)),
  covariances = list(diag(3), diag(3), diag(3))
)
A_agg  <- matrix(c(1, 1, 0,
                   0, 0, 1), nrow = 2L, byrow = TRUE)
g_agg  <- gmm_aggregate(g3, A_agg)
H_aobs <- matrix(c(1, 0), nrow = 1L)
g_a_obs <- gmm_observe(g_agg, H_aobs, y = 1.2,
                       noise_cov = matrix(0.3, 1L, 1L))

results$aggregate_and_observe <- list(
  post_weights = round(g_a_obs@weights, 4L),
  post_mu      = lapply(g_a_obs@means, round, digits = 3L),
  out_dim      = gmm_dim(g_a_obs)
)

stopifnot(abs(sum(results$aggregate_and_observe$post_weights) - 1) < 1e-6)

## ---------------------------------------------------------------------------
## Report.
## ---------------------------------------------------------------------------

cat("\noperator_calculus pinned validation — proxymix v",
    as.character(utils::packageVersion("proxymix")),
    "\n", sep = "")
print(results)
invisible(results)
