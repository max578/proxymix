## Regime (i) of Hoek and Elliott (2024).
##
## Closed-form moment-matching fit.
##   N == 1 : exact moment match, mu = E[x], Sigma = Cov(x).
##   N > 1  : moment-seeded mixture along the principal direction (see
##            `init_moment_seed()`). No iterative refinement is applied;
##            the result is a deterministic seed that the user can then
##            feed into regime (ii) or (iii) for further refinement.

#' Closed-form moment-matching fit
#'
#' Implements regime (i) of Hoek and Elliott (2024). When `N == 1`, this is
#' the exact moment match: `mu` is the target mean and `Sigma` is the target
#' covariance. When `N > 1`, the function returns the deterministic
#' moment-seed of [init_moment_seed()] wrapped as a [gmm_fit], without
#' iterative refinement — useful as a starting point for the iterative
#' regimes.
#'
#' Either the `target` must carry an `n` by `p` `samples` matrix, or its
#' `metadata` slot must contain pre-computed `moments` of the form
#' `list(mean = <p-vec>, cov = <p-by-p>)`.
#'
#' @param target A [gmm_target].
#' @param N Number of components. `N >= 2` returns a moment-seeded mixture
#'   without iterative refinement.
#' @param ridge_eps Ridge added to the empirical covariance for numerical
#'   stability.
#' @param canonicalise Logical. If `TRUE` (the default), the fitted
#'   mixture is post-processed by [gmm_canonicalise()] so that components
#'   are sorted by descending weight and (as a tiebreaker) by descending
#'   `||mu||`. Set `FALSE` to retain the raw component order.
#'
#' @returns A [gmm_fit] with `regime = "moment"`.
#' @family fitting
#' @export
#' @examples
#' x <- matrix(stats::rnorm(200), ncol = 2)
#' tgt <- gmm_target_from_samples(x)
#' fit_moment_match(tgt, N = 1L)
fit_moment_match <- function(target, N = 1L, ridge_eps = 1e-6,
                             canonicalise = TRUE) {
  if (!S7::S7_inherits(target, gmm_target)) {
    cli::cli_abort("`target` must be a {.cls gmm_target} object.")
  }
  N <- as.integer(N)
  if (N < 1L) cli::cli_abort("`N` must be a positive integer.")
  p <- target@n_dim

  if (!is.null(target@samples)) {
    mu_hat <- colMeans(target@samples)
    S_hat <- ridge(stats::cov(target@samples), epsilon = ridge_eps)
    n_used <- nrow(target@samples)
  } else if (!is.null(target@metadata$moments)) {
    mu_hat <- as.numeric(target@metadata$moments$mean)
    S_hat <- target@metadata$moments$cov
    n_used <- NA_integer_
  } else {
    cli::cli_abort(c(
      "Cannot moment-match without samples or pre-computed moments.",
      "i" = "Either attach `samples` to the target or place {.code moments = list(mean = ., cov = .)} in {.code target@metadata}."
    ))
  }

  if (N == 1L) {
    g <- gmm_fit(
      weights = 1,
      means = list(mu_hat),
      covariances = list(S_hat),
      target = target,
      regime = "moment",
      diagnostics = list(
        n_used = n_used,
        empirical_mean = mu_hat,
        empirical_cov = S_hat
      ),
      converged = TRUE,
      iterations = 0L,
      call = match.call(),
      name = sprintf("moment_match[N=1] on %s", target@name)
    )
    return(if (isTRUE(canonicalise)) gmm_canonicalise(g) else g)
  }

  ## N > 1: spread N component means along PC1.
  eig <- eigen(S_hat, symmetric = TRUE)
  pc1 <- eig$vectors[, 1L]
  sd1 <- sqrt(eig$values[1L])
  offsets <- seq(-1, 1, length.out = N) * 1.5 * sd1
  means <- lapply(offsets, function(d) as.numeric(mu_hat + d * pc1))
  covs <- replicate(N, S_hat, simplify = FALSE)
  g <- gmm_fit(
    weights = rep(1 / N, N),
    means = means,
    covariances = covs,
    target = target,
    regime = "moment",
    diagnostics = list(
      n_used = n_used,
      empirical_mean = mu_hat,
      empirical_cov = S_hat,
      note = "moment-seeded mixture along PC1; no iterative refinement"
    ),
    converged = TRUE,
    iterations = 0L,
    call = match.call(),
    name = sprintf("moment_match[N=%d] on %s", N, target@name)
  )
  if (isTRUE(canonicalise)) gmm_canonicalise(g) else g
}
