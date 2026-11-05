## Regime (i) of Hoek and Elliott (2024).
##
## Closed-form moment-matching fit.
##   N == 1 : exact moment match, mu = E[x], Sigma = Cov(x).
##   N > 1  : equal-weight mixture spread along the principal direction, with
##            within-component covariance reduced so that the mixture's mean
##            and covariance equal the target's exactly. Moments do not
##            determine a shape, so the result is a seed for regime (ii) or
##            (iii) and says so.

#' Closed-form moment-matching fit
#'
#' Implements regime (i) of Hoek and Elliott (2024). When `N == 1`, this is
#' the exact moment match: `mu` is the target mean and `Sigma` is the target
#' covariance. When `N > 1`, the function returns an equal-weight mixture
#' spread along the first principal direction whose mean and covariance also
#' equal the target's exactly. Half of the variance along that direction lies
#' between components, which keeps two components from forming two peaks.
#' Moments do not determine the shape of a distribution, so this mixture adds
#' no information to the `N = 1` fit: it is a starting point for the iterative
#' regimes, and a message says so.
#'
#' Either the `target` must carry an `n` by `p` `samples` matrix, or its
#' `metadata` slot must contain pre-computed `moments` of the form
#' `list(mean = <p-vec>, cov = <p-by-p>)`.
#'
#' @param target A [gmm_target].
#' @param N Number of components. `N >= 2` returns a mixture with the target's
#'   mean and covariance, without iterative refinement.
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
    if (nrow(target@samples) < 1L) {
      cli::cli_abort("`target@samples` has no rows; regime {.val moment} needs at least one sample.")
    }
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
      name = sprintf("moment_match[N=1] on %s", target@name),
      metadata = list(quality = .moment_quality())
    )
    return(if (isTRUE(canonicalise)) gmm_canonicalise(g) else g)
  }

  # N > 1: spread the means along PC1, keeping the moments exact -------------
  eig <- eigen(S_hat, symmetric = TRUE)
  pc1 <- eig$vectors[, 1L]
  var_between <- 0.5 * eig$values[1L]
  offsets <- seq(-1, 1, length.out = N)
  offsets <- offsets * sqrt(var_between / mean(offsets^2))
  means <- lapply(offsets, function(d) as.numeric(mu_hat + d * pc1))
  S_within <- symmetrise(S_hat - var_between * tcrossprod(pc1))
  covs <- replicate(N, S_within, simplify = FALSE)
  cli::cli_inform(c(
    "Moments do not determine a shape: this {N}-component mixture adds nothing to {.code N = 1}.",
    "i" = "It has the target's mean and covariance. Use it as {.arg init} for regime {.val sample} or {.val kld}."
  ))
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
      note = "equal-weight mixture along PC1 with the target's mean and covariance; no iterative refinement"
    ),
    converged = TRUE,
    iterations = 0L,
    call = match.call(),
    name = sprintf("moment_match[N=%d] on %s", N, target@name),
    metadata = list(quality = .moment_quality())
  )
  if (isTRUE(canonicalise)) gmm_canonicalise(g) else g
}

## The certificate of a closed-form moment match: deterministic, always
## converged, with the sampling-diagnostic fields not applicable.
.moment_quality <- function() {
  list(
    regime = "moment",
    converged = TRUE,
    degenerate = FALSE,
    ess = NA_real_,
    ess_relative = NA_real_,
    min_component_ess = NA_real_,
    max_weight = NA_real_,
    support_fraction = NA_real_,
    kld_final = NA_real_,
    kld_approx = NA_real_,
    validation_gap = NA_real_,
    heldout_kld = NA_real_,
    heldout_kld_se = NA_real_
  )
}
