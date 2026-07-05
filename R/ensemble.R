## The proxy's own error budget. Everything downstream of a fit conditions
## on the fitted mixture as if it were exact; the ensemble makes the fit's
## sampling variability a first-class object. In regime (iii) this is
## unusually cheap: the fit is a functional of the weighted sample
## {(x_i, W_i)} and the expensive part -- the target evaluations -- is
## cached in the fit, so a Bayesian-bootstrap replicate (Rubin, 1981)
## re-weights the SAME evaluations and refits by a warm-started weighted
## EM with zero new target evaluations.

## The bootstrap-replicate engine: plain weighted EM on fixed observation
## weights, warm-started, with no empty-component reseeding (a replicate
## that loses a component simply carries it at its floor weight).
.em_weighted <- function(x, W, init, max_iter = 25L, tol = 1e-5,
                         ridge_eps = 1e-8) {
  n <- nrow(x)
  N <- length(init@weights)
  weights <- init@weights
  means <- init@means
  covs <- init@covariances
  prev <- -Inf
  for (it in seq_len(max_iter)) {
    log_unnorm <- gmm_log_unnorm(x, weights, means, covs)
    log_g <- logsumexp_rows(log_unnorm)
    obj <- sum(W * log_g)
    if (is.finite(prev) &&
          abs(obj - prev) / (abs(prev) + 1e-12) < tol) {
      break
    }
    prev <- obj
    resp <- exp(log_unnorm - log_g)
    W_resp <- resp * W
    Nk <- pmax(colSums(W_resp), 1e-12)
    weights <- as.numeric(Nk / sum(Nk))
    for (k in seq_len(N)) {
      mu_new <- as.numeric(colSums(W_resp[, k] * x) / Nk[k])
      diff <- x - matrix(mu_new, nrow = n, ncol = ncol(x), byrow = TRUE)
      S_new <- crossprod(diff * sqrt(W_resp[, k])) / Nk[k]
      means[[k]] <- mu_new
      covs[[k]] <- symmetrise(ridge(S_new, ridge_eps))
    }
  }
  gmm(weights = weights, means = means, covariances = covs,
      name = "ensemble_member")
}

#' Bootstrap ensemble of a fitted proxy
#'
#' Quantifies the sampling variability of the fitted mixture itself by a
#' Bayesian (weighted) bootstrap: each replicate re-weights the fit's own
#' observations with Dirichlet(1, ..., 1) weights and refits by a
#' warm-started weighted EM. In regime `"kld"` the observations are the
#' fit's cached importance draws with their self-normalised weights, so a
#' replicate costs **zero new target evaluations**; in regimes `"sample"`
#' and `"moment"` the observations are the target's samples. Summaries of
#' any functional of the proxy are then read off the ensemble with
#' [proxy_functional_ci()] -- functional-space intervals sidestep the
#' label-switching that makes parameter-space intervals incoherent.
#'
#' @param fit A [gmm_fit] whose fitting inputs are recoverable (a
#'   regime-`"kld"` fit carries its importance draws; regimes `"sample"`
#'   and `"moment"` require the target to carry its samples).
#' @param B Number of bootstrap replicates.
#' @param max_iter,tol Convergence controls for the per-replicate
#'   warm-started weighted EM.
#' @param seed Optional integer seed for the replicate weights.
#'
#' @returns A list of class `gmm_ensemble`: `fit` (the base fit),
#'   `members` (a length-`B` list of [gmm]), `B`, and `regime`.
#' @family diagnostics
#' @references Rubin, D. B. (1981) The Bayesian bootstrap. *The Annals of
#'   Statistics* 9(1), 130--134.
#' @export
#' @examples
#' fit <- fit_proxymix(banana_target(), N = 2L, regime = "kld",
#'                     is_size = 1500L, max_iter = 20L, seed = 1L)
#' ens <- gmm_fit_ensemble(fit, B = 30L, seed = 2L)
#' proxy_functional_ci(ens, gmm_mean)
gmm_fit_ensemble <- function(fit, B = 200L, max_iter = 25L, tol = 1e-5,
                             seed = NULL) {
  if (!S7::S7_inherits(fit, gmm_fit)) {
    cli::cli_abort("`fit` must be a {.cls gmm_fit}.")
  }
  B <- as.integer(B)
  if (length(B) != 1L || is.na(B) || B < 10L) {
    cli::cli_abort("`B` must be a single integer of at least 10.")
  }
  if (identical(fit@regime, "kld")) {
    x <- fit@diagnostics$is_sample
    log_W <- fit@diagnostics$is_log_weights
    if (is.null(x) || is.null(log_W)) {
      cli::cli_abort("the fit does not carry its importance draws; refit before building an ensemble.")
    }
    W0 <- exp(log_W)
    W0[!is.finite(W0)] <- 0
    W0 <- W0 / sum(W0)
  } else {
    tgt <- fit@target
    if (is.null(tgt) || is.null(tgt@samples)) {
      cli::cli_abort("regime {.val {fit@regime}} needs the target's samples to bootstrap.")
    }
    x <- tgt@samples
    W0 <- rep(1 / nrow(x), nrow(x))
  }
  n <- nrow(x)
  init <- gmm(weights = fit@weights, means = fit@means,
              covariances = fit@covariances)

  run <- function() {
    lapply(seq_len(B), function(b) {
      xi <- stats::rexp(n)
      Wb <- xi * W0
      Wb <- Wb / sum(Wb)
      .em_weighted(x, Wb, init, max_iter = max_iter, tol = tol)
    })
  }
  members <- if (is.null(seed)) run() else withr::with_seed(seed, run())

  structure(
    list(fit = fit, members = members, B = B, regime = fit@regime),
    class = "gmm_ensemble"
  )
}

#' @export
print.gmm_ensemble <- function(x, ...) {
  cat(sprintf("<gmm_ensemble>: %d bootstrap replicates of a %s-regime fit (K = %d, p = %d)\n",
              x$B, x$regime, gmm_n_components(x$fit), gmm_dim(x$fit)))
  invisible(x)
}

#' Percentile interval for any functional of a fitted proxy
#'
#' Applies a functional to every member of a bootstrap ensemble and
#' returns the base fit's point estimate with percentile confidence
#' limits. The functional may return a scalar or a fixed-length numeric
#' vector -- the moments ([gmm_mean()], [gmm_cov()]), a tail probability
#' via [pgmm()] on a marginal, an entropy, a conditional mean, or any
#' composition of the operator calculus.
#'
#' @param ensemble A `gmm_ensemble` from [gmm_fit_ensemble()].
#' @param fn A function mapping a [gmm] to a numeric scalar or vector.
#' @param level Confidence level. Default `0.9`.
#' @param ... Forwarded to `fn`.
#'
#' @returns A data frame with one row per element of `fn`'s value:
#'   `term`, `estimate` (the base fit's value), `conf.low`, `conf.high`.
#' @family diagnostics
#' @export
#' @examples
#' fit <- fit_proxymix(banana_target(), N = 2L, regime = "kld",
#'                     is_size = 1500L, max_iter = 20L, seed = 1L)
#' ens <- gmm_fit_ensemble(fit, B = 30L, seed = 2L)
#' proxy_functional_ci(ens, function(g) gmm_mean(g)[1L])
proxy_functional_ci <- function(ensemble, fn, level = 0.9, ...) {
  if (!inherits(ensemble, "gmm_ensemble")) {
    cli::cli_abort("`ensemble` must be a {.cls gmm_ensemble} from {.fn gmm_fit_ensemble}.")
  }
  if (!is.function(fn)) {
    cli::cli_abort("`fn` must be a function of a {.cls gmm}.")
  }
  if (!is.numeric(level) || length(level) != 1L || level <= 0 || level >= 1) {
    cli::cli_abort("`level` must be a single number strictly inside (0, 1).")
  }
  est <- as.numeric(fn(ensemble$fit, ...))
  vals <- vapply(ensemble$members,
                 function(g) as.numeric(fn(g, ...)),
                 numeric(length(est)))
  vals <- matrix(vals, nrow = length(est))
  alpha <- (1 - level) / 2
  lo <- apply(vals, 1L, stats::quantile, probs = alpha, names = FALSE)
  hi <- apply(vals, 1L, stats::quantile, probs = 1 - alpha, names = FALSE)
  nms <- names(fn(ensemble$fit, ...))
  data.frame(
    term = if (!is.null(nms)) nms else paste0("f", seq_along(est)),
    estimate = est,
    conf.low = lo,
    conf.high = hi
  )
}
