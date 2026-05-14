## Initialisations for the EM-style fitters.
##
## Each `init_*()` returns a `gmm` carrying the initial weights, means, and
## covariances. The `multi_start_best_of()` helper drives multiple inits
## through a fitter and keeps the best.

#' Random initialisation
#'
#' Generates an `N`-component random initialisation by perturbing isotropic
#' means around a centre. Useful as one starting point in a multi-start
#' best-of strategy.
#'
#' @param N Number of components.
#' @param p Ambient dimension.
#' @param centre Length-`p` numeric vector — the location around which means
#'   are drawn.
#' @param scale Standard deviation of the mean perturbation.
#' @param sigma_diag Diagonal value used for the initial component
#'   covariances.
#' @param seed Optional integer seed.
#'
#' @returns A [gmm] of `N` components in dimension `p`.
#' @family init
#' @export
#' @examples
#' init_random(N = 3L, p = 2L, seed = 1L)
init_random <- function(N = 1L, p = 2L,
                        centre = rep(0, p),
                        scale = 1,
                        sigma_diag = 1,
                        seed = NULL) {
  N <- as.integer(N)
  p <- as.integer(p)
  draw <- function() {
    means <- lapply(seq_len(N), function(k) {
      as.numeric(centre + stats::rnorm(p, mean = 0, sd = scale))
    })
    covs <- replicate(N, diag(sigma_diag, p), simplify = FALSE)
    gmm(
      weights = rep(1 / N, N),
      means = means,
      covariances = covs,
      name = "init_random"
    )
  }
  if (is.null(seed)) draw() else withr::with_seed(seed, draw())
}

#' k-means initialisation
#'
#' Runs `stats::kmeans()` on the supplied samples and uses the resulting
#' cluster centres and within-cluster covariances to seed an EM-style
#' fitter.
#'
#' @param samples An `n` by `p` numeric matrix of samples from (or close to)
#'   the target.
#' @param N Number of components.
#' @param ridge_eps Ridge added to each cluster covariance for numerical
#'   stability when a cluster has fewer than two points.
#' @param nstart `stats::kmeans` `nstart` argument.
#'
#' @returns A [gmm] of `N` components in dimension `ncol(samples)`.
#' @family init
#' @export
#' @examples
#' x <- matrix(stats::rnorm(200), ncol = 2)
#' init_kmeans(x, N = 3L)
init_kmeans <- function(samples, N = 2L, ridge_eps = 1e-6, nstart = 10L) {
  if (!is.matrix(samples) || !is.numeric(samples)) {
    cli::cli_abort("`samples` must be a numeric matrix.")
  }
  N <- as.integer(N)
  p <- ncol(samples)
  km <- stats::kmeans(samples, centers = N, nstart = nstart)
  weights <- as.numeric(km$size / sum(km$size))
  means <- lapply(seq_len(N), function(k) as.numeric(km$centers[k, ]))
  covs <- lapply(seq_len(N), function(k) {
    idx <- which(km$cluster == k)
    if (length(idx) >= 2L) {
      ridge(stats::cov(samples[idx, , drop = FALSE]), epsilon = ridge_eps)
    } else if (length(idx) == 1L) {
      ridge(diag(1, p), epsilon = ridge_eps)
    } else {
      ridge(diag(1, p), epsilon = ridge_eps)
    }
  })
  gmm(
    weights = weights,
    means = means,
    covariances = covs,
    name = "init_kmeans"
  )
}

#' Moment-seed initialisation
#'
#' Computes the global mean and covariance of the supplied samples and
#' spreads `N` components along the leading principal direction. Useful as
#' a deterministic starting point that survives multi-modal targets better
#' than a single-Gaussian fit.
#'
#' @param samples An `n` by `p` numeric matrix.
#' @param N Number of components.
#' @param spread Multiplier on the principal-direction standard deviation
#'   used to place the component means symmetrically about the global mean.
#'
#' @returns A [gmm] of `N` components in dimension `ncol(samples)`.
#' @family init
#' @export
#' @examples
#' x <- matrix(stats::rnorm(200), ncol = 2)
#' init_moment_seed(x, N = 3L)
init_moment_seed <- function(samples, N = 2L, spread = 1.5) {
  if (!is.matrix(samples) || !is.numeric(samples)) {
    cli::cli_abort("`samples` must be a numeric matrix.")
  }
  N <- as.integer(N)
  p <- ncol(samples)
  mu_global <- colMeans(samples)
  S_global <- ridge(stats::cov(samples), epsilon = 1e-6)
  eig <- eigen(S_global, symmetric = TRUE)
  pc1 <- eig$vectors[, 1L]
  sd1 <- sqrt(eig$values[1L])
  offsets <- if (N == 1L) {
    list(rep(0, p))
  } else {
    seq(-1, 1, length.out = N) * spread * sd1
  }
  means <- lapply(offsets, function(d) as.numeric(mu_global + d * pc1))
  covs <- replicate(N, S_global, simplify = FALSE)
  gmm(
    weights = rep(1 / N, N),
    means = means,
    covariances = covs,
    name = "init_moment_seed"
  )
}

#' Warm-start initialisation from an existing fit
#'
#' Returns the input as-is. Provided as a name so that the multi-start
#' driver can include "warm starts" by symbolic name.
#'
#' @param g A [gmm] (or [gmm_fit]).
#'
#' @returns The input `g`, validated.
#' @family init
#' @export
#' @examples
#' g <- gmm(weights = 1, means = list(c(0, 0)),
#'          covariances = list(diag(2)))
#' init_warm_start(g)
init_warm_start <- function(g) {
  if (!S7::S7_inherits(g, gmm)) {
    cli::cli_abort("`g` must be a {.cls gmm} object.")
  }
  g
}

#' Multi-start best-of wrapper
#'
#' Runs the supplied fitter from each of several initialisations and
#' returns the fit with the best score, following Karlis and Xekalaki
#' (2003)'s recommendation.
#'
#' @param fit_fn A function with signature `function(init, ...)` returning
#'   a [gmm_fit].
#' @param inits A list of [gmm] initialisations.
#' @param score_fn A function `function(fit)` returning a numeric score —
#'   *larger is better* (typically the final log-target evidence).
#' @param ... Additional arguments forwarded to `fit_fn`.
#'
#' @returns The [gmm_fit] with the largest `score_fn(fit)`.
#' @family init
#' @export
#' @examples
#' x <- matrix(stats::rnorm(200), ncol = 2)
#' tgt <- gmm_target_from_samples(x)
#' inits <- list(init_random(2L, 2L, seed = 1L),
#'               init_moment_seed(x, N = 2L))
#' best <- multi_start_best_of(
#'   fit_fn   = function(init, ...) fit_em_samples(tgt, init = init, ...),
#'   inits    = inits,
#'   score_fn = function(fit) fit@diagnostics$loglik_final,
#'   max_iter = 25L
#' )
#' best@diagnostics$loglik_final
multi_start_best_of <- function(fit_fn, inits, score_fn, ...) {
  if (!is.function(fit_fn)) cli::cli_abort("`fit_fn` must be a function.")
  if (!is.function(score_fn)) cli::cli_abort("`score_fn` must be a function.")
  if (!is.list(inits) || length(inits) == 0L) {
    cli::cli_abort("`inits` must be a non-empty list of {.cls gmm} initialisations.")
  }
  fits <- lapply(inits, function(init) {
    tryCatch(fit_fn(init = init, ...), error = function(e) e)
  })
  ok <- vapply(fits, function(f) S7::S7_inherits(f, gmm_fit), logical(1L))
  if (!any(ok)) {
    cli::cli_abort("all multi-start fits failed.")
  }
  fits <- fits[ok]
  scores <- vapply(fits, function(f) {
    s <- tryCatch(score_fn(f), error = function(e) NA_real_)
    if (length(s) != 1L || !is.finite(s)) NA_real_ else as.numeric(s)
  }, FUN.VALUE = numeric(1L))
  if (all(is.na(scores))) {
    cli::cli_abort("all multi-start fits scored as `NA`.")
  }
  best_idx <- which.max(scores)
  fits[[best_idx]]
}
