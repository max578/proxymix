## Read-only accessors and tidiers, so user code never reaches into the
## S7 property layout (`g@weights`, `fit@diagnostics$...`): property paths
## are the package's private layout, and freezing them into user scripts
## would block refactors after release.

#' Component parameters of a Gaussian mixture
#'
#' Read-only accessors for the component weights, means, and covariances.
#'
#' @param g A [gmm] (or [gmm_fit]).
#'
#' @returns `gmm_weights()` returns a numeric vector of length `K`;
#'   `gmm_means()` and `gmm_covariances()` return length-`K` lists of
#'   length-`p` numeric vectors and `p` by `p` matrices respectively.
#' @family classes
#' @seealso [gmm_mean()] / [gmm_cov()] for the moments of the mixture as a
#'   whole, [gmm_fit_quality()] for the fit-quality certificate.
#' @export
#' @examples
#' g <- gmm(weights = c(0.3, 0.7), means = list(-1, 2),
#'          covariances = list(matrix(1), matrix(0.5)))
#' gmm_weights(g)
#' gmm_means(g)
gmm_weights <- function(g) {
  if (!S7::S7_inherits(g, gmm)) {
    cli::cli_abort("`g` must be a {.cls gmm} object.")
  }
  g@weights
}

#' @rdname gmm_weights
#' @export
gmm_means <- function(g) {
  if (!S7::S7_inherits(g, gmm)) {
    cli::cli_abort("`g` must be a {.cls gmm} object.")
  }
  g@means
}

#' @rdname gmm_weights
#' @export
gmm_covariances <- function(g) {
  if (!S7::S7_inherits(g, gmm)) {
    cli::cli_abort("`g` must be a {.cls gmm} object.")
  }
  g@covariances
}

## Broom-style tidiers, registered against the `generics` generics in
## `.onLoad()` when that package is installed (kept in Suggests).

#' Tidy a Gaussian mixture into a component table
#'
#' A broom-style `tidy()` method: one row per component, with the weight,
#' the mean coordinates (`mean_1`, ...), and the marginal variances
#' (`var_1`, ...). Available as `generics::tidy(g)` (or `broom::tidy(g)`)
#' when the `generics` package is installed.
#'
#' @param x A [gmm] (or [gmm_fit]).
#' @param ... Ignored, for generic compatibility.
#'
#' @returns A data frame with `K` rows.
#' @family classes
#' @name tidy.gmm
#' @examplesIf requireNamespace("generics", quietly = TRUE)
#' g <- gmm(weights = c(0.3, 0.7), means = list(c(-1, 0), c(2, 1)),
#'          covariances = list(diag(2), 0.5 * diag(2)))
#' generics::tidy(g)
NULL

.tidy_gmm <- function(x, ...) {
  K <- gmm_n_components(x)
  p <- gmm_dim(x)
  out <- data.frame(component = seq_len(K), weight = x@weights)
  mean_mat <- do.call(rbind, x@means)
  var_mat <- do.call(rbind, lapply(x@covariances, diag))
  colnames(mean_mat) <- paste0("mean_", seq_len(p))
  colnames(var_mat) <- paste0("var_", seq_len(p))
  cbind(out, mean_mat, var_mat)
}

#' Glance at a fitted Gaussian-mixture proxy
#'
#' A broom-style `glance()` method: a one-row summary of a [gmm_fit] with
#' the regime, the component count and dimension, convergence, iteration
#' count, and the regime's headline fit statistics. Available as
#' `generics::glance(fit)` when the `generics` package is installed.
#'
#' @param x A [gmm_fit].
#' @param ... Ignored, for generic compatibility.
#'
#' @returns A one-row data frame.
#' @family classes
#' @name glance.gmm_fit
#' @examplesIf requireNamespace("generics", quietly = TRUE)
#' fit <- fit_proxymix(banana_target(), N = 2L, regime = "kld",
#'                     is_size = 1000L, max_iter = 10L, seed = 1L)
#' generics::glance(fit)
NULL

.glance_gmm_fit <- function(x, ...) {
  d <- x@diagnostics
  data.frame(
    regime = x@regime,
    n_components = gmm_n_components(x),
    dim = gmm_dim(x),
    converged = x@converged,
    iterations = x@iterations,
    ess = d$ess %||% NA_real_,
    kld_final = d$kld_final %||% NA_real_,
    loglik_final = d$loglik_final %||% NA_real_,
    bic = d$bic %||% NA_real_
  )
}
