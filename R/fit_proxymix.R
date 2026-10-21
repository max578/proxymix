## The top-level fitting verb.

#' Fit a Gaussian-mixture proxy to a target distribution
#'
#' The unified front door of `proxymix`. Picks a fitting regime (or honours
#' an explicit choice) and dispatches to the corresponding regime-specific
#' fitter:
#'
#' * `"moment"` - closed-form moment matching ([fit_moment_match()]).
#' * `"sample"` - classical EM on i.i.d. samples ([fit_em_samples()]).
#' * `"kld"` - importance-sampled KLD-EM ([fit_kld_em()]) for a target
#'   that can be evaluated but not sampled.
#' * `"cf"` - a fit to the target's characteristic function ([from_cf()]),
#'   for a target known only through it ([cf_target()]).
#'
#' With `regime = "auto"` the choice is made from the shape of the supplied
#' target:
#'
#' * `N == 1` and the target carries samples or moments: `"moment"`.
#' * `N >= 2` and the target carries samples: `"sample"`.
#' * The target carries `log_density` only (no samples): `"kld"`.
#' * The target carries `cf` only (no samples, no `log_density`): `"cf"`.
#'
#' @param target A [gmm_target].
#' @param N Number of components.
#' @param regime One of `"auto"`, `"moment"`, `"sample"`, `"kld"`, `"cf"`.
#' @param ... Additional arguments forwarded to the regime-specific fitter.
#'   The most useful pass-throughs are `canonicalise` (whether to apply
#'   [gmm_canonicalise()] to the returned fit; default `TRUE`),
#'   `validation_size` and `validation_proposal` (held-out IS validation
#'   for regime `"kld"`), `max_iter`, `tol`, `n_starts`, and `seed`. For
#'   regime `"cf"` the arguments of [from_cf()] apply, among them `n_freq`,
#'   `n_start`, `seed` and `bracket`; its fit is always returned in canonical
#'   order, and `canonicalise` is ignored.
#'
#' @returns A [gmm_fit].
#' @family interface
#' @export
#' @examples
#' ## auto: samples + N=2 -> classical EM.
#' x <- matrix(stats::rnorm(200), ncol = 2)
#' tgt_s <- gmm_target_from_samples(x)
#' fit_proxymix(tgt_s, N = 2L, max_iter = 25L)
#'
#' ## explicit "kld" on a log-density-only target.
#' fit_proxymix(banana_target(), N = 3L, regime = "kld",
#'              is_size = 1000L, max_iter = 20L, seed = 1L)
#'
#' ## auto on a target known only through its characteristic function:
#' ## a Gamma(3) variable plus N(0, 0.25) noise -> "cf".
#' tgt_cf <- cf_target(
#'   function(t) (1 - 1i * t[, 1])^(-3) * exp(-0.125 * t[, 1]^2),
#'   n_dim = 1L
#' )
#' fit_proxymix(tgt_cf, N = 3L, n_freq = 1024L, seed = 1L)
fit_proxymix <- function(target,
                         N = 1L,
                         regime = c("auto", "moment", "sample", "kld", "cf"),
                         ...) {
  if (!S7::S7_inherits(target, gmm_target)) {
    cli::cli_abort("`target` must be a {.cls gmm_target} object.")
  }
  regime <- rlang::arg_match(regime)
  N <- as.integer(N)

  if (regime == "auto") {
    have_samples <- !is.null(target@samples)
    have_log_f   <- !is.null(target@log_density)
    have_moments <- !is.null(target@metadata$moments)
    have_cf      <- !is.null(target@cf)
    regime <- if (N == 1L && (have_samples || have_moments)) {
      "moment"
    } else if (N >= 2L && have_samples) {
      "sample"
    } else if (have_log_f) {
      "kld"
    } else if (have_cf) {
      "cf"
    } else {
      cli::cli_abort(c(
        "Cannot choose a regime automatically for the supplied target.",
        "i" = "Attach `samples`, `log_density` or `cf` to the target, or pass an explicit `regime`."
      ))
    }
    cli::cli_inform("Auto-selected regime: {.val {regime}}.")
  }
  if (regime == "cf" && is.null(target@cf)) {
    cli::cli_abort(c(
      "Regime {.val cf} needs a target with a characteristic function.",
      "i" = "Build the target with {.fn cf_target}."
    ))
  }

  switch(regime,
    moment = fit_moment_match(target, N = N, ...),
    sample = fit_em_samples(target, N = N, ...),
    kld    = fit_kld_em(target, N = N, ...),
    cf     = .fit_cf_target(target, N = N, ...)
  )
}

.fit_cf_target <- function(target, N, ...) {
  args <- list(...)
  args$canonicalise <- NULL
  fit <- do.call(from_cf, c(list(target@cf, target@n_dim, N), args))
  fit@target <- target
  fit
}
