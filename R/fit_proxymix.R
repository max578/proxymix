## The top-level fitting verb.

#' Fit a Gaussian-mixture proxy to a target density
#'
#' The unified front door of `proxymix`. Picks a fitting regime (or honours
#' an explicit choice) and dispatches to the corresponding regime-specific
#' fitter:
#'
#' * `"moment"` - closed-form moment matching ([fit_moment_match()]).
#' * `"sample"` - classical EM on i.i.d. samples ([fit_em_samples()]).
#' * `"kld"` - importance-sampled KLD-EM ([fit_kld_em()]) — the wedge.
#'
#' With `regime = "auto"` the choice is made from the shape of the supplied
#' target:
#'
#' * `N == 1` and the target carries samples or moments: `"moment"`.
#' * `N >= 2` and the target carries samples: `"sample"`.
#' * The target carries `log_density` only (no samples): `"kld"`.
#'
#' @param target A [gmm_target].
#' @param N Number of components.
#' @param regime One of `"auto"`, `"moment"`, `"sample"`, `"kld"`.
#' @param ... Additional arguments forwarded to the regime-specific fitter.
#'   The most useful pass-throughs are `canonicalise` (whether to apply
#'   [gmm_canonicalise()] to the returned fit; default `TRUE`),
#'   `validation_size` and `validation_proposal` (held-out IS validation
#'   for regime `"kld"`), `max_iter`, `tol`, `n_starts`, and `seed`.
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
fit_proxymix <- function(target,
                         N = 1L,
                         regime = c("auto", "moment", "sample", "kld"),
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
    regime <- if (N == 1L && (have_samples || have_moments)) {
      "moment"
    } else if (N >= 2L && have_samples) {
      "sample"
    } else if (have_log_f) {
      "kld"
    } else {
      cli::cli_abort(c(
        "Cannot choose a regime automatically for the supplied target.",
        "i" = "Attach `samples` or `log_density` to the target, or pass an explicit `regime`."
      ))
    }
    cli::cli_inform("Auto-selected regime: {.val {regime}}.")
  }

  switch(regime,
    moment = fit_moment_match(target, N = N, ...),
    sample = fit_em_samples(target, N = N, ...),
    kld    = fit_kld_em(target, N = N, ...)
  )
}
