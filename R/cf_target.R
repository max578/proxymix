## Constructor for a target known through its characteristic function.

#' A target given by its characteristic function
#'
#' Builds a [gmm_target] that carries only a characteristic function
#' \eqn{\varphi(t) = E \exp(i t^\top X)}. [fit_proxymix()] fits such a target
#' with regime `"cf"`, which calls [from_cf()].
#'
#' The function is checked before the target is built: it must return one
#' finite value for each row of its argument, equal 1 at zero, have modulus
#' at most 1, satisfy \eqn{\varphi(-t) = \overline{\varphi(t)}}, and fall
#' to one half along each coordinate axis between frequencies `1e-4` and
#' `1e4`. The checks use random frequencies, and the state of the random
#' number generator is left unchanged.
#'
#' @param cf A function that takes an `n` by `n_dim` matrix of frequencies,
#'   one per row, and returns the characteristic function at each row as a
#'   complex vector of length `n`, as in [from_cf()].
#' @param n_dim Dimension of the target.
#' @param name Human-readable name.
#' @param metadata Optional list of additional descriptors.
#'
#' @returns A [gmm_target] with `cf` set and no `log_density` or
#'   `samples`.
#' @family targets
#' @seealso [fit_proxymix()], [from_cf()], [from_cf_elliptical()],
#'   [from_cf_sources()], [gmm_cf_distance()].
#' @export
#' @examples
#' \dontshow{if (requireNamespace("RhpcBLASctl", quietly = TRUE)) {
#'   blas_threads <- RhpcBLASctl::blas_get_num_procs()
#'   RhpcBLASctl::blas_set_num_threads(1L)
#' }}
#' ## The sum of a Gamma(3) variable and independent N(0, 0.25) noise.
#' tgt <- cf_target(
#'   function(t) (1 - 1i * t[, 1])^(-3) * exp(-0.125 * t[, 1]^2),
#'   n_dim = 1L
#' )
#' tgt
#' fit <- fit_proxymix(tgt, N = 3L, n_freq = 1024L, seed = 1L)
#' gmm_fit_quality(fit)$cf_upper
#' \dontshow{if (exists("blas_threads")) {
#'   RhpcBLASctl::blas_set_num_threads(blas_threads)
#'   rm(blas_threads)
#' }}
cf_target <- function(cf, n_dim, name = "cf_target", metadata = list()) {
  if (!is.function(cf)) {
    cli::cli_abort("{.arg cf} must be a function of a matrix of frequencies.")
  }
  if (!is.numeric(n_dim) || length(n_dim) != 1L || is.na(n_dim) ||
      n_dim < 1 || n_dim != round(n_dim)) {
    cli::cli_abort("{.arg n_dim} must be a positive integer scalar.")
  }
  p <- as.integer(n_dim)
  phi <- .cf_checked(cf, p)
  withr::with_preserve_seed({
    th <- .cf_axis_scales(phi, p)
    .check_cf_properties(phi, p, th)
  })
  gmm_target(n_dim = p, cf = cf, name = name, metadata = metadata)
}
