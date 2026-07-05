## The normalising constant of the target, estimated with the fitted proxy
## as the importance proposal. In regime (iii) the user hands over an
## unnormalised density; once a proxy is fitted, log Z is estimable at the
## cost of n proxy draws and n target evaluations.

#' Estimate the target's normalising constant from a fitted proxy
#'
#' Importance-sampling estimate of \eqn{Z = \int f(x)\, dx} for the fit's
#' target \eqn{f}, using the fitted mixture \eqn{\hat g} as the proposal:
#' \deqn{\widehat{Z} = \frac{1}{n} \sum_{i=1}^n
#'   \frac{f(x_i)}{\hat g(x_i)}, \qquad x_i \sim \hat g,}
#' computed in the log domain. For a Bayesian posterior handed over as
#' `likelihood x prior`, \eqn{\log Z} is the log marginal likelihood, so a
#' fitted proxy doubles as a model-comparison device.
#'
#' The estimator is exact in expectation for any proposal that dominates
#' \eqn{f}, and its Monte Carlo error is driven by how well \eqn{\hat g}
#' matches \eqn{f} -- which is precisely what the fit optimised. The
#' variance is finite only when \eqn{\hat g} has tails at least as heavy
#' as \eqn{f}; a right-tail diagnostic is returned (the effective sample
#' size and the share of the estimate carried by the largest ten percent
#' of weights, which sits near 0.10 for a well-matched proxy), and a
#' classed warning (`proxymix_heavy_tail`) is raised when it indicates an
#' untrustworthy tail. Results also report the delta-method standard error
#' of \eqn{\log \widehat{Z}}.
#'
#' When the target declares itself normalised (`normalised = TRUE`), the
#' true value is \eqn{\log Z = 0} and the function still estimates it --
#' a useful end-to-end diagnostic of the fit.
#'
#' @param fit A [gmm_fit] whose target carries a `log_density`.
#' @param n Number of evidence draws from the fitted proxy.
#' @param seed Optional integer seed for the evidence draw.
#'
#' @returns A list of class `proxymix_evidence` with elements `log_z`,
#'   `se_log_z`, `n`, `ess`, `max_weight_share`, `top_decile_share`, and
#'   `flagged` (the heavy-tail indicator).
#' @family diagnostics
#' @export
#' @examples
#' ## An unnormalised target with a known constant: log f = log N(., 0, I) + 3.
#' tgt <- gmm_target(
#'   n_dim = 2L,
#'   log_density = function(x) {
#'     if (is.null(dim(x))) x <- matrix(x, ncol = 2L)
#'     -0.5 * rowSums(x^2) - log(2 * pi) + 3
#'   },
#'   normalised = FALSE, name = "shifted_gaussian"
#' )
#' fit <- fit_kld_em(tgt, N = 1L, is_size = 2000L, max_iter = 40L, seed = 1L)
#' ev <- gmm_evidence(fit, n = 2000L, seed = 2L)
#' ev$log_z   # close to 3
gmm_evidence <- function(fit, n = 4000L, seed = NULL) {
  if (!S7::S7_inherits(fit, gmm_fit)) {
    cli::cli_abort("`fit` must be a {.cls gmm_fit}.")
  }
  tgt <- fit@target
  if (is.null(tgt) || is.null(tgt@log_density)) {
    cli::cli_abort("`fit@target@log_density` must be supplied to estimate the normalising constant.")
  }
  n <- as.integer(n)
  if (length(n) != 1L || is.na(n) || n < 10L) {
    cli::cli_abort("`n` must be a single integer of at least 10.")
  }
  x <- if (is.null(seed)) rgmm(n, fit) else withr::with_seed(seed, rgmm(n, fit))
  log_f <- tgt@log_density(x)
  log_g <- dgmm(x, fit, log = TRUE)
  log_w <- log_f - log_g
  finite <- is.finite(log_w)
  if (!any(finite)) {
    cli::cli_abort("every evidence draw received a non-finite weight; the proxy does not overlap the target.")
  }
  lw <- log_w[finite]
  m <- max(lw)
  a <- exp(lw - m)
  ## log Z-hat = m + log(mean(a)); delta-method SE of log Z-hat.
  log_z <- m + log(mean(a))
  se_log_z <- stats::sd(a) / (sqrt(length(a)) * mean(a))
  W <- a / sum(a)
  ess <- 1 / sum(W^2)
  max_share <- max(W)
  ## Tail diagnostic: the share of the estimate carried by the largest 10%
  ## of weights. Uniform weights give 0.10; a well-matched proxy sits near
  ## that; a proxy with too-light tails concentrates the estimate in a few
  ## draws long before the Kish ESS collapses (the ESS is a weak detector
  ## of infinite-variance importance sampling).
  top_decile_share <- sum(sort(W, decreasing = TRUE)[
    seq_len(ceiling(0.1 * length(W)))])
  flagged <- (ess / length(a)) < 0.05 || top_decile_share > 0.4
  if (flagged) {
    cli::cli_warn(c(
      "The evidence weights are heavy-tailed (relative ESS = {signif(ess / length(a), 2)}, top-decile share = {signif(top_decile_share, 2)}).",
      "i" = "The proxy's tails may be lighter than the target's; the estimate (and especially its standard error) may be unreliable."
    ), class = "proxymix_heavy_tail")
  }
  structure(
    list(log_z = log_z, se_log_z = se_log_z, n = length(a),
         ess = ess, max_weight_share = max_share,
         top_decile_share = top_decile_share, flagged = flagged),
    class = "proxymix_evidence"
  )
}

#' @export
print.proxymix_evidence <- function(x, ...) {
  cat(sprintf("log Z = %.4f  (SE %.4f, n = %d, ESS = %.0f%s)\n",
              x$log_z, x$se_log_z, x$n, x$ess,
              if (x$flagged) ", heavy-tailed: interpret with care" else ""))
  invisible(x)
}
