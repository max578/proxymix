## End-of-sample instability testing on a Gaussian-state-space filter.
##
## The filter's one-step predictive density gives, at each step, a standardised
## innovation z_t = S_t^{-1/2} (y_t - C x_{t|t-1}); under a correct, stable model
## these are iid standard normal. To ask whether the LAST m observations are
## consistent with a model fitted on the rest -- the regime where m is tiny (even
## m = 1) and m < the parameter count, so Chow / structural-break tests are
## undefined -- score the last m standardised innovations and calibrate the score.
## Two calibrations are offered: a parametric chi-square (exact when the
## innovations are Gaussian) and a distribution-free Andrews (2003) subsampling
## P-test (the score's rank among the in-sample m-blocks), which stays valid when
## the innovations are non-Gaussian. This is the end-of-sample / "next
## observation" instability test on the proxymix operator calculus.
##
## Reference: Andrews, D. W. K. (2003). End-of-Sample Instability Tests.
## Econometrica, 71(6), 1661-1694.

# ---------------------------------------------------------------------------
# Internal: standardised one-step innovations from the filter recursion.
# Runs the operator calculus (gmm_affine predict, gmm_observe update) and
# returns, per step, the squared standardised innovation z_t' z_t, which is
# chi-square on ncol(Y) degrees of freedom under the null.
# ---------------------------------------------------------------------------
.eos_innovations <- function(prior, dynamics, measurement, Y) {
  n <- nrow(Y)
  p <- gmm_dim(prior)
  z2 <- numeric(n)
  g <- prior
  for (t in seq_len(n)) {
    dyn <- .resolve_dynamics(dynamics, t, p)
    meas <- .resolve_measurement(measurement, t, p)
    if (S7::S7_inherits(dyn$Q, gmm) || S7::S7_inherits(meas$R, gmm)) {
      cli::cli_abort(c(
        "`gmm_eos_test` requires Gaussian noise: `Q` and `R` must be covariance matrices, not {.cls gmm} mixtures.",
        "i" = "Both calibrations are defined for standardised Gaussian innovations; Gaussian-sum noise makes the innovation non-Gaussian."
      ))
    }
    pred <- gmm_affine(g, dyn$A, b = dyn$b, noise_cov = dyn$Q, ridge_eps = 0)
    mu_pred <- as.numeric(pred@means[[1L]])
    p_pred <- pred@covariances[[1L]]
    e <- as.numeric(Y[t, ]) - (as.numeric(meas$C %*% mu_pred) + meas$d)
    s_mat <- meas$C %*% p_pred %*% t(meas$C) + meas$R
    z2[t] <- drop(crossprod(e, solve(s_mat, e)))
    g <- gmm_observe(pred, A = meas$C, y = as.numeric(Y[t, ]),
                     noise_cov = meas$R, b = meas$d, ridge_eps = 0)
  }
  z2
}

#' End-of-sample instability test on a Gaussian state-space filter
#'
#' Tests whether the last `m` observations of a series are consistent with a
#' linear-Gaussian state-space model fitted on the rest, in the regime where `m`
#' is small (even `m = 1`) and smaller than the parameter count -- where ordinary
#' structural-break tests (Chow, sup-Wald) are undefined because the post-break
#' parameters cannot be estimated. The statistic is the sum of the last `m`
#' squared standardised one-step innovations from the filter; a break inflates it.
#'
#' Two calibrations are offered. `method = "chisq"` refers the statistic to a
#' chi-square distribution on `m * ncol(y)` degrees of freedom, which is exact
#' when the standardised innovations are Gaussian. `method = "andrews"` is the
#' distribution-free Andrews (2003) subsampling P-test: the statistic's rank among
#' the in-sample overlapping `m`-blocks of the same innovations gives the p-value,
#' so it stays calibrated when the innovations are non-Gaussian (heavy-tailed
#' observation noise, say). The model is supplied exactly as for [gmm_filter()].
#'
#' Two finite-sample cautions. The chi-square calibration is exact when the
#' model is given; with parameters estimated on a short series it over-rejects
#' (size 0.068 at `n = 30` against a nominal 0.05 in the validation study,
#' settling to 0.044 by `n = 120`), so prefer `method = "andrews"` when the
#' model is estimated on little data. The subsampling p-value has a floor of
#' `1 / (n - 2m + 2)`, so it can reject at level 0.05 only when `n > 2m + 18`.
#'
#' @param prior A single-component [gmm] giving the state prior (the test is
#'   defined for a fitted linear-Gaussian model; multi-component priors are not
#'   yet supported).
#' @param dynamics A list with `A` (the state-transition matrix), `Q` (the
#'   process-noise covariance) and an optional offset `b`, or a function
#'   `function(t)` returning such a list, exactly as for [gmm_filter()].
#'   Gaussian-sum (mixture) process noise is rejected: the calibrations
#'   below are defined for Gaussian innovations only.
#' @param measurement A list with `C` (the observation matrix), `R` (the
#'   observation-noise covariance) and an optional offset `d`, or a function
#'   `function(t)` returning such a list, exactly as for [gmm_filter()].
#'   Gaussian-sum measurement noise is rejected for the same reason.
#' @param y A numeric vector or an `n x d` matrix of observations.
#' @param m Integer; the number of end-of-sample observations to test, `1 <= m <
#'   nrow(y)`. The tiny-`m` regime (`m = 1, 2, 3`) is the point of the test.
#' @param method Either `"chisq"` (parametric) or `"andrews"` (distribution-free
#'   subsampling P-test).
#' @param alpha The nominal level used to set `reject`.
#'
#' @returns An object of class `gmm_eos_test`: a list with the `statistic`, the
#'   `p_value`, the logical `reject`, the `method`, `m`, `alpha`, and the
#'   in-sample block statistics used for the subsampling calibration.
#'
#' @references Andrews, D. W. K. (2003). End-of-Sample Instability Tests.
#'   \emph{Econometrica}, 71(6), 1661--1694.
#' @seealso [gmm_filter()]
#' @examples
#' prior <- gmm(weights = 1, means = list(0), covariances = list(matrix(10)))
#' dyn   <- list(A = matrix(1), Q = matrix(0.04))
#' meas  <- list(C = matrix(1), R = matrix(1))
#' set.seed(1)
#' y <- c(rnorm(119), 6)                # a stable series with a jump at the end
#' gmm_eos_test(prior, dyn, meas, y, m = 1L, method = "andrews")
#' @export
gmm_eos_test <- function(prior, dynamics, measurement, y, m = 1L,
                         method = c("chisq", "andrews"), alpha = 0.05) {
  method <- match.arg(method)
  if (!S7::S7_inherits(prior, gmm)) {
    cli::cli_abort("`prior` must be a {.cls gmm}.")
  }
  if (gmm_n_components(prior) != 1L) {
    cli::cli_abort(c(
      "`prior` must be a single-component {.cls gmm}.",
      i = "End-of-sample testing is defined for a fitted linear-Gaussian model; \\
           the Gaussian-sum extension is not yet implemented."))
  }
  if (!is.numeric(alpha) || length(alpha) != 1L || alpha <= 0 || alpha >= 1) {
    cli::cli_abort("`alpha` must be a single number in (0, 1).")
  }
  ## Model specs are validated and normalised by the same resolvers as
  ## `gmm_filter()`, so offsets (`b`, `d`) and function-valued (time-varying)
  ## specs are honoured rather than silently dropped.
  meas1 <- .resolve_measurement(measurement, 1L, gmm_dim(prior))
  Y <- if (is.matrix(y)) y else matrix(as.numeric(y), ncol = 1L)
  if (ncol(Y) != meas1$m) {
    cli::cli_abort("`y` must have {meas1$m} column(s) to match `measurement$C` (got {ncol(Y)}).")
  }
  if (!all(is.finite(Y))) cli::cli_abort("`y` must be finite.")
  n <- nrow(Y)
  m <- as.integer(m)
  if (length(m) != 1L || is.na(m) || m < 1L || m >= n) {
    cli::cli_abort("`m` must be a single integer with 1 <= m < nrow(y).")
  }

  z2 <- .eos_innovations(prior, dynamics, measurement, Y)
  d_obs <- ncol(Y)
  out_idx <- (n - m + 1L):n
  statistic <- sum(z2[out_idx])
  in_sample <- z2[seq_len(n - m)]

  if (method == "chisq") {
    p_value <- stats::pchisq(statistic, df = m * d_obs, lower.tail = FALSE)
    blocks <- NULL
  } else {
    n_block <- length(in_sample) - m + 1L
    if (n_block < 1L) {
      cli::cli_abort("too few in-sample observations for the subsampling \\
                      calibration; use `method = \"chisq\"` or a smaller `m`.")
    }
    blocks <- vapply(seq_len(n_block),
                     function(i) sum(in_sample[i:(i + m - 1L)]), numeric(1L))
    p_value <- (1 + sum(blocks >= statistic)) / (1 + n_block)
  }

  structure(
    list(statistic = statistic, p_value = p_value, reject = p_value < alpha,
         method = method, m = m, alpha = alpha, df = m * d_obs,
         in_sample_blocks = blocks),
    class = "gmm_eos_test")
}

#' @export
print.gmm_eos_test <- function(x, ...) {
  cat(sprintf("End-of-sample instability test (%s), m = %d\n", x$method, x$m))
  cat(sprintf("  statistic = %.4f   p-value = %.4f   %s at alpha = %.2f\n",
              x$statistic, x$p_value,
              if (x$reject) "REJECT (instability)" else "no rejection", x$alpha))
  invisible(x)
}
