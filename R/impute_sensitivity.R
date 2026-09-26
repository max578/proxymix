## Missing-not-at-random sensitivity analysis.
##
## Sweeps the selection slope over a grid, imputes and pools at each value, and
## returns the pooled estimate with its interval, the observed-data log-likelihood
## and the fit's convergence. The grid value zero is missing at random.

#' Missing-not-at-random sensitivity analysis for a coordinate mean
#'
#' Sweeps the missing-not-at-random sensitivity slope `beta` over a grid and, at
#' each value, multiply-imputes `coord` under the selection model
#' \eqn{P(\text{missing}\mid y) = g(\alpha + \beta y)} and pools its mean by
#' Rubin's rules. The result traces how the estimate and its confidence interval
#' move as the assumed dependence of missingness on the unobserved value
#' strengthens, so an analyst can read off the value of `beta` at which a
#' conclusion would change. `beta = 0` is the missing-at-random anchor.
#'
#' The mixture is refitted at each slope. Under this selection model the observed
#' data carry information about the slope, but only through the assumed shape of
#' the outcome distribution; a different shape can fit the observed data equally
#' well at a different slope. The `loglik` column shows how strongly the data
#' favour one slope over another under the mixture. Treat a large difference as a
#' consequence of the model's shape, not as evidence about the missingness
#' mechanism itself, and report the whole curve.
#'
#' @param data A numeric matrix or data frame with `NA` in `coord` only (its
#'   other columns must be fully observed).
#' @param coord Name or index of the coordinate the mechanism acts on.
#' @param beta_grid Numeric vector of sensitivity slopes. Positive values make
#'   larger unobserved values more likely to be missing.
#' @param link Selection link, `"logit"` (the default) or `"probit"`.
#' @param N,m,seed,... Passed to [gmm_impute()]; a single `seed` makes the whole
#'   sweep reproducible and keeps the curve smooth across the grid.
#' @param max_iter Maximum EM iterations per fit, passed to [gmm_impute()].
#'   Default `500L`; fits near missing at random can need more than the
#'   `gmm_impute()` default of 100.
#'
#' @returns A data frame with one row per grid value: `beta`, `estimate`,
#'   `std.error`, `conf.low`, `conf.high`, `fmi`, `loglik` (the observed-data
#'   log-likelihood of the selection model at the fitted mixture and intercept)
#'   and `converged` (whether the point fit converged). A warning names any slope
#'   whose fit did not converge.
#' @family imputation
#' @seealso [gmm_impute()], [mnar()].
#' @export
#' @examples
#' set.seed(1)
#' x1 <- rnorm(300)
#' y <- x1 + rnorm(300)
#' y[runif(300) < plogis(-0.4 + 0.8 * y)] <- NA      # MNAR on y
#' dat <- data.frame(x1 = x1, y = y)
#' proxy_mnar_sensitivity(dat, "y", beta_grid = c(0, 0.5, 1), m = 5L, seed = 1L)
proxy_mnar_sensitivity <- function(data, coord, beta_grid = seq(0, 1, by = 0.25),
                                   link = c("logit", "probit"), N = NULL,
                                   m = 20L, seed = NULL, max_iter = 500L, ...) {
  link <- match.arg(link)
  if (!is.numeric(beta_grid) || !length(beta_grid) || any(!is.finite(beta_grid))) {
    cli::cli_abort("`beta_grid` must be a vector of finite numbers.")
  }
  var_names <- if (is.data.frame(data)) names(data) else colnames(data)
  if (is.null(var_names)) var_names <- paste0("V", seq_len(ncol(as.matrix(data))))
  col_name <- if (is.character(coord)) {
    if (!coord %in% var_names) cli::cli_abort("`coord` is not a column of `data`.")
    coord
  } else {
    var_names[as.integer(coord)]
  }

  X <- as.matrix(data)
  cj <- match(col_name, var_names)

  rows <- lapply(beta_grid, function(b) {
    imp <- gmm_impute(data, N = N, m = m,
                      mechanism = mnar(coord, beta = b, link = link),
                      seed = seed, max_iter = max_iter, ...)
    pooled <- proxy_pool(imp, col_name, method = "rubin")
    ll <- .mnar_obs_loglik(imp@point_fit, X, cj, imp@diagnostics$mnar_alpha,
                           b, link)
    data.frame(beta = b, estimate = pooled$estimate, std.error = pooled$std.error,
               conf.low = pooled$conf.low, conf.high = pooled$conf.high,
               fmi = pooled$fmi, loglik = ll,
               converged = isTRUE(imp@diagnostics$converged))
  })
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  if (!all(out$converged)) {
    cli::cli_warn(c(
      "The fit did not converge at {sum(!out$converged)} slope value{?s}: {out$beta[!out$converged]}.",
      "i" = "Increase {.arg max_iter}; estimates at those slopes are not final."
    ))
  }
  out
}
