## Missing-not-at-random sensitivity analysis.
##
## A missing-not-at-random departure is not identified from the observed data, so
## the responsible analysis is not a single estimate but a curve: how does the
## conclusion move as the assumed strength of the mechanism varies? This sweeps
## the sensitivity slope over a grid, imputes and pools at each value, and returns
## the pooled estimate with its interval -- the input to a tipping-point or
## tornado display. The grid value zero is missing at random.

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
#' The slope is a sensitivity parameter, not an estimate: the data do not identify
#' it. Report the curve, not a single point.
#'
#' @param data A numeric matrix or data frame with `NA` in `coord` only (its
#'   other columns must be fully observed).
#' @param coord Name or index of the coordinate the mechanism acts on.
#' @param beta_grid Numeric vector of sensitivity slopes. Positive values make
#'   larger unobserved values more likely to be missing.
#' @param link Selection link, `"logit"` (the default) or `"probit"`.
#' @param N,m,seed,... Passed to [gmm_impute()]; a single `seed` makes the whole
#'   sweep reproducible and keeps the curve smooth across the grid.
#'
#' @returns A data frame with one row per grid value: `beta`, `estimate`,
#'   `std.error`, `conf.low`, `conf.high`, `fmi`.
#' @family imputation
#' @seealso [gmm_impute()], [mnar()].
#' @export
#' @examples
#' set.seed(1)
#' x1 <- rnorm(300)
#' y <- x1 + rnorm(300)
#' y[runif(300) < plogis(-0.4 + 0.8 * y)] <- NA      # MNAR on y
#' dat <- data.frame(x1 = x1, y = y)
#' proxy_mnar_sensitivity(dat, "y", beta_grid = c(0, 0.4, 0.8, 1.2), m = 10L, seed = 1L)
proxy_mnar_sensitivity <- function(data, coord, beta_grid = seq(0, 1, by = 0.25),
                                   link = c("logit", "probit"), N = NULL,
                                   m = 20L, seed = NULL, ...) {
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

  rows <- lapply(beta_grid, function(b) {
    imp <- gmm_impute(data, N = N, m = m,
                      mechanism = mnar(coord, beta = b, link = link),
                      seed = seed, ...)
    pooled <- proxy_pool(imp, col_name, method = "rubin")
    data.frame(beta = b, estimate = pooled$estimate, std.error = pooled$std.error,
               conf.low = pooled$conf.low, conf.high = pooled$conf.high,
               fmi = pooled$fmi)
  })
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}
