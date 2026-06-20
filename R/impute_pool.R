## Pooling a column mean across the completed datasets.
##
## proxymix pools the one estimand whose between-imputation variance it knows
## in closed form: a column mean. Because the per-row imputation law is a
## Gaussian-mixture conditional, the imputation variance of the mean is exact
## -- the m = infinity limit of Rubin's between-imputation term, with no
## Monte-Carlo noise and an imputation / parameter variance split. For
## regression coefficients and other model estimands, hand the completions to
## mice with `as_mids()` and pool with `mice::pool()`; there is no reason to
## re-implement Rubin's rules here.

# ---------------------------------------------------------------------------
# Rubin's rules (for the column-mean comparison only)
# ---------------------------------------------------------------------------

## Combine per-imputation point estimates Q and variances U (length m).
.rubin_combine <- function(Q, U, m) {
  qbar <- mean(Q)
  ubar <- mean(U)
  B <- if (m > 1L) stats::var(Q) else 0
  total <- ubar + (1 + 1 / m) * B
  df <- if (B <= 0) Inf else {
    r <- (1 + 1 / m) * B / ubar
    (m - 1) * (1 + 1 / r)^2
  }
  fmi <- if (total > 0) (1 + 1 / m) * B / total else 0
  list(estimate = qbar, var = total, df = df, fmi = fmi)
}

## Per-row conditional mean and variance of column `col` given the observed
## coordinates of each row, under mixture `g`. Observed rows contribute the
## value with zero variance.
.col_cond_moments <- function(g, X, obs, col) {
  n <- nrow(X)
  K <- length(g@weights)
  cmean <- numeric(n)
  cvar <- numeric(n)
  obs_col <- obs[, col]
  cmean[obs_col] <- X[obs_col, col]
  miss <- which(!obs_col)
  if (!length(miss)) return(list(mean = cmean, var = cvar))
  pr <- .pattern_rows(obs[miss, , drop = FALSE])
  for (key in names(pr)) {
    local <- pr[[key]]
    rows <- miss[local]
    o <- .pattern_observed(key)
    nr <- length(rows)
    logw <- matrix(0, nr, K)
    mk <- matrix(0, nr, K)
    vk <- numeric(K)
    for (k in seq_len(K)) {
      if (length(o)) {
        Xo <- X[rows, o, drop = FALSE]
        Sooi <- chol2inv(chol(g@covariances[[k]][o, o, drop = FALSE]))
        reg <- g@covariances[[k]][col, o, drop = FALSE] %*% Sooi
        mk[, k] <- as.numeric(g@means[[k]][col] +
          sweep(Xo, 2L, g@means[[k]][o]) %*% t(reg))
        vk[k] <- g@covariances[[k]][col, col] -
          as.numeric(reg %*% g@covariances[[k]][o, col, drop = FALSE])
        logw[, k] <- log(g@weights[k]) +
          .dmvn_rows(Xo, g@means[[k]][o], g@covariances[[k]][o, o, drop = FALSE])
      } else {
        mk[, k] <- g@means[[k]][col]
        vk[k] <- g@covariances[[k]][col, col]
        logw[, k] <- log(g@weights[k])
      }
    }
    gamma <- exp(logw - .logsumexp_rows(logw))
    mu_row <- rowSums(gamma * mk)
    var_row <- rowSums(gamma * (sweep(mk^2, 2L, vk, "+"))) - mu_row^2
    cmean[rows] <- mu_row
    cvar[rows] <- pmax(var_row, 0)
  }
  list(mean = cmean, var = cvar)
}

## Analytic (m = infinity) pooling of a column mean. Returns the Rubin
## quantities plus the imputation / parameter variance split.
.analytic_mean <- function(object, col) {
  X <- object@data
  obs <- object@observed
  n <- nrow(X)
  m <- object@m
  miss <- !obs[, col]
  per <- vapply(object@fits, function(g) {
    cm <- .col_cond_moments(g, X, obs, col)
    qb <- mean(cm$mean)
    bimp <- sum(cm$var[miss]) / n^2
    ss <- sum((cm$mean - qb)^2) + sum(cm$var[miss])
    c(qb, bimp, ss / ((n - 1) * n))
  }, numeric(3))
  qbar <- mean(per[1, ])
  bimp <- mean(per[2, ])
  ubar <- mean(per[3, ])
  vparam <- (1 + 1 / m) * stats::var(per[1, ])
  total <- ubar + bimp + vparam
  miss_var <- bimp + vparam
  df <- if (miss_var <= 0) Inf else (m - 1) * (total / vparam)^2
  list(
    estimate = qbar, var = total, df = df,
    fmi = if (total > 0) miss_var / total else 0,
    components = c(complete = ubar, imputation = bimp, parameter = vparam)
  )
}

.pool_row <- function(term, comb) {
  se <- sqrt(comb$var)
  tcrit <- stats::qt(0.975, comb$df)
  data.frame(
    term = term, estimate = comb$estimate, std.error = se,
    statistic = comb$estimate / se, df = comb$df,
    conf.low = comb$estimate - tcrit * se,
    conf.high = comb$estimate + tcrit * se,
    fmi = comb$fmi, stringsAsFactors = FALSE
  )
}

# ---------------------------------------------------------------------------
# proxy_pool
# ---------------------------------------------------------------------------

#' Pool a column mean across imputations
#'
#' Pools the mean of one column over the `m` completed datasets in a
#' [gmm_imputation], returning the estimate with a standard error, degrees of
#' freedom, confidence interval, and fraction of missing information. The
#' default `method = "analytic"` computes the between-imputation variance in
#' closed form -- the exact \eqn{m \to \infty} limit, with no Monte-Carlo
#' noise -- from the mixture conditional, and splits the total variance into
#' complete-data, imputation, and parameter parts. `method = "rubin"` instead
#' applies Rubin's rules to the drawn completions (useful as a check).
#'
#' For a regression or any other model estimand, do not pool here: convert the
#' imputations to a `mice` object with [as_mids()] and pool with
#' [mice::pool()], which is the established workflow and reports the same
#' diagnostics.
#'
#' @param object A [gmm_imputation].
#' @param column Name of a single numeric column whose mean is pooled.
#' @param method `"analytic"` (the default) for the closed-form pooling, or
#'   `"rubin"` for Rubin's rules over the drawn completions.
#'
#' @returns A one-row data frame: `term`, `estimate`, `std.error`,
#'   `statistic`, `df`, `conf.low`, `conf.high`, `fmi`.
#' @family imputation
#' @seealso [proxy_fmi()], [as_mids()] to pool models with [mice::pool()].
#' @export
#' @examples
#' set.seed(1)
#' x1 <- rnorm(150); x2 <- x1 + rnorm(150)
#' x2[runif(150) < 0.3] <- NA
#' imp <- gmm_impute(cbind(x1, x2), N = 1L, m = 10L, seed = 1L)
#' proxy_pool(imp, "x2")                       # analytic column mean
proxy_pool <- function(object, column, method = c("analytic", "rubin")) {
  if (!S7::S7_inherits(object, gmm_imputation)) {
    cli::cli_abort("`object` must be a {.cls gmm_imputation}.")
  }
  method <- match.arg(method)
  if (!is.character(column) || length(column) != 1L ||
      !column %in% object@var_names) {
    cli::cli_abort(c(
      "`column` must name a single column of the data.",
      "i" = "For a model estimand, pool with mice: {.code mice::pool(with(as_mids(imp), lm(...)))}."
    ))
  }
  col <- match(column, object@var_names)
  if (method == "analytic") {
    out <- .pool_row(column, .analytic_mean(object, col))
    attr(out, "method") <- "analytic"
    return(out)
  }
  QU <- vapply(object@completions, function(cm)
    c(mean(cm[, col]), stats::var(cm[, col]) / nrow(cm)), numeric(2))
  out <- .pool_row(column, .rubin_combine(QU[1, ], QU[2, ], object@m))
  attr(out, "method") <- "rubin"
  out
}

#' Fraction of missing information for a column mean
#'
#' The share of a column mean's total variance attributable to the missing
#' data, read from [proxy_pool()].
#'
#' @inheritParams proxy_pool
#'
#' @returns A named numeric scalar.
#' @family imputation
#' @seealso [proxy_pool()].
#' @export
#' @examples
#' set.seed(1)
#' x1 <- rnorm(150); x2 <- x1 + rnorm(150); x2[runif(150) < 0.3] <- NA
#' imp <- gmm_impute(cbind(x1, x2), N = 1L, m = 10L, seed = 1L)
#' proxy_fmi(imp, "x2")
proxy_fmi <- function(object, column, method = c("analytic", "rubin")) {
  pooled <- proxy_pool(object, column, method = match.arg(method))
  stats::setNames(pooled$fmi, pooled$term)
}
