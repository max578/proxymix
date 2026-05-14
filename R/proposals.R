## Importance-sampling proposals for regime (iii).

#' Uniform-on-a-box proposal
#'
#' Builds an [is_proposal] that samples uniformly on the hyperrectangle
#' \eqn{[\text{lower}_1, \text{upper}_1] \times \cdots \times [\text{lower}_p, \text{upper}_p]}
#' and reports the (constant) log-density on that box. Outside the box the
#' log-density is `-Inf`.
#'
#' @param n_dim Ambient dimension `p`.
#' @param lower Length-`p` numeric vector of lower bounds (recycled from a
#'   single value).
#' @param upper Length-`p` numeric vector of upper bounds (recycled from a
#'   single value).
#'
#' @returns An [is_proposal] object.
#' @family proposals
#' @export
#' @examples
#' q <- is_uniform(n_dim = 2L, lower = -5, upper = 5)
#' q
#' q@sample(3L)
is_uniform <- function(n_dim, lower = -1, upper = 1) {
  n_dim <- as.integer(n_dim)
  lower <- as.numeric(rep_len(lower, n_dim))
  upper <- as.numeric(rep_len(upper, n_dim))
  if (any(upper <= lower)) {
    cli::cli_abort("`upper` must be strictly greater than `lower` in each coordinate.")
  }
  log_vol <- sum(log(upper - lower))
  sampler <- function(n) {
    n <- as.integer(n)
    matrix(
      stats::runif(n * n_dim,
                   min = rep(lower, each = n),
                   max = rep(upper, each = n)),
      nrow = n, ncol = n_dim
    )
  }
  log_dens <- function(x) {
    if (is.null(dim(x))) x <- matrix(x, nrow = 1L)
    inside <- rowSums(
      (x >= matrix(lower, nrow = nrow(x), ncol = n_dim, byrow = TRUE)) &
        (x <= matrix(upper, nrow = nrow(x), ncol = n_dim, byrow = TRUE))
    ) == n_dim
    out <- rep(-Inf, nrow(x))
    out[inside] <- -log_vol
    out
  }
  is_proposal(
    n_dim = n_dim,
    sample = sampler,
    log_density = log_dens,
    name = "is_uniform",
    metadata = list(lower = lower, upper = upper)
  )
}

#' Multivariate-normal proposal
#'
#' Builds an [is_proposal] using a multivariate-normal `N(mean, cov)`
#' density and sampler.
#'
#' @param n_dim Ambient dimension `p`.
#' @param mean Length-`p` numeric mean vector. Defaults to the zero vector.
#' @param cov A `p`-by-`p` symmetric positive-definite covariance matrix.
#'   Defaults to the identity.
#'
#' @returns An [is_proposal] object.
#' @family proposals
#' @export
#' @examples
#' q <- is_mvn(n_dim = 2L, mean = c(0, 0), cov = 4 * diag(2))
#' q
is_mvn <- function(n_dim, mean = rep(0, n_dim), cov = diag(n_dim)) {
  n_dim <- as.integer(n_dim)
  mean <- as.numeric(mean)
  if (length(mean) != n_dim) {
    cli::cli_abort("`mean` must have length {n_dim}.")
  }
  if (!is.matrix(cov) || nrow(cov) != n_dim || ncol(cov) != n_dim) {
    cli::cli_abort("`cov` must be a {n_dim} x {n_dim} matrix.")
  }
  sampler <- function(n) {
    mvnfast::rmvn(as.integer(n), mu = mean, sigma = cov)
  }
  log_dens <- function(x) {
    if (is.null(dim(x))) x <- matrix(x, nrow = 1L)
    mvnfast::dmvn(x, mu = mean, sigma = cov, log = TRUE)
  }
  is_proposal(
    n_dim = n_dim,
    sample = sampler,
    log_density = log_dens,
    name = "is_mvn",
    metadata = list(mean = mean, cov = cov)
  )
}

#' Multivariate-t proposal
#'
#' Builds an [is_proposal] using a multivariate-Student-t density and
#' sampler with `df` degrees of freedom, location `mean`, and scale matrix
#' `sigma`. Heavier tails than `is_mvn()`, so often a safer importance
#' proposal at moderate dimensions.
#'
#' @param n_dim Ambient dimension `p`.
#' @param mean Length-`p` numeric location vector. Defaults to the zero
#'   vector.
#' @param sigma A `p`-by-`p` symmetric positive-definite scale matrix.
#'   Defaults to the identity.
#' @param df Degrees of freedom (`df > 2` recommended for finite variance).
#'
#' @returns An [is_proposal] object.
#' @family proposals
#' @export
#' @examples
#' q <- is_mvt(n_dim = 2L, df = 5)
#' q
is_mvt <- function(n_dim, mean = rep(0, n_dim), sigma = diag(n_dim), df = 5) {
  n_dim <- as.integer(n_dim)
  mean <- as.numeric(mean)
  if (length(mean) != n_dim) {
    cli::cli_abort("`mean` must have length {n_dim}.")
  }
  if (!is.matrix(sigma) || nrow(sigma) != n_dim || ncol(sigma) != n_dim) {
    cli::cli_abort("`sigma` must be a {n_dim} x {n_dim} matrix.")
  }
  if (length(df) != 1L || !is.numeric(df) || df <= 0) {
    cli::cli_abort("`df` must be a positive scalar.")
  }
  df <- as.numeric(df)
  sampler <- function(n) {
    mvnfast::rmvt(as.integer(n), mu = mean, sigma = sigma, df = df)
  }
  log_dens <- function(x) {
    if (is.null(dim(x))) x <- matrix(x, nrow = 1L)
    mvnfast::dmvt(x, mu = mean, sigma = sigma, df = df, log = TRUE)
  }
  is_proposal(
    n_dim = n_dim,
    sample = sampler,
    log_density = log_dens,
    name = "is_mvt",
    metadata = list(mean = mean, sigma = sigma, df = df)
  )
}
