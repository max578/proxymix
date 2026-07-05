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

## ---------------------------------------------------------------------------
## Support-aware proposal auto-selection (regime-(iii) robustness layer).
## ---------------------------------------------------------------------------

## Normalise a `gmm_target@support` declaration into a pair of length-`p`
## bound vectors, or NULL when the support is unconstrained. Recycles
## length-1 bounds and validates lower < upper (defence in depth alongside
## the S7 validator).
.normalise_support <- function(support, n_dim) {
  if (is.null(support)) {
    return(NULL)
  }
  lo <- rep_len(as.numeric(support$lower), n_dim)
  hi <- rep_len(as.numeric(support$upper), n_dim)
  if (any(hi <= lo)) {
    cli::cli_abort("`support$upper` must exceed `support$lower` in every coordinate.")
  }
  list(lower = lo, upper = hi)
}

## Select an importance proposal matched to a target's declared support.
##
## Returns NULL when the support is unconstrained (every coordinate spans the
## whole line), signalling the caller to keep the heavy-tailed default
## proposal -- a uniform must never truncate an unbounded target. For a
## bounded or one-sided support it returns an `is_uniform` proposal matched to
## the support:
##   * coordinates finite on both ends use the declared box, inset inward by
##     `inset` of the box width so a density that blows up or vanishes exactly
##     at the boundary cannot poison the importance weights;
##   * coordinates with an infinite end derive a finite working bound from the
##     target's samples (sample range extended by `tail_margin` of its span).
## When an infinite end cannot be resolved (no samples), it aborts with
## guidance rather than silently falling back to the crashing default.
.support_matched_proposal <- function(target, inset = 1e-3, tail_margin = 0.1) {
  sup <- .normalise_support(target@support, target@n_dim)
  if (is.null(sup)) {
    return(NULL)
  }
  lo <- sup$lower
  hi <- sup$upper
  p <- target@n_dim

  ## Fully unbounded (declared, but every coordinate spans the line): defer to
  ## the heavy-tailed default rather than truncate.
  if (all(is.infinite(lo) & lo < 0) && all(is.infinite(hi) & hi > 0)) {
    return(NULL)
  }

  needs_data <- any(is.infinite(lo) | is.infinite(hi))
  smp <- target@samples
  if (needs_data && is.null(smp)) {
    cli::cli_abort(c(
      "Cannot build a support-matched proposal for a one-sided or unbounded coordinate without samples.",
      "i" = "Attach `samples` to the target, declare a finite `support` bound, or pass an explicit {.arg proposal}."
    ))
  }

  lo_work <- numeric(p)
  hi_work <- numeric(p)
  for (j in seq_len(p)) {
    lo_j <- lo[j]
    hi_j <- hi[j]
    if (is.finite(lo_j) && is.finite(hi_j)) {
      delta <- inset * (hi_j - lo_j)
      lo_work[j] <- lo_j + delta
      hi_work[j] <- hi_j - delta
      next
    }
    col <- smp[, j]
    span <- diff(range(col))
    if (!is.finite(span) || span <= 0) {
      span <- max(abs(col), 1)
    }
    margin <- tail_margin * span
    lo_work[j] <- if (is.finite(lo_j)) lo_j + inset * span else min(col) - margin
    hi_work[j] <- if (is.finite(hi_j)) hi_j - inset * span else max(col) + margin
  }

  q <- is_uniform(n_dim = p, lower = lo_work, upper = hi_work)
  q@name <- "is_uniform[support-matched]"
  q@metadata <- c(q@metadata,
                  list(support_lower = lo, support_upper = hi,
                       data_derived = needs_data))
  q
}

# ---------------------------------------------------------------------------
# Preferred constructor names
# ---------------------------------------------------------------------------

#' Preferred names for the importance-proposal constructors
#'
#' `proposal_uniform()`, `proposal_mvn()`, and `proposal_mvt()` are the
#' preferred names of [is_uniform()], [is_mvn()], and [is_mvt()]: the
#' historical `is_*` prefix reads as a logical predicate, which these
#' constructors are not. The `is_*` names remain available as aliases and
#' are not scheduled for removal.
#'
#' @inheritParams is_uniform
#' @returns An [is_proposal] object.
#' @family proposals
#' @export
proposal_uniform <- is_uniform

#' @rdname proposal_uniform
#' @export
proposal_mvn <- is_mvn

#' @rdname proposal_uniform
#' @export
proposal_mvt <- is_mvt
