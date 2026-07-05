## Bounded Gaussian-sum filtering over an observation series.
##
## The operator calculus already is a filter: `gmm_affine` is the predict step
## (affine-Gaussian propagation), `gmm_observe` is the update step, and run as a
## loop they are the Kalman filter at one component and the Gaussian-sum filter
## at many. `gmm_filter()` is the loop: it alternates predict, update and an
## optional `gmm_reduce()` over a series of observations and returns the filtered
## mixtures.
##
## Non-Gaussian process or measurement noise is supplied as a Gaussian sum -- a
## `gmm` in place of a covariance matrix. Each such mixture multiplies the
## component count by its number of components every step, so a long horizon is
## only runnable with reduction enabled; `k_max` caps the count after each step.
## The channels stay affine per component and the arithmetic stays closed form.

# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

## Mixture mean and covariance of a gmm (independent of the filter recursion).
##   mu  = sum_k w_k mu_k
##   Sig = sum_k w_k (Sig_k + mu_k mu_k') - mu mu'.
.gmm_moments <- function(g) {
  w <- g@weights
  mu <- Reduce(`+`, Map(function(wk, mk) wk * mk, w, g@means))
  m2 <- Reduce(`+`, Map(function(wk, mk, Sk) wk * (Sk + tcrossprod(mk)),
                        w, g@means, g@covariances))
  list(mean = as.numeric(mu), cov = symmetrise(m2 - tcrossprod(mu)))
}

## Concatenate a list of gmms (all in the same dimension) into one mixture,
## scaling block `b`'s weights by `exp(log_block_scale[b] - max)` before
## renormalising. Blocks with a non-finite scale (zero evidence) are dropped.
.stack_gmms <- function(parts, log_block_scale, name) {
  keep <- is.finite(log_block_scale)
  parts <- parts[keep]
  log_block_scale <- log_block_scale[keep]
  scale <- exp(log_block_scale - max(log_block_scale))
  weights <- vector("list", length(parts))
  means <- list()
  covs <- list()
  for (b in seq_along(parts)) {
    g <- parts[[b]]
    weights[[b]] <- g@weights * scale[b]
    means <- c(means, g@means)
    covs <- c(covs, g@covariances)
  }
  w <- unlist(weights)
  gmm(weights = w / sum(w), means = means, covariances = covs, name = name)
}

## Validate and normalise the dynamics list at a single step.
.resolve_dynamics <- function(dynamics, t, p) {
  dyn <- if (is.function(dynamics)) dynamics(t) else dynamics
  if (!is.list(dyn) || is.null(dyn$A)) {
    cli::cli_abort("`dynamics` must be a list with at least an `A` matrix (or a function returning one).")
  }
  A <- dyn$A
  if (!is.matrix(A) || !is.numeric(A) || nrow(A) != p || ncol(A) != p) {
    cli::cli_abort("`dynamics$A` must be a {p} x {p} numeric matrix (state-to-state).")
  }
  b <- if (is.null(dyn$b)) 0 else dyn$b
  Q <- dyn$Q
  if (!is.null(Q) && !is.matrix(Q) && !S7::S7_inherits(Q, gmm)) {
    cli::cli_abort("`dynamics$Q` must be a {p} x {p} matrix, a {.cls gmm} on R^{p} (mixture noise), or `NULL`.")
  }
  if (S7::S7_inherits(Q, gmm) && gmm_dim(Q) != p) {
    cli::cli_abort("`dynamics$Q` as a {.cls gmm} must live in R^{p} (got R^{gmm_dim(Q)}).")
  }
  list(A = A, b = b, Q = Q)
}

## Validate and normalise the measurement list at a single step.
.resolve_measurement <- function(measurement, t, p) {
  meas <- if (is.function(measurement)) measurement(t) else measurement
  if (!is.list(meas) || is.null(meas$C) || is.null(meas$R)) {
    cli::cli_abort("`measurement` must be a list with a `C` matrix and an `R` noise (or a function returning them).")
  }
  C <- meas$C
  if (!is.matrix(C) || !is.numeric(C) || ncol(C) != p) {
    cli::cli_abort("`measurement$C` must be a numeric matrix with {p} columns (observation-from-state).")
  }
  m <- nrow(C)
  R <- meas$R
  if (S7::S7_inherits(R, gmm)) {
    if (gmm_dim(R) != m) {
      cli::cli_abort("`measurement$R` as a {.cls gmm} must live in R^{m} to match `C` (got R^{gmm_dim(R)}).")
    }
  } else if (!is.matrix(R) || !is.numeric(R) || nrow(R) != m || ncol(R) != m) {
    cli::cli_abort("`measurement$R` must be a {m} x {m} numeric matrix or a {.cls gmm} on R^{m}.")
  }
  d <- if (is.null(meas$d)) 0 else meas$d
  list(C = C, R = R, d = d, m = m)
}

## Predict step: push the belief through y = A x + b + w. Gaussian process
## noise (`Q` a matrix or NULL) keeps the component count; mixture process
## noise (`Q` a gmm) multiplies it by the noise component count.
.filter_predict <- function(g, A, b, Q, ridge_eps) {
  if (S7::S7_inherits(Q, gmm)) {
    parts <- lapply(seq_len(gmm_n_components(Q)), function(j) {
      gmm_affine(g, A, b = b + Q@means[[j]], noise_cov = Q@covariances[[j]],
                 ridge_eps = ridge_eps)
    })
    .stack_gmms(parts, log(Q@weights), name = "predict")
  } else {
    out <- gmm_affine(g, A, b = b, noise_cov = Q, ridge_eps = ridge_eps)
    out@name <- "predict"
    out
  }
}

## Update step: condition the belief on y = C x + d + v. Gaussian measurement
## noise (`R` a matrix) keeps the component count; mixture measurement noise
## (`R` a gmm) multiplies it. Returns the posterior and the step's log marginal
## evidence log p(y | y_{1:t-1}).
.filter_update <- function(g, C, d, R, y, ridge_eps) {
  if (S7::S7_inherits(R, gmm)) {
    J <- gmm_n_components(R)
    parts <- vector("list", J)
    log_scale <- numeric(J)
    for (j in seq_len(J)) {
      gj <- gmm_observe(g, A = C, y = y, noise_cov = R@covariances[[j]],
                        b = d + R@means[[j]], ridge_eps = ridge_eps)
      parts[[j]] <- gj
      log_scale[j] <- if (isTRUE(gj@metadata$gmm_observe_no_update)) {
        -Inf
      } else {
        log(R@weights[j]) + gj@metadata$log_marginal_evidence
      }
    }
    if (all(!is.finite(log_scale))) {
      return(list(g = g, log_evidence = -Inf, updated = FALSE))
    }
    out <- .stack_gmms(parts, log_scale, name = "update")
    list(g = out, log_evidence = logsumexp_rows(matrix(log_scale, nrow = 1L)),
         updated = TRUE)
  } else {
    gj <- gmm_observe(g, A = C, y = y, noise_cov = R, b = d,
                      ridge_eps = ridge_eps)
    if (isTRUE(gj@metadata$gmm_observe_no_update)) {
      list(g = gj, log_evidence = -Inf, updated = FALSE)
    } else {
      list(g = gj, log_evidence = gj@metadata$log_marginal_evidence,
           updated = TRUE)
    }
  }
}

## Coerce the observation series to a list of length-m numeric vectors.
.resolve_observations <- function(y, m) {
  if (is.matrix(y)) {
    if (ncol(y) != m) {
      cli::cli_abort("`y` must have {m} column(s) to match `measurement$C` (got {ncol(y)}).")
    }
    n <- nrow(y)
    y_list <- lapply(seq_len(n), function(t) as.numeric(y[t, ]))
  } else if (is.list(y)) {
    n <- length(y)
    y_list <- lapply(y, as.numeric)
  } else if (is.numeric(y) && m == 1L) {
    n <- length(y)
    y_list <- as.list(as.numeric(y))
  } else {
    cli::cli_abort(c(
      "`y` must be an n-by-{m} numeric matrix or a length-n list of length-{m} vectors.",
      "i" = "When `m` is 1 a plain numeric vector of length n is also accepted."
    ))
  }
  if (n < 1L) {
    cli::cli_abort("`y` must contain at least one observation.")
  }
  for (t in seq_len(n)) {
    yt <- y_list[[t]]
    if (length(yt) != m || anyNA(yt) || !all(is.finite(yt))) {
      cli::cli_abort("observation {t} must be a finite length-{m} vector.")
    }
  }
  y_list
}

# ---------------------------------------------------------------------------
# gmm_filter
# ---------------------------------------------------------------------------

#' Bounded Gaussian-sum filtering over an observation series
#'
#' Runs a Gaussian-sum filter over a series of `n` observations by alternating
#' the predict operator [gmm_affine()], the update operator [gmm_observe()] and
#' an optional reduction [gmm_reduce()]. At one component this is the classical
#' Kalman filter; at several components, with Gaussian-mixture process or
#' measurement noise, it is the Gaussian-sum filter of Alspach and Sorenson
#' (1972). The filter is a pure composition of the existing closed-form
#' operators; nothing here leaves the affine-Gaussian world.
#'
#' At step \eqn{t} the belief is propagated through the linear-Gaussian dynamics
#' \eqn{x_t = A x_{t-1} + b + w}, then conditioned on the observation
#' \eqn{y_t = C x_t + d + v}, then (if `k_max` is set) reduced back to at most
#' `k_max` components. With Gaussian noise the component count is constant and
#' the recursion is the Kalman filter. Non-Gaussian noise is supplied as a
#' Gaussian sum -- a [gmm] in place of the covariance matrix `Q` or `R`. Each
#' such mixture multiplies the component count by its own component count every
#' step, so a long horizon is only runnable with reduction enabled. Reduction is
#' moment-preserving (see [gmm_reduce()]), so the filtered mean and covariance
#' are unaffected to the moment-matching order while the count stays bounded.
#'
#' @param prior A [gmm] (or [gmm_fit]) in \eqn{\mathbb{R}^p}: the initial belief.
#' @param dynamics A list `list(A =, b =, Q =)` describing the predict step, or
#'   a function `function(t)` returning such a list for time-varying dynamics.
#'   `A` is the `p`-by-`p` state-transition matrix; `b` is an optional length-`p`
#'   offset (default `0`); `Q` is the process noise -- a `p`-by-`p` covariance
#'   matrix, a [gmm] on \eqn{\mathbb{R}^p} (a Gaussian-sum process noise), or
#'   `NULL` (a deterministic predict).
#' @param measurement A list `list(C =, R =, d =)` describing the update step,
#'   or a function `function(t)` returning such a list. `C` is the `m`-by-`p`
#'   observation matrix; `R` is the measurement noise -- an `m`-by-`m`
#'   covariance matrix or a [gmm] on \eqn{\mathbb{R}^m} (a Gaussian-sum
#'   measurement noise); `d` is an optional length-`m` offset (default `0`).
#' @param y The observations: an `n`-by-`m` numeric matrix, a length-`n` list of
#'   length-`m` numeric vectors, or (when `m = 1`) a numeric vector of length
#'   `n`.
#' @param k_max Optional component cap. `NULL` (the default) disables reduction
#'   and runs the exact Kalman / Gaussian-sum filter; a positive integer caps
#'   the component count after each step via [gmm_reduce()].
#' @param reduce The reduction method passed to [gmm_reduce()] when `k_max` is
#'   set: `"merge"` (the default, a deterministic moment-preserving merge) or
#'   `"anneal"`.
#' @param ridge_eps Tiny ridge passed to [gmm_affine()] and [gmm_observe()] for
#'   numerical hygiene. Set to `0` for an exact (ridge-free) recursion.
#'
#' @returns A list with elements
#'   \describe{
#'     \item{`filtered`}{a length-`n` list of the filtered [gmm] beliefs, one
#'       per step (after predict, update and any reduction).}
#'     \item{`mean`}{an `n`-by-`p` matrix of the filtered mixture means.}
#'     \item{`cov`}{a length-`n` list of the `p`-by-`p` filtered mixture
#'       covariances.}
#'     \item{`summary`}{a data frame with one row per step: the step index, the
#'       per-coordinate filtered mean (`mean_1`, ...) and standard deviation
#'       (`sd_1`, ...), the component count `n_components`, and the step log
#'       marginal evidence `log_evidence` (whose sum over steps is the log
#'       likelihood of the series).}
#'   }
#' @family operators
#' @references Alspach, D. L. and Sorenson, H. W. (1972) Nonlinear Bayesian
#'   estimation using Gaussian sum approximations. *IEEE Transactions on
#'   Automatic Control* 17(4), 439--448. \doi{10.1109/TAC.1972.1100034}
#' @export
#' @examples
#' ## A one-dimensional local-level filter (random walk observed in noise).
#' prior <- gmm(weights = 1, means = list(0), covariances = list(matrix(1)))
#' truth <- cumsum(stats::rnorm(20, sd = 0.3))
#' y <- truth + stats::rnorm(20, sd = 0.5)
#' out <- gmm_filter(
#'   prior,
#'   dynamics    = list(A = matrix(1), Q = matrix(0.09)),
#'   measurement = list(C = matrix(1), R = matrix(0.25)),
#'   y = y
#' )
#' head(out$summary)
#'
#' ## Heavy-tailed process noise as a two-component Gaussian sum, capped at
#' ## four components per step.
#' q_noise <- gmm(weights = c(0.9, 0.1),
#'                means = list(0, 0),
#'                covariances = list(matrix(0.05), matrix(1)))
#' out2 <- gmm_filter(
#'   prior,
#'   dynamics    = list(A = matrix(1), Q = q_noise),
#'   measurement = list(C = matrix(1), R = matrix(0.25)),
#'   y = y, k_max = 4L
#' )
#' max(out2$summary$n_components)
gmm_filter <- function(prior, dynamics, measurement, y,
                       k_max = NULL, reduce = c("merge", "anneal"),
                       ridge_eps = 1e-8) {
  if (!S7::S7_inherits(prior, gmm)) {
    cli::cli_abort("`prior` must be a {.cls gmm} object.")
  }
  .check_quality(prior, "gmm_filter")
  reduce <- rlang::arg_match(reduce)
  if (!is.null(k_max)) {
    k_max <- as.integer(k_max)
    if (length(k_max) != 1L || is.na(k_max) || k_max < 1L) {
      cli::cli_abort("`k_max` must be `NULL` or a single positive integer.")
    }
  }
  p <- gmm_dim(prior)

  ## The observation dimension is fixed across the series; read it from the
  ## measurement at the first step.
  m <- .resolve_measurement(measurement, 1L, p)$m
  y_list <- .resolve_observations(y, m)
  n <- length(y_list)

  filtered <- vector("list", n)
  mean_mat <- matrix(NA_real_, nrow = n, ncol = p)
  sd_mat <- matrix(NA_real_, nrow = n, ncol = p)
  cov_list <- vector("list", n)
  ncomp <- integer(n)
  log_ev <- numeric(n)

  belief <- prior
  for (t in seq_len(n)) {
    dyn <- .resolve_dynamics(dynamics, t, p)
    meas <- .resolve_measurement(measurement, t, p)
    if (meas$m != m) {
      cli::cli_abort("the observation dimension must stay constant across the series (step {t} gives {meas$m}, expected {m}).")
    }

    belief <- .filter_predict(belief, dyn$A, dyn$b, dyn$Q, ridge_eps)
    upd <- .filter_update(belief, meas$C, meas$d, meas$R, y_list[[t]], ridge_eps)
    belief <- upd$g
    if (!is.null(k_max) && gmm_n_components(belief) > k_max) {
      belief <- gmm_reduce(belief, k_max, method = reduce)
    }

    mom <- .gmm_moments(belief)
    filtered[[t]] <- belief
    mean_mat[t, ] <- mom$mean
    sd_mat[t, ] <- sqrt(pmax(diag(mom$cov), 0))
    cov_list[[t]] <- mom$cov
    ncomp[t] <- gmm_n_components(belief)
    log_ev[t] <- upd$log_evidence
  }

  summary_df <- data.frame(step = seq_len(n))
  summary_df[paste0("mean_", seq_len(p))] <- mean_mat
  summary_df[paste0("sd_", seq_len(p))] <- sd_mat
  summary_df$n_components <- ncomp
  summary_df$log_evidence <- log_ev

  list(filtered = filtered, mean = mean_mat, cov = cov_list,
       summary = summary_df)
}
