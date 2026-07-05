## Affine-Gaussian operator calculus on Gaussian mixtures.
##
## Added in v0.3.0. Closed-form pushforward, conditioning, aggregation,
## and missing-coordinate updates. Mathematically: the finite-mixture
## analogue of the Kalman update applied component-wise with weight
## reweighting by the per-component marginal evidence.
##
## Affine channels and Gaussian noise only -- no silent non-affine fallbacks.

## Validation helpers --------------------------------------------------------

.validate_A <- function(A, p, m_arg = NULL) {
  if (!is.matrix(A) || !is.numeric(A)) {
    cli::cli_abort("`A` must be a numeric matrix.")
  }
  if (anyNA(A) || !all(is.finite(A))) {
    cli::cli_abort("`A` must be finite (no NA / Inf entries).")
  }
  if (ncol(A) != p) {
    cli::cli_abort("`A` must have {p} columns (got {ncol(A)}).")
  }
  if (!is.null(m_arg) && nrow(A) != m_arg) {
    cli::cli_abort("`A` must have {m_arg} rows (got {nrow(A)}).")
  }
  invisible(A)
}

.validate_b <- function(b, m) {
  if (length(b) == 1L) {
    b <- rep_len(as.numeric(b), m)
  }
  if (!is.numeric(b) || length(b) != m) {
    cli::cli_abort("`b` must be a numeric scalar or length-{m} vector.")
  }
  if (anyNA(b) || !all(is.finite(b))) {
    cli::cli_abort("`b` must be finite.")
  }
  as.numeric(b)
}

.validate_noise_cov <- function(noise_cov, m, allow_null = TRUE) {
  if (is.null(noise_cov)) {
    if (!allow_null) {
      cli::cli_abort("`noise_cov` is required for this operator (got `NULL`).")
    }
    return(matrix(0, nrow = m, ncol = m))
  }
  if (!is.matrix(noise_cov) || !is.numeric(noise_cov)) {
    cli::cli_abort("`noise_cov` must be a numeric matrix.")
  }
  if (nrow(noise_cov) != m || ncol(noise_cov) != m) {
    cli::cli_abort("`noise_cov` must be a {m} x {m} matrix.")
  }
  noise_cov <- symmetrise(noise_cov)
  if (anyNA(noise_cov) || !all(is.finite(noise_cov))) {
    cli::cli_abort("`noise_cov` must be finite.")
  }
  noise_cov
}

## ---------------------------------------------------------------------------
## gmm_affine — pushforward through y = A x + b + epsilon
## ---------------------------------------------------------------------------

#' Affine pushforward of a Gaussian mixture
#'
#' Returns the (closed-form) distribution of \eqn{Y = A X + b + \epsilon}
#' when \eqn{X \sim g} is a Gaussian mixture and \eqn{\epsilon \sim
#' \mathcal{N}(0, R)} is independent additive Gaussian noise.
#'
#' For each component `k`, the pushed-forward parameters are
#' \deqn{\mu'_k = A \mu_k + b, \qquad \Sigma'_k = A \Sigma_k A^\top + R,}
#' and the mixture weights are unchanged. This is the finite-mixture
#' analogue of a Kalman-style predict step.
#'
#' The channel is required to be **affine** in `x` and the noise is
#' required to be Gaussian. Non-linear channels are not closed form and
#' are not silently approximated: push samples through the map instead
#' (`rgmm()` then the transform) and refit with [fit_em_samples()] when a
#' mixture of the image is needed.
#'
#' @param g A [gmm] (or [gmm_fit]) in `R^p`.
#' @param A An `m` by `p` numeric matrix.
#' @param b Numeric scalar or length-`m` vector. Default `0`.
#' @param noise_cov `m` by `m` SPD numeric matrix, or `NULL` (treated as
#'   the zero matrix — a deterministic channel).
#' @param ridge_eps Tiny ridge added to the output covariances for
#'   numerical hygiene. Set to `0` to disable.
#'
#' @returns A [gmm] in `R^m` with the same number of components and the
#'   same weights as `g`.
#' @family operators
#' @export
#' @examples
#' g <- gmm(weights = c(0.5, 0.5),
#'          means = list(c(-1, 0), c(1, 0)),
#'          covariances = list(diag(2), diag(2)))
#' A <- matrix(c(1, 0, 0, 1, 1, -1), nrow = 3L, byrow = TRUE)
#' gmm_affine(g, A, b = c(0, 0, 0), noise_cov = 0.01 * diag(3))
gmm_affine <- function(g, A, b = 0, noise_cov = NULL, ridge_eps = 1e-6) {
  if (!S7::S7_inherits(g, gmm)) {
    cli::cli_abort("`g` must be a {.cls gmm} object.")
  }
  p <- gmm_dim(g)
  A <- .validate_A(A, p)
  m <- nrow(A)
  b <- .validate_b(b, m)
  R <- .validate_noise_cov(noise_cov, m, allow_null = TRUE)

  K <- gmm_n_components(g)
  out_means <- vector("list", K)
  out_covs  <- vector("list", K)
  for (k in seq_len(K)) {
    out_means[[k]] <- as.numeric(A %*% g@means[[k]] + b)
    S_k <- A %*% g@covariances[[k]] %*% t(A) + R
    out_covs[[k]] <- if (ridge_eps > 0) {
      ridge(symmetrise(S_k), ridge_eps)
    } else {
      symmetrise(S_k)
    }
  }

  gmm(
    weights     = g@weights,
    means       = out_means,
    covariances = out_covs,
    name        = .op_name("affine", g),
    metadata    = .op_result_meta(g, "affine",
                                  list(input_dim = p, output_dim = m))
  )
}

## ---------------------------------------------------------------------------
## gmm_observe — Kalman update on a Gaussian mixture
## ---------------------------------------------------------------------------

#' Bayesian update of a Gaussian mixture on a noisy linear observation
#'
#' Conditions a Gaussian mixture `g` on a single noisy linear observation
#' \eqn{y = A x + b + \epsilon}, \eqn{\epsilon \sim \mathcal{N}(0, R)}.
#' Per component, applies the Kalman gain
#' \deqn{K_k = \Sigma_k A^\top S_k^{-1}, \qquad S_k = A \Sigma_k A^\top + R,}
#' and updates
#' \deqn{\mu'_k = \mu_k + K_k (y - A \mu_k - b), \qquad
#'       \Sigma'_k = (I - K_k A) \Sigma_k.}
#' Component weights are multiplied by the marginal evidence
#' \eqn{\pi_k \mathcal{N}(y; A \mu_k + b, S_k)} and renormalised. This is
#' the finite-mixture analogue of a Kalman update step.
#'
#' If the marginal evidence vanishes at every component (e.g. `y` is
#' many standard deviations from every component), the function issues a
#' warning and returns `g` unchanged with `metadata$gmm_observe_no_update =
#' TRUE`.
#'
#' @param g A [gmm] (or [gmm_fit]) in `R^p`.
#' @param A An `m` by `p` numeric matrix.
#' @param y A length-`m` numeric vector (the observation).
#' @param noise_cov An `m` by `m` SPD numeric matrix (the observation
#'   noise covariance `R`). Required.
#' @param b Numeric scalar or length-`m` vector. Default `0`.
#' @param ridge_eps Tiny ridge added to updated covariances for numerical
#'   hygiene. Set to `0` to disable.
#'
#' @returns A [gmm] in `R^p` with the same number of components and the
#'   reweighted component weights.
#' @family operators
#' @export
#' @examples
#' g <- gmm(weights = c(0.5, 0.5),
#'          means = list(c(-1, 0), c(1, 0)),
#'          covariances = list(diag(2), diag(2)))
#' A <- matrix(c(1, 0), nrow = 1L)
#' gmm_observe(g, A = A, y = 0.8, noise_cov = matrix(0.25, 1, 1))
gmm_observe <- function(g, A, y, noise_cov, b = 0, ridge_eps = 1e-6) {
  if (!S7::S7_inherits(g, gmm)) {
    cli::cli_abort("`g` must be a {.cls gmm} object.")
  }
  p <- gmm_dim(g)
  A <- .validate_A(A, p)
  m <- nrow(A)
  b <- .validate_b(b, m)
  R <- .validate_noise_cov(noise_cov, m, allow_null = FALSE)

  if (!is.numeric(y) || length(y) != m) {
    cli::cli_abort("`y` must be a length-{m} numeric vector.")
  }
  if (anyNA(y) || !all(is.finite(y))) {
    cli::cli_abort("`y` must be finite.")
  }

  K <- gmm_n_components(g)
  innov <- y - b
  new_means <- vector("list", K)
  new_covs  <- vector("list", K)
  log_evidence <- numeric(K)
  for (k in seq_len(K)) {
    mu_k <- g@means[[k]]
    S_k  <- g@covariances[[k]]
    Sx   <- A %*% S_k %*% t(A) + R
    Sx   <- symmetrise(Sx)
    L_S  <- tryCatch(chol(Sx), error = function(e) NULL)
    if (is.null(L_S)) {
      Sx   <- ridge(Sx, max(ridge_eps, 1e-8))
      L_S  <- chol(Sx)
    }
    Sx_inv <- chol2inv(L_S)
    gain   <- S_k %*% t(A) %*% Sx_inv
    resid  <- innov - as.numeric(A %*% mu_k)
    new_means[[k]] <- as.numeric(mu_k + gain %*% resid)
    upd_cov <- S_k - gain %*% A %*% S_k
    upd_cov <- symmetrise(upd_cov)
    if (ridge_eps > 0) upd_cov <- ridge(upd_cov, ridge_eps)
    new_covs[[k]] <- upd_cov
    log_evidence[k] <- mvnfast::dmvn(
      matrix(innov, nrow = 1L),
      mu = as.numeric(A %*% mu_k), sigma = Sx, log = TRUE
    )
  }

  log_new <- log(g@weights) + log_evidence
  mx <- max(log_new)
  ## Catastrophic case: every component sees the observation at numerical
  ## density zero (log_evidence below ~-700 underflows exp() in double).
  catastrophic <- !is.finite(mx) || max(log_evidence) < -700
  if (catastrophic) {
    cli::cli_warn(c(
      "Marginal evidence is (numerically) zero at every component; observation is far from the prior.",
      "i" = "Returning the prior mixture unchanged (`metadata$gmm_observe_no_update = TRUE`)."
    ))
    return(gmm(
      weights = g@weights,
      means = g@means,
      covariances = g@covariances,
      name = sprintf("%s [no update]", .op_name("observe", g)),
      metadata = .op_result_meta(g, "observe",
                                 list(gmm_observe_no_update = TRUE))
    ))
  }
  new_weights_unnormalised <- exp(log_new - mx)
  log_marginal_evidence <- mx + log(sum(new_weights_unnormalised))
  new_weights <- new_weights_unnormalised / sum(new_weights_unnormalised)

  gmm(
    weights     = new_weights,
    means       = new_means,
    covariances = new_covs,
    name        = .op_name("observe", g),
    metadata    = .op_result_meta(g, "observe", list(
      gmm_observe_no_update = FALSE,
      log_marginal_evidence = log_marginal_evidence
    ))
  )
}

## ---------------------------------------------------------------------------
## gmm_aggregate — affine pushforward through an aggregation matrix
## ---------------------------------------------------------------------------

#' Aggregation pushforward of a Gaussian mixture
#'
#' A named alias for [gmm_affine()] when `A` is a (row-wise) aggregation
#' matrix — e.g. a block-sum, block-average, or unequal-weight aggregation
#' used in downscaling pipelines. The mathematics is identical to
#' [gmm_affine()]; the alias gives the public API a clearer hook for
#' aggregation-specific diagnostics in later releases.
#'
#' @param g A [gmm] (or [gmm_fit]) in `R^p`.
#' @param A An `m` by `p` numeric matrix.
#' @param noise_cov Optional `m` by `m` SPD numeric matrix. Default
#'   `NULL` (deterministic aggregation).
#' @param ridge_eps Tiny ridge added to the output covariances for
#'   numerical hygiene.
#'
#' @returns A [gmm] in `R^m`.
#' @family operators
#' @export
#' @examples
#' g <- gmm(weights = c(0.5, 0.5),
#'          means = list(c(-1, 0, 1), c(1, 0, -1)),
#'          covariances = list(diag(3), diag(3)))
#' # Sum coordinates 1 and 2 into a single aggregate; pass coord 3 through.
#' A <- matrix(c(1, 1, 0,
#'               0, 0, 1), nrow = 2L, byrow = TRUE)
#' gmm_aggregate(g, A)
gmm_aggregate <- function(g, A, noise_cov = NULL, ridge_eps = 1e-6) {
  out <- gmm_affine(g, A, b = 0, noise_cov = noise_cov, ridge_eps = ridge_eps)
  out@name <- .op_name("aggregate", g)
  out@metadata$provenance[length(out@metadata$provenance)] <- "aggregate"
  out
}

## ---------------------------------------------------------------------------
## gmm_missing — Schur conditioning on a fully observed coordinate subset
## ---------------------------------------------------------------------------

#' Condition a Gaussian mixture on the exact values of some coordinates
#'
#' A structured wrapper around [gmm_conditionalise()] for the common case
#' where the observed coordinates are specified by integer index rather
#' than `NA`-padded vector. Equivalent to [gmm_observe()] with a
#' selection matrix `A` and zero noise covariance, but routed through
#' the Schur-complement path for efficiency.
#'
#' @param g A [gmm] (or [gmm_fit]) in `R^p`.
#' @param observed Integer vector of indices in `seq_len(p)`. The
#'   coordinates to condition on (fully observed).
#' @param values Numeric vector of length `length(observed)`. The
#'   observed values, in the same order as `observed`.
#'
#' @returns A [gmm] in `R^(p - length(observed))`.
#' @family operators
#' @export
#' @examples
#' g <- gmm(weights = c(0.4, 0.6),
#'          means = list(c(-1, 0), c(1, 0)),
#'          covariances = list(diag(2), diag(2)))
#' # Condition coord 2 on the value 0.5; keep coord 1.
#' gmm_missing(g, observed = 2L, values = 0.5)
gmm_missing <- function(g, observed, values) {
  if (!S7::S7_inherits(g, gmm)) {
    cli::cli_abort("`g` must be a {.cls gmm} object.")
  }
  p <- gmm_dim(g)
  observed <- as.integer(observed)
  if (anyNA(observed) || any(observed < 1L) || any(observed > p)) {
    cli::cli_abort("`observed` must be integer indices in 1..{p}.")
  }
  if (anyDuplicated(observed)) {
    cli::cli_abort("`observed` must be unique.")
  }
  if (length(observed) >= p) {
    cli::cli_abort("`observed` must leave at least one coordinate free.")
  }
  if (!is.numeric(values) || length(values) != length(observed)) {
    cli::cli_abort("`values` must be a numeric vector of length {length(observed)}.")
  }
  if (anyNA(values) || !all(is.finite(values))) {
    cli::cli_abort("`values` must be finite.")
  }

  given <- rep(NA_real_, p)
  given[observed] <- as.numeric(values)
  out <- gmm_conditionalise(g, given = given)
  out@name <- .op_name("missing", g)
  out@metadata$provenance[length(out@metadata$provenance)] <- "missing"
  out@metadata$observed_idx <- observed
  out@metadata$values <- as.numeric(values)
  out
}
