## uplift_internal.R -- Closed-form machinery shared by the decision layer.
##
## The decision module reads three rungs of Pearl's ladder off one fitted
## Gaussian mixture over (outcome, treatment, covariates): seeing via
## component conditioning, doing via an X-only regime gate, imagining via
## per-component abduction. The helpers here are the closed-form pieces those
## reads share -- a component-wise Schur conditional, the abduction
## responsibilities of an observed unit, the mixture mean, and the per-unit
## overlap (positivity) score. None are exported.
##
## All formulae are the finite-mixture analogues used in `gmm_ops.R` /
## `operator_calculus.R`; the within-component conditioning matches
## `gmm_conditionalise()` exactly so the K = 1 reduction is the OLS prediction.

## Component Schur conditional ----------------------------------------------

## Within a single Gaussian component (`mu`, `S`), the conditional law of the
## `free` coordinates given the `cond_idx` coordinates fixed at `cond_val`.
## Returns the conditional mean and covariance of the free block. This is the
## same Schur-complement update applied component-wise by
## `gmm_conditionalise()`; isolating it lets the do-operator gate and condition
## on different coordinate sets.
.schur_cond_component <- function(mu, S, free, cond_idx, cond_val,
                                  ridge_eps = 1e-8) {
  mu_a <- mu[free]
  mu_b <- mu[cond_idx]
  S_aa <- S[free, free, drop = FALSE]
  S_ab <- S[free, cond_idx, drop = FALSE]
  S_bb <- S[cond_idx, cond_idx, drop = FALSE]
  L_bb <- tryCatch(chol(S_bb), error = function(e) NULL)
  if (is.null(L_bb)) {
    S_bb <- ridge(S_bb, max(ridge_eps, 1e-10))
    L_bb <- chol(S_bb)
  }
  S_bb_inv <- chol2inv(L_bb)
  gain <- S_ab %*% S_bb_inv
  mean_free <- as.numeric(mu_a + gain %*% (cond_val - mu_b))
  cov_free <- symmetrise(S_aa - gain %*% t(S_ab))
  list(mean = mean_free, cov = cov_free)
}

## Abduction responsibilities -----------------------------------------------

## Posterior over the latent component given an observed unit `obs_val` on the
## coordinates `obs_idx` -- the abduction step. Returns a length-K vector of
## responsibilities `pi_k(evidence)` that sum to one, computed in the log
## domain for stability.
.abduction_responsibilities <- function(g, obs_idx, obs_val) {
  K <- gmm_n_components(g)
  log_resp <- numeric(K)
  has_obs <- length(obs_idx) > 0L
  for (k in seq_len(K)) {
    log_resp[k] <- log(g@weights[k])
    if (has_obs) {
      mu_b <- g@means[[k]][obs_idx]
      S_bb <- g@covariances[[k]][obs_idx, obs_idx, drop = FALSE]
      log_resp[k] <- log_resp[k] +
        mvnfast::dmvn(matrix(obs_val, nrow = 1L),
                      mu = mu_b, sigma = S_bb, log = TRUE)
    }
  }
  mx <- max(log_resp)
  if (!is.finite(mx)) {
    cli::cli_abort(c(
      "Evidence has (numerically) zero density under every component.",
      "i" = "The observed unit is far outside the fitted mixture's support."
    ))
  }
  w <- exp(log_resp - mx)
  w / sum(w)
}

## Batch-safe responsibilities ----------------------------------------------

## Like `.abduction_responsibilities()` but never aborts: when the evidence
## underflows every component (the unit is outside the fitted support), it
## falls back to the prior weights and marks `attr(., "covered") <- FALSE`.
## Used by the batch scoring verbs, where an out-of-support unit must be
## flagged, not error.
.responsibilities_safe <- function(g, obs_idx, obs_val) {
  K <- gmm_n_components(g)
  log_resp <- log(g@weights)
  if (length(obs_idx) > 0L) {
    for (k in seq_len(K)) {
      mu_b <- g@means[[k]][obs_idx]
      S_bb <- g@covariances[[k]][obs_idx, obs_idx, drop = FALSE]
      log_resp[k] <- log_resp[k] +
        mvnfast::dmvn(matrix(obs_val, nrow = 1L),
                      mu = mu_b, sigma = S_bb, log = TRUE)
    }
  }
  mx <- max(log_resp)
  if (!is.finite(mx)) {
    out <- g@weights
    attr(out, "covered") <- FALSE
    return(out)
  }
  w <- exp(log_resp - mx)
  out <- w / sum(w)
  attr(out, "covered") <- TRUE
  out
}

## Vectorised responsibilities over a matrix of units -----------------------

## Batch form of `.responsibilities_safe()`: `obs_mat` is `n` by
## `length(obs_idx)`. Returns an `n` by `K` matrix of responsibilities, one row
## per unit, with `attr(., "covered")` a length-`n` logical. Rows that underflow
## every component fall back to the prior weights and are marked not covered.
## `mvnfast::dmvn` is vectorised over rows, so this is the O(K) serving path.
.responsibilities_batch <- function(g, obs_idx, obs_mat) {
  n <- nrow(obs_mat)
  K <- gmm_n_components(g)
  log_r <- matrix(0, nrow = n, ncol = K)
  log_w <- log(g@weights)
  for (k in seq_len(K)) {
    mu_b <- g@means[[k]][obs_idx]
    S_bb <- g@covariances[[k]][obs_idx, obs_idx, drop = FALSE]
    log_r[, k] <- log_w[k] +
      mvnfast::dmvn(obs_mat, mu = mu_b, sigma = S_bb, log = TRUE)
  }
  mx <- apply(log_r, 1L, max)
  covered <- is.finite(mx)
  r <- matrix(g@weights, nrow = n, ncol = K, byrow = TRUE)
  if (any(covered)) {
    w <- exp(log_r[covered, , drop = FALSE] - mx[covered])
    r[covered, ] <- w / rowSums(w)
  }
  attr(r, "covered") <- covered
  r
}

## Vectorised positivity / mass coverage over a matrix of units. Returns the
## length-`n` minimum-over-arms coverage probability.
.coverage_batch <- function(g, z_idx, arms, X) {
  n <- nrow(X)
  df <- length(z_idx)
  K <- gmm_n_components(g)
  cover <- matrix(Inf, nrow = n, ncol = length(arms))
  for (a in seq_along(arms)) {
    Z <- cbind(arms[a], X)
    best <- rep(Inf, n)
    for (k in seq_len(K)) {
      Pz <- chol2inv(chol(g@covariances[[k]][z_idx, z_idx, drop = FALSE]))
      D <- sweep(Z, 2L, g@means[[k]][z_idx])
      d2 <- rowSums((D %*% Pz) * D)
      best <- pmin(best, d2)
    }
    cover[, a] <- stats::pchisq(best, df = df, lower.tail = FALSE)
  }
  apply(cover, 1L, min)
}

## Per-unit positivity / mass coverage --------------------------------------

## Coverage of a (treatment, covariate) configuration by the fitted joint. For
## each arm value in `arms`, the squared Mahalanobis distance to the nearest
## component centre over the (treatment, covariate) block is converted to an
## upper-tail chi-square probability (1 at a centre, -> 0 in the tail). The
## returned coverage is the minimum across arms: CATE needs both arms covered.
.unit_coverage <- function(g, z_idx, arms, x) {
  df <- length(z_idx)
  K <- gmm_n_components(g)
  cover_arm <- vapply(arms, function(t_val) {
    z <- c(t_val, x)
    d2 <- vapply(seq_len(K), function(k) {
      mu <- g@means[[k]][z_idx]
      S <- g@covariances[[k]][z_idx, z_idx, drop = FALSE]
      dv <- z - mu
      as.numeric(crossprod(dv, chol2inv(chol(S)) %*% dv))
    }, numeric(1L))
    stats::pchisq(min(d2), df = df, lower.tail = FALSE)
  }, numeric(1L))
  min(cover_arm)
}

