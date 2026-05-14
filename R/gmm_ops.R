## Closed-form operators on Gaussian mixtures.

#' Canonicalise the component ordering of a Gaussian mixture
#'
#' Returns a new [gmm] (or [gmm_fit]) with the components permuted into a
#' canonical order: weight descending, then `||mu||` descending as a
#' tiebreaker. The mixture distribution is unchanged — only the bookkeeping
#' order is — but the canonical ordering removes the EM label-switching
#' nuisance from snapshot tests, cross-run comparisons, and printed
#' summaries.
#'
#' Applied automatically by the regime-specific fitters
#' ([fit_moment_match()], [fit_em_samples()], [fit_kld_em()]) and by the
#' top-level dispatcher [fit_proxymix()] when `canonicalise = TRUE`
#' (the default).
#'
#' @param g A [gmm] (or [gmm_fit]) object.
#'
#' @returns A [gmm] (or [gmm_fit]) of the same subclass as `g`, with the
#'   components permuted into canonical order.
#' @family ops
#' @export
#' @examples
#' g <- gmm(weights = c(0.1, 0.6, 0.3),
#'          means = list(c(0, 0), c(3, 0), c(-1, 1)),
#'          covariances = list(diag(2), diag(2), diag(2)))
#' gmm_canonicalise(g)
gmm_canonicalise <- function(g) {
  if (!S7::S7_inherits(g, gmm)) {
    cli::cli_abort("`g` must be a {.cls gmm} object.")
  }
  K <- gmm_n_components(g)
  if (K <= 1L) return(g)
  norms <- vapply(g@means, function(mu) sqrt(sum(mu^2)), numeric(1L))
  ord <- order(g@weights, norms, decreasing = TRUE)
  if (identical(ord, seq_len(K))) return(g)
  out_weights <- g@weights[ord]
  out_means <- g@means[ord]
  out_covs <- g@covariances[ord]
  if (S7::S7_inherits(g, gmm_fit)) {
    gmm_fit(
      weights = out_weights,
      means = out_means,
      covariances = out_covs,
      target = g@target,
      regime = g@regime,
      diagnostics = g@diagnostics,
      converged = g@converged,
      iterations = g@iterations,
      call = g@call,
      name = g@name,
      metadata = g@metadata
    )
  } else {
    gmm(
      weights = out_weights,
      means = out_means,
      covariances = out_covs,
      name = g@name,
      metadata = g@metadata
    )
  }
}

#' Density of a Gaussian mixture
#'
#' Evaluates the density (or log-density) of a Gaussian mixture at one or
#' more points.
#'
#' @param x A numeric matrix with one observation per row, or a length-`p`
#'   numeric vector (treated as a single observation).
#' @param g A [gmm] (or [gmm_fit]) object.
#' @param log Logical. If `TRUE`, return log-densities.
#'
#' @returns A numeric vector of length `nrow(x)`.
#' @family ops
#' @export
#' @examples
#' g <- gmm(weights = c(0.5, 0.5),
#'          means = list(c(-1, 0), c(1, 0)),
#'          covariances = list(diag(2), diag(2)))
#' dgmm(c(0, 0), g)
#' dgmm(c(0, 0), g, log = TRUE)
dgmm <- function(x, g, log = FALSE) {
  if (!S7::S7_inherits(g, gmm)) {
    cli::cli_abort("`g` must be a {.cls gmm} object.")
  }
  p <- gmm_dim(g)
  x <- as_input_matrix(x, p)
  parts <- gmm_log_unnorm(x, g@weights, g@means, g@covariances)
  ldens <- logsumexp_rows(parts)
  if (isTRUE(log)) ldens else exp(ldens)
}

#' Sample from a Gaussian mixture
#'
#' Draws `n` independent samples from a Gaussian mixture.
#'
#' @param n Number of samples (positive integer scalar).
#' @param g A [gmm] (or [gmm_fit]) object.
#'
#' @returns A numeric matrix of dimension `n` by `p`.
#' @family ops
#' @export
#' @examples
#' g <- gmm(weights = c(0.5, 0.5),
#'          means = list(c(-1, 0), c(1, 0)),
#'          covariances = list(diag(2), diag(2)))
#' x <- rgmm(50L, g)
#' dim(x)
rgmm <- function(n, g) {
  if (!S7::S7_inherits(g, gmm)) {
    cli::cli_abort("`g` must be a {.cls gmm} object.")
  }
  if (length(n) != 1L || !is.numeric(n) || n < 0L) {
    cli::cli_abort("`n` must be a non-negative integer scalar.")
  }
  n <- as.integer(n)
  p <- gmm_dim(g)
  K <- gmm_n_components(g)
  if (n == 0L) return(matrix(0, nrow = 0L, ncol = p))
  z <- sample.int(K, size = n, replace = TRUE, prob = g@weights)
  out <- matrix(0, nrow = n, ncol = p)
  for (k in seq_len(K)) {
    idx <- which(z == k)
    if (length(idx) > 0L) {
      out[idx, ] <- mvnfast::rmvn(length(idx),
                                  mu = g@means[[k]],
                                  sigma = g@covariances[[k]])
    }
  }
  out
}

#' Marginal of a Gaussian mixture
#'
#' Computes the marginal distribution of a Gaussian mixture over a subset of
#' coordinates. The marginal of a Gaussian mixture is itself a Gaussian
#' mixture with the same weights.
#'
#' @param g A [gmm] (or [gmm_fit]) object.
#' @param keep Integer vector of coordinate indices to retain (in `1..p`).
#'
#' @returns A [gmm] object in dimension `length(keep)`.
#' @family ops
#' @export
#' @examples
#' g <- gmm(weights = c(0.5, 0.5),
#'          means = list(c(-1, 0, 2), c(1, 0, -2)),
#'          covariances = list(diag(3), diag(3)))
#' gmm_marginalise(g, keep = c(1L, 3L))
gmm_marginalise <- function(g, keep) {
  if (!S7::S7_inherits(g, gmm)) {
    cli::cli_abort("`g` must be a {.cls gmm} object.")
  }
  p <- gmm_dim(g)
  if (!is.numeric(keep) || any(keep < 1L) || any(keep > p) || anyDuplicated(keep)) {
    cli::cli_abort("`keep` must be an integer vector of unique indices in {.val 1}:{.val {p}}.")
  }
  keep <- as.integer(keep)
  K <- gmm_n_components(g)
  new_means <- vector("list", K)
  new_covs <- vector("list", K)
  for (k in seq_len(K)) {
    new_means[[k]] <- g@means[[k]][keep]
    new_covs[[k]] <- g@covariances[[k]][keep, keep, drop = FALSE]
  }
  gmm(
    weights = g@weights,
    means = new_means,
    covariances = new_covs,
    name = sprintf("marginal[%s] of %s",
                   paste(keep, collapse = ","),
                   g@name)
  )
}

#' Conditional of a Gaussian mixture
#'
#' Computes the conditional distribution of a Gaussian mixture given fixed
#' values of a subset of coordinates, by the Schur-complement formula
#' applied component-wise and re-weighted by the marginal evidence
#' \eqn{p(\textit{x}_b)} of each component.
#'
#' @param g A [gmm] (or [gmm_fit]) object.
#' @param given A length-`p` numeric vector. Coordinates to *condition on*
#'   take their numeric value; coordinates left *free* are `NA`.
#'
#' @returns A [gmm] object in dimension equal to the number of free
#'   coordinates.
#' @family ops
#' @export
#' @examples
#' g <- gmm(weights = c(0.5, 0.5),
#'          means = list(c(-1, 0), c(1, 0)),
#'          covariances = list(diag(2), diag(2)))
#' gmm_conditionalise(g, given = c(NA, 0.5))
gmm_conditionalise <- function(g, given) {
  if (!S7::S7_inherits(g, gmm)) {
    cli::cli_abort("`g` must be a {.cls gmm} object.")
  }
  p <- gmm_dim(g)
  if (!is.numeric(given) || length(given) != p) {
    cli::cli_abort("`given` must be a numeric vector of length {p}.")
  }
  free <- which(is.na(given))
  fixed <- which(!is.na(given))
  if (length(free) == 0L) {
    cli::cli_abort("`given` has no free (NA) coordinates; result would be zero-dimensional.")
  }
  if (length(fixed) == 0L) {
    return(g) ## nothing to condition on
  }
  xb <- given[fixed]
  K <- gmm_n_components(g)
  new_means <- vector("list", K)
  new_covs <- vector("list", K)
  log_evidence <- numeric(K)
  for (k in seq_len(K)) {
    mu <- g@means[[k]]
    S <- g@covariances[[k]]
    mu_a <- mu[free]
    mu_b <- mu[fixed]
    S_aa <- S[free, free, drop = FALSE]
    S_ab <- S[free, fixed, drop = FALSE]
    S_bb <- S[fixed, fixed, drop = FALSE]
    L_bb <- chol(S_bb)
    S_bb_inv <- chol2inv(L_bb)
    new_means[[k]] <- as.numeric(mu_a + S_ab %*% S_bb_inv %*% (xb - mu_b))
    new_covs[[k]] <- symmetrise(S_aa - S_ab %*% S_bb_inv %*% t(S_ab))
    log_evidence[k] <- mvnfast::dmvn(matrix(xb, nrow = 1L),
                                     mu = mu_b, sigma = S_bb, log = TRUE)
  }
  log_new_weights <- log(g@weights) + log_evidence
  log_new_weights <- log_new_weights - max(log_new_weights)
  new_weights <- exp(log_new_weights)
  new_weights <- new_weights / sum(new_weights)
  gmm(
    weights = new_weights,
    means = new_means,
    covariances = new_covs,
    name = sprintf("conditional of %s", g@name)
  )
}

#' Kullback-Leibler divergence between two Gaussian mixtures
#'
#' Estimates `KL(p || q)` between two Gaussian mixtures by Monte Carlo and,
#' optionally, evaluates the Hershey--Olsen variational approximation as a
#' deterministic sanity check.
#'
#' The Monte Carlo estimator draws `n_mc` samples from `p` and returns the
#' empirical mean of `log p(x) - log q(x)`, together with a Monte Carlo
#' standard error.
#'
#' The variational approximation is
#' \deqn{\widehat{D}_{\mathrm{var}}(p \Vert q) =
#'   \sum_a \pi_a \log\!\left(\frac{\sum_{a'} \pi_{a'} \, e^{-\mathrm{KL}(p_a \Vert p_{a'})}}{\sum_b \omega_b \, e^{-\mathrm{KL}(p_a \Vert q_b)}}\right),}
#' which is exact when `p == q` and tends to be a usable lower bound when
#' the components of `p` and `q` are well-separated. The closed-form
#' Gaussian--Gaussian KL `KL(p_a || q_b)` is used internally.
#'
#' @param p,q Two [gmm] (or [gmm_fit]) objects of the same ambient
#'   dimension.
#' @param n_mc Number of Monte Carlo samples drawn from `p`.
#' @param variational If `TRUE`, also return the Hershey--Olsen variational
#'   approximation.
#'
#' @returns A list with components
#'   * `mc` - the Monte Carlo estimate of `KL(p || q)`,
#'   * `mc_se` - its Monte Carlo standard error,
#'   * `variational` - the variational approximation (`NA` if
#'     `variational = FALSE`),
#'   * `n_mc` - the number of Monte Carlo samples used.
#' @family ops
#' @export
#' @examples
#' p <- gmm(weights = c(0.5, 0.5),
#'          means = list(c(-1, 0), c(1, 0)),
#'          covariances = list(diag(2), diag(2)))
#' q <- gmm(weights = 1,
#'          means = list(c(0, 0)),
#'          covariances = list(diag(2) * 2))
#' gmm_kld(p, q, n_mc = 500L)
gmm_kld <- function(p, q, n_mc = 5000L, variational = TRUE) {
  if (!S7::S7_inherits(p, gmm) || !S7::S7_inherits(q, gmm)) {
    cli::cli_abort("`p` and `q` must both be {.cls gmm} objects.")
  }
  if (gmm_dim(p) != gmm_dim(q)) {
    cli::cli_abort("`p` and `q` must have the same ambient dimension.")
  }
  n_mc <- as.integer(n_mc)
  x <- rgmm(n_mc, p)
  log_p <- dgmm(x, p, log = TRUE)
  log_q <- dgmm(x, q, log = TRUE)
  d <- log_p - log_q
  finite <- is.finite(d)
  mc <- mean(d[finite])
  mc_se <- stats::sd(d[finite]) / sqrt(sum(finite))

  var_est <- NA_real_
  if (isTRUE(variational)) {
    Kp <- gmm_n_components(p)
    Kq <- gmm_n_components(q)
    pi_p <- p@weights
    pi_q <- q@weights
    kl_pp <- matrix(0, Kp, Kp)
    for (a in seq_len(Kp)) {
      for (b in seq_len(Kp)) {
        kl_pp[a, b] <- kl_gauss(p@means[[a]], p@covariances[[a]],
                                p@means[[b]], p@covariances[[b]])
      }
    }
    kl_pq <- matrix(0, Kp, Kq)
    for (a in seq_len(Kp)) {
      for (b in seq_len(Kq)) {
        kl_pq[a, b] <- kl_gauss(p@means[[a]], p@covariances[[a]],
                                q@means[[b]], q@covariances[[b]])
      }
    }
    num <- log(rowSums(matrix(rep(pi_p, each = Kp), nrow = Kp) * exp(-kl_pp)))
    den <- log(rowSums(matrix(rep(pi_q, each = Kp), nrow = Kp) * exp(-kl_pq)))
    var_est <- sum(pi_p * (num - den))
  }

  list(mc = mc, mc_se = mc_se, variational = var_est, n_mc = n_mc)
}
