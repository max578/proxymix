## Completing the closed-form mixture algebra: moments, pointwise product,
## convolution, mixing, and the one-dimensional distribution/quantile
## functions. Every operation here is exact (no Monte Carlo) and follows
## the operator metadata policy (`.op_result_meta()`): the source fit's
## quality certificate travels, and each operation appends to the
## provenance chain.

#' Mean and covariance of a Gaussian mixture
#'
#' The exact first two moments of a mixture:
#' \deqn{\mu = \sum_k w_k \mu_k, \qquad
#'   \Sigma = \sum_k w_k \left(\Sigma_k + \mu_k \mu_k^\top\right)
#'   - \mu \mu^\top.}
#'
#' @param g A [gmm] (or [gmm_fit]).
#'
#' @returns `gmm_mean()` returns a length-`p` numeric vector; `gmm_cov()`
#'   returns a `p` by `p` numeric matrix.
#' @family ops
#' @export
#' @examples
#' g <- gmm(weights = c(0.3, 0.7), means = list(c(-1, 0), c(1, 2)),
#'          covariances = list(diag(2), 0.5 * diag(2)))
#' gmm_mean(g)
#' gmm_cov(g)
gmm_mean <- function(g) {
  if (!S7::S7_inherits(g, gmm)) {
    cli::cli_abort("`g` must be a {.cls gmm} object.")
  }
  .gmm_moments(g)$mean
}

#' @rdname gmm_mean
#' @export
gmm_cov <- function(g) {
  if (!S7::S7_inherits(g, gmm)) {
    cli::cli_abort("`g` must be a {.cls gmm} object.")
  }
  .gmm_moments(g)$cov
}

#' Pointwise product of two Gaussian mixtures
#'
#' The normalised pointwise product of two mixture densities -- the
#' conjugate Bayes update when one mixture plays the prior and the other
#' the likelihood. The product of two Gaussian mixtures is again a
#' Gaussian mixture with \eqn{K_1 K_2} components:
#' \deqn{p(x)\, q(x) \propto \sum_{ij} w_i v_j z_{ij}\,
#'   \mathcal{N}\!\left(x;\ \mu_{ij},\ \Sigma_{ij}\right),}
#' with \eqn{\Sigma_{ij} = (\Sigma_i^{-1} + S_j^{-1})^{-1}},
#' \eqn{\mu_{ij} = \Sigma_{ij} (\Sigma_i^{-1} \mu_i + S_j^{-1} m_j)}, and
#' \eqn{z_{ij} = \mathcal{N}(\mu_i;\ m_j,\ \Sigma_i + S_j)} evaluated in
#' the log domain. The overall normaliser
#' \eqn{\int p\, q = \sum_{ij} w_i v_j z_{ij}} is returned in the result's
#' metadata as `log_integral` (it is the marginal evidence of the update).
#'
#' Because the component count multiplies, chained products grow
#' geometrically; pass `reduce` to cap the count via [gmm_reduce()] after
#' the product (the standard assumed-density trick), or manage the budget
#' yourself.
#'
#' @param g1,g2 Two [gmm] (or [gmm_fit]) objects of the same ambient
#'   dimension.
#' @param reduce Optional positive integer: reduce the product to at most
#'   this many components via [gmm_reduce()]. Default `NULL` (no
#'   reduction).
#'
#' @returns A [gmm] with `K1 * K2` components (or at most `reduce`), with
#'   `metadata$log_integral` recording \eqn{\log \int p\, q}.
#' @family operators
#' @export
#' @examples
#' prior <- gmm(weights = c(0.5, 0.5), means = list(-2, 2),
#'              covariances = list(matrix(1), matrix(1)))
#' lik <- gmm(weights = 1, means = list(0.5), covariances = list(matrix(0.5)))
#' post <- gmm_product(prior, lik)
#' post@weights
gmm_product <- function(g1, g2, reduce = NULL) {
  if (!S7::S7_inherits(g1, gmm) || !S7::S7_inherits(g2, gmm)) {
    cli::cli_abort("`g1` and `g2` must both be {.cls gmm} objects.")
  }
  p <- gmm_dim(g1)
  if (gmm_dim(g2) != p) {
    cli::cli_abort("`g1` and `g2` must have the same ambient dimension.")
  }
  K1 <- gmm_n_components(g1)
  K2 <- gmm_n_components(g2)
  if (is.null(reduce) && K1 * K2 > 100L) {
    cli::cli_warn(c(
      "The product has {K1 * K2} components ({K1} x {K2}).",
      "i" = "Pass {.code reduce =} to cap the count via {.fn gmm_reduce}."
    ), class = "proxymix_k_growth")
  }

  n_out <- K1 * K2
  out_means <- vector("list", n_out)
  out_covs <- vector("list", n_out)
  log_w <- numeric(n_out)
  idx <- 0L
  for (i in seq_len(K1)) {
    mi <- g1@means[[i]]
    Si <- g1@covariances[[i]]
    for (j in seq_len(K2)) {
      mj <- g2@means[[j]]
      Sj <- g2@covariances[[j]]
      Ssum <- Si + Sj
      ## (Si^{-1} + Sj^{-1})^{-1} = Si (Si + Sj)^{-1} Sj, and the paired
      ## mean uses the same single solve of (Si + Sj).
      sol <- solve(Ssum)
      Sij <- symmetrise(Si %*% sol %*% Sj)
      mu_ij <- as.numeric(Sj %*% sol %*% mi + Si %*% sol %*% mj)
      idx <- idx + 1L
      out_means[[idx]] <- mu_ij
      out_covs[[idx]] <- Sij
      log_w[idx] <- log(g1@weights[i]) + log(g2@weights[j]) +
        mvnfast::dmvn(matrix(mi, nrow = 1L), mu = mj, sigma = Ssum,
                      log = TRUE)
    }
  }
  log_integral <- as.numeric(logsumexp_rows(matrix(log_w, nrow = 1L)))
  w <- exp(log_w - log_integral)

  out <- gmm(
    weights = w / sum(w),
    means = out_means,
    covariances = out_covs,
    name = sprintf("product(%s, %s)",
                   if (nchar(g1@name) > 20L) "..." else g1@name,
                   if (nchar(g2@name) > 20L) "..." else g2@name),
    metadata = .op_result_meta(g1, "product",
                               list(log_integral = log_integral))
  )
  if (!is.null(reduce)) {
    li <- out@metadata$log_integral
    out <- gmm_reduce(out, k_max = as.integer(reduce))
    out@metadata$log_integral <- li
  }
  out
}

#' Convolution of two independent Gaussian mixtures
#'
#' The exact distribution of \eqn{X + Y} for independent \eqn{X \sim g_1}
#' and \eqn{Y \sim g_2}: a Gaussian mixture with \eqn{K_1 K_2} components,
#' \deqn{g_1 * g_2 = \sum_{ij} w_i v_j\,
#'   \mathcal{N}\!\left(\mu_i + m_j,\ \Sigma_i + S_j\right).}
#'
#' For the affine special case \eqn{X + c} with a constant `c`, or
#' \eqn{A X + \epsilon} with Gaussian \eqn{\epsilon}, use [gmm_affine()];
#' the convolution operator is the general mixture-plus-mixture case that
#' `gmm_affine()` cannot express.
#'
#' @param g1,g2 Two [gmm] (or [gmm_fit]) objects of the same ambient
#'   dimension.
#'
#' @returns A [gmm] with `K1 * K2` components.
#' @family operators
#' @export
#' @examples
#' g1 <- gmm(weights = c(0.5, 0.5), means = list(-1, 1),
#'           covariances = list(matrix(0.5), matrix(0.5)))
#' g2 <- gmm(weights = 1, means = list(2), covariances = list(matrix(1)))
#' gmm_convolve(g1, g2)
gmm_convolve <- function(g1, g2) {
  if (!S7::S7_inherits(g1, gmm) || !S7::S7_inherits(g2, gmm)) {
    cli::cli_abort("`g1` and `g2` must both be {.cls gmm} objects.")
  }
  p <- gmm_dim(g1)
  if (gmm_dim(g2) != p) {
    cli::cli_abort("`g1` and `g2` must have the same ambient dimension.")
  }
  K1 <- gmm_n_components(g1)
  K2 <- gmm_n_components(g2)
  n_out <- K1 * K2
  out_means <- vector("list", n_out)
  out_covs <- vector("list", n_out)
  w <- numeric(n_out)
  idx <- 0L
  for (i in seq_len(K1)) {
    for (j in seq_len(K2)) {
      idx <- idx + 1L
      out_means[[idx]] <- as.numeric(g1@means[[i]] + g2@means[[j]])
      out_covs[[idx]] <- symmetrise(g1@covariances[[i]] + g2@covariances[[j]])
      w[idx] <- g1@weights[i] * g2@weights[j]
    }
  }
  gmm(
    weights = w / sum(w),
    means = out_means,
    covariances = out_covs,
    name = sprintf("convolve(%s, %s)",
                   if (nchar(g1@name) > 20L) "..." else g1@name,
                   if (nchar(g2@name) > 20L) "..." else g2@name),
    metadata = .op_result_meta(g1, "convolve")
  )
}

#' Mix Gaussian mixtures into one mixture
#'
#' Flattens a list of Gaussian mixtures into a single mixture whose
#' density is the weighted average \eqn{\sum_b \alpha_b\, g_b(x)} -- model
#' averaging, prior pooling, or the mixture-of-mixtures construction.
#'
#' @param gmms A non-empty list of [gmm] (or [gmm_fit]) objects sharing
#'   one ambient dimension.
#' @param weights Optional non-negative mixing weights, one per element of
#'   `gmms` (normalised internally). Default: equal weights.
#'
#' @returns A [gmm] with `sum(K_b)` components.
#' @family operators
#' @export
#' @examples
#' g1 <- gmm(weights = 1, means = list(-1), covariances = list(matrix(1)))
#' g2 <- gmm(weights = 1, means = list(2), covariances = list(matrix(0.5)))
#' gmm_mix(list(g1, g2), weights = c(0.7, 0.3))
gmm_mix <- function(gmms, weights = NULL) {
  if (!is.list(gmms) || length(gmms) < 1L ||
        !all(vapply(gmms, function(g) S7::S7_inherits(g, gmm), logical(1L)))) {
    cli::cli_abort("`gmms` must be a non-empty list of {.cls gmm} objects.")
  }
  B <- length(gmms)
  p <- gmm_dim(gmms[[1L]])
  dims <- vapply(gmms, gmm_dim, integer(1L))
  if (!all(dims == p)) {
    cli::cli_abort("all elements of `gmms` must share one ambient dimension (got {unique(dims)}).")
  }
  if (is.null(weights)) weights <- rep(1 / B, B)
  if (!is.numeric(weights) || length(weights) != B || anyNA(weights) ||
        any(weights < 0) || sum(weights) <= 0) {
    cli::cli_abort("`weights` must be {B} non-negative numbers with a positive sum.")
  }
  alpha <- weights / sum(weights)
  .stack_gmms(gmms, log(alpha), name = sprintf("mix[%d]", B))
}

#' Distribution and quantile functions of a one-dimensional mixture
#'
#' The exact cumulative distribution function
#' \eqn{F(x) = \sum_k w_k \Phi\!\left((x - \mu_k) / \sigma_k\right)} of a
#' one-dimensional Gaussian mixture, and its inverse by monotone
#' root-finding. Together with [dgmm()] and [rgmm()] these complete the
#' usual d/p/q/r quartet for the one-dimensional case; for tail
#' probabilities of a multivariate mixture, marginalise first
#' ([gmm_marginalise()]) or push through the relevant linear functional
#' ([gmm_affine()]).
#'
#' @param q Numeric vector of quantiles.
#' @param p Numeric vector of probabilities in `(0, 1)`.
#' @param g A one-dimensional [gmm] (or [gmm_fit]).
#' @param lower.tail Logical; if `TRUE` (default), probabilities are
#'   \eqn{P(X \le x)}.
#'
#' @returns A numeric vector the length of the first argument.
#' @family ops
#' @export
#' @examples
#' g <- gmm(weights = c(0.4, 0.6), means = list(-2, 1),
#'          covariances = list(matrix(0.5), matrix(1)))
#' pgmm(c(-2, 0, 2), g)
#' qgmm(c(0.1, 0.5, 0.9), g)
pgmm <- function(q, g, lower.tail = TRUE) {
  .check_gmm_1d(g)
  sds <- sqrt(vapply(g@covariances, function(S) S[1L, 1L], numeric(1L)))
  mus <- vapply(g@means, function(m) m[1L], numeric(1L))
  out <- vapply(q, function(qi) {
    sum(g@weights * stats::pnorm(qi, mean = mus, sd = sds))
  }, numeric(1L))
  if (isTRUE(lower.tail)) out else 1 - out
}

#' @rdname pgmm
#' @export
qgmm <- function(p, g) {
  .check_gmm_1d(g)
  if (!is.numeric(p) || anyNA(p) || any(p <= 0) || any(p >= 1)) {
    cli::cli_abort("`p` must be probabilities strictly inside (0, 1).")
  }
  sds <- sqrt(vapply(g@covariances, function(S) S[1L, 1L], numeric(1L)))
  mus <- vapply(g@means, function(m) m[1L], numeric(1L))
  vapply(p, function(pi) {
    ## Bracket from the extreme per-component normal quantiles: F is a
    ## weighted average of the component CDFs, so its inverse lies inside.
    lo <- min(stats::qnorm(pi, mean = mus, sd = sds))
    hi <- max(stats::qnorm(pi, mean = mus, sd = sds))
    if (hi - lo < .Machine$double.eps) return(lo)
    stats::uniroot(function(x) pgmm(x, g) - pi, lower = lo, upper = hi,
                   tol = 1e-10, extendInt = "upX")$root
  }, numeric(1L))
}

## Guard shared by the one-dimensional distribution functions.
.check_gmm_1d <- function(g) {
  if (!S7::S7_inherits(g, gmm)) {
    cli::cli_abort("`g` must be a {.cls gmm} object.")
  }
  if (gmm_dim(g) != 1L) {
    cli::cli_abort(c(
      "`g` must be one-dimensional (got p = {gmm_dim(g)}).",
      "i" = "Marginalise first ({.fn gmm_marginalise}) or push through a linear functional ({.fn gmm_affine})."
    ))
  }
  invisible(NULL)
}
