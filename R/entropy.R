## Closed-form entropy and divergence diagnostics for Gaussian mixtures.
##
## Spine: the integral of a product of two Gaussian densities is itself a
## Gaussian density evaluation,
##   integral N(x; a, A) N(x; b, B) dx = N(a; b, A + B),
## so any functional that is quadratic in a Gaussian mixture reduces to a finite
## sum of Gaussian-density evaluations. The order-2 (quadratic) Renyi entropy
## and the Cauchy-Schwarz divergence are therefore closed-form; Shannon entropy
## (the log of a sum) is estimated by Monte Carlo and bracketed by an analytic
## upper bound.

## Scalar log-sum-exp, robust to -Inf entries (zero-weight components).
.lse <- function(v) {
  v <- v[is.finite(v)]
  if (length(v) == 0L) {
    return(-Inf)
  }
  m <- max(v)
  m + log(sum(exp(v - m)))
}

## log of the cross-information potential
##   V(p, q) = sum_{i,j} w^p_i w^q_j N(mu^p_i; mu^q_j, Sigma^p_i + Sigma^q_j),
## evaluated on the log scale for numerical stability.
.log_cross_information_potential <- function(p, q) {
  wp <- p@weights
  wq <- q@weights
  terms <- vector("numeric", length(wp) * length(wq))
  idx <- 1L
  for (i in seq_along(wp)) {
    mu_i <- matrix(p@means[[i]], nrow = 1L)
    cov_i <- p@covariances[[i]]
    for (j in seq_along(wq)) {
      ld <- mvnfast::dmvn(mu_i,
                          mu = q@means[[j]],
                          sigma = cov_i + q@covariances[[j]],
                          log = TRUE)
      terms[idx] <- log(wp[i]) + log(wq[j]) + ld
      idx <- idx + 1L
    }
  }
  .lse(terms)
}

## Analytic upper bound on Shannon differential entropy of a mixture:
##   H(g) <= sum_k w_k H(N_k) + H(w),  H(N_k) = 1/2 log((2 pi e)^p |Sigma_k|),
##   H(w) = -sum_k w_k log w_k.
.mixture_entropy_upper_bound <- function(g) {
  p <- gmm_dim(g)
  w <- g@weights
  comp <- vapply(seq_along(w), function(k) {
    ld <- as.numeric(determinant(g@covariances[[k]], logarithm = TRUE)$modulus)
    0.5 * (p * (1 + log(2 * pi)) + ld)
  }, numeric(1))
  pos <- w > 0
  hw <- -sum(w[pos] * log(w[pos]))
  sum(w * comp) + hw
}

#' Differential entropy of a Gaussian mixture
#'
#' Computes the differential entropy of a Gaussian mixture. The quadratic
#' (order-2) Renyi entropy \eqn{H_2(g) = -\log \int g(x)^2 \, dx} is available in
#' closed form, because \eqn{\int g^2} is a finite sum of Gaussian-density
#' evaluations. Shannon entropy has no closed form for a mixture (the integrand
#' carries the logarithm of a sum) and is estimated by Monte Carlo, reported
#' with its standard error and an analytic upper bound that brackets it from
#' above.
#'
#' @param g A [gmm] (or [gmm_fit]) object.
#' @param order `"renyi2"` (closed-form quadratic Renyi entropy, the default) or
#'   `"shannon"` (Monte-Carlo estimate with an analytic upper bound).
#' @param n_mc Number of Monte Carlo samples for `order = "shannon"`.
#' @param seed Optional integer seed for the Monte Carlo draw.
#'
#' @returns For `order = "renyi2"`, a numeric scalar. For `order = "shannon"`, a
#'   list with components `mc` (the estimate), `mc_se` (its standard error),
#'   `upper_bound` (the analytic upper bound), and `n_mc`.
#' @family diagnostics
#' @seealso [gmm_divergence()], [gmm_kld()]
#' @export
#' @examples
#' g <- gmm(weights = c(0.5, 0.5),
#'          means = list(c(-2, 0), c(2, 0)),
#'          covariances = list(diag(2), diag(2)))
#' gmm_entropy(g)
#' gmm_entropy(g, order = "shannon", n_mc = 2000L, seed = 1L)
gmm_entropy <- function(g, order = c("renyi2", "shannon"),
                        n_mc = 5000L, seed = NULL) {
  if (!S7::S7_inherits(g, gmm)) {
    cli::cli_abort("`g` must be a {.cls gmm} object.")
  }
  order <- match.arg(order)
  if (order == "renyi2") {
    return(-.log_cross_information_potential(g, g))
  }
  n_mc <- as.integer(n_mc)
  draw <- function() rgmm(n_mc, g)
  x <- if (is.null(seed)) draw() else withr::with_seed(seed, draw())
  lg <- dgmm(x, g, log = TRUE)
  finite <- is.finite(lg)
  list(
    mc = -mean(lg[finite]),
    mc_se = stats::sd(lg[finite]) / sqrt(sum(finite)),
    upper_bound = .mixture_entropy_upper_bound(g),
    n_mc = n_mc
  )
}

#' Divergence between two Gaussian mixtures
#'
#' Computes a divergence between two Gaussian mixtures of the same ambient
#' dimension. The Cauchy-Schwarz divergence
#' \deqn{D_{\mathrm{CS}}(p, q) = \tfrac{1}{2}\log V(p, p)
#'   + \tfrac{1}{2}\log V(q, q) - \log V(p, q),}
#' with \eqn{V(p, q) = \int p(x) q(x)\, dx}, is closed-form, symmetric,
#' non-negative, and zero exactly when \eqn{p \propto q}. The `"kl"` option
#' delegates to [gmm_kld()], a Monte-Carlo estimate of the asymmetric
#' Kullback-Leibler divergence \eqn{\mathrm{KL}(p \Vert q)}.
#'
#' @param p,q Two [gmm] (or [gmm_fit]) objects of the same ambient dimension.
#' @param type `"cs"` (closed-form symmetric Cauchy-Schwarz divergence, the
#'   default) or `"kl"` (delegates to [gmm_kld()]).
#' @param n_mc Number of Monte Carlo samples used when `type = "kl"`.
#'
#' @returns For `type = "cs"`, a non-negative numeric scalar. For `type = "kl"`,
#'   the list returned by [gmm_kld()].
#' @family ops
#' @seealso [gmm_entropy()], [gmm_kld()]
#' @export
#' @examples
#' p <- gmm(weights = c(0.5, 0.5),
#'          means = list(c(-1, 0), c(1, 0)),
#'          covariances = list(diag(2), diag(2)))
#' q <- gmm(weights = 1, means = list(c(0, 0)),
#'          covariances = list(diag(2) * 2))
#' gmm_divergence(p, q)
#' gmm_divergence(p, p)
gmm_divergence <- function(p, q, type = c("cs", "kl"), n_mc = 5000L) {
  if (!S7::S7_inherits(p, gmm) || !S7::S7_inherits(q, gmm)) {
    cli::cli_abort("`p` and `q` must both be {.cls gmm} objects.")
  }
  if (gmm_dim(p) != gmm_dim(q)) {
    cli::cli_abort("`p` and `q` must have the same ambient dimension.")
  }
  type <- match.arg(type)
  if (type == "kl") {
    return(gmm_kld(p, q, n_mc = n_mc))
  }
  d <- 0.5 * .log_cross_information_potential(p, p) +
    0.5 * .log_cross_information_potential(q, q) -
    .log_cross_information_potential(p, q)
  ## D_CS >= 0 by the Cauchy-Schwarz inequality; clean floating-point noise.
  if (d < 0 && d > -1e-9) {
    d <- 0
  }
  d
}
