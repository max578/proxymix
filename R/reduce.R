## Mixture reduction: collapse a Gaussian mixture to fewer components.
##
## A repeated, noisy operation (a Gaussian-sum filter, a product of mixtures, a
## fit with redundant components) leaves a mixture with more components than the
## structure warrants. `gmm_reduce()` collapses it to a budget by a greedy,
## moment-preserving pairwise merge: at each step the pair of components whose
## merge costs the least is replaced by the single Gaussian that preserves their
## combined weight, mean and covariance, until the budget is met.
##
## Two merge costs are offered. The Kullback-Leibler bound of Runnalls (2007) is
## the textbook choice; the Cauchy-Schwarz cost reuses the package's closed-form
## Gaussian-product identity (the same one behind `gmm_divergence`). Because each
## merge is moment-preserving, the whole reduction preserves the mixture's global
## mean and covariance exactly, and reducing to a single component returns the
## moment-matched Gaussian.

# ---------------------------------------------------------------------------
# Moment-preserving merge and the Gaussian-product helper
# ---------------------------------------------------------------------------

## Log-determinant of a symmetric positive-definite matrix via its Cholesky
## factor, with a determinant() fallback for the rare non-PD input.
.logdet <- function(S) {
  L <- tryCatch(chol(S), error = function(e) NULL)
  if (!is.null(L)) {
    2 * sum(log(diag(L)))
  } else {
    as.numeric(determinant(S, logarithm = TRUE)$modulus)
  }
}

## Moment-preserving merge of two weighted Gaussian components into one,
## preserving combined weight, mean and (second-moment) covariance:
##   w   = w_i + w_j
##   mu  = (w_i mu_i + w_j mu_j) / w
##   Sig = (w_i Sig_i + w_j Sig_j) / w
##         + (w_i w_j / w^2) (mu_i - mu_j)(mu_i - mu_j)'.
.merge_two <- function(wi, mi, Si, wj, mj, Sj) {
  w <- wi + wj
  mu <- (wi * mi + wj * mj) / w
  diff <- mi - mj
  Sig <- (wi * Si + wj * Sj) / w +
    (wi * wj / w^2) * tcrossprod(diff)
  list(w = w, mu = as.numeric(mu), Sigma = symmetrise(Sig))
}

## Gaussian-product integral N(a; b, A + B) = integral N(x;a,A) N(x;b,B) dx,
## evaluated as a single multivariate-normal density.
.gauss_product <- function(a, b, AB) {
  mvnfast::dmvn(matrix(a, nrow = 1L), mu = b, sigma = AB)
}

# ---------------------------------------------------------------------------
# Pairwise merge costs
# ---------------------------------------------------------------------------

## Runnalls (2007) Kullback-Leibler upper-bound cost of merging components i and
## j, given the merged covariance Sij:
##   B(i,j) = 1/2 [ (w_i + w_j) logdet(Sij) - w_i logdet(Si) - w_j logdet(Sj) ].
## Non-negative; zero for identical components.
.runnalls_cost <- function(wi, Si, wj, Sj, Sij) {
  0.5 * ((wi + wj) * .logdet(Sij) - wi * .logdet(Si) - wj * .logdet(Sj))
}

## Cauchy-Schwarz merge cost: the closed-form Cauchy-Schwarz divergence between
## the (normalised) two-component sub-mixture and its single-Gaussian merge,
## scaled by the merged mass so heavier merges cost more. Reuses the Gaussian-
## product identity, so it is closed form.
.cs_cost <- function(wi, mi, Si, wj, mj, Sj, mij, Sij) {
  w <- wi + wj
  p1 <- wi / w
  p2 <- wj / w
  v_pp <- p1^2 * .gauss_product(mi, mi, 2 * Si) +
    2 * p1 * p2 * .gauss_product(mi, mj, Si + Sj) +
    p2^2 * .gauss_product(mj, mj, 2 * Sj)
  v_qq <- .gauss_product(mij, mij, 2 * Sij)
  v_pq <- p1 * .gauss_product(mi, mij, Si + Sij) +
    p2 * .gauss_product(mj, mij, Sj + Sij)
  d_cs <- 0.5 * log(v_pp) + 0.5 * log(v_qq) - log(v_pq)
  w * max(d_cs, 0)
}

# ---------------------------------------------------------------------------
# gmm_reduce
# ---------------------------------------------------------------------------

#' Reduce a Gaussian mixture to fewer components
#'
#' Collapses a Gaussian mixture to at most `k_max` components by a greedy,
#' moment-preserving pairwise merge. At each step the cheapest pair of components
#' is replaced by the single Gaussian that preserves their combined weight, mean
#' and covariance, until the component budget is met. Because every merge is
#' moment-preserving, the reduced mixture has the **same global mean and
#' covariance** as the original, and reducing all the way to a single component
#' returns the moment-matched Gaussian.
#'
#' The merge cost decides which pair is collapsed first. `cost = "kl"` is the
#' Kullback-Leibler upper bound of Runnalls (2007), the standard mixture-
#' reduction criterion; `cost = "cs"` is the closed-form Cauchy-Schwarz
#' divergence between the two-component sub-mixture and its merge (the same
#' Gaussian-product identity that underlies [gmm_divergence()]), scaled by the
#' merged mass. Both are non-negative and zero for identical components.
#'
#' Reduction is the closing operation of a Gaussian-sum filter: repeated
#' [gmm_observe()] / [gmm_affine()] steps under a Gaussian-mixture noise or
#' dynamics model multiply the component count, and `gmm_reduce()` returns it to
#' a fixed budget while preserving the first two moments.
#'
#' @param g A [gmm] (or [gmm_fit]).
#' @param k_max Maximum number of components to keep (a positive integer). When
#'   `k_max` is at least the current component count the mixture is returned
#'   unchanged.
#' @param cost The pairwise merge cost: `"kl"` (the Runnalls Kullback-Leibler
#'   bound, the default) or `"cs"` (the Cauchy-Schwarz divergence).
#' @param ridge_eps Ridge added to each merged covariance. The moment-preserving
#'   merge is positive-definite by construction, so the default `0` keeps the
#'   global moments exact; set a small positive value for extra numerical
#'   headroom in very long reduction chains (at the cost of a tiny moment drift).
#'
#' @returns A [gmm] with at most `k_max` components.
#' @family operators
#' @references Runnalls, A. R. (2007) Kullback-Leibler approach to Gaussian
#'   mixture reduction. *IEEE Transactions on Aerospace and Electronic Systems*
#'   43(3), 989--999. \doi{10.1109/TAES.2007.4383588}
#' @export
#' @examples
#' ## A six-component mixture with three near-duplicate pairs.
#' g <- gmm(
#'   weights = rep(1 / 6, 6),
#'   means = list(c(-4, 0), c(-4, 0.1), c(4, 0), c(4.1, 0), c(0, 5), c(0, 5.1)),
#'   covariances = rep(list(diag(2)), 6)
#' )
#' gmm_reduce(g, k_max = 3L)
gmm_reduce <- function(g, k_max, cost = c("kl", "cs"), ridge_eps = 0) {
  if (!S7::S7_inherits(g, gmm)) {
    cli::cli_abort("`g` must be a {.cls gmm} object.")
  }
  k_max <- as.integer(k_max)
  if (length(k_max) != 1L || is.na(k_max) || k_max < 1L) {
    cli::cli_abort("`k_max` must be a single positive integer.")
  }
  cost <- rlang::arg_match(cost)

  w <- g@weights
  m <- g@means
  S <- g@covariances
  K <- length(w)

  if (k_max >= K) {
    return(gmm(weights = w, means = m, covariances = S,
               name = sprintf("%s (<= %d comp)", g@name, k_max)))
  }

  while (K > k_max) {
    best <- Inf
    bi <- NA_integer_
    bj <- NA_integer_
    best_merge <- NULL
    for (i in seq_len(K - 1L)) {
      for (j in (i + 1L):K) {
        merged <- .merge_two(w[i], m[[i]], S[[i]], w[j], m[[j]], S[[j]])
        c_ij <- if (cost == "kl") {
          .runnalls_cost(w[i], S[[i]], w[j], S[[j]], merged$Sigma)
        } else {
          .cs_cost(w[i], m[[i]], S[[i]], w[j], m[[j]], S[[j]],
                   merged$mu, merged$Sigma)
        }
        if (c_ij < best) {
          best <- c_ij
          bi <- i
          bj <- j
          best_merge <- merged
        }
      }
    }
    w[bi] <- best_merge$w
    m[[bi]] <- best_merge$mu
    S[[bi]] <- symmetrise(ridge(best_merge$Sigma, ridge_eps))
    w <- w[-bj]
    m <- m[-bj]
    S <- S[-bj]
    K <- K - 1L
  }

  w <- w / sum(w)
  gmm(weights = w, means = m, covariances = S,
      name = sprintf("%s (reduced to %d comp)", g@name, K))
}
