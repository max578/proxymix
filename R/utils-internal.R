## Internal utilities. Not exported.

## Stable row-wise log-sum-exp on an n-by-K matrix.
## Callers must supply an actual matrix; we never guess the orientation.
logsumexp_rows <- function(m) {
  if (!is.matrix(m)) {
    cli::cli_abort("`m` must be a matrix; got {.cls {class(m)[1L]}}.")
  }
  n <- nrow(m)
  if (n == 0L) return(numeric(0L))
  m_safe <- m
  m_safe[is.nan(m_safe)] <- -Inf
  mx <- vapply(seq_len(n), function(i) max(m_safe[i, ]), numeric(1L))
  finite <- is.finite(mx)
  out <- rep(-Inf, n)
  if (any(finite)) {
    out[finite] <- mx[finite] +
      log(rowSums(exp(m[finite, , drop = FALSE] - mx[finite])))
  }
  out
}

## Unnormalised per-observation log-mixture-component evaluations.
## Returns an n-by-K matrix where entry [i, k] is
##   log(weights[k]) + log N(x[i, ] | means[[k]], covariances[[k]]).
## Always a matrix, including when n == 1 or K == 1.
gmm_log_unnorm <- function(x, weights, means, covariances) {
  K <- length(weights)
  n <- nrow(x)
  out <- matrix(0, nrow = n, ncol = K)
  log_w <- log(weights)
  for (k in seq_len(K)) {
    out[, k] <- log_w[k] +
      mvnfast::dmvn(x, mu = means[[k]], sigma = covariances[[k]], log = TRUE)
  }
  out
}

## Closed-form Gaussian-to-Gaussian KL.
##   KL(N(mu_a, S_a) || N(mu_b, S_b))
##   = 0.5 * ( tr(S_b^{-1} S_a) + (mu_b - mu_a)' S_b^{-1} (mu_b - mu_a)
##             - p + log det(S_b) - log det(S_a) ).
kl_gauss <- function(mu_a, S_a, mu_b, S_b) {
  p <- length(mu_a)
  L_b <- chol(S_b)
  S_b_inv <- chol2inv(L_b)
  log_det_S_b <- 2 * sum(log(diag(L_b)))
  L_a <- tryCatch(chol(S_a), error = function(e) NULL)
  log_det_S_a <- if (!is.null(L_a)) {
    2 * sum(log(diag(L_a)))
  } else {
    as.numeric(determinant(S_a, logarithm = TRUE)$modulus)
  }
  diff <- mu_b - mu_a
  0.5 * (sum(diag(S_b_inv %*% S_a)) +
           sum(diff * as.numeric(S_b_inv %*% diff)) -
           p + log_det_S_b - log_det_S_a)
}

## Ridge a matrix to make it positive-definite-ish: S + epsilon * I.
ridge <- function(S, epsilon = 1e-6) {
  S + diag(epsilon, nrow = nrow(S))
}

## Test a covariance for positive definiteness by attempting chol().
is_pd <- function(S) {
  inherits(tryCatch(chol(S), error = function(e) e), "matrix") ||
    isTRUE(tryCatch({ chol(S); TRUE }, error = function(e) FALSE))
}

## Force a covariance to be symmetric (numerical hygiene after EM updates).
symmetrise <- function(S) 0.5 * (S + t(S))

## Coerce x to a matrix with `p` columns. Vectors of length p become a single
## row; vectors of length n*p become an n-by-p matrix only if n*p divides
## cleanly *and* the user passed a matrix-like coercion intent — to avoid
## that ambiguity we require either a matrix or a length-p vector.
as_input_matrix <- function(x, p) {
  if (is.matrix(x)) {
    if (ncol(x) != p) {
      cli::cli_abort("`x` must have {p} column(s), not {ncol(x)}.")
    }
    return(x)
  }
  if (is.numeric(x)) {
    if (length(x) != p) {
      cli::cli_abort("vector `x` must have length {p}, not {length(x)}.")
    }
    return(matrix(x, nrow = 1L))
  }
  cli::cli_abort("`x` must be a numeric matrix or vector.")
}
