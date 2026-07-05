## Internal utilities. Not exported.

## Stable row-wise log-sum-exp on an n-by-K matrix.
## Callers must supply an actual matrix; we never guess the orientation.
## NaN entries are treated as -Inf (an excluded term) everywhere, not just
## in the row-max; a row whose maximum is +Inf returns +Inf (the sum of a
## divergent term is divergent, not zero).
logsumexp_rows <- function(m) {
  if (!is.matrix(m)) {
    cli::cli_abort("`m` must be a matrix; got {.cls {class(m)[1L]}}.")
  }
  n <- nrow(m)
  if (n == 0L) return(numeric(0L))
  m_safe <- m
  m_safe[is.nan(m_safe)] <- -Inf
  ## Row max as a K-pass pmax reduce: this sits on every E-step's hot path,
  ## where the per-row vapply loop dominated at large n.
  mx <- m_safe[, 1L]
  for (k in seq_len(ncol(m_safe))[-1L]) mx <- pmax(mx, m_safe[, k])
  out <- rep(-Inf, n)
  out[is.infinite(mx) & mx > 0] <- Inf
  finite <- is.finite(mx)
  if (any(finite)) {
    out[finite] <- mx[finite] +
      log(rowSums(exp(m_safe[finite, , drop = FALSE] - mx[finite])))
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
## `epsilon` is absolute. Callers that regularise across EM iterations pass
## a DATA-scaled epsilon (see `.data_scaled_eps()`), so the floor is
## invariant to the data's units but constant across iterations -- a ridge
## relative to the component's own diagonal would shrink together with a
## collapsing component and stop flooring exactly when needed.
ridge <- function(S, epsilon = 1e-6) {
  S + diag(epsilon, nrow = nrow(S))
}

## Scale a user-facing relative ridge epsilon by the data's covariance
## scale (mean diagonal), so the same default regularises identically at
## data scale 1e-8 and 1e+8 while remaining a fixed floor within a fit.
.data_scaled_eps <- function(epsilon, data_cov_diag_mean) {
  s <- data_cov_diag_mean
  if (!is.finite(s) || s <= 0) s <- 1
  epsilon * s
}

## Force a covariance to be symmetric (numerical hygiene after EM updates).
symmetrise <- function(S) 0.5 * (S + t(S))

## Single metadata policy for the closed-form operators: every operator
## result carries (a) the source fit's quality certificate, unchanged, and
## (b) a `provenance` character vector with one entry per operation, plus
## any operator-specific extras. Nothing else is merged forward, so stale
## keys from earlier operations cannot masquerade as current.
.op_result_meta <- function(g, op, extra = list()) {
  meta <- list()
  if (!is.null(g@metadata$quality)) meta$quality <- g@metadata$quality
  meta$provenance <- c(g@metadata$provenance %||% character(0L), op)
  utils::modifyList(meta, extra)
}

## Operator result names cap their nesting: a chained pipeline previously
## grew names like observe(affine(observe(...))) without bound.
.op_name <- function(op, g) {
  nm <- g@name
  if (is.null(nm) || !nzchar(nm) || nchar(nm) > 40L) nm <- "..."
  sprintf("%s(%s)", op, nm)
}

## One-shot advisory when a downstream verb consumes a fit whose quality
## certificate is flagged. Plain mixtures (no certificate) pass silently.
.check_quality <- function(g, verb) {
  q <- g@metadata$quality
  if (is.null(q)) return(invisible(NULL))
  low_ess <- is.finite(q$ess_relative %||% NA_real_) && q$ess_relative < 0.05
  bad <- isTRUE(q$degenerate) || isFALSE(q$converged) || low_ess
  if (bad) {
    rlang::inform(
      c(sprintf(
        "`%s` received a fit whose quality certificate is flagged (converged = %s, degenerate = %s, relative ESS = %s).",
        verb, format(q$converged), format(isTRUE(q$degenerate)),
        format(signif(q$ess_relative %||% NA_real_, 2L))),
        "i" = "Downstream quantities condition on this proxy; consider refitting before relying on them."),
      class = "proxymix_low_quality",
      .frequency = "once",
      .frequency_id = paste0("proxymix_low_quality_", verb)
    )
  }
  invisible(NULL)
}

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
