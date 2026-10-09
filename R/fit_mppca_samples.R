## Mixture of probabilistic principal component analysers (Tipping and
## Bishop, 1999), fitted to samples by EM.

#' Mixture of probabilistic PCA fitted to samples
#'
#' Fits a Gaussian mixture whose component covariances are restricted to
#' \eqn{\Sigma_k = W_k W_k^\top + \sigma_k^2 I}, where \eqn{W_k} has `q`
#' columns (Tipping and Bishop, 1999). Each component varies mainly in a
#' `q`-dimensional subspace, with equal noise variance in every other
#' direction. With `p` variables a component has \eqn{pq - q(q-1)/2 + 1}
#' covariance parameters instead of \eqn{p(p+1)/2}, so the model can be fitted
#' when a component holds fewer points than there are variables, where a
#' full-covariance fit ([fit_em_samples()]) is singular.
#'
#' The fit is by EM. Given the responsibilities, the update of
#' \eqn{(W_k, \sigma_k^2)} is the closed-form maximum-likelihood solution of
#' probabilistic PCA on the responsibility-weighted covariance \eqn{S_k}: the
#' columns of \eqn{W_k} are the leading `q` eigenvectors of \eqn{S_k}, scaled
#' by \eqn{(\lambda_j - \sigma_k^2)^{1/2}}, and \eqn{\sigma_k^2} is the mean
#' of the remaining \eqn{p - q} eigenvalues. Starts are k-means partitions.
#'
#' The result is an ordinary [gmm_fit], so every operator of the package
#' (conditioning, marginals, imputation, sampling) applies to it.
#'
#' @param target A [gmm_target] carrying an `n` by `p` `samples` matrix.
#' @param N Number of mixture components.
#' @param q Dimension of each component's principal subspace, an integer with
#'   `1 <= q < p`. Compare values of `q` (and `N`) with the BIC in the
#'   diagnostics.
#' @param n_starts Number of k-means starts. The fit with the largest final
#'   log-likelihood is returned.
#' @param max_iter Maximum number of EM iterations per start.
#' @param tol Relative log-likelihood convergence tolerance.
#' @param seed Optional integer seed for the k-means starts.
#' @param canonicalise Logical. If `TRUE` (the default), components are sorted
#'   by [gmm_canonicalise()].
#'
#' @returns A [gmm_fit] with `regime = "sample"`. The diagnostics carry the
#'   log-likelihood trace and final value, the number of free parameters, BIC
#'   and AIC. `metadata$mppca` holds `q`, the loading matrices `loadings` and
#'   the noise variances `noise_var`, in component order.
#' @references
#' Tipping, M. E. and Bishop, C. M. (1999). Mixtures of probabilistic
#' principal component analyzers. *Neural Computation*, 11(2), 443-482.
#' \doi{10.1162/089976699300016728}
#' @family fitting
#' @export
#' @examples
#' set.seed(1)
#' p <- 10
#' w1 <- qr.Q(qr(matrix(rnorm(p * 2), p)))
#' w2 <- qr.Q(qr(matrix(rnorm(p * 2), p)))
#' x <- rbind(matrix(rnorm(60 * 2), 60) %*% t(3 * w1) + matrix(rnorm(60 * p, sd = 0.3), 60),
#'            matrix(rnorm(60 * 2), 60) %*% t(3 * w2) + matrix(rnorm(60 * p, sd = 0.3), 60) + 4)
#' fit <- fit_mppca_samples(gmm_target_from_samples(x), N = 2L, q = 2L, seed = 1L)
#' fit@metadata$mppca$noise_var
#' fit@diagnostics$bic
fit_mppca_samples <- function(target, N = 2L, q = 1L,
                              n_starts = 5L,
                              max_iter = 500L,
                              tol = 1e-8,
                              seed = NULL,
                              canonicalise = TRUE) {
  if (!S7::S7_inherits(target, gmm_target)) {
    cli::cli_abort("`target` must be a {.cls gmm_target} object.")
  }
  if (is.null(target@samples)) {
    cli::cli_abort("`target@samples` must be supplied.")
  }
  samples <- target@samples
  n <- nrow(samples)
  p <- ncol(samples)
  N <- as.integer(N)
  q <- as.integer(q)
  if (length(N) != 1L || is.na(N) || N < 1L) cli::cli_abort("`N` must be a positive integer.")
  if (length(q) != 1L || is.na(q) || q < 1L || q >= p) {
    cli::cli_abort("`q` must be an integer with 1 <= q < p (here p = {p}).")
  }
  if (n < N) cli::cli_abort("`target@samples` has {n} row{?s}, fewer than `N` = {N}.")

  starts <- lapply(seq_len(max(1L, as.integer(n_starts))), function(i) {
    km <- if (is.null(seed)) {
      stats::kmeans(samples, centers = N, nstart = 1L)
    } else {
      withr::with_seed(as.integer(seed) + i - 1L,
                       stats::kmeans(samples, centers = N, nstart = 1L))
    }
    km$cluster
  })
  fits <- lapply(starts, function(cl) {
    tryCatch(mppca_one_run(samples, cl, N, q, max_iter, tol, target, match.call()),
             error = function(e) e)
  })
  ok <- vapply(fits, function(f) S7::S7_inherits(f, gmm_fit), logical(1L))
  if (!any(ok)) cli::cli_abort("all EM starts failed.")
  fits <- fits[ok]
  scores <- vapply(fits, function(f) f@diagnostics$loglik_final, numeric(1L))
  best <- fits[[which.max(scores)]]
  if (!isTRUE(canonicalise)) return(best)
  ## Keep the loadings in the same order as the canonicalised components.
  out <- gmm_canonicalise(best)
  ord <- vapply(out@means, function(m) {
    which.min(vapply(best@means, function(b) sum((b - m)^2), numeric(1L)))
  }, integer(1L))
  out@metadata$mppca$loadings <- best@metadata$mppca$loadings[ord]
  out@metadata$mppca$noise_var <- best@metadata$mppca$noise_var[ord]
  out
}

## Single EM run from a hard partition.
mppca_one_run <- function(samples, cluster, N, q, max_iter, tol, target, call) {
  n <- nrow(samples)
  p <- ncol(samples)
  resp <- outer(cluster, seq_len(N), "==") * 1
  var_floor <- 1e-10 * mean(apply(samples, 2L, stats::var))
  loglik_trace <- numeric(0L)
  converged <- FALSE
  it <- 0L
  for (it in seq_len(max_iter)) {
    Nk <- colSums(resp)
    if (any(Nk < 1e-12)) cli::cli_abort("a component lost all its responsibility.")
    weights <- as.numeric(Nk / n)
    means <- vector("list", N)
    covs <- vector("list", N)
    loadings <- vector("list", N)
    noise_var <- numeric(N)
    for (k in seq_len(N)) {
      mu <- as.numeric(colSums(resp[, k] * samples) / Nk[k])
      diff <- samples - matrix(mu, nrow = n, ncol = p, byrow = TRUE)
      S <- crossprod(diff * sqrt(resp[, k])) / Nk[k]
      e <- eigen(symmetrise(S), symmetric = TRUE)
      s2 <- max(mean(e$values[(q + 1L):p]), var_floor)
      W <- e$vectors[, seq_len(q), drop = FALSE] %*%
        diag(sqrt(pmax(e$values[seq_len(q)] - s2, 0)), q)
      means[[k]] <- mu
      loadings[[k]] <- W
      noise_var[k] <- s2
      covs[[k]] <- symmetrise(tcrossprod(W) + diag(s2, p))
    }
    log_unnorm <- gmm_log_unnorm(samples, weights, means, covs)
    ll_row <- logsumexp_rows(log_unnorm)
    loglik_trace <- c(loglik_trace, sum(ll_row))
    resp <- exp(log_unnorm - ll_row)
    if (it > 1L) {
      delta <- abs(loglik_trace[it] - loglik_trace[it - 1L]) /
        (abs(loglik_trace[it - 1L]) + 1e-12)
      if (delta < tol) {
        converged <- TRUE
        break
      }
    }
  }

  ## Free parameters: N - 1 weights, N * p means, and per component
  ## p * q - q * (q - 1) / 2 loadings (up to rotation) plus one noise variance.
  n_params <- (N - 1L) + N * p + N * (p * q - q * (q - 1L) / 2 + 1)
  loglik_final <- loglik_trace[length(loglik_trace)]
  gmm_fit(
    weights = weights,
    means = means,
    covariances = covs,
    target = target,
    regime = "sample",
    diagnostics = list(
      loglik_trace = loglik_trace,
      loglik_final = loglik_final,
      n_used = n,
      n_params = n_params,
      bic = -2 * loglik_final + n_params * log(n),
      aic = -2 * loglik_final + 2 * n_params
    ),
    converged = converged,
    iterations = as.integer(it),
    call = call,
    name = sprintf("mppca_samples[N=%d, q=%d] on %s", N, q, target@name),
    metadata = list(mppca = list(q = q, loadings = loadings, noise_var = noise_var))
  )
}
