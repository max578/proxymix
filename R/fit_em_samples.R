## Regime (ii) of Hoek and Elliott (2024).
##
## Classical expectation-maximisation on independent samples.

#' Classical EM fit on samples
#'
#' Implements regime (ii) of Hoek and Elliott (2024). Runs the textbook
#' expectation-maximisation algorithm for Gaussian mixtures on the supplied
#' samples, with diagonal ridge regularisation for numerical stability,
#' optional multi-start, and monotone-log-likelihood checking.
#'
#' @param target A [gmm_target] carrying an `n` by `p` `samples` matrix.
#' @param N Number of mixture components.
#' @param init A [gmm] initialisation, or `NULL` to use [init_kmeans()].
#' @param max_iter Maximum number of EM iterations.
#' @param tol Relative-log-likelihood convergence tolerance.
#' @param ridge_eps Ridge added to each component covariance at every M-step.
#' @param n_starts Number of multi-start initialisations (only when `init`
#'   is `NULL`). The best fit by final log-likelihood is returned.
#' @param anneal Logical. If `TRUE`, a deterministic-annealing warm-start
#'   (see [gmm_anneal_path()]) replaces the multi-start: the components are
#'   annealed from a high temperature down to one, and the resulting
#'   parameters seed a single final (cold) EM polish. This attacks the
#'   local-optima sensitivity of cold EM at the cost of the schedule length.
#'   Defaults to `FALSE` (cold best-of-`n_starts`).
#' @param temp_schedule Optional numeric vector of descending temperatures for
#'   the annealing warm-start. `NULL` (the default) uses a geometric schedule
#'   from `10` down to `1` in covariance-whitened units. Ignored when
#'   `anneal = FALSE`.
#' @param seed Optional integer seed for the annealing perturbations (the
#'   warm-start is deterministic given a seed). Ignored when `anneal = FALSE`.
#' @param canonicalise Logical. If `TRUE` (the default), the fitted
#'   mixture is post-processed by [gmm_canonicalise()] so that components
#'   are sorted by descending weight and (as a tiebreaker) by descending
#'   `||mu||`.
#'
#' @returns A [gmm_fit] with `regime = "sample"`. When `anneal = TRUE` the
#'   diagnostics list also carries `annealed = TRUE` and the `temp_schedule`
#'   used.
#' @family fitting
#' @export
#' @examples
#' x <- matrix(stats::rnorm(200), ncol = 2)
#' tgt <- gmm_target_from_samples(x)
#' fit <- fit_em_samples(tgt, N = 2L, max_iter = 30L, n_starts = 2L)
#' fit@diagnostics$loglik_final
fit_em_samples <- function(target, N = 2L,
                           init = NULL,
                           max_iter = 100L,
                           tol = 1e-6,
                           ridge_eps = 1e-6,
                           n_starts = 5L,
                           anneal = FALSE,
                           temp_schedule = NULL,
                           seed = NULL,
                           canonicalise = TRUE) {
  if (!S7::S7_inherits(target, gmm_target)) {
    cli::cli_abort("`target` must be a {.cls gmm_target} object.")
  }
  if (is.null(target@samples)) {
    cli::cli_abort("regime {.val sample} requires `target@samples` to be supplied.")
  }
  N <- as.integer(N)
  if (N < 1L) cli::cli_abort("`N` must be a positive integer.")
  samples <- target@samples
  n <- nrow(samples)
  p <- ncol(samples)

  if (isTRUE(anneal)) {
    warm <- .anneal_em_warmstart(samples, rw = rep(1 / n, n), N = N,
                                 temp_schedule = temp_schedule,
                                 ridge_eps = ridge_eps, seed = seed)
    fit <- em_samples_one_run(samples, warm, target,
                              max_iter, tol, ridge_eps,
                              call = match.call())
    fit@diagnostics[["annealed"]] <- TRUE
    fit@diagnostics[["temp_schedule"]] <- warm@metadata$temp_schedule
    return(if (isTRUE(canonicalise)) gmm_canonicalise(fit) else fit)
  }

  if (!is.null(init)) {
    fit <- em_samples_one_run(samples, init, target,
                              max_iter, tol, ridge_eps,
                              call = match.call())
    return(if (isTRUE(canonicalise)) gmm_canonicalise(fit) else fit)
  }

  ## Multi-start: kmeans + (n_starts - 1) random restarts. A supplied `seed`
  ## drives both the kmeans start and the restart draws, so the multi-start
  ## is reproducible end-to-end; without one, the restarts keep their
  ## historical fixed seeds (identical across calls, documented) and kmeans
  ## uses the ambient stream.
  inits <- vector("list", max(1L, n_starts))
  inits[[1L]] <- if (is.null(seed)) {
    init_kmeans(samples, N = N, ridge_eps = ridge_eps)
  } else {
    withr::with_seed(as.integer(seed),
                     init_kmeans(samples, N = N, ridge_eps = ridge_eps))
  }
  if (n_starts > 1L) {
    mu_global <- colMeans(samples)
    sd_global <- sqrt(diag(stats::cov(samples)))
    for (i in seq_len(n_starts - 1L)) {
      inits[[1L + i]] <- init_random(
        N = N, p = p, centre = mu_global,
        scale = mean(sd_global), sigma_diag = mean(sd_global^2),
        seed = if (is.null(seed)) 1000L + i else as.integer(seed) + i
      )
    }
  }

  fits <- lapply(inits, function(this_init) {
    tryCatch(
      em_samples_one_run(samples, this_init, target,
                         max_iter, tol, ridge_eps,
                         call = match.call()),
      error = function(e) e
    )
  })
  ok <- vapply(fits, function(f) S7::S7_inherits(f, gmm_fit), logical(1L))
  if (!any(ok)) {
    cli::cli_abort("all EM starts failed.")
  }
  fits <- fits[ok]
  scores <- vapply(fits, function(f) f@diagnostics$loglik_final,
                   FUN.VALUE = numeric(1L))
  best <- fits[[which.max(scores)]]
  if (isTRUE(canonicalise)) gmm_canonicalise(best) else best
}

## Single EM run, given a fully-formed initialisation.
em_samples_one_run <- function(samples, init, target,
                               max_iter, tol, ridge_eps, call) {
  N <- gmm_n_components(init)
  n <- nrow(samples)
  p <- ncol(samples)

  weights <- init@weights
  means <- init@means
  ## Data-scaled ridge: invariant to the data's units, constant within the
  ## fit (see `.data_scaled_eps()`).
  ridge_eps <- .data_scaled_eps(ridge_eps, mean(diag(stats::cov(samples))))
  covs <- lapply(init@covariances, function(S) ridge(S, ridge_eps))

  loglik_trace <- numeric(0L)
  converged <- FALSE
  it <- 0L
  n_reseeds <- 0L

  for (it in seq_len(max_iter)) {
    log_resp_unnorm <- gmm_log_unnorm(samples, weights, means, covs)
    ll_row <- logsumexp_rows(log_resp_unnorm)
    loglik <- sum(ll_row)
    loglik_trace <- c(loglik_trace, loglik)

    if (it > 1L) {
      delta <- abs(loglik_trace[it] - loglik_trace[it - 1L]) /
        (abs(loglik_trace[it - 1L]) + 1e-12)
      if (delta < tol) {
        converged <- TRUE
        break
      }
    }

    log_resp <- log_resp_unnorm - ll_row
    resp <- exp(log_resp)
    Nk <- colSums(resp)
    empty <- Nk < 1e-12
    ## A dead component must get its WEIGHT reset too: `Nk / n` gives an
    ## exact zero, `log(0) = -Inf` zeroes its responsibilities forever, and
    ## a reseeded mean would be dead on arrival (the reseed then re-fires
    ## every iteration). Granting it one observation's worth of mass lets
    ## the next E-step actually reach it.
    Nk_adj <- Nk
    Nk_adj[empty] <- 1
    n_reseeds <- n_reseeds + sum(empty)
    weights <- as.numeric(Nk_adj / sum(Nk_adj))
    for (k in seq_len(N)) {
      if (empty[k]) {
        ## Re-seed at the data point the current mixture explains worst
        ## (the lowest per-row log-likelihood), at data scale.
        worst <- which.min(ll_row)
        means[[k]] <- as.numeric(samples[worst, ])
        covs[[k]] <- ridge(stats::cov(samples), ridge_eps)
        next
      }
      mu_new <- as.numeric(colSums(resp[, k] * samples) / Nk[k])
      diff <- samples - matrix(mu_new, nrow = n, ncol = p, byrow = TRUE)
      S_new <- crossprod(diff * sqrt(resp[, k])) / Nk[k]
      means[[k]] <- mu_new
      covs[[k]] <- symmetrise(ridge(S_new, ridge_eps))
    }
  }

  ## The EM objective is monotone under exact updates. An empty-component
  ## reseed is a deliberate restart-like intervention that legitimately dips
  ## the objective before recovering, so the invariant is only asserted for
  ## reseed-free runs; reseed counts are recorded in the diagnostics either
  ## way.
  if (n_reseeds == 0L && length(loglik_trace) > 2L) {
    worst_drop <- min(diff(loglik_trace))
    if (worst_drop < -1e-6 * max(1, abs(loglik_trace[1L]))) {
      cli::cli_warn(
        c("The EM log-likelihood decreased by {signif(-worst_drop, 3)} during fitting.",
          "i" = "This usually signals ridge regularisation interacting with a near-singular component."),
        class = "proxymix_nonmonotone"
      )
    }
  }

  ## Number of free parameters for BIC/AIC:
  ##   N - 1 (weights) + N * p (means) + N * p * (p + 1) / 2 (covariances).
  n_params <- (N - 1L) + N * p + N * p * (p + 1L) / 2L
  loglik_final <- if (length(loglik_trace) > 0L) loglik_trace[length(loglik_trace)] else NA_real_
  bic <- -2 * loglik_final + n_params * log(n)
  aic <- -2 * loglik_final + 2 * n_params

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
      n_reseeds = n_reseeds,
      bic = bic,
      aic = aic
    ),
    converged = converged,
    iterations = as.integer(it),
    call = call,
    name = sprintf("em_samples[N=%d] on %s", N, target@name),
    metadata = list(quality = list(
      regime = "sample",
      converged = converged,
      degenerate = FALSE,
      ess = NA_real_,
      ess_relative = NA_real_,
      min_component_ess = NA_real_,
      max_weight = NA_real_,
      support_fraction = NA_real_,
      kld_final = NA_real_,
      validation_gap = NA_real_
    ))
  )
}
