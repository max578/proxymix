## Regime (iii) of Hoek and Elliott (2024): the wedge.
##
## Importance-sampled KLD-EM minimises KL(f || g_theta) when the target f
## can be evaluated point-wise but not (cheaply) sampled.

#' Importance-sampled KLD-EM fit (the wedge)
#'
#' Implements regime (iii) of Hoek and Elliott (2024). Minimises
#' `KL(f || g_theta)` where `f` is supplied as an evaluable log-density on
#' the target, via expectation-maximisation against importance-sampled
#' draws from a user-chosen proposal `q`.
#'
#' The Monte Carlo draws from `q` are computed once at the start; the
#' resulting self-normalised importance-sampling weights are reused at every
#' EM iteration. Adaptive importance sampling — redrawing each round — is
#' a Tier-3 deferral.
#'
#' Since v0.1.1 the function also draws an *independent* validation IS
#' sample when `validation_size > 0` and reports its own KLD estimate,
#' effective sample size, and largest weight share. This lets users tell
#' the difference between in-sample EM overfit to one particular IS draw
#' and a fit that generalises across independent IS draws.
#'
#' When the target's `normalised` property is `FALSE` or `NA`, the
#' importance-sampled `kld_final` and `kld_trace` measure
#' \eqn{\widehat{KL}(f \Vert g) - \log Z(f)} rather than the absolute
#' divergence. The fit's diagnostics list records this via
#' `kld_is_shifted = TRUE` and a `kld_shift_explanation` string. When the
#' target also supplies a finite `log_normalizer`, a corrected absolute
#' estimate is reported as `kld_final_absolute`.
#'
#' @param target A [gmm_target] with a non-NULL `log_density`.
#' @param N Number of mixture components.
#' @param proposal An [is_proposal]. When `NULL` (the default) the proposal is
#'   chosen automatically: a support-matched [is_uniform()] when the target
#'   declares a bounded or one-sided `support`, otherwise a multivariate-t with
#'   `df = 5` in `target@n_dim` dimensions. The automatic choice is announced
#'   with a one-line message so it is never silent.
#' @param is_size Number of importance-sampling draws used for fitting.
#' @param init A [gmm] initialisation, or `NULL` to use a kmeans pass on
#'   the importance-resampled draws.
#' @param max_iter Maximum number of EM iterations.
#' @param tol Convergence tolerance on the relative change in the
#'   importance-sampled KLD estimate.
#' @param ridge_eps Ridge added to each component covariance at every
#'   M-step.
#' @param min_ess Minimum effective sample size below which a warning is
#'   issued.
#' @param seed Optional integer seed used when drawing the fitting IS
#'   sample.
#' @param validation_size Number of independent importance-sampling draws
#'   to use for held-out validation. Default `0L` (no validation split).
#'   Set to a positive integer (typically `is_size / 2` to `is_size`) to
#'   enable validation diagnostics.
#' @param validation_proposal Optional [is_proposal] for the validation
#'   sample. Defaults to the same proposal used for fitting.
#' @param validation_seed Optional integer seed used when drawing the
#'   validation sample. Defaults to `seed + 1L` when `seed` is supplied,
#'   `NULL` otherwise.
#' @param support_warn Logical. If `TRUE` (the default), issue a warning
#'   when more than 5% of IS draws receive non-finite weights (typically
#'   because the proposal does not dominate the target's support).
#' @param canonicalise Logical. If `TRUE` (the default), the fitted
#'   mixture is post-processed by [gmm_canonicalise()].
#'
#' @returns A [gmm_fit] with `regime = "kld"`. The diagnostics list
#'   contains, among others, `kld_trace`, `kld_final`,
#'   `kld_is_shifted`, `kld_final_absolute` (when computable), `ess`,
#'   `ess_relative` (`ess / is_size`), `max_weight`, `support_fraction`,
#'   `mc_se_kld`, `validation_kld`, `validation_ess`, and
#'   `validation_max_weight`.
#' @family fitting
#' @export
#' @examples
#' tgt <- banana_target()
#' q <- is_mvt(n_dim = 2L, mean = c(0, 0),
#'             sigma = 4 * diag(2), df = 5)
#' fit <- fit_kld_em(tgt, N = 3L, proposal = q,
#'                   is_size = 1500L, max_iter = 25L, seed = 1L,
#'                   validation_size = 1500L)
#' fit@diagnostics$kld_final
#' fit@diagnostics$validation_kld
fit_kld_em <- function(target,
                       N = 3L,
                       proposal = NULL,
                       is_size = 5000L,
                       init = NULL,
                       max_iter = 100L,
                       tol = 1e-5,
                       ridge_eps = 1e-6,
                       min_ess = 50,
                       seed = NULL,
                       validation_size = 0L,
                       validation_proposal = NULL,
                       validation_seed = NULL,
                       support_warn = TRUE,
                       canonicalise = TRUE) {
  if (!S7::S7_inherits(target, gmm_target)) {
    cli::cli_abort("`target` must be a {.cls gmm_target} object.")
  }
  if (is.null(target@log_density)) {
    cli::cli_abort("regime {.val kld} requires `target@log_density` to be supplied.")
  }
  N <- as.integer(N)
  is_size <- as.integer(is_size)
  validation_size <- as.integer(validation_size)
  p <- target@n_dim

  if (is.null(proposal)) {
    proposal <- .support_matched_proposal(target)
    if (is.null(proposal)) {
      proposal <- is_mvt(n_dim = p,
                         mean = if (!is.null(target@samples)) colMeans(target@samples) else rep(0, p),
                         sigma = if (!is.null(target@samples)) ridge(stats::cov(target@samples), 1e-3) else diag(p),
                         df = 5)
    } else {
      cli::cli_inform(c(
        "Auto-selected a support-matched proposal ({.val {proposal@name}}) for the declared target support.",
        "i" = "Pass an explicit {.arg proposal} to override."
      ))
    }
  }
  if (!S7::S7_inherits(proposal, is_proposal)) {
    cli::cli_abort("`proposal` must be a {.cls is_proposal} object.")
  }
  if (proposal@n_dim != p) {
    cli::cli_abort("`proposal@n_dim` ({proposal@n_dim}) must equal `target@n_dim` ({p}).")
  }

  ## ---- Draw fitting IS sample and compute self-normalised weights ----
  draws <- draw_is_weights(target, proposal, is_size, seed)
  x <- draws$x
  log_f <- draws$log_f
  log_W <- draws$log_W
  W <- draws$W
  ess <- draws$ess
  max_weight <- draws$max_weight
  support_fraction <- draws$support_fraction

  if (ess < min_ess) {
    cli::cli_warn(c(
      "Effective sample size is low: ESS = {round(ess, 1)} out of {is_size}.",
      "i" = "Consider a heavier-tailed proposal or more IS draws."
    ))
  }
  if (isTRUE(support_warn) && support_fraction < 0.95) {
    cli::cli_warn(c(
      "Importance proposal does not dominate target support: only {round(100 * support_fraction, 1)}% of IS draws received finite weight.",
      "i" = "Verify {.code log_density(x) - log q(x)} is finite on the target's support; consider a wider or heavier-tailed proposal."
    ))
  }

  ## ---- Initialisation ----
  if (is.null(init)) {
    ## Bootstrap-resample x by IS weights to obtain a pseudo-sample with
    ## approximately target distribution, then kmeans on it.
    idx <- sample.int(is_size, size = min(is_size, 5L * N * 50L),
                      replace = TRUE, prob = W)
    pseudo <- x[idx, , drop = FALSE]
    init <- tryCatch(init_kmeans(pseudo, N = N, ridge_eps = ridge_eps),
                     error = function(e) NULL)
    if (is.null(init)) {
      init <- init_random(N = N, p = p,
                          centre = if (!is.null(target@samples)) colMeans(target@samples) else rep(0, p),
                          scale = 1, sigma_diag = 1, seed = 42L)
    }
  }
  if (!S7::S7_inherits(init, gmm)) {
    cli::cli_abort("`init` must be a {.cls gmm} initialisation.")
  }

  weights <- init@weights
  means <- init@means
  covs <- lapply(init@covariances, function(S) ridge(S, ridge_eps))

  ## ---- EM iterations ----
  kld_trace <- numeric(0L)
  weighted_obj_trace <- numeric(0L)
  converged <- FALSE
  it <- 0L
  for (it in seq_len(max_iter)) {
    log_resp_unnorm <- gmm_log_unnorm(x, weights, means, covs)
    log_g <- logsumexp_rows(log_resp_unnorm)
    log_resp <- log_resp_unnorm - log_g
    resp <- exp(log_resp)

    ## Two paired bookkeeping quantities, both evaluated under fixed W:
    ##   * the IS-estimated KLD,
    ##     KL(f || g) ~ sum_n W_n (log f(x_n) - log g(x_n));
    ##   * the IS-weighted Q-objective minimised by EM,
    ##     Q(theta) = sum_n W_n log g(x_n).
    kld <- sum(W * (log_f - log_g))
    weighted_obj <- sum(W * log_g)
    kld_trace <- c(kld_trace, kld)
    weighted_obj_trace <- c(weighted_obj_trace, weighted_obj)

    if (it > 1L) {
      delta <- abs(kld_trace[it] - kld_trace[it - 1L]) /
        (abs(kld_trace[it - 1L]) + 1e-12)
      if (delta < tol) {
        converged <- TRUE
        break
      }
    }

    ## M-step.
    W_resp <- resp * W
    Nk_W <- colSums(W_resp)
    weights <- pmax(as.numeric(Nk_W), 1e-300)
    weights <- weights / sum(weights)
    for (k in seq_len(N)) {
      if (Nk_W[k] < 1e-12) {
        ## Empty component: re-seed at IS-weighted random sample.
        idx <- sample.int(is_size, size = 1L, prob = W)
        means[[k]] <- as.numeric(x[idx, ])
        covs[[k]] <- ridge(diag(1, p), ridge_eps)
        next
      }
      mu_new <- as.numeric(colSums(W_resp[, k] * x) / Nk_W[k])
      diff <- x - matrix(mu_new, nrow = is_size, ncol = p, byrow = TRUE)
      S_new <- crossprod(diff * sqrt(W_resp[, k])) / Nk_W[k]
      means[[k]] <- mu_new
      covs[[k]] <- symmetrise(ridge(S_new, ridge_eps))
    }
  }

  ## ---- Final fitting-sample diagnostics ----
  log_resp_unnorm <- gmm_log_unnorm(x, weights, means, covs)
  log_g <- logsumexp_rows(log_resp_unnorm)
  d_n <- log_f - log_g
  finite_d <- is.finite(d_n) & is.finite(W) & W > 0
  kld_final <- sum(W[finite_d] * d_n[finite_d])
  mc_se_kld <- sqrt(sum(W[finite_d]^2 * (d_n[finite_d] - kld_final)^2))

  norm_status <- target@normalised
  kld_is_shifted <- !isTRUE(norm_status)
  kld_shift_explanation <- if (kld_is_shifted) {
    "target@normalised is not TRUE; reported KLD is f-evaluated up to an additive log Z(f) offset."
  } else {
    NA_character_
  }
  kld_final_absolute <- if (!kld_is_shifted) {
    kld_final
  } else if (is.finite(target@log_normalizer)) {
    kld_final + target@log_normalizer
  } else {
    NA_real_
  }

  ## ---- Validation-sample diagnostics ----
  validation_diag <- if (validation_size > 0L) {
    vp <- validation_proposal %||% proposal
    if (!S7::S7_inherits(vp, is_proposal)) {
      cli::cli_abort("`validation_proposal` must be a {.cls is_proposal} object.")
    }
    if (vp@n_dim != p) {
      cli::cli_abort("`validation_proposal@n_dim` must equal {p}.")
    }
    v_seed <- validation_seed %||% (if (!is.null(seed)) as.integer(seed) + 1L else NULL)
    vdraws <- draw_is_weights(target, vp, validation_size, v_seed)
    log_g_v <- dgmm_log_via_components(vdraws$x, weights, means, covs)
    d_v <- vdraws$log_f - log_g_v
    finite_v <- is.finite(d_v) & is.finite(vdraws$W) & vdraws$W > 0
    val_kld <- sum(vdraws$W[finite_v] * d_v[finite_v])
    list(
      validation_kld = val_kld,
      validation_kld_absolute = if (!kld_is_shifted) val_kld
        else if (is.finite(target@log_normalizer)) val_kld + target@log_normalizer
        else NA_real_,
      validation_ess = vdraws$ess,
      validation_ess_relative = vdraws$ess / validation_size,
      validation_max_weight = vdraws$max_weight,
      validation_support_fraction = vdraws$support_fraction,
      validation_size = validation_size,
      validation_proposal_name = vp@name
    )
  } else {
    list(
      validation_kld = NA_real_,
      validation_kld_absolute = NA_real_,
      validation_ess = NA_real_,
      validation_ess_relative = NA_real_,
      validation_max_weight = NA_real_,
      validation_support_fraction = NA_real_,
      validation_size = 0L,
      validation_proposal_name = NA_character_
    )
  }

  diagnostics <- c(
    list(
      kld_trace = kld_trace,
      kld_final = kld_final,
      kld_final_absolute = kld_final_absolute,
      kld_is_shifted = kld_is_shifted,
      kld_shift_explanation = kld_shift_explanation,
      weighted_obj_trace = weighted_obj_trace,
      mc_se_kld = mc_se_kld,
      ess = ess,
      ess_relative = ess / is_size,
      max_weight = max_weight,
      support_fraction = support_fraction,
      is_size = is_size,
      is_sample = x,
      is_log_weights = log_W,
      proposal_name = proposal@name
    ),
    validation_diag
  )

  fit <- gmm_fit(
    weights = weights,
    means = means,
    covariances = covs,
    target = target,
    regime = "kld",
    diagnostics = diagnostics,
    converged = converged,
    iterations = as.integer(it),
    call = match.call(),
    name = sprintf("kld_em[N=%d] on %s", N, target@name)
  )
  if (isTRUE(canonicalise)) gmm_canonicalise(fit) else fit
}

## Internal: draw an IS sample, evaluate target/proposal log-densities,
## form self-normalised weights, and report the headline diagnostics.
draw_is_weights <- function(target, proposal, n, seed = NULL) {
  draw <- function() proposal@sample(n)
  x <- if (is.null(seed)) draw() else withr::with_seed(seed, draw())
  log_f <- target@log_density(x)
  log_q <- proposal@log_density(x)
  log_w <- log_f - log_q
  finite <- is.finite(log_w)
  support_fraction <- mean(finite)
  if (!any(finite)) {
    cli::cli_abort(c(
      "No importance-sampling draws received finite weight.",
      "i" = "`log_density(x) - log q(x)` is non-finite for every draw; check proposal support and target log-density."
    ))
  }
  log_w_max <- max(log_w[finite])
  log_W <- log_w - log_w_max
  log_W[!is.finite(log_W)] <- -Inf
  log_W <- log_W - log(sum(exp(log_W)))
  W <- exp(log_W)
  W[!is.finite(W)] <- 0
  W <- W / sum(W)
  ess <- 1 / sum(W^2)
  max_weight <- max(W)
  list(
    x = x,
    log_f = log_f,
    log_q = log_q,
    log_W = log_W,
    W = W,
    ess = ess,
    max_weight = max_weight,
    support_fraction = support_fraction
  )
}

## Internal: log g(x) for a mixture given as weights/means/covs, where
## g is the current EM iterate (kept out of the S7 layer to avoid
## constructing a temporary gmm_fit each iteration).
dgmm_log_via_components <- function(x, weights, means, covariances) {
  parts <- gmm_log_unnorm(x, weights, means, covariances)
  logsumexp_rows(parts)
}
