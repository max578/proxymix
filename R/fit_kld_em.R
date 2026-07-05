## Regime (iii) of Hoek and Elliott (2024).
##
## Importance-sampled KLD-EM minimises KL(f || g_theta) when the target f
## can be evaluated point-wise but not (cheaply) sampled.

#' Importance-sampled KLD-EM fit (regime iii)
#'
#' Implements regime (iii) of Hoek and Elliott (2024). Minimises
#' `KL(f || g_theta)` where `f` is supplied as an evaluable log-density on
#' the target, via expectation-maximisation against importance-sampled
#' draws from a user-chosen proposal `q`.
#'
#' With `adapt = "none"` (the default) the Monte Carlo draws from `q` are
#' computed once at the start and the resulting self-normalised
#' importance-sampling weights are reused at every EM iteration. With
#' `adapt = "pmc"` the proposal is refreshed every `refresh_every`
#' iterations with a defensive mixture built from the current iterate --
#' the population-Monte-Carlo scheme: the fitted mixture (covariances
#' inflated by `inflate`) carries `1 - defensive_gamma` of the proposal
#' mass and the original proposal `q` keeps `defensive_gamma` as a
#' heavy-tailed anchor, a fresh IS batch is drawn, and EM continues on the
#' refreshed weights. Because the refreshed proposal tracks the target,
#' the effective sample size recovers from a poor initial proposal and the
#' usable dimension range extends well beyond what a fixed proposal
#' reaches; the per-batch ESS trace is reported as
#' `diagnostics$ess_history`. While a batch is degenerate (its effective
#' sample size is below `min_ess`), the refresh fires every iteration
#' with an escalating covariance inflation floored at a growing fraction
#' of the batch's sample covariance, so a collapsed iterate walks back
#' out toward the target instead of freezing; and convergence is only
#' accepted on an adapted batch, so a run that stabilises on the original
#' proposal's draw is refreshed at least once before it is allowed to
#' stop. The scheme is the mixture population-Monte-Carlo idea of Cappé
#' et al. (2008) with the defensive-mixture safeguard of Owen and Zhou
#' (2000); it re-draws rather than recycles batches (compare the adaptive
#' multiple importance sampling of Cornuet et al., 2012).
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
#'   importance-weighted EM objective `Q(theta) = sum_n W_n log g(x_n)`.
#'   `Q` is invariant to the target's normalising constant, so the
#'   stopping rule behaves identically for normalised and unnormalised
#'   targets (the importance-sampled KLD estimate carries an additive
#'   `-log Z(f)` offset and is therefore never used for stopping).
#' @param ridge_eps Ridge added to each component covariance at every
#'   M-step.
#' @param min_ess Minimum effective sample size below which the fit is
#'   flagged as degenerate: a classed warning (`proxymix_low_ess`) is
#'   issued (or, with `on_low_ess = "abort"`, a classed error
#'   `proxymix_degenerate_fit`), the fit's `converged` flag is forced to
#'   `FALSE`, and `degenerate = TRUE` is recorded in the diagnostics and
#'   the quality certificate.
#' @param on_low_ess What to do when the effective sample size falls below
#'   `min_ess`: `"warn"` (the default) flags and continues, `"abort"`
#'   refuses to return a degenerate fit.
#' @param seed Optional integer seed. When supplied, the fit is
#'   reproducible end-to-end: the fitting IS draw, the initialisation
#'   resample and kmeans pass, and any empty-component reseed draws are
#'   all derived from it. When `NULL`, those draws consume the ambient
#'   random-number stream.
#' @param validation_size Number of independent importance-sampling draws
#'   to use for held-out validation. The default `NULL` uses
#'   `ceiling(is_size / 4)`, so the overfit-vs-generalise diagnostic
#'   (`validation_kld` and the certificate's `validation_gap`) exists by
#'   default; set `0L` to disable the validation split.
#' @param validation_proposal Optional [is_proposal] for the validation
#'   sample. Defaults to the same proposal used for fitting.
#' @param validation_seed Optional integer seed used when drawing the
#'   validation sample. Defaults to `seed + 1L` when `seed` is supplied,
#'   `NULL` otherwise.
#' @param support_warn Logical. If `TRUE` (the default), issue a warning
#'   when more than 5% of IS draws receive non-finite weights (typically
#'   because the proposal does not dominate the target's support).
#' @param adapt Proposal adaptation: `"none"` (the default; one fixed IS
#'   draw, the historical behaviour) or `"pmc"` (population-Monte-Carlo
#'   refresh of the proposal from the current iterate; see Details).
#' @param refresh_every With `adapt = "pmc"`, refresh the proposal after
#'   this many EM iterations on the current batch. Default `5L`.
#' @param defensive_gamma With `adapt = "pmc"`, the mass kept on the
#'   original proposal as a heavy-tailed defensive anchor at every
#'   refresh (bounds the importance-weight variance). Default `0.15`.
#' @param inflate With `adapt = "pmc"`, the factor inflating the current
#'   iterate's covariances inside the refreshed proposal. Default `1.5`.
#' @param anneal Logical. If `TRUE`, a deterministic-annealing warm-start
#'   (see [gmm_anneal_path()]) replaces the kmeans initialisation: components
#'   are annealed from a high temperature down to one on the importance-weighted
#'   draws, and the resulting parameters seed the (unchanged) cold KLD-EM loop.
#'   This attacks the local-optima sensitivity of cold EM. Defaults to `FALSE`.
#' @param temp_schedule Optional numeric vector of descending temperatures for
#'   the annealing warm-start. `NULL` (the default) uses a geometric schedule
#'   from `10` down to `1` in covariance-whitened units. Ignored when
#'   `anneal = FALSE`.
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
#' @references Cappé, O., Douc, R., Guillin, A., Marin, J.-M. and
#'   Robert, C. P. (2008) Adaptive importance sampling in general mixture
#'   classes. *Statistics and Computing* 18, 447--459.
#'   \doi{10.1007/s11222-008-9059-x}
#'
#'   Cornuet, J.-M., Marin, J.-M., Mira, A. and Robert, C. P. (2012)
#'   Adaptive multiple importance sampling. *Scandinavian Journal of
#'   Statistics* 39, 798--812. \doi{10.1111/j.1467-9469.2011.00756.x}
#'
#'   Owen, A. and Zhou, Y. (2000) Safe and effective importance sampling.
#'   *Journal of the American Statistical Association* 95(449), 135--143.
#'   \doi{10.1080/01621459.2000.10473909}
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
                       on_low_ess = c("warn", "abort"),
                       seed = NULL,
                       validation_size = NULL,
                       validation_proposal = NULL,
                       validation_seed = NULL,
                       support_warn = TRUE,
                       adapt = c("none", "pmc"),
                       refresh_every = 5L,
                       defensive_gamma = 0.15,
                       inflate = 1.5,
                       anneal = FALSE,
                       temp_schedule = NULL,
                       canonicalise = TRUE) {
  if (!S7::S7_inherits(target, gmm_target)) {
    cli::cli_abort("`target` must be a {.cls gmm_target} object.")
  }
  if (is.null(target@log_density)) {
    cli::cli_abort("regime {.val kld} requires `target@log_density` to be supplied.")
  }
  on_low_ess <- rlang::arg_match(on_low_ess)
  adapt <- rlang::arg_match(adapt)
  refresh_every <- as.integer(refresh_every)
  if (length(refresh_every) != 1L || is.na(refresh_every) || refresh_every < 1L) {
    cli::cli_abort("`refresh_every` must be a single positive integer.")
  }
  if (!is.numeric(defensive_gamma) || length(defensive_gamma) != 1L ||
        defensive_gamma <= 0 || defensive_gamma >= 1) {
    cli::cli_abort("`defensive_gamma` must be a single number strictly inside (0, 1).")
  }
  if (!is.numeric(inflate) || length(inflate) != 1L || inflate < 1) {
    cli::cli_abort("`inflate` must be a single number of at least 1.")
  }
  N <- as.integer(N)
  is_size <- as.integer(is_size)
  validation_size <- if (is.null(validation_size)) {
    as.integer(ceiling(is_size / 4))
  } else {
    as.integer(validation_size)
  }
  p <- target@n_dim

  ## The dimension guard lives at the core fitter, not only in the wrapper
  ## entry points: importance sampling loses effective sample size sharply
  ## with dimension, and a direct call deserves the same disclosure.
  if (p > 10L) {
    cli::cli_warn(c(
      "Importance-sampled KLD-EM in {p} dimensions: expect severe effective-sample-size loss.",
      "i" = "Inspect {.code ess_summary(fit)}; consider fitting a lower-dimensional sub-vector and composing with the operator calculus."
    ), class = "proxymix_high_dimension")
  } else if (p > 5L) {
    cli::cli_inform(c(
      "Importance-sampled KLD-EM is well-characterised for {.code p <= 5} (got p = {p}).",
      "i" = "Inspect {.code ess_summary(fit)} for effective-sample-size loss."
    ))
  }

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

  flag_low_ess <- function(ess_now) {
    low_ess_msg <- c(
      "Effective sample size is low: ESS = {round(ess_now, 1)} out of {is_size}; the fit is flagged as degenerate.",
      "i" = "Consider a heavier-tailed proposal, more IS draws, or {.code adapt = \"pmc\"}."
    )
    if (on_low_ess == "abort") {
      cli::cli_abort(low_ess_msg, class = "proxymix_degenerate_fit")
    }
    cli::cli_warn(low_ess_msg, class = "proxymix_low_ess")
  }
  degenerate <- ess < min_ess
  ## Under PMC adaptation a poor initial draw is expected to recover, so
  ## degeneracy is judged on the FINAL batch after the loop instead.
  if (degenerate && adapt == "none") {
    flag_low_ess(ess)
  }
  if (isTRUE(support_warn) && support_fraction < 0.95) {
    cli::cli_warn(c(
      "Importance proposal does not dominate target support: only {round(100 * support_fraction, 1)}% of IS draws received finite weight.",
      "i" = "Verify {.code log_density(x) - log q(x)} is finite on the target's support; consider a wider or heavier-tailed proposal."
    ), class = "proxymix_support")
  }

  ## IS-weighted data moments, used to reseed an empty component at data
  ## scale rather than at the arbitrary unit scale.
  mu_w <- colSums(W * x)
  S_w <- crossprod((x - matrix(mu_w, nrow = is_size, ncol = p, byrow = TRUE)) *
                     sqrt(W))
  ## Data-scaled ridge: invariant to the target's units, constant within
  ## the fit (see `.data_scaled_eps()`).
  ridge_eps <- .data_scaled_eps(ridge_eps, mean(diag(S_w)))

  ## ---- Deterministic-annealing warm-start (optional) ----
  anneal_schedule_used <- NULL
  if (isTRUE(anneal)) {
    init <- .anneal_em_warmstart(
      x, rw = W, N = N, temp_schedule = temp_schedule,
      ridge_eps = ridge_eps,
      seed = if (!is.null(seed)) as.integer(seed) + 2L else NULL
    )
    anneal_schedule_used <- init@metadata$temp_schedule
  }

  ## ---- Initialisation ----
  if (is.null(init)) {
    ## Bootstrap-resample x by IS weights to obtain a pseudo-sample with
    ## approximately target distribution, then kmeans on it. Derived from
    ## `seed` when one is supplied, so the whole fit is reproducible.
    build_init <- function() {
      idx <- sample.int(is_size, size = min(is_size, 5L * N * 50L),
                        replace = TRUE, prob = W)
      pseudo <- x[idx, , drop = FALSE]
      tryCatch(init_kmeans(pseudo, N = N, ridge_eps = ridge_eps),
               error = function(e) NULL)
    }
    init <- if (is.null(seed)) {
      build_init()
    } else {
      withr::with_seed(as.integer(seed) + 3L, build_init())
    }
    if (is.null(init)) {
      init <- init_random(N = N, p = p,
                          centre = if (!is.null(target@samples)) colMeans(target@samples) else rep(0, p),
                          scale = 1, sigma_diag = 1,
                          seed = if (is.null(seed)) 42L else as.integer(seed) + 4L)
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
  n_reseeds <- 0L
  n_refresh <- 0L
  ess_history <- ess
  batch_start <- 0L
  esc <- 1
  force_refresh <- FALSE
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

    ## Convergence is judged on the importance-weighted EM objective
    ## Q(theta) = sum_n W_n log g(x_n), never on the KLD trace: for an
    ## unnormalised target the KLD trace carries the additive -log Z(f)
    ## offset, so a relative-change rule on it would depend on the
    ## target's arbitrary normalisation. Q is offset-free by construction
    ## (the self-normalised W and log g never touch the constant). Both
    ## compared iterations must belong to the same IS batch: across a PMC
    ## refresh, Q jumps because the sample changed, not the parameters.
    ## Under PMC, convergence cannot fire while the current batch is
    ## degenerate -- stabilising on a handful of effective draws is not
    ## convergence, and the refreshes need iterations to walk the proposal
    ## toward the target.
    if (it > 1L && (it - 1L) > batch_start &&
          !(adapt == "pmc" && ess < min_ess)) {
      delta <- abs(weighted_obj_trace[it] - weighted_obj_trace[it - 1L]) /
        (abs(weighted_obj_trace[it - 1L]) + 1e-12)
      if (delta < tol) {
        if (adapt == "pmc" && batch_start == 0L) {
          ## The parameters stabilised on the ORIGINAL proposal's batch.
          ## Adapt before accepting: the refreshed proposal (built from
          ## the stabilised iterate) is usually far more efficient, and
          ## the final diagnostics should come from an adapted batch.
          force_refresh <- TRUE
        } else {
          converged <- TRUE
          break
        }
      }
    }

    ## M-step.
    W_resp <- resp * W
    Nk_W <- colSums(W_resp)
    empty <- Nk_W < 1e-12
    ## A dead component must get its WEIGHT reset too, not just its mean:
    ## a zero (or 1e-300-floored) weight zeroes its responsibilities on the
    ## next E-step, so a reseeded mean is dead on arrival and the reseed
    ## re-fires every iteration. Grant it one draw's worth of mass.
    Nk_adj <- Nk_W
    Nk_adj[empty] <- 1 / is_size
    n_reseeds <- n_reseeds + sum(empty)
    weights <- as.numeric(Nk_adj / sum(Nk_adj))
    for (k in seq_len(N)) {
      if (empty[k]) {
        ## Re-seed at an IS-weighted random draw, at data scale.
        idx <- if (is.null(seed)) {
          sample.int(is_size, size = 1L, prob = W)
        } else {
          withr::with_seed(as.integer(seed) + 100L + it,
                           sample.int(is_size, size = 1L, prob = W))
        }
        means[[k]] <- as.numeric(x[idx, ])
        covs[[k]] <- ridge(S_w, max(ridge_eps, 1e-6))
        next
      }
      mu_new <- as.numeric(colSums(W_resp[, k] * x) / Nk_W[k])
      diff <- x - matrix(mu_new, nrow = is_size, ncol = p, byrow = TRUE)
      S_new <- crossprod(diff * sqrt(W_resp[, k])) / Nk_W[k]
      means[[k]] <- mu_new
      covs[[k]] <- symmetrise(ridge(S_new, ridge_eps))
    }

    ## PMC refresh: rebuild the proposal from the current iterate (a
    ## defensive mixture keeping `defensive_gamma` mass on the original
    ## proposal as a heavy-tailed anchor), draw a fresh IS batch, and
    ## continue EM on the refreshed weights. While the current batch is
    ## degenerate the refresh fires EVERY iteration with an escalating
    ## covariance inflation, so a proposal far from the target walks
    ## toward it (each refresh re-centres on the current iterate) instead
    ## of freezing after a single hop; the escalation resets as soon as a
    ## healthy batch is drawn.
    if (adapt == "pmc" && it < max_iter &&
          (force_refresh || it %% refresh_every == 0L || ess < min_ess)) {
      force_refresh <- FALSE
      ## On a degenerate batch the iterate has typically collapsed onto a
      ## handful of draws, so multiplying its (near-zero) covariance can
      ## never re-acquire the target. The escalation therefore also floors
      ## the proposal covariance at a growing fraction of the current
      ## batch's sample covariance -- the batch is dominated by the wide
      ## anchor, so the floor interpolates from "trust the iterate" up to
      ## "cover what the anchor covers", while the component MEANS stay at
      ## the iterate (the effective draws' best guess of the target).
      prop_covs <- if (ess < min_ess) {
        S_floor <- diag(diag(stats::cov(x)), nrow = p) * (esc / 64)
        lapply(covs, function(S) (inflate * esc) * S + S_floor)
      } else {
        lapply(covs, function(S) inflate * S)
      }
      g_infl <- gmm(
        weights = weights, means = means,
        covariances = prop_covs
      )
      draws <- .draw_pmc_refresh(
        target, g_infl, proposal, defensive_gamma, is_size,
        seed = if (is.null(seed)) NULL else as.integer(seed) + 500L + it
      )
      x <- draws$x
      log_f <- draws$log_f
      log_W <- draws$log_W
      W <- draws$W
      ess <- draws$ess
      max_weight <- draws$max_weight
      support_fraction <- draws$support_fraction
      ess_history <- c(ess_history, ess)
      n_refresh <- n_refresh + 1L
      batch_start <- it
      esc <- if (ess < min_ess) min(esc * 2, 64) else 1
    }
  }

  ## Under PMC, degeneracy is judged on the final batch (a poor initial
  ## proposal is expected to recover through the refreshes).
  if (adapt == "pmc") {
    degenerate <- ess < min_ess
    if (degenerate) {
      flag_low_ess(ess)
    }
  }

  ## The IS-weighted EM objective is monotone under exact updates. An
  ## empty-component reseed is a deliberate restart-like intervention that
  ## legitimately dips the objective, so the invariant is only asserted for
  ## reseed-free runs.
  if (n_reseeds == 0L && length(weighted_obj_trace) > 2L) {
    worst_drop <- min(diff(weighted_obj_trace))
    if (worst_drop < -1e-6 * max(1, abs(weighted_obj_trace[1L]))) {
      cli::cli_warn(
        c("The importance-weighted EM objective decreased by {signif(-worst_drop, 3)} during fitting.",
          "i" = "This usually signals ridge regularisation interacting with a near-singular component."),
        class = "proxymix_nonmonotone"
      )
    }
  }

  ## ---- Final fitting-sample diagnostics ----
  log_resp_unnorm <- gmm_log_unnorm(x, weights, means, covs)
  log_g <- logsumexp_rows(log_resp_unnorm)
  d_n <- log_f - log_g
  finite_d <- is.finite(d_n) & is.finite(W) & W > 0
  kld_final <- sum(W[finite_d] * d_n[finite_d])
  mc_se_kld <- sqrt(sum(W[finite_d]^2 * (d_n[finite_d] - kld_final)^2))

  ## Per-component Kish ESS at the final responsibilities: a tail component
  ## can sit on a handful of effective draws while the global ESS looks
  ## healthy.
  W_resp_final <- exp(log_resp_unnorm - log_g) * W
  per_component_ess <- colSums(W_resp_final)^2 /
    pmax(colSums(W_resp_final^2), 1e-300)

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
    val_mc_se <- sqrt(sum(vdraws$W[finite_v]^2 *
                            (d_v[finite_v] - val_kld)^2))
    list(
      validation_kld = val_kld,
      validation_mc_se = val_mc_se,
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
      validation_mc_se = NA_real_,
      validation_kld_absolute = NA_real_,
      validation_ess = NA_real_,
      validation_ess_relative = NA_real_,
      validation_max_weight = NA_real_,
      validation_support_fraction = NA_real_,
      validation_size = 0L,
      validation_proposal_name = NA_character_
    )
  }

  ## A degenerate (ESS-collapsed) fit is never reported as converged: the
  ## EM may have stabilised, but on a handful of effective draws.
  converged <- converged && !degenerate

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
      per_component_ess = per_component_ess,
      max_weight = max_weight,
      support_fraction = support_fraction,
      degenerate = degenerate,
      n_reseeds = n_reseeds,
      adapt = adapt,
      n_refresh = n_refresh,
      ess_history = ess_history,
      n_target_evals = is_size * (1L + n_refresh) + validation_size,
      is_size = is_size,
      is_sample = x,
      is_log_weights = log_W,
      proposal_name = proposal@name,
      annealed = isTRUE(anneal),
      temp_schedule = anneal_schedule_used
    ),
    validation_diag
  )

  quality <- list(
    regime = "kld",
    converged = converged,
    degenerate = degenerate,
    ess = ess,
    ess_relative = ess / is_size,
    min_component_ess = min(per_component_ess),
    max_weight = max_weight,
    support_fraction = support_fraction,
    kld_final = kld_final,
    validation_gap = if (validation_size > 0L) {
      validation_diag$validation_kld - kld_final
    } else {
      NA_real_
    }
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
    name = sprintf("kld_em[N=%d] on %s", N, target@name),
    metadata = list(quality = quality)
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
  c(list(x = x, log_f = log_f, log_q = log_q),
    .finalise_is_weights(log_f, log_q))
}

## Internal: self-normalise log importance weights and report the headline
## diagnostics. Shared by the initial draw and the PMC refreshes.
.finalise_is_weights <- function(log_f, log_q) {
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
  list(
    log_W = log_W,
    W = W,
    ess = 1 / sum(W^2),
    max_weight = max(W),
    support_fraction = support_fraction
  )
}

## Internal: a PMC refresh draw from the defensive proposal
##   q_t = (1 - gamma) * g_infl + gamma * q0,
## where g_infl is the current EM iterate with inflated covariances and q0
## is the original (heavy-tailed anchor) proposal. The mixture density of
## the two parts is evaluated exactly, so the weights are unbiased for the
## defensive proposal actually sampled from.
.draw_pmc_refresh <- function(target, g_infl, q0, gamma, n, seed = NULL) {
  draw <- function() {
    from_anchor <- stats::runif(n) < gamma
    n_anchor <- sum(from_anchor)
    x <- matrix(0, nrow = n, ncol = gmm_dim(g_infl))
    if (n_anchor > 0L) {
      x[from_anchor, ] <- q0@sample(n_anchor)
    }
    if (n_anchor < n) {
      x[!from_anchor, ] <- rgmm(n - n_anchor, g_infl)
    }
    x
  }
  x <- if (is.null(seed)) draw() else withr::with_seed(seed, draw())
  log_f <- target@log_density(x)
  log_q <- logsumexp_rows(cbind(
    log1p(-gamma) + dgmm(x, g_infl, log = TRUE),
    log(gamma) + q0@log_density(x)
  ))
  c(list(x = x, log_f = log_f, log_q = log_q),
    .finalise_is_weights(log_f, log_q))
}

## Internal: log g(x) for a mixture given as weights/means/covs, where
## g is the current EM iterate (kept out of the S7 layer to avoid
## constructing a temporary gmm_fit each iteration).
dgmm_log_via_components <- function(x, weights, means, covariances) {
  parts <- gmm_log_unnorm(x, weights, means, covariances)
  logsumexp_rows(parts)
}
