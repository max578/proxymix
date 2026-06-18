## Diagnostics for fitted proxies.

#' Per-iteration KLD trace of a fit
#'
#' Returns the per-iteration estimate of `KL(f || g_theta)` produced during
#' a regime-(iii) fit, or `NA` for regimes that do not estimate the KLD
#' internally.
#'
#' @param fit A [gmm_fit].
#'
#' @returns Numeric vector (or `NA_real_`).
#' @family diagnostics
#' @export
#' @examples
#' fit <- fit_proxymix(banana_target(), N = 2L, regime = "kld",
#'                     is_size = 1000L, max_iter = 15L, seed = 1L)
#' kld_trace(fit)
kld_trace <- function(fit) {
  if (!S7::S7_inherits(fit, gmm_fit)) {
    cli::cli_abort("`fit` must be a {.cls gmm_fit}.")
  }
  fit@diagnostics$kld_trace %||% NA_real_
}

#' Effective sample size of the importance-sampling weights
#'
#' Returns the effective sample size (`1 / sum(W^2)`) of the
#' self-normalised importance weights used by a regime-(iii) fit. `NA` for
#' regimes that do not use importance sampling.
#'
#' @param fit A [gmm_fit].
#'
#' @returns Numeric scalar (or `NA_real_`).
#' @family diagnostics
#' @export
#' @examples
#' fit <- fit_proxymix(banana_target(), N = 2L, regime = "kld",
#'                     is_size = 1000L, max_iter = 15L, seed = 1L)
#' ess_trace(fit)
ess_trace <- function(fit) {
  if (!S7::S7_inherits(fit, gmm_fit)) {
    cli::cli_abort("`fit` must be a {.cls gmm_fit}.")
  }
  fit@diagnostics$ess %||% NA_real_
}

#' Monte-Carlo Hellinger distance between a fit and its target
#'
#' Estimates the squared Hellinger distance `H^2(f, g) = 1 - integral
#' sqrt(f(x) g(x)) dx` by importance sampling against the proposal stored
#' in the fit (for regime `"kld"`) or by sampling from the fit itself (for
#' regime `"sample"`). The target's `log_density` must be supplied **and
#' normalised**; otherwise the Monte Carlo integral is biased by the
#' missing \eqn{\sqrt{Z(f)}}. When the target's `normalised` property is
#' not `TRUE`, a warning is issued and the returned value is flagged.
#'
#' @param fit A [gmm_fit] whose target carries a `log_density`.
#' @param n_mc Number of Monte Carlo samples.
#' @param seed Optional integer seed.
#'
#' @returns A list with components
#'   * `h2` - estimate of `H^2(f, g)`,
#'   * `se` - Monte Carlo standard error,
#'   * `n_mc` - sample size used.
#' @family diagnostics
#' @export
#' @examples
#' fit <- fit_proxymix(banana_target(), N = 3L, regime = "kld",
#'                     is_size = 2000L, max_iter = 25L, seed = 1L)
#' hellinger_mc(fit, n_mc = 1000L, seed = 1L)
hellinger_mc <- function(fit, n_mc = 5000L, seed = NULL) {
  if (!S7::S7_inherits(fit, gmm_fit)) {
    cli::cli_abort("`fit` must be a {.cls gmm_fit}.")
  }
  tgt <- fit@target
  if (is.null(tgt) || is.null(tgt@log_density)) {
    cli::cli_abort("`fit@target@log_density` must be supplied for a Hellinger estimate.")
  }
  trustworthy <- isTRUE(tgt@normalised)
  if (!trustworthy) {
    cli::cli_warn(c(
      "Target is not declared {.code normalised = TRUE}; Hellinger estimate is biased.",
      "i" = "Set {.code normalised = TRUE} on the target (or supply {.code log_normalizer}) for a meaningful H^2."
    ))
  }

  draw <- function() {
    if (fit@regime == "kld" && !is.null(fit@diagnostics$is_sample)) {
      ## Reuse the IS sample from the fit.
      x <- fit@diagnostics$is_sample
      log_w <- fit@diagnostics$is_log_weights
      list(x = x, log_w = log_w, source = "kld_is")
    } else {
      ## Sample from the fit and importance-sample target.
      n <- as.integer(n_mc)
      x <- rgmm(n, fit)
      list(x = x, log_w = rep(0, n), source = "gmm_proposal")
    }
  }
  d <- if (is.null(seed)) draw() else withr::with_seed(seed, draw())
  x <- d$x
  log_f <- tgt@log_density(x)
  log_g <- dgmm(x, fit, log = TRUE)

  if (d$source == "kld_is") {
    ## Use IS weights normalised to sum to 1.
    log_W <- d$log_w
    log_W <- log_W - max(log_W[is.finite(log_W)])
    log_W <- log_W - log(sum(exp(log_W)))
    W <- exp(log_W)
    ## h^2 = 1 - E_q [ sqrt(f g) / q ] but we are sampling from q via x and
    ## weighting with f/q already absorbed: instead express
    ##   integral sqrt(f g) dx = E_W[sqrt(g/f)] when sampling from p~f.
    ## We use the self-normalised IS form: integral sqrt(f g) dx
    ##   approximately sum_n (W_n / f(x_n)) * sqrt(f(x_n) * g(x_n))
    ##   = sum_n W_n * sqrt(g(x_n) / f(x_n))
    ratio <- exp(0.5 * (log_g - log_f))
    finite <- is.finite(ratio) & is.finite(W)
    integral_est <- sum(W[finite] * ratio[finite])
    se <- stats::sd(ratio[finite]) / sqrt(sum(finite))
    n_used <- sum(finite)
  } else {
    ## Sampling from g_theta directly. integral sqrt(f g) dx
    ##   = E_g [ sqrt(f / g) ] (sample x ~ g).
    ratio <- exp(0.5 * (log_f - log_g))
    finite <- is.finite(ratio)
    integral_est <- mean(ratio[finite])
    se <- stats::sd(ratio[finite]) / sqrt(sum(finite))
    n_used <- sum(finite)
  }
  h2 <- 1 - integral_est
  list(h2 = h2, se = se, n_mc = n_used, trustworthy = trustworthy)
}

#' Summary of importance-sampling diagnostics
#'
#' Convenience accessor returning the headline IS-quality numbers for a
#' regime-(iii) fit: effective sample size and its ratio to `is_size`,
#' the largest self-normalised weight, the support fraction (proportion
#' of draws that received a finite weight), and the Monte-Carlo standard
#' error of the final KLD estimate. Returns `NA` fields for regimes that
#' do not use importance sampling.
#'
#' Validation-side numbers (`validation_*`) are populated only when the
#' fit was called with `validation_size > 0`.
#'
#' @param fit A [gmm_fit].
#'
#' @returns A list of numeric scalars (or `NA`s where not applicable).
#' @family diagnostics
#' @export
#' @examples
#' fit <- fit_proxymix(banana_target(), N = 3L, regime = "kld",
#'                     is_size = 1500L, max_iter = 20L, seed = 1L,
#'                     validation_size = 1500L)
#' ess_summary(fit)
ess_summary <- function(fit) {
  if (!S7::S7_inherits(fit, gmm_fit)) {
    cli::cli_abort("`fit` must be a {.cls gmm_fit}.")
  }
  d <- fit@diagnostics
  list(
    is_size                   = d$is_size %||% NA_integer_,
    ess                       = d$ess %||% NA_real_,
    ess_relative              = d$ess_relative %||% NA_real_,
    max_weight                = d$max_weight %||% NA_real_,
    support_fraction          = d$support_fraction %||% NA_real_,
    mc_se_kld                 = d$mc_se_kld %||% NA_real_,
    validation_size           = d$validation_size %||% NA_integer_,
    validation_ess            = d$validation_ess %||% NA_real_,
    validation_ess_relative   = d$validation_ess_relative %||% NA_real_,
    validation_max_weight     = d$validation_max_weight %||% NA_real_,
    validation_support_fraction = d$validation_support_fraction %||% NA_real_,
    validation_kld            = d$validation_kld %||% NA_real_
  )
}

#' Information criteria: BIC, AIC, and ICL
#'
#' Returns the Bayesian and Akaike information criteria of a regime-(ii) fit,
#' together with the integrated completed likelihood (ICL). All three are
#' computed against the *empirical* log-likelihood of the samples used to fit
#' the model and are reported on the same scale (smaller is better). They are
#' `NA` for regimes that do not have an empirical likelihood (`"moment"`,
#' `"kld"`).
#'
#' The ICL of Biernacki, Celeux and Govaert (2000) adds to the BIC twice the
#' entropy of the fitted classification,
#' \eqn{\mathrm{ICL} = \mathrm{BIC} + 2 E_N}, where
#' \eqn{E_N = -\sum_{i,k} \gamma_{ik} \log \gamma_{ik} \ge 0} is the entropy of
#' the responsibilities \eqn{\gamma_{ik}}. It therefore penalises mixtures whose
#' components overlap (uncertain assignments), and favours well-separated
#' clustering solutions over the merely best-fitting ones. Because
#' \eqn{E_N \ge 0}, the ICL is never smaller than the BIC, and the two coincide
#' for a single component (\eqn{K = 1}), where every responsibility is one. The
#' classification entropy itself is returned as `classification_entropy`.
#'
#' @param fit A [gmm_fit].
#'
#' @returns A list with `bic`, `aic`, `icl`, `classification_entropy`, and
#'   `n_params`.
#' @family diagnostics
#' @references Biernacki, C., Celeux, G. and Govaert, G. (2000) Assessing a
#'   mixture model for clustering with the integrated completed likelihood.
#'   *IEEE Transactions on Pattern Analysis and Machine Intelligence* 22(7),
#'   719--725. \doi{10.1109/34.865189}
#' @export
#' @examples
#' x <- matrix(stats::rnorm(200), ncol = 2)
#' tgt <- gmm_target_from_samples(x)
#' fit <- fit_proxymix(tgt, N = 2L, regime = "sample", max_iter = 25L)
#' bic_aic(fit)
bic_aic <- function(fit) {
  if (!S7::S7_inherits(fit, gmm_fit)) {
    cli::cli_abort("`fit` must be a {.cls gmm_fit}.")
  }
  bic <- fit@diagnostics$bic %||% NA_real_
  aic <- fit@diagnostics$aic %||% NA_real_
  n_params <- fit@diagnostics$n_params %||% NA_integer_

  ## ICL = BIC + 2 E_N, with the classification entropy
  ##   E_N = -sum_{i,k} gamma_ik log gamma_ik   (0 log 0 = 0).
  ## Defined for the sample regime, where the responsibilities are over the
  ## fitted samples; NA elsewhere, matching the BIC / AIC policy.
  en <- NA_real_
  icl <- NA_real_
  if (identical(fit@regime, "sample") && !is.null(fit@target) &&
        !is.null(fit@target@samples) && is.finite(bic)) {
    samples <- fit@target@samples
    log_unnorm <- gmm_log_unnorm(samples, fit@weights, fit@means,
                                 fit@covariances)
    log_resp <- log_unnorm - logsumexp_rows(log_unnorm)
    resp <- exp(log_resp)
    en <- -sum(resp * log_resp, na.rm = TRUE)
    icl <- bic + 2 * en
  }

  list(
    bic = bic,
    aic = aic,
    icl = icl,
    classification_entropy = en,
    n_params = n_params
  )
}

`%||%` <- function(a, b) if (is.null(a)) b else a
