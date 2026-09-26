## Compress a Gaussian-kernel KDE, itself an n-component Gaussian mixture,
## into an N-component proxy by regime (iii) against the KDE as a normalised
## target.

## Per-coordinate bandwidth selection.
##   "silverman"  : (4 / (p + 2))^(1 / (p + 4)) * n^(-1 / (p + 4)) * sd
##   "scott"      : n^(-1 / (p + 4)) * sd
##   numeric(1)   : absolute scalar h applied to every coordinate
##   numeric(p)   : per-coordinate absolute h
choose_bandwidth <- function(samples, bandwidth) {
  n <- nrow(samples)
  p <- ncol(samples)
  if (is.character(bandwidth)) {
    if (length(bandwidth) != 1L) {
      cli::cli_abort("`bandwidth` must be length 1 if a string.")
    }
    bandwidth <- rlang::arg_match0(
      bandwidth, c("silverman", "scott"), arg_nm = "bandwidth"
    )
    sigma <- apply(samples, 2L, stats::sd)
    if (any(sigma <= 0)) {
      cli::cli_abort(c(
        "Cannot derive a rule-of-thumb bandwidth: at least one coordinate has zero variance.",
        "i" = "Supply a numeric `bandwidth` explicitly."
      ))
    }
    factor <- switch(bandwidth,
      silverman = (4 / (p + 2))^(1 / (p + 4)) * n^(-1 / (p + 4)),
      scott     =                              n^(-1 / (p + 4))
    )
    return(factor * sigma)
  }
  if (is.numeric(bandwidth)) {
    if (length(bandwidth) == 1L) {
      if (bandwidth <= 0 || !is.finite(bandwidth)) {
        cli::cli_abort("`bandwidth` must be a positive finite scalar.")
      }
      return(rep_len(as.numeric(bandwidth), p))
    }
    if (length(bandwidth) == p) {
      if (any(bandwidth <= 0) || any(!is.finite(bandwidth))) {
        cli::cli_abort("`bandwidth` entries must be positive and finite.")
      }
      return(as.numeric(bandwidth))
    }
    cli::cli_abort("Numeric `bandwidth` must be length 1 or {p} (got {length(bandwidth)}).")
  }
  cli::cli_abort("`bandwidth` must be a character string or numeric vector.")
}

## Vectorised KDE log-density with a diagonal Gaussian kernel.
## Returns function(x) producing length-nrow(x) numeric. Uses chunking
## to keep peak memory bounded regardless of how many query rows the
## fitter passes in.
make_kde_log_density <- function(samples, h, chunk = 256L) {
  n <- nrow(samples)
  p <- ncol(samples)
  inv_h <- 1 / h
  log_norm_const <- -log(n) - 0.5 * p * log(2 * pi) - sum(log(h))
  function(x) {
    if (is.null(dim(x))) x <- matrix(x, nrow = 1L)
    m <- nrow(x)
    if (m == 0L) return(numeric(0L))
    out <- numeric(m)
    for (start in seq.int(1L, m, by = chunk)) {
      end <- min(start + chunk - 1L, m)
      idx <- seq.int(start, end)
      mi <- length(idx)
      logw <- matrix(0, nrow = mi, ncol = n)
      for (d in seq_len(p)) {
        diff <- outer(x[idx, d], samples[, d], FUN = "-") * inv_h[d]
        logw <- logw - 0.5 * diff * diff
      }
      out[idx] <- log_norm_const + logsumexp_rows(logw)
    }
    out
  }
}

#' Compile a kernel-density estimate into a Gaussian-mixture proxy
#'
#' Fits an `N`-component Gaussian-mixture proxy to a (Gaussian, diagonal-
#' bandwidth) kernel-density estimate over `samples`, via regime (iii)
#' KLD-EM.
#'
#' A Gaussian-kernel KDE on `n` points is already an `n`-component Gaussian
#' mixture with equal weights, so it can be marginalised, conditioned and
#' sampled in closed form as it stands. `from_kde()` is a **compression**:
#' it replaces the `n` components with the closest `N`-component mixture in
#' the Kullback-Leibler sense, with `N` much smaller than `n`. Density
#' evaluation, storage and the operator calculus then cost in proportion to
#' `N` rather than `n` (products and convolutions of two mixtures have as
#' many components as the product of their counts). Bias inherited from the
#' KDE is reproduced in the proxy; the bandwidth controls the bias-variance
#' trade-off.
#'
#' Dimensional scope. The dimensional guard is `p <= 5` (recommended),
#' `p <= 10` (allowed with warning), `p > 10` (rejected). Regime-(iii)
#' KLD-EM is driven by importance sampling, whose effective sample size
#' collapses sharply in high dimensions.
#'
#' @param samples An `n` by `p` numeric matrix of points. `n >= 5`,
#'   `p <= 10`.
#' @param N Number of mixture components in the proxy.
#' @param bandwidth Either `"silverman"`, `"scott"`, a positive numeric
#'   scalar (absolute bandwidth applied to every coordinate), or a
#'   length-`p` positive numeric vector of per-coordinate absolute
#'   bandwidths. Default `"silverman"`.
#' @param proposal Optional [is_proposal]. Default is a multivariate-t
#'   centred at `colMeans(samples)`, scale = ridge(cov(samples)) +
#'   diag(h^2), `df = 5`.
#' @param is_size Importance-sample size for fitting. The default `NULL`
#'   uses `10000 + 2500 * N`: the fitted components overfit a fixed draw, and
#'   larger `N` needs more draws to keep the held-out divergence close to its
#'   large-sample value.
#' @param max_iter Maximum EM iterations. Forwarded to [fit_kld_em()].
#' @param tol Convergence tolerance. Forwarded to [fit_kld_em()].
#' @param ridge_eps Ridge added to each component covariance at every
#'   M-step. Forwarded to [fit_kld_em()].
#' @param min_ess Minimum effective sample size below which a warning is
#'   issued. Forwarded to [fit_kld_em()].
#' @param seed Optional integer seed for the fitting IS draw.
#' @param validation_size Held-out IS sample size. The default `NULL` adds
#'   draws until `validation_kld` is precise, as described in
#'   [fit_kld_em()]; `0L` disables it. Forwarded to [fit_kld_em()].
#' @param validation_proposal Optional [is_proposal] for the held-out
#'   sample. Forwarded to [fit_kld_em()].
#' @param validation_seed Optional integer seed for the held-out IS draw.
#'   Forwarded to [fit_kld_em()].
#' @param support_warn Logical. Forwarded to [fit_kld_em()].
#' @param canonicalise Logical. If `TRUE`, the fitted mixture is
#'   post-processed by [gmm_canonicalise()]. Forwarded to [fit_kld_em()].
#'
#' @returns A [gmm_fit] with `regime = "kld"` and metadata recording the
#'   KDE inputs (`kde_samples_n`, `bandwidth`, `bandwidth_method`). Its
#'   in-sample `kld_final` reads low, often below zero; the held-out
#'   `validation_kld` (see [ess_summary()]) measures the compression error.
#' @family fitting
#' @export
#' @examples
#' set.seed(1L)
#' x <- rbind(
#'   mvnfast::rmvn(120L, mu = c(-2, 0), sigma = diag(2)),
#'   mvnfast::rmvn(120L, mu = c( 2, 0), sigma = diag(2))
#' )
#' fit <- from_kde(x, N = 2L, is_size = 2000L, max_iter = 40L, seed = 1L)
#' fit
#' ess_summary(fit)
from_kde <- function(samples, N = 3L,
                     bandwidth = "silverman",
                     proposal = NULL, is_size = NULL,
                     max_iter = 100L, tol = 1e-5,
                     ridge_eps = 1e-6, min_ess = 50L,
                     seed = NULL,
                     validation_size = NULL,
                     validation_proposal = NULL,
                     validation_seed = NULL,
                     support_warn = TRUE,
                     canonicalise = TRUE) {
  if (!is.matrix(samples) || !is.numeric(samples)) {
    cli::cli_abort("`samples` must be a numeric matrix.")
  }
  if (anyNA(samples)) {
    cli::cli_abort("`samples` must not contain `NA` or `NaN`.")
  }
  n <- nrow(samples)
  p <- ncol(samples)
  if (n < 5L) {
    cli::cli_abort("`samples` must have at least 5 rows (got {n}).")
  }
  if (p < 1L) {
    cli::cli_abort("`samples` must have at least 1 column.")
  }
  if (p > 10L) {
    cli::cli_abort(c(
      "{.fn from_kde} currently supports {.code p <= 10} (got p = {p}).",
      "i" = "Higher ambient dimensions await the v0.3 affine-Gaussian operator calculus.",
      "i" = "See {.file vignettes/from_kde.Rmd} for the design rationale."
    ))
  }
  if (p > 5L) {
    cli::cli_warn(c(
      "{.fn from_kde} is well-tested for {.code p <= 5} (got p = {p}).",
      "i" = "Expect importance-sampling ESS to fall sharply beyond p = 5; inspect {.code ess_summary(fit)}."
    ))
  }
  N <- as.integer(N)
  if (length(N) != 1L || is.na(N) || N < 1L) {
    cli::cli_abort("`N` must be a positive integer scalar.")
  }
  if (is.null(is_size)) {
    is_size <- 10000L + 2500L * N
  }

  h <- choose_bandwidth(samples, bandwidth)
  kde_log_density <- make_kde_log_density(samples, h)

  tgt <- gmm_target(
    n_dim = p,
    log_density = kde_log_density,
    samples = NULL,
    normalised = TRUE,
    log_normalizer = 0,
    name = "from_kde",
    metadata = list(
      source            = "from_kde",
      kde_samples_n     = n,
      bandwidth         = h,
      bandwidth_method  = if (is.character(bandwidth)) bandwidth else "user"
    )
  )

  if (is.null(proposal)) {
    centre  <- colMeans(samples)
    scale_m <- ridge(stats::cov(samples), 1e-3) + diag(h * h, nrow = p)
    proposal <- is_mvt(n_dim = p, mean = centre, sigma = scale_m, df = 5)
  }

  fit <- fit_kld_em(
    target = tgt, N = N, proposal = proposal,
    is_size = is_size, max_iter = max_iter, tol = tol,
    ridge_eps = ridge_eps, min_ess = min_ess, seed = seed,
    validation_size = validation_size,
    validation_proposal = validation_proposal,
    validation_seed = validation_seed,
    support_warn = support_warn,
    canonicalise = canonicalise
  )

  fit@metadata <- modifyList(fit@metadata, list(
    from_kde = list(
      kde_samples_n     = n,
      bandwidth         = h,
      bandwidth_method  = if (is.character(bandwidth)) bandwidth else "user"
    )
  ))
  fit
}
