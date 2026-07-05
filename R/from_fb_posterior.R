## from_fb_posterior.R -- the posterior-producer seam.
##
## Compress an external Bayesian posterior into a closed-form Gaussian-mixture
## proxy. The producer side (a package exporting an `fb_log_posterior`
## generic) is addressed only through a small, documented producer interface
## (`fb_log_posterior_spec()`); the path degrades gracefully behind a
## capability probe (`fb_producer_available()`), and a synthetic mock producer
## (`mock_fb_posterior()`) makes the path testable with no producer package
## installed. proxymix never `Imports:` a producer package -- the seam is a
## soft contract and `R CMD check` is clean with none installed.

# ---------------------------------------------------------------------------
# Producer interface -- the seam an external posterior producer must satisfy
# ---------------------------------------------------------------------------

#' The `fb_log_posterior` producer interface
#'
#' Materialises the posterior-producer contract as a small, validated
#' specification object that [from_fb_posterior()] consumes. It is the single
#' seam between proxymix (the consumer) and an external Bayesian fitting
#' package (the producer): any object that satisfies the interface described
#' here can be compressed into a Gaussian-mixture proxy without proxymix ever
#' depending on the producer package at runtime.
#'
#' A conforming producer is one of:
#'
#' * a **bare callable** `function(theta)` that takes a numeric matrix with
#'   rows indexing independent parameter draws and columns indexing parameters,
#'   and returns a length-`nrow(theta)` numeric vector of
#'   `log p(theta | data) + const` (the unnormalised log-posterior). The
#'   callable must be vectorised, side-effect free, and domain-safe (return
#'   `-Inf` outside support rather than erroring). Supply `parameter_names`
#'   (or attach them as `attr(callable, "parameter_names")`);
#' * a **fitted model object** whose class registers an
#'   `fb_log_posterior()` method *in its own package*. When such a producer
#'   is installed, `fb_log_posterior(fit)` is expected to return an object
#'   carrying the callable plus the metadata fields below, which this
#'   constructor normalises;
#' * an already-built **`fb_log_posterior_spec`**, returned unchanged.
#'
#' The metadata fields the producer should supply, in addition to the callable,
#' are:
#'
#' \describe{
#'   \item{`parameter_names`}{Character vector naming the parameters; its
#'     length fixes the proxy's ambient dimension.}
#'   \item{`log_normalizer`}{The log marginal likelihood `log Z`, when the
#'     producer can supply it, or `NA_real_` otherwise. A finite value lets
#'     downstream diagnostics report an absolute, rather than shifted, KLD.}
#'   \item{`support_lower`, `support_upper`}{Optional length-`n_dim` numeric
#'     bounds on each parameter's support (e.g. a variance parameter is bounded
#'     below by zero). Used to centre and scale the default importance proposal;
#'     `NA` entries are treated as unbounded.}
#'   \item{`draws`}{Optional `n` by `n_dim` numeric matrix of posterior draws
#'     from the producer, used only to seed the default importance proposal's
#'     location and scale. Never required -- the contract is the
#'     *log-density*, not the draws.}
#' }
#'
#' Where the contract leaves a detail open, the choice made here (vectorised
#' matrix input, `-Inf`-outside-support, optional `support_*` / `draws`
#' proposal seeds) is recorded so both sides of the seam agree.
#'
#' @param producer A bare callable, a fitted model object from a producer
#'   package, or an `fb_log_posterior_spec`.
#' @param parameter_names Character vector of parameter names. Required for the
#'   bare-callable form unless attached as
#'   `attr(producer, "parameter_names")`.
#' @param log_normalizer Numeric scalar `log Z`, or `NA_real_` (the default)
#'   when unknown.
#' @param support_lower,support_upper Optional length-`n_dim` numeric support
#'   bounds (`NA` for an unbounded coordinate). Default `NULL` (all unbounded).
#' @param draws Optional `n` by `n_dim` numeric matrix of posterior draws used
#'   only to seed the default proposal. Default `NULL`.
#' @param name Optional human-readable name. Defaults to `"fb_posterior"`.
#'
#' @returns An `fb_log_posterior_spec`: a classed list with elements
#'   `log_density` (the validated, vectorised callable), `n_dim`,
#'   `parameter_names`, `log_normalizer`, `support_lower`, `support_upper`,
#'   `draws`, and `name`.
#' @family interop
#' @seealso [from_fb_posterior()] for the consumer verb, [mock_fb_posterior()]
#'   for a synthetic producer, [fb_producer_available()] for the capability
#'   probe.
#' @export
#' @examples
#' ## Build a spec from a bare unnormalised log-posterior (a 2-D Gaussian).
#' log_post <- function(theta) -0.5 * rowSums(theta^2)
#' spec <- fb_log_posterior_spec(
#'   log_post,
#'   parameter_names = c("mu", "log_sigma"),
#'   log_normalizer = -log(2 * pi)
#' )
#' spec
fb_log_posterior_spec <- function(producer,
                                  parameter_names = NULL,
                                  log_normalizer = NA_real_,
                                  support_lower = NULL,
                                  support_upper = NULL,
                                  draws = NULL,
                                  name = NULL) {
  if (inherits(producer, "fb_log_posterior_spec")) {
    return(producer)
  }

  ## A fitted model object dispatches to the producer package's own
  ## `fb_log_posterior()` method. proxymix does not implement that method --
  ## it lives in the producer -- so when no producer is installed this branch
  ## raises the seam error rather than guessing.
  if (!is.function(producer)) {
    .fb_abort_no_producer(producer)
  }

  spec <- .fb_spec_from_callable(
    log_density     = producer,
    parameter_names = parameter_names,
    log_normalizer  = log_normalizer,
    support_lower   = support_lower,
    support_upper   = support_upper,
    draws           = draws,
    name            = name
  )
  spec
}

#' @export
print.fb_log_posterior_spec <- function(x, ...) {
  cat(sprintf("<fb_log_posterior_spec>: \"%s\" in p = %d dimensions\n",
              x$name, x$n_dim))
  cat(sprintf("  parameters     : %s\n",
              paste(x$parameter_names, collapse = ", ")))
  cat(sprintf("  log Z          : %s\n",
              if (is.na(x$log_normalizer)) "unknown"
              else format(x$log_normalizer, digits = 6L)))
  bounded <- !is.null(x$support_lower) || !is.null(x$support_upper)
  cat(sprintf("  support bounds : %s\n",
              if (bounded) "supplied" else "unbounded"))
  cat(sprintf("  proposal seed  : %s\n",
              if (!is.null(x$draws))
                sprintf("%d draws", nrow(x$draws))
              else "<none>"))
  invisible(x)
}

# ---------------------------------------------------------------------------
# Capability probe -- degrade gracefully when no real producer is present
# ---------------------------------------------------------------------------

#' Is an external posterior producer available?
#'
#' Capability probe for the posterior-producer seam. Returns `TRUE` only when
#' the package named by `getOption("proxymix.producer_package")` is installed
#' and exports the `fb_log_posterior` generic. When no producer package is
#' named or installed, the probe returns `FALSE`, and [from_fb_posterior()]
#' still works against a bare callable or the synthetic [mock_fb_posterior()]
#' -- the probe lets a caller branch on whether a fitted model object can be
#' compressed directly.
#'
#' The contract is intentionally a soft seam: proxymix never `Imports:` a
#' producer package and `R CMD check` passes with none installed. This probe
#' is the runtime expression of that: it never errors, never loads the
#' producer as a side effect of a negative answer, and degrades to `FALSE`
#' rather than failing.
#'
#' @returns A logical scalar: `TRUE` when an installed producer package
#'   exports `fb_log_posterior`, `FALSE` otherwise.
#' @family interop
#' @seealso [from_fb_posterior()], [fb_log_posterior_spec()].
#' @export
#' @examples
#' ## FALSE on a machine with no producer package configured.
#' fb_producer_available()
fb_producer_available <- function() {
  ## The producer package is a soft, optional peer, deliberately absent from
  ## DESCRIPTION (proxymix must `R CMD check` clean with none installed, and
  ## the producer must never reverse-depend on proxymix). Its name is read
  ## from an option so no unreleased package is hard-coded here.
  pkg <- getOption("proxymix.producer_package", "")
  if (!is.character(pkg) || length(pkg) != 1L || !nzchar(pkg)) {
    return(FALSE)
  }
  if (!requireNamespace(pkg, quietly = TRUE)) {
    return(FALSE)
  }
  exported <- tryCatch(
    getNamespaceExports(pkg),
    error = function(e) character(0L)
  )
  "fb_log_posterior" %in% exported
}

# ---------------------------------------------------------------------------
# Synthetic mock producer -- a known target for tests and examples
# ---------------------------------------------------------------------------

#' A synthetic posterior producer for testing the seam
#'
#' Returns a conforming `fb_log_posterior_spec` whose log-density is a known
#' target, so the consumer path can be exercised end-to-end with no producer
#' package installed. It is the stand-in for a fitted model object: the spec
#' it returns is byte-for-byte the shape [from_fb_posterior()] expects from a
#' real producer.
#'
#' Two shapes are provided. `"gaussian"` is an isotropic `n_dim`-dimensional
#' Gaussian with a known normaliser, against which proxy recovery and absolute
#' KLD can be checked exactly. `"banana"` is the warped 2-D Gaussian of
#' [banana_target()], a non-Gaussian target that genuinely needs more than one
#' mixture component -- the realistic compression case.
#'
#' @param shape One of `"gaussian"` (default) or `"banana"`.
#' @param n_dim Ambient dimension for `shape = "gaussian"`. Ignored (fixed at
#'   2) for `shape = "banana"`.
#' @param sd Standard deviation of the isotropic Gaussian (`shape =
#'   "gaussian"`). Default `1`.
#' @param unnormalised If `TRUE`, the returned log-density carries an arbitrary
#'   additive constant and `log_normalizer` is set so the absolute KLD remains
#'   recoverable. If `FALSE` (the default), the log-density is exactly
#'   normalised and `log_normalizer` is `0`.
#'
#' @returns An `fb_log_posterior_spec`.
#' @family interop
#' @seealso [from_fb_posterior()], [fb_log_posterior_spec()].
#' @export
#' @examples
#' spec <- mock_fb_posterior(shape = "gaussian", n_dim = 2L)
#' spec
#' ## A non-Gaussian target needing several components.
#' mock_fb_posterior(shape = "banana")
mock_fb_posterior <- function(shape = c("gaussian", "banana"),
                              n_dim = 2L,
                              sd = 1,
                              unnormalised = FALSE) {
  shape <- rlang::arg_match(shape)
  n_dim <- as.integer(n_dim)
  if (length(n_dim) != 1L || is.na(n_dim) || n_dim < 1L) {
    cli::cli_abort("`n_dim` must be a positive integer scalar.")
  }

  # Gaussian shape -------------------------------------------------------------
  if (shape == "gaussian") {
    if (length(sd) != 1L || !is.numeric(sd) || !is.finite(sd) || sd <= 0) {
      cli::cli_abort("`sd` must be a positive finite scalar.")
    }
    ## A fixed additive offset makes the evaluated density integrate to
    ## exp(offset), i.e. log Z(evaluated) = offset. The contract's
    ## `log_normalizer` is the correction *added* to the log-density to
    ## normalise it, so it is -offset (see `fit_kld_em()`'s shifted-KLD
    ## correction). For the normalised case both are zero.
    offset <- if (isTRUE(unnormalised)) 17.3 else 0
    log_norm_const <- -n_dim * (0.5 * log(2 * pi) + log(sd))
    log_density <- function(theta) {
      if (is.null(dim(theta))) theta <- matrix(theta, nrow = 1L)
      quad <- rowSums(theta^2) / (2 * sd^2)
      log_norm_const - quad + offset
    }
    return(fb_log_posterior_spec(
      log_density,
      parameter_names = paste0("theta", seq_len(n_dim)),
      log_normalizer = -offset,
      name = sprintf("mock_fb_gaussian[sd=%g]", sd)
    ))
  }

  # Banana shape ---------------------------------------------------------------
  offset <- if (isTRUE(unnormalised)) -9.1 else 0
  log_density <- function(theta) {
    if (is.null(dim(theta))) theta <- matrix(theta, nrow = 1L)
    z1 <- theta[, 1L]
    z2 <- theta[, 2L] - 0.5 * (z1^2 - 1)
    -0.5 * (z1^2 + z2^2) - log(2 * pi) + offset
  }
  fb_log_posterior_spec(
    log_density,
    parameter_names = c("z1", "z2"),
    log_normalizer = -offset,
    name = "mock_fb_banana"
  )
}

# ---------------------------------------------------------------------------
# Consumer verb -- compress an external posterior into a proxy
# ---------------------------------------------------------------------------

#' Compile an external Bayesian posterior into a Gaussian-mixture proxy
#'
#' The consumer entry point of the posterior-producer seam. Takes a producer
#' -- a fitted model object from a producer package, a bare unnormalised
#' log-posterior callable, or a pre-built [fb_log_posterior_spec()] -- and
#' returns a closed-form `N`-component Gaussian-mixture proxy that compresses
#' it, via importance-sampled KLD-EM ([fit_kld_em()], regime (iii)).
#'
#' This generalises proxymix's input source: where [from_kde()] compresses a
#' kernel-density estimate, this compresses a Bayesian posterior addressed
#' only through its (unnormalised) log-density. The fitting machinery is
#' unchanged; only the source differs. The proxy is closed-form
#' marginalisable, conditionable, and samplable -- a few-component object
#' replacing thousands of MCMC draws -- and, when the producer supplies a
#' finite `log_normalizer`, carries an absolute (not shifted) KLD against the
#' posterior.
#'
#' The fitted-object path activates when an installed producer package
#' exports the `fb_log_posterior` generic (see [fb_producer_available()]).
#' Otherwise a fitted model object raises an informative seam error; the
#' bare-callable and [mock_fb_posterior()] paths work with no producer
#' installed. proxymix never depends on a producer package at runtime.
#'
#' Dimensional scope. As with [from_kde()], importance-sampled KLD-EM loses
#' effective sample size sharply in high dimensions. The guard is
#' `n_dim <= 5` (recommended), `n_dim <= 10` (allowed with a warning),
#' `n_dim > 10` (rejected). Inspect [ess_summary()] on the result.
#'
#' @param producer A fitted model object from a producer package, a bare
#'   callable satisfying the [fb_log_posterior_spec()] contract, or an
#'   `fb_log_posterior_spec`.
#' @param N Number of mixture components in the proxy.
#' @param parameter_names Character vector of parameter names, forwarded to
#'   [fb_log_posterior_spec()] for the bare-callable form. Ignored when
#'   `producer` is already a spec.
#' @param log_normalizer Numeric scalar `log Z`, forwarded to
#'   [fb_log_posterior_spec()]. Default `NA_real_`.
#' @param proposal Optional [is_proposal] for the importance-sampling draws.
#'   The default is a multivariate-t (`df = 5`) centred and scaled from the
#'   spec's `draws` when supplied, from its `support_*` bounds otherwise, and
#'   on the origin with a wide isotropic scale as a last resort.
#' @param is_size Number of importance-sampling draws used for fitting. Default
#'   `5000L`.
#' @param max_iter Maximum EM iterations. Forwarded to [fit_kld_em()].
#' @param tol Convergence tolerance. Forwarded to [fit_kld_em()].
#' @param ridge_eps Ridge added to each component covariance at every M-step.
#'   Forwarded to [fit_kld_em()].
#' @param min_ess Minimum effective sample size below which a warning is issued.
#'   Forwarded to [fit_kld_em()].
#' @param seed Optional integer seed for the fitting IS draw.
#' @param validation_size Held-out IS sample size. Forwarded to [fit_kld_em()].
#' @param validation_proposal Optional [is_proposal] for the held-out sample.
#'   Forwarded to [fit_kld_em()].
#' @param validation_seed Optional integer seed for the held-out IS draw.
#'   Forwarded to [fit_kld_em()].
#' @param support_warn Logical. Forwarded to [fit_kld_em()].
#' @param canonicalise Logical. If `TRUE`, the fitted mixture is post-processed
#'   by [gmm_canonicalise()]. Forwarded to [fit_kld_em()].
#'
#' @returns A [gmm_fit] with `regime = "kld"` and metadata recording the
#'   source (`from_fb_posterior` with the producer's `name`, `n_dim`,
#'   `parameter_names`, and `log_normalizer`).
#' @family fitting
#' @family interop
#' @seealso [fb_log_posterior_spec()] for the producer interface,
#'   [mock_fb_posterior()] for a synthetic producer,
#'   [fb_producer_available()] for the capability probe.
#' @export
#' @examples
#' ## Compress a synthetic posterior (no producer package needed).
#' spec <- mock_fb_posterior(shape = "banana")
#' fit <- from_fb_posterior(spec, N = 3L, is_size = 2000L,
#'                          max_iter = 40L, seed = 1L)
#' fit
#' ess_summary(fit)
from_fb_posterior <- function(producer,
                              N = 3L,
                              parameter_names = NULL,
                              log_normalizer = NA_real_,
                              proposal = NULL,
                              is_size = 5000L,
                              max_iter = 100L,
                              tol = 1e-5,
                              ridge_eps = 1e-6,
                              min_ess = 50L,
                              seed = NULL,
                              validation_size = 0L,
                              validation_proposal = NULL,
                              validation_seed = NULL,
                              support_warn = TRUE,
                              canonicalise = TRUE) {
  spec <- fb_log_posterior_spec(
    producer,
    parameter_names = parameter_names,
    log_normalizer  = log_normalizer
  )
  p <- spec$n_dim

  # Dimensional guard ----------------------------------------------------------
  if (p > 10L) {
    cli::cli_abort(c(
      "{.fn from_fb_posterior} currently supports {.code n_dim <= 10} (got n_dim = {p}).",
      "i" = "Importance-sampled KLD-EM loses effective sample size sharply beyond p = 10.",
      "i" = "Marginalise the posterior to a lower-dimensional sub-vector before compressing."
    ))
  }
  if (p > 5L) {
    cli::cli_warn(c(
      "{.fn from_fb_posterior} is well-tested for {.code n_dim <= 5} (got n_dim = {p}).",
      "i" = "Expect importance-sampling ESS to fall beyond p = 5; inspect {.code ess_summary(fit)}."
    ))
  }
  N <- as.integer(N)
  if (length(N) != 1L || is.na(N) || N < 1L) {
    cli::cli_abort("`N` must be a positive integer scalar.")
  }

  # Build the gmm_target from the producer's log-density -----------------------
  normalised <- isTRUE(is.finite(spec$log_normalizer) &&
                         spec$log_normalizer == 0)
  tgt <- gmm_target(
    n_dim          = p,
    log_density    = spec$log_density,
    samples        = spec$draws,
    normalised     = if (normalised) TRUE else NA,
    log_normalizer = spec$log_normalizer,
    name           = spec$name,
    metadata       = list(
      source          = "from_fb_posterior",
      parameter_names = spec$parameter_names
    )
  )

  # Default proposal: seed location and scale from the spec --------------------
  if (is.null(proposal)) {
    proposal <- .fb_default_proposal(spec)
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

  fit@metadata <- utils::modifyList(fit@metadata, list(
    from_fb_posterior = list(
      producer_name   = spec$name,
      n_dim           = p,
      parameter_names = spec$parameter_names,
      log_normalizer  = spec$log_normalizer
    )
  ))
  fit
}

# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

## Normalise a bare callable into the canonical spec list. Reuses the same
## vectorisation / domain-safety probe contract as
## `gmm_target_from_posterior.function()`, so the two interop entry points
## agree on what a conforming callable is.
.fb_spec_from_callable <- function(log_density,
                                   parameter_names,
                                   log_normalizer,
                                   support_lower,
                                   support_upper,
                                   draws,
                                   name) {
  if (is.null(parameter_names)) {
    parameter_names <- attr(log_density, "parameter_names")
  }
  if (is.null(parameter_names)) {
    cli::cli_abort(c(
      "`parameter_names` must be supplied (or attached as {.code attr(producer, 'parameter_names')}).",
      "i" = "Example: {.code from_fb_posterior(my_log_post, parameter_names = c('mu', 'sigma'))}."
    ))
  }
  if (!is.character(parameter_names) || length(parameter_names) < 1L ||
      anyNA(parameter_names) || any(!nzchar(parameter_names))) {
    cli::cli_abort("`parameter_names` must be a non-empty character vector of distinct, non-empty strings.")
  }
  if (anyDuplicated(parameter_names)) {
    cli::cli_abort("`parameter_names` must be unique.")
  }
  n_dim <- length(parameter_names)

  .fb_probe_callable(log_density, n_dim, parameter_names)

  log_z <- if (is.null(log_normalizer) || length(log_normalizer) == 0L) {
    NA_real_
  } else {
    as.numeric(log_normalizer)
  }
  if (length(log_z) != 1L) {
    cli::cli_abort("`log_normalizer` must be a length-1 numeric (or `NA_real_`).")
  }

  support_lower <- .fb_check_support(support_lower, n_dim, "support_lower")
  support_upper <- .fb_check_support(support_upper, n_dim, "support_upper")
  if (!is.null(support_lower) && !is.null(support_upper) &&
      any(support_upper <= support_lower, na.rm = TRUE)) {
    cli::cli_abort("`support_upper` must exceed `support_lower` in every bounded coordinate.")
  }

  draws <- .fb_check_draws(draws, n_dim)

  stored <- log_density
  attr(stored, "parameter_names") <- parameter_names

  structure(
    list(
      log_density     = stored,
      n_dim           = n_dim,
      parameter_names = parameter_names,
      log_normalizer  = log_z,
      support_lower   = support_lower,
      support_upper   = support_upper,
      draws           = draws,
      name            = name %||% "fb_posterior"
    ),
    class = "fb_log_posterior_spec"
  )
}

## Probe a candidate log-density: a pair of zero rows is the cheapest
## non-trivial input. Reject scalar-style or erroring callables here, at
## construction, rather than deep inside the EM loop.
.fb_probe_callable <- function(log_density, n_dim, parameter_names) {
  probe <- matrix(0, nrow = 2L, ncol = n_dim)
  colnames(probe) <- parameter_names
  result <- tryCatch(log_density(probe), error = function(e) e)
  if (inherits(result, "error")) {
    cli::cli_abort(c(
      "The supplied log-posterior raised an error on a probe call.",
      "i" = "Contract: {.code function(theta_matrix)} returning a length-{.code nrow(theta)} numeric, vectorised, domain-safe.",
      "x" = conditionMessage(result)
    ))
  }
  if (!is.numeric(result)) {
    cli::cli_abort(c(
      "The supplied log-posterior returned a {.cls {class(result)[1L]}}, expected numeric.",
      "i" = "Vectorisation contract: {.code function(theta_matrix)} -> {.cls numeric}."
    ))
  }
  if (length(result) != 2L) {
    cli::cli_abort(c(
      "The supplied log-posterior returned length {length(result)}, expected 2 (one per row).",
      "i" = "Vectorisation contract: {.code function(theta_matrix)} returns one value per row."
    ))
  }
  invisible(NULL)
}

## Validate an optional support-bound vector.
.fb_check_support <- function(bound, n_dim, arg_nm) {
  if (is.null(bound)) {
    return(NULL)
  }
  if (!is.numeric(bound) || length(bound) != n_dim) {
    cli::cli_abort("`{arg_nm}` must be a length-{n_dim} numeric vector (or `NULL`).")
  }
  as.numeric(bound)
}

## Validate the optional proposal-seed draws matrix.
.fb_check_draws <- function(draws, n_dim) {
  if (is.null(draws)) {
    return(NULL)
  }
  if (!is.matrix(draws) || !is.numeric(draws)) {
    cli::cli_abort("`draws` must be a numeric matrix (or `NULL`).")
  }
  if (ncol(draws) != n_dim) {
    cli::cli_abort("`draws` must have {n_dim} columns, not {ncol(draws)}.")
  }
  if (anyNA(draws)) {
    cli::cli_abort("`draws` must not contain `NA` or `NaN`.")
  }
  draws
}

## Default importance proposal for a spec: a heavy-tailed multivariate-t,
## located and scaled from whatever the producer supplied -- draws first
## (most informative), then the midpoint / width of any finite support
## bounds, then a wide isotropic default on the origin.
.fb_default_proposal <- function(spec) {
  p <- spec$n_dim
  centre <- rep(0, p)
  scale_diag <- rep(9, p)

  if (!is.null(spec$draws) && nrow(spec$draws) >= 2L) {
    centre <- colMeans(spec$draws)
    scale_m <- ridge(stats::cov(spec$draws), 1e-3)
    return(is_mvt(n_dim = p, mean = centre, sigma = scale_m, df = 5))
  }

  if (!is.null(spec$support_lower) || !is.null(spec$support_upper)) {
    lo <- spec$support_lower %||% rep(NA_real_, p)
    hi <- spec$support_upper %||% rep(NA_real_, p)
    for (j in seq_len(p)) {
      if (is.finite(lo[j]) && is.finite(hi[j])) {
        centre[j] <- 0.5 * (lo[j] + hi[j])
        ## A t with this scale covers the bounded box generously.
        scale_diag[j] <- max(((hi[j] - lo[j]) / 4)^2, 1e-3)
      }
    }
  }

  is_mvt(n_dim = p, mean = centre, sigma = diag(scale_diag, nrow = p), df = 5)
}

## Seam error for a producer that is not a callable and not a recognised
## spec -- the fitted-model-object path when no producer package is present.
.fb_abort_no_producer <- function(producer) {
  cls <- class(producer)[1L]
  if (cls == "function") {
    return(invisible(NULL))
  }
  cli::cli_abort(c(
    "No posterior producer is available for an object of class {.cls {cls}}.",
    "i" = "No installed package supplies an {.fn fb_log_posterior} method for this class.",
    "i" = "Pass a bare callable via {.code from_fb_posterior(my_log_post, parameter_names = ...)}",
    "i" = "or a synthetic producer from {.fn mock_fb_posterior} for testing.",
    "i" = "Probe availability with {.fn fb_producer_available}."
  ))
}
