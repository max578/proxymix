## S7 class hierarchy for proxymix.
##
## Hierarchy:
##   gmm           : a Gaussian mixture (weights, means, covariances).
##   gmm_fit       : extends gmm; adds the target it was fitted to and
##                   the iteration diagnostics.
##   gmm_target    : a target density on R^p, defined by samples and/or
##                   an evaluable log-density.
##   is_proposal   : an importance-sampling proposal — pair of (sampler,
##                   log-density) on R^p.

# ---------------------------------------------------------------------------
# gmm
# ---------------------------------------------------------------------------

#' A Gaussian mixture
#'
#' Lightweight S7 class representing an `N`-component multivariate Gaussian
#' mixture on \eqn{\mathbb{R}^p}. Use `gmm()` to construct, [dgmm()] / [rgmm()]
#' to evaluate or sample, and [gmm_marginalise()] / [gmm_conditionalise()]
#' for closed-form operations.
#'
#' @param weights Numeric vector of length `K`, non-negative, summing to one.
#' @param means List of length `K`, each element a length-`p` numeric vector.
#' @param covariances List of length `K`, each element a `p`-by-`p`
#'   symmetric positive-definite numeric matrix.
#' @param name Optional human-readable name.
#' @param metadata Optional list of arbitrary metadata (regime tags,
#'   diagnostic snapshots, etc.).
#'
#' @returns An S7 object inheriting from `gmm`.
#' @family classes
#' @export
#' @examples
#' g <- gmm(
#'   weights = c(0.4, 0.6),
#'   means = list(c(-1, 0), c(1, 0)),
#'   covariances = list(diag(2), diag(2))
#' )
#' g
gmm <- S7::new_class(
  name = "gmm",
  package = "proxymix",
  properties = list(
    weights = S7::class_double,
    means = S7::class_list,
    covariances = S7::class_list,
    name = S7::new_property(
      class = S7::class_character,
      default = "gmm"
    ),
    metadata = S7::new_property(
      class = S7::class_list,
      default = list()
    )
  ),
  validator = function(self) {
    K <- length(self@weights)
    if (K < 1L) {
      return("`weights` must have at least one element")
    }
    if (length(self@means) != K) {
      return("`means` must have the same length as `weights`")
    }
    if (length(self@covariances) != K) {
      return("`covariances` must have the same length as `weights`")
    }
    if (any(!is.finite(self@weights)) || any(self@weights < 0)) {
      return("`weights` must be finite and non-negative")
    }
    if (abs(sum(self@weights) - 1) > 1e-6) {
      return("`weights` must sum to 1")
    }
    p <- length(self@means[[1L]])
    for (k in seq_len(K)) {
      mu <- self@means[[k]]
      if (length(mu) != p || !is.numeric(mu)) {
        return(sprintf("component %d: `means[[k]]` must be a length-%d numeric vector", k, p))
      }
      S <- self@covariances[[k]]
      if (!is.matrix(S) || nrow(S) != p || ncol(S) != p) {
        return(sprintf("component %d: `covariances[[k]]` must be a %dx%d matrix", k, p, p))
      }
    }
    NULL
  }
)

#' Dimension of a Gaussian mixture
#'
#' Convenience accessor returning the ambient dimension \eqn{p}.
#'
#' @param x A [gmm] (or [gmm_fit]) object.
#'
#' @returns Integer scalar.
#' @family classes
#' @export
#' @examples
#' g <- gmm(weights = 1, means = list(c(0, 0)),
#'          covariances = list(diag(2)))
#' gmm_dim(g)
gmm_dim <- function(x) {
  if (!S7::S7_inherits(x, gmm)) {
    cli::cli_abort("`x` must be a {.cls gmm} object.")
  }
  length(x@means[[1L]])
}

#' Number of components in a Gaussian mixture
#'
#' @param x A [gmm] (or [gmm_fit]) object.
#'
#' @returns Integer scalar.
#' @family classes
#' @export
#' @examples
#' g <- gmm(weights = c(0.5, 0.5), means = list(c(0, 0), c(1, 1)),
#'          covariances = list(diag(2), diag(2)))
#' gmm_n_components(g)
gmm_n_components <- function(x) {
  if (!S7::S7_inherits(x, gmm)) {
    cli::cli_abort("`x` must be a {.cls gmm} object.")
  }
  length(x@weights)
}

# ---------------------------------------------------------------------------
# gmm_fit (extends gmm)
# ---------------------------------------------------------------------------

#' A fitted Gaussian-mixture proxy
#'
#' A `gmm_fit` is the result of [fit_proxymix()] (or one of the regime-specific
#' fitters). It inherits the mixture parameters of [gmm] and adds a record of
#' the target it was fitted to, the regime used, and the iteration
#' diagnostics.
#'
#' @param weights,means,covariances,name,metadata See [gmm].
#' @param target The [gmm_target] the mixture was fitted to.
#' @param regime One of `"moment"`, `"sample"`, `"kld"`.
#' @param diagnostics A list of regime-specific diagnostics
#'   (see [kld_trace()], [ess_trace()]).
#' @param converged Logical scalar.
#' @param iterations Integer scalar.
#' @param call The matched call.
#'
#' @returns An S7 object inheriting from `gmm_fit` (and `gmm`).
#' @family classes
#' @export
#' @examples
#' samples <- matrix(stats::rnorm(200), ncol = 2)
#' tgt <- gmm_target_from_samples(samples)
#' fit <- fit_proxymix(tgt, N = 2L, regime = "sample", max_iter = 25L)
#' inherits(fit, "proxymix::gmm_fit")
gmm_fit <- S7::new_class(
  name = "gmm_fit",
  parent = gmm,
  package = "proxymix",
  properties = list(
    target = S7::class_any,
    regime = S7::new_property(
      class = S7::class_character,
      default = NA_character_
    ),
    diagnostics = S7::new_property(
      class = S7::class_list,
      default = list()
    ),
    converged = S7::new_property(
      class = S7::class_logical,
      default = NA
    ),
    iterations = S7::new_property(
      class = S7::class_integer,
      default = NA_integer_
    ),
    call = S7::class_any
  ),
  validator = function(self) {
    if (length(self@regime) != 1L) {
      return("`regime` must be a single string")
    }
    if (!is.na(self@regime) &&
        !self@regime %in% c("moment", "sample", "kld")) {
      return("`regime` must be one of \"moment\", \"sample\" or \"kld\"")
    }
    NULL
  }
)

# ---------------------------------------------------------------------------
# gmm_target
# ---------------------------------------------------------------------------

#' A target density on R^p
#'
#' An S7 representation of a target density that `proxymix` is asked to
#' approximate. A target may carry an evaluable `log_density`, a matrix of
#' i.i.d. `samples`, or both. Each of the three fitting regimes consumes a
#' different subset:
#'
#' * regime `"moment"` needs samples (or both moments via metadata);
#' * regime `"sample"` needs samples;
#' * regime `"kld"` needs the log-density.
#'
#' Use [gmm_target()] or [gmm_target_from_samples()] to construct.
#'
#' Importance-sampled KLD-EM (regime `"kld"`) only requires `log_density`
#' to be specified up to an unknown additive constant — the self-normalised
#' weights are invariant to scaling. The package's *diagnostics* downstream,
#' however, do depend on normalisation: an importance-sampled KLD estimate
#' against an unnormalised log-density measures
#' \eqn{\widehat{KL}(f \Vert g) - \log Z(f)} rather than
#' \eqn{\widehat{KL}(f \Vert g)}, and a squared-Hellinger Monte Carlo
#' estimate is only meaningful when both densities integrate to one.
#' Declare the target's normalisation explicitly via `normalised` (and,
#' where possible, supply `log_normalizer`) so that the package can label
#' shifted KLDs as shifted and refuse misleading Hellinger reports.
#'
#' @param n_dim Integer scalar — the ambient dimension `p` (the property is
#'   called `n_dim` rather than `dim` because S7 reserves `dim` as an
#'   attribute name).
#' @param log_density Optional function: `function(x)` taking a numeric matrix
#'   `n` by `p` and returning a length-`n` numeric vector of `log f(x)`.
#' @param samples Optional `n` by `p` numeric matrix of i.i.d. samples
#'   from the target.
#' @param normalised Logical scalar declaring whether `log_density`
#'   integrates to one. `TRUE`, `FALSE`, or `NA` (unknown). Defaults to
#'   `NA`. Downstream diagnostics treat `NA` and `FALSE` identically and
#'   label any KLD estimate as shifted.
#' @param log_normalizer Numeric scalar `log Z(f)` of the supplied
#'   `log_density`, if known. Default `NA_real_`. When `normalised = FALSE`
#'   and `log_normalizer` is finite, downstream diagnostics can correct
#'   shifted KLD estimates by `+ log_normalizer`.
#' @param name Human-readable name.
#' @param metadata Optional list of additional descriptors.
#'
#' @returns An S7 object of class `gmm_target`.
#' @family classes
#' @export
#' @examples
#' tgt <- banana_target()
#' tgt
gmm_target <- S7::new_class(
  name = "gmm_target",
  package = "proxymix",
  properties = list(
    n_dim = S7::class_integer,
    log_density = S7::new_property(
      class = S7::class_any,
      default = NULL
    ),
    samples = S7::new_property(
      class = S7::class_any,
      default = NULL
    ),
    normalised = S7::new_property(
      class = S7::class_logical,
      default = NA
    ),
    log_normalizer = S7::new_property(
      class = S7::class_double,
      default = NA_real_
    ),
    name = S7::new_property(
      class = S7::class_character,
      default = "gmm_target"
    ),
    metadata = S7::new_property(
      class = S7::class_list,
      default = list()
    )
  ),
  validator = function(self) {
    if (length(self@n_dim) != 1L || self@n_dim < 1L) {
      return("`n_dim` must be a positive integer of length 1")
    }
    if (is.null(self@log_density) && is.null(self@samples)) {
      return("at least one of `log_density` or `samples` must be supplied")
    }
    if (!is.null(self@log_density) && !is.function(self@log_density)) {
      return("`log_density` must be a function (or NULL)")
    }
    if (!is.null(self@samples)) {
      if (!is.matrix(self@samples) || !is.numeric(self@samples)) {
        return("`samples` must be a numeric matrix (or NULL)")
      }
      if (ncol(self@samples) != self@n_dim) {
        return(sprintf("`samples` must have %d columns, not %d",
                       self@n_dim, ncol(self@samples)))
      }
    }
    if (length(self@normalised) != 1L) {
      return("`normalised` must be a length-1 logical")
    }
    if (length(self@log_normalizer) != 1L) {
      return("`log_normalizer` must be a length-1 numeric")
    }
    NULL
  }
)

# ---------------------------------------------------------------------------
# is_proposal
# ---------------------------------------------------------------------------

#' An importance-sampling proposal
#'
#' An `is_proposal` packages a sampler with its corresponding log-density.
#' Pass one to [fit_kld_em()] (or [fit_proxymix()] with `regime = "kld"`)
#' to plug an alternative proposal into the regime-(iii) loop.
#'
#' @param n_dim Integer scalar — the ambient dimension `p` (the property is
#'   called `n_dim` rather than `dim` because S7 reserves `dim` as an
#'   attribute name).
#' @param sample Function: `function(n)` returning an `n` by `p` numeric
#'   matrix of independent draws.
#' @param log_density Function: `function(x)` taking a numeric matrix `n` by
#'   `p` and returning a length-`n` numeric vector of `log q(x)`.
#' @param name Human-readable name.
#' @param metadata Optional list of additional descriptors.
#'
#' @returns An S7 object of class `is_proposal`.
#' @family classes
#' @export
#' @examples
#' q <- is_mvn(n_dim = 2L, mean = c(0, 0), cov = diag(2))
#' q
is_proposal <- S7::new_class(
  name = "is_proposal",
  package = "proxymix",
  properties = list(
    n_dim = S7::class_integer,
    sample = S7::class_any,
    log_density = S7::class_any,
    name = S7::new_property(
      class = S7::class_character,
      default = "is_proposal"
    ),
    metadata = S7::new_property(
      class = S7::class_list,
      default = list()
    )
  ),
  validator = function(self) {
    if (length(self@n_dim) != 1L || self@n_dim < 1L) {
      return("`dim` must be a positive integer of length 1")
    }
    if (!is.function(self@sample)) {
      return("`sample` must be a function of one argument `n`")
    }
    if (!is.function(self@log_density)) {
      return("`log_density` must be a function of one argument `x`")
    }
    NULL
  }
)

# ---------------------------------------------------------------------------
# print methods
# ---------------------------------------------------------------------------

#' @export
S7::method(print, gmm) <- function(x, ...) {
  cat(sprintf("<%s>: K = %d components in p = %d dimensions\n",
              x@name, gmm_n_components(x), gmm_dim(x)))
  K <- gmm_n_components(x)
  show_max <- min(K, 5L)
  for (k in seq_len(show_max)) {
    cat(sprintf("  [%d] w = %.4f, |mu| = %.4f, tr(Sigma) = %.4f\n",
                k, x@weights[k], sqrt(sum(x@means[[k]]^2)),
                sum(diag(x@covariances[[k]]))))
  }
  if (K > show_max) {
    cat(sprintf("  ... %d more components\n", K - show_max))
  }
  invisible(x)
}

#' @export
S7::method(print, gmm_fit) <- function(x, ...) {
  cat(sprintf("<gmm_fit>: regime = \"%s\", K = %d, p = %d\n",
              x@regime, gmm_n_components(x), gmm_dim(x)))
  cat(sprintf("  target     : %s\n",
              if (!is.null(x@target)) x@target@name else "<none>"))
  cat(sprintf("  iterations : %s\n",
              if (is.na(x@iterations)) "NA" else as.character(x@iterations)))
  cat(sprintf("  converged  : %s\n", as.character(x@converged)))
  K <- gmm_n_components(x)
  show_max <- min(K, 5L)
  for (k in seq_len(show_max)) {
    cat(sprintf("  [%d] w = %.4f, |mu| = %.4f, tr(Sigma) = %.4f\n",
                k, x@weights[k], sqrt(sum(x@means[[k]]^2)),
                sum(diag(x@covariances[[k]]))))
  }
  if (K > show_max) {
    cat(sprintf("  ... %d more components\n", K - show_max))
  }
  invisible(x)
}

#' @export
S7::method(print, gmm_target) <- function(x, ...) {
  cat(sprintf("<gmm_target>: \"%s\" in p = %d dimensions\n",
              x@name, x@n_dim))
  cat(sprintf("  log_density : %s\n",
              if (!is.null(x@log_density)) "supplied" else "<absent>"))
  cat(sprintf("  samples     : %s\n",
              if (!is.null(x@samples))
                sprintf("%d x %d matrix", nrow(x@samples), ncol(x@samples))
              else "<absent>"))
  norm_str <- if (is.na(x@normalised)) {
    "unknown"
  } else if (isTRUE(x@normalised)) {
    "TRUE"
  } else {
    "FALSE"
  }
  cat(sprintf("  normalised  : %s\n", norm_str))
  if (!is.na(x@log_normalizer)) {
    cat(sprintf("  log Z(f)    : %s\n",
                format(x@log_normalizer, digits = 6L)))
  }
  invisible(x)
}

#' @export
S7::method(print, is_proposal) <- function(x, ...) {
  cat(sprintf("<is_proposal>: \"%s\" in p = %d dimensions\n",
              x@name, x@n_dim))
  invisible(x)
}
