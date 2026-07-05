## Convenience constructor for compiling an (unnormalised) Bayesian
## posterior into a proxymix `gmm_target`.
##
## proxymix never depends on any Bayesian fitting package (Stan, pymc,
## brms, ...). The generic dispatches on the supplied object: a raw
## callable goes through `gmm_target_from_posterior.function`; a model
## object goes through a class-specific method that the relevant Bayesian
## package is expected to register.

#' Compile an unnormalised Bayesian posterior into a `gmm_target`
#'
#' Generic S3 constructor that turns a Bayesian posterior — represented
#' either by a fitted model object (e.g. from `brms` or `Stan`) or by a
#' bare callable — into a `gmm_target` suitable for regime (iii) of
#' [fit_proxymix()] / [fit_kld_em()].
#'
#' The contract for the underlying callable is:
#'
#' * **Vectorised**: accepts a numeric matrix with rows indexing
#'   independent parameter draws and columns indexing parameters; returns
#'   a length-`nrow(theta)` numeric vector of `log p(theta | data) + const`.
#' * **Unnormalised is fine**: the marginal likelihood `log Z` is not
#'   required. Where the source package can supply it, pass
#'   `log_normalizer`.
#' * **Side-effect free**: no plotting, no mutable state. Pure function.
#' * **Domain-safe**: returns `-Inf` outside support rather than raising
#'   an error.
#'
#' The default method errors with a hint pointing the user at either
#' (a) a Bayesian package that registers a method, or (b) the
#' `function` method below.
#'
#' @param model One of:
#'   * a function — a bare callable satisfying the contract above;
#'   * a fitted model object whose class registers a
#'     `gmm_target_from_posterior.<class>` method in its own package.
#' @param ... Forwarded to method-specific implementations.
#'
#' @returns A [gmm_target] with `normalised = FALSE` and the user-supplied
#'   `log_normalizer` (or `NA_real_`).
#' @family interop
#' @export
#' @examples
#' # A trivial unnormalised log-posterior: a 2D banana centred near (1, 0).
#' log_post <- function(theta) {
#'   x <- theta[, 1L]
#'   y <- theta[, 2L]
#'   -0.5 * (x^2 + (y - 0.1 * x^2 + 1)^2)
#' }
#' tgt <- gmm_target_from_posterior(
#'   log_post,
#'   parameter_names = c("x", "y")
#' )
#' tgt
gmm_target_from_posterior <- function(model, ...) {
  UseMethod("gmm_target_from_posterior")
}

#' @rdname gmm_target_from_posterior
#' @export
gmm_target_from_posterior.default <- function(model, ...) {
  cli::cli_abort(c(
    "No {.fn gmm_target_from_posterior} method registered for class {.cls {class(model)[1L]}}.",
    "i" = "Either load a Bayesian package that registers a method for this class,",
    "i" = "or pass a bare callable via {.code gmm_target_from_posterior(my_log_post, parameter_names = c('a', 'b'))}.",
    "i" = "See {.help [gmm_target_from_posterior](proxymix::gmm_target_from_posterior)} for the contract."
  ))
}

#' @rdname gmm_target_from_posterior
#' @param parameter_names Character vector of parameter names. Required
#'   for the `function` method (or attached as
#'   `attr(model, "parameter_names")`). The length determines `n_dim`.
#' @param log_normalizer Numeric scalar `log Z` of the posterior, if known.
#'   `NA_real_` (the default) otherwise; downstream diagnostics will label
#'   any KLD estimate as shifted.
#' @param name Optional human-readable target name. Defaults to
#'   `"posterior"`.
#' @export
gmm_target_from_posterior.function <- function(model, ...,
                                               parameter_names = NULL,
                                               log_normalizer = NA_real_,
                                               name = NULL) {
  if (is.null(parameter_names)) {
    parameter_names <- attr(model, "parameter_names")
  }
  if (is.null(parameter_names)) {
    cli::cli_abort(c(
      "`parameter_names` must be supplied (or attached as {.code attr(model, 'parameter_names')}).",
      "i" = "Example: {.code gmm_target_from_posterior(my_log_post, parameter_names = c('mu', 'sigma'))}."
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

  ## Vectorisation probe. A pair of zero rows is the cheapest non-trivial
  ## input; reject the callable here if it is scalar-style.
  probe <- matrix(0, nrow = 2L, ncol = n_dim)
  colnames(probe) <- parameter_names
  result <- tryCatch(model(probe), error = function(e) e)
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

  log_post <- model
  attr(log_post, "parameter_names") <- parameter_names

  log_norm <- if (is.null(log_normalizer) || length(log_normalizer) == 0L) {
    NA_real_
  } else {
    as.numeric(log_normalizer)
  }
  if (length(log_norm) != 1L) {
    cli::cli_abort("`log_normalizer` must be a length-1 numeric (or `NA_real_`).")
  }

  gmm_target(
    n_dim          = n_dim,
    log_density    = log_post,
    samples        = NULL,
    normalised     = FALSE,
    log_normalizer = log_norm,
    name           = name %||% "posterior",
    metadata       = list(
      source          = "from_posterior",
      parameter_names = parameter_names,
      model_class     = class(model)[1L]
    )
  )
}
