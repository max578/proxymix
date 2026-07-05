## Missingness-mechanism gates for gmm_impute().
##
## Imputation in proxymix is conditioning the fitted mixture on whatever was
## actually observed. A mechanism gate says what "observed" means for a missing
## entry: nothing beyond the other coordinates (missing at random), a known
## interval (censoring, a detection limit), or a value-dependent propensity
## (missing not at random). Each gate is a small specification object that
## gmm_impute() dispatches on; mar() is the default and reproduces the earlier
## behaviour, so a bare gmm_impute(data) is unchanged.

# ---------------------------------------------------------------------------
# constructors
# ---------------------------------------------------------------------------

#' Missingness mechanisms for multiple imputation
#'
#' Specifications passed to the `mechanism` argument of [gmm_impute()]. They
#' describe how an entry came to be missing, which sets the conditional the
#' missing value is drawn from.
#'
#' `mar()` is missing at random: the probability that an entry is missing may
#' depend on the observed entries but not on the missing value, so the imputation
#' conditional is the plain mixture conditional. This is the default.
#'
#' `censored()` is a known interval. A missing entry of `coord` is known only to
#' lie in `[lower, upper]` -- a detection limit (`upper = LOD`), a ceiling
#' (`lower = cap`), or interval censoring. The imputation conditional is the
#' mixture conditional truncated to that interval, which proxymix evaluates in
#' closed form, so the imputations respect the bound instead of substituting a
#' constant such as half the detection limit.
#'
#' `mnar()` is missing not at random through a selection model: an entry of
#' `coord` is missing with probability \eqn{g(\alpha + \beta\, y)} in its own
#' unobserved value \eqn{y}, where \eqn{g} is the logistic or normal link. The
#' slope `beta` is the sensitivity parameter and is supplied, not estimated --
#' missing-not-at-random departures are not identified from the observed data, so
#' the appropriate use is to posit `beta`, propagate it, and report how conclusions
#' move with it (see [proxy_mnar_sensitivity()]). The intercept is calibrated to
#' the observed missingness rate. `beta = 0` is missing at random.
#'
#' @param coord Name or index of the single coordinate the mechanism acts on.
#' @param beta Sensitivity slope of the log-odds (or probit score) of being
#'   missing in the missing value itself. Positive `beta` makes larger values
#'   more likely to be missing.
#' @param link Selection link, `"logit"` (the default) or `"probit"`.
#' @param lower,upper Bounds of the interval a censored missing entry is known to
#'   lie in. At least one must be finite.
#'
#' @returns A `proxymix_gate` object for [gmm_impute()].
#' @name mechanism
#' @family imputation
#' @seealso [gmm_impute()], [proxy_mnar_sensitivity()].
#' @examples
#' mar()
#' censored("y", upper = 0.5)             # a lower detection limit at 0.5
#' mnar("y", beta = 0.8)                  # larger y more likely missing
NULL

#' @rdname mechanism
#' @export
mar <- function() {
  structure(list(type = "mar"), class = c("proxymix_mar", "proxymix_gate"))
}

#' @rdname mechanism
#' @export
mnar <- function(coord, beta, link = c("logit", "probit")) {
  link <- match.arg(link)
  .check_coord(coord)
  if (!is.numeric(beta) || length(beta) != 1L || !is.finite(beta)) {
    cli::cli_abort("`beta` must be a single finite number.")
  }
  structure(list(type = "mnar", coord = coord, beta = beta, link = link),
            class = c("proxymix_mnar", "proxymix_gate"))
}

#' @rdname mechanism
#' @export
censored <- function(coord, lower = -Inf, upper = Inf) {
  .check_coord(coord)
  if (!is.numeric(lower) || length(lower) != 1L ||
      !is.numeric(upper) || length(upper) != 1L) {
    cli::cli_abort("`lower` and `upper` must be single numbers.")
  }
  if (!(lower < upper)) cli::cli_abort("`lower` must be less than `upper`.")
  if (!is.finite(lower) && !is.finite(upper)) {
    cli::cli_abort("At least one of `lower`, `upper` must be finite.")
  }
  structure(list(type = "censored", coord = coord, lower = lower, upper = upper),
            class = c("proxymix_censored", "proxymix_gate"))
}

## A coordinate reference is a single column name or a positive integer index.
.check_coord <- function(coord) {
  ok <- (is.character(coord) && length(coord) == 1L) ||
    (is.numeric(coord) && length(coord) == 1L && coord >= 1L)
  if (!ok) {
    cli::cli_abort("`coord` must be a single column name or positive index.")
  }
  invisible(coord)
}

## Coerce whatever the user passed as `mechanism` to a gate object. A bare
## "mar" string stays valid for backwards compatibility.
.as_gate <- function(mechanism) {
  if (inherits(mechanism, "proxymix_gate")) return(mechanism)
  if (is.character(mechanism) && length(mechanism) == 1L) {
    return(switch(match.arg(mechanism, c("mar")), mar = mar()))
  }
  cli::cli_abort(c(
    "`mechanism` must be a gate object or {.val mar}.",
    "i" = "Use {.fn mar}, {.fn mnar}, or {.fn censored}."
  ))
}

## Resolve a gate's `coord` to a 1-based column index against `var_names`.
.gate_coord_index <- function(gate, var_names) {
  if (is.null(gate$coord)) return(NA_integer_)
  if (is.character(gate$coord)) {
    j <- match(gate$coord, var_names)
    if (is.na(j)) cli::cli_abort("`coord` {.val {gate$coord}} is not a column.")
  } else {
    j <- as.integer(gate$coord)
    if (j > length(var_names)) cli::cli_abort("`coord` index out of range.")
  }
  j
}

# ---------------------------------------------------------------------------
# print
# ---------------------------------------------------------------------------

#' @export
print.proxymix_gate <- function(x, ...) {
  msg <- switch(x$type,
    mar = "missing at random",
    censored = sprintf("censored on %s to [%g, %g]", x$coord, x$lower, x$upper),
    mnar = sprintf("missing not at random on %s (%s, beta = %g)",
                   x$coord, x$link, x$beta))
  cat(sprintf("<proxymix gate>: %s\n", msg))
  invisible(x)
}
