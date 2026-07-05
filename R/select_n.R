## Principled component-count selection across the regimes. Regime (ii)
## has an empirical likelihood, so BIC/ICL apply; regime (iii) does not
## (the information criteria are deliberately NA there), so the selector
## scores candidates on a held-out importance draw and applies the
## one-standard-error rule.

#' Select the number of mixture components
#'
#' Fits every candidate component count and chooses one by a
#' regime-appropriate criterion. With samples (regime ii) the choice is
#' the smallest BIC. With an evaluable-only target (regime iii) each
#' candidate is scored by its held-out validation KLD -- an independent
#' importance draw the fit never trained on -- and the choice follows the
#' one-standard-error rule: the smallest `N` whose validation score is
#' within one Monte Carlo standard error of the best score. The scored
#' table is returned alongside the choice, and callers who prefer to
#' choose by eye can ignore the recommendation.
#'
#' @param target A [gmm_target].
#' @param candidates Integer vector of component counts to try.
#' @param regime `"auto"` (default: `"sample"` when the target carries
#'   samples, `"kld"` otherwise), `"sample"`, or `"kld"`.
#' @param seed Optional integer seed shared across candidates (paired
#'   fits).
#' @param ... Forwarded to the regime's fitter ([fit_em_samples()] or
#'   [fit_kld_em()]).
#'
#' @returns A list of class `proxymix_selection`: `best_n` (the chosen
#'   count), `best_fit` (its [gmm_fit]), and `table` (one row per
#'   candidate with the criterion values).
#' @family fitting
#' @export
#' @examples
#' sel <- select_N(banana_target(), candidates = 1:3,
#'                 is_size = 1500L, max_iter = 20L, seed = 1L)
#' sel$table
#' sel$best_n
select_N <- function(target, candidates = 1:6,
                     regime = c("auto", "sample", "kld"),
                     seed = NULL, ...) {
  if (!S7::S7_inherits(target, gmm_target)) {
    cli::cli_abort("`target` must be a {.cls gmm_target}.")
  }
  regime <- rlang::arg_match(regime)
  candidates <- sort(unique(as.integer(candidates)))
  if (anyNA(candidates) || any(candidates < 1L)) {
    cli::cli_abort("`candidates` must be positive integers.")
  }
  if (regime == "auto") {
    regime <- if (!is.null(target@samples)) "sample" else "kld"
  }

  if (regime == "sample") {
    if (is.null(target@samples)) {
      cli::cli_abort("regime {.val sample} requires the target to carry samples.")
    }
    fits <- lapply(candidates, function(N) {
      fit_em_samples(target, N = N, seed = seed, ...)
    })
    bics <- vapply(fits, function(f) f@diagnostics$bic, numeric(1L))
    icls <- vapply(fits, function(f) bic_aic(f)$icl, numeric(1L))
    best <- which.min(bics)
    tab <- data.frame(N = candidates, bic = bics, icl = icls,
                      chosen = seq_along(candidates) == best)
  } else {
    if (is.null(target@log_density)) {
      cli::cli_abort("regime {.val kld} requires the target to carry a log-density.")
    }
    fits <- lapply(candidates, function(N) {
      fit_kld_em(target, N = N, seed = seed, ...)
    })
    scores <- vapply(fits, function(f) f@diagnostics$validation_kld,
                     numeric(1L))
    ses <- vapply(fits, function(f) {
      f@diagnostics$validation_mc_se %||% NA_real_
    }, numeric(1L))
    if (anyNA(scores)) {
      cli::cli_abort(c(
        "regime {.val kld} selection needs a validation split.",
        "i" = "Do not pass {.code validation_size = 0} through {.code ...}."
      ))
    }
    i_min <- which.min(scores)
    band <- scores[i_min] + (if (is.finite(ses[i_min])) ses[i_min] else 0)
    ## One-standard-error rule: the smallest N whose held-out score is
    ## within one MC standard error of the best.
    best <- which(scores <= band)[1L]
    tab <- data.frame(N = candidates, validation_kld = scores,
                      validation_mc_se = ses,
                      chosen = seq_along(candidates) == best)
  }

  structure(
    list(best_n = candidates[best], best_fit = fits[[best]],
         regime = regime, table = tab),
    class = "proxymix_selection"
  )
}

#' @export
print.proxymix_selection <- function(x, ...) {
  cat(sprintf("Component-count selection (regime \"%s\"): N = %d\n",
              x$regime, x$best_n))
  print(x$table, row.names = FALSE, digits = 4L)
  invisible(x)
}
