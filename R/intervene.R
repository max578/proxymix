## intervene.R -- The do-operator and the K-atom counterfactual.
##
## Promotes the operators `gmm_intervene()` and
## `gmm_counterfactual()` to first-class functions. Both read a fitted joint
## Gaussian mixture
## as a latent-class structural causal model: each component is an affine
## structural equation and the prior `pi_k` is the regime prevalence. The
## do-operator severs the incoming edges of the intervened coordinates -- it
## sets them inside each component but does not let them re-weight the regime
## gate. The counterfactual adds abduction (recover the regime posterior and
## the structural residual of an observed unit) before the same surgery.
##
## What is identified: interventional and counterfactual *means*. What is not:
## the individual counterfactual *law* (its spread / tail probabilities), which
## depends on an unidentified cross-world coupling. `gmm_counterfactual()`
## returns the mean and refuses the rest -- only the mean is identified.

## Validation ---------------------------------------------------------------

.validate_coord_vector <- function(v, p, arg) {
  if (is.logical(v) && all(is.na(v))) {
    v <- as.double(v)
  }
  if (!is.numeric(v) || length(v) != p) {
    cli::cli_abort("`{arg}` must be a numeric vector of length {p} (NA marks non-participating coordinates).")
  }
  active <- which(!is.na(v))
  if (length(active) > 0L && (anyNA(v[active]) || !all(is.finite(v[active])))) {
    cli::cli_abort("the active entries of `{arg}` must be finite.")
  }
  active
}

## gmm_intervene ------------------------------------------------------------

#' Interventional law of a Gaussian mixture (the do-operator)
#'
#' Reads a fitted Gaussian mixture as a latent-class structural causal model
#' and returns the interventional distribution of the free coordinates under
#' `do()` of some coordinates, optionally conditioning on others.
#'
#' Intervened (`do`) coordinates are *set* inside every component but do not
#' re-weight the regime gate -- this is the graph surgery that distinguishes
#' \eqn{p(\cdot \mid do(T = t))} from \eqn{p(\cdot \mid T = t)}. Conditioned
#' (`given`) coordinates re-weight the gate in the usual Bayesian way. Writing
#' the component prior as \eqn{\pi_k}, the within-component conditional mean as
#' \eqn{\mu_k}, and the given-coordinate evidence as \eqn{e_k}, the returned
#' mixture has weights \eqn{\pi_k(given) \propto \pi_k\, e_k} and per-component
#' parameters from the Schur conditional on the union of the `do` and `given`
#' coordinates.
#'
#' For a joint fit over \eqn{(Y, T, X)}, `gmm_intervene(fit, do = T = 1,
#' given = X = x)` is the do-response \eqn{p(Y \mid do(T = 1), X = x)}; its mean
#' is \eqn{\sum_k \pi_k(x)\, \mu_k^{y \mid 1, x}}, the latent-confounder-mode
#' interventional mean of \link{proxy_cate}.
#'
#' @param g A [gmm] (or [gmm_fit]) in `R^p`.
#' @param do A length-`p` numeric vector. Coordinates to intervene on take
#'   their `do`-value; coordinates not intervened on are `NA`.
#' @param given A length-`p` numeric vector, or `NULL` (the default, meaning
#'   no conditioning). Coordinates to condition on take their value;
#'   coordinates not conditioned on are `NA`. A coordinate may not appear in
#'   both `do` and `given`.
#' @param ridge_eps Tiny ridge added to the conditional covariances for
#'   numerical hygiene. Set to `0` to disable.
#'
#' @returns A [gmm] over the free coordinates (those `NA` in both `do` and
#'   `given`), with weights re-gated by the `given` evidence only.
#' @family operators
#' @family decision
#' @seealso [gmm_counterfactual()], [proxy_cate()]
#' @export
#' @examples
#' ## Joint (Y, T, X): set T = 1 while conditioning on X = 0.3.
#' g <- gmm(weights = c(0.5, 0.5),
#'          means = list(c(0, 0, 0), c(2, 1, 1)),
#'          covariances = list(diag(3), diag(3)))
#' gmm_intervene(g, do = c(NA, 1, NA), given = c(NA, NA, 0.3))
gmm_intervene <- function(g, do, given = NULL, ridge_eps = 1e-6) {
  if (!S7::S7_inherits(g, gmm)) {
    cli::cli_abort("`g` must be a {.cls gmm} object.")
  }
  p <- gmm_dim(g)
  do_idx <- .validate_coord_vector(do, p, "do")
  if (is.null(given)) {
    given <- rep(NA_real_, p)
  }
  given_idx <- .validate_coord_vector(given, p, "given")

  if (length(do_idx) == 0L) {
    cli::cli_abort("`do` must intervene on at least one coordinate.")
  }
  if (length(intersect(do_idx, given_idx)) > 0L) {
    cli::cli_abort("a coordinate may not appear in both `do` and `given`.")
  }
  free <- setdiff(seq_len(p), c(do_idx, given_idx))
  if (length(free) == 0L) {
    cli::cli_abort("`do` and `given` together leave no free coordinate.")
  }

  # Regime gate: re-weight by the `given` evidence only -----------------------
  # The do-coordinates are deliberately excluded -- that omission is the
  # do-operator's graph surgery.
  K <- gmm_n_components(g)
  if (length(given_idx) > 0L) {
    gate <- .abduction_responsibilities(g, given_idx, given[given_idx])
  } else {
    gate <- g@weights
  }

  # Per-component Schur conditional on (do union given) -----------------------
  cond_idx <- c(do_idx, given_idx)
  cond_val <- c(do[do_idx], given[given_idx])
  new_means <- vector("list", K)
  new_covs <- vector("list", K)
  for (k in seq_len(K)) {
    sc <- .schur_cond_component(g@means[[k]], g@covariances[[k]],
                                free, cond_idx, cond_val)
    new_covs[[k]] <- if (ridge_eps > 0) ridge(sc$cov, ridge_eps) else sc$cov
    new_means[[k]] <- sc$mean
  }

  gmm(
    weights = gate,
    means = new_means,
    covariances = new_covs,
    name = sprintf("intervene(%s)", g@name),
    metadata = list(
      operator = "gmm_intervene",
      do_idx = do_idx,
      given_idx = given_idx,
      free_idx = free
    )
  )
}

## gmm_counterfactual -------------------------------------------------------

#' Counterfactual law of one unit (abduction, action, prediction)
#'
#' Computes the per-unit counterfactual of a `query` coordinate (the outcome)
#' for an observed unit, under a `do()` intervention -- Pearl's third rung on
#' the latent-class structural causal model read off a fitted Gaussian
#' mixture.
#'
#' The three steps are closed form. *Abduction* recovers the regime posterior
#' \eqn{\pi_k(\text{evidence})} of the observed unit and, within each
#' component, its structural residual. *Action* sets the `do` coordinates.
#' *Prediction* re-evaluates the query coordinate. For a binary treatment the
#' result is a discrete law on `K` atoms,
#' \deqn{Y_{t'} \mid (y, t, x) = \sum_k \pi_k(y, t, x)\,
#'       \delta\!\big(y + \beta_k^{T}(t' - t)\big),}
#' where \eqn{\beta_k^T} is component `k`'s within-class treatment slope. Only
#' the **mean** of this law is identified; its spread reflects regime
#' uncertainty, not the (unidentified) cross-world coupling, so the variance
#' and tail accessors refuse to answer (see [gmm_cf_variance()]).
#'
#' @param g A [gmm] (or [gmm_fit]) in `R^p`.
#' @param evidence A length-`p` numeric vector of the observed unit.
#'   Unobserved coordinates are `NA`; the `query` coordinate must be observed.
#' @param do A length-`p` numeric vector. Intervened coordinates take their
#'   counterfactual value; the rest are `NA`. Intervened coordinates must be
#'   observed in `evidence`.
#' @param query A single integer coordinate index in `seq_len(p)` -- the
#'   outcome whose counterfactual is sought.
#' @param ridge_eps Tiny ridge added to the conditioning covariances for
#'   numerical hygiene.
#'
#' @returns A `gmm_counterfactual_law` object carrying the `K` atoms, their
#'   abduction weights, and the identified counterfactual `mean`.
#' @family operators
#' @family decision
#' @seealso [gmm_intervene()], [gmm_cf_variance()], [proxy_retrospective_uplift()]
#' @export
#' @examples
#' ## Observed unit (y = 1.2, t = 0, x = 0.5); imagine t = 1.
#' g <- gmm(weights = c(0.6, 0.4),
#'          means = list(c(0, 0, 0), c(2, 1, 1)),
#'          covariances = list(diag(3), diag(3)))
#' cf <- gmm_counterfactual(g, evidence = c(1.2, 0, 0.5),
#'                          do = c(NA, 1, NA), query = 1L)
#' cf@mean
gmm_counterfactual <- function(g, evidence, do, query, ridge_eps = 1e-6) {
  if (!S7::S7_inherits(g, gmm)) {
    cli::cli_abort("`g` must be a {.cls gmm} object.")
  }
  p <- gmm_dim(g)
  obs_idx <- .validate_coord_vector(evidence, p, "evidence")
  do_idx <- .validate_coord_vector(do, p, "do")
  query <- as.integer(query)
  if (length(query) != 1L || is.na(query) || query < 1L || query > p) {
    cli::cli_abort("`query` must be a single coordinate index in 1..{p}.")
  }
  if (length(do_idx) == 0L) {
    cli::cli_abort("`do` must intervene on at least one coordinate.")
  }
  if (!query %in% obs_idx) {
    cli::cli_abort("the `query` coordinate must be observed in `evidence`.")
  }
  if (!all(do_idx %in% obs_idx)) {
    cli::cli_abort("every intervened coordinate must be observed in `evidence`.")
  }
  if (query %in% do_idx) {
    cli::cli_abort("`query` may not also be an intervened coordinate.")
  }

  # Abduction: regime posterior of the observed unit ---------------------------
  resp <- .abduction_responsibilities(g, obs_idx, evidence[obs_idx])

  # Factual vs counterfactual conditioning sets --------------------------------
  # Both condition the query coordinate on the other observed coordinates; the
  # counterfactual set overwrites the do-coordinates with their do-values.
  cond_idx <- setdiff(obs_idx, query)
  cond_factual <- evidence[cond_idx]
  cond_cf <- cond_factual
  cond_cf[match(do_idx, cond_idx)] <- do[do_idx]

  # The counterfactual returns means only, so the within-component conditioning
  # uses the raw covariance (ridging the input would bias the mean); `ridge_eps`
  # is the singular-block fallback inside `.schur_cond_component()`.
  K <- gmm_n_components(g)
  atoms <- numeric(K)
  for (k in seq_len(K)) {
    mu <- g@means[[k]]
    S <- g@covariances[[k]]
    mu_fac <- .schur_cond_component(mu, S, query, cond_idx, cond_factual,
                                    ridge_eps = ridge_eps)$mean
    mu_cf <- .schur_cond_component(mu, S, query, cond_idx, cond_cf,
                                   ridge_eps = ridge_eps)$mean
    atoms[k] <- evidence[query] + (mu_cf - mu_fac)
  }
  cf_mean <- sum(resp * atoms)

  gmm_counterfactual_law(
    atoms = atoms,
    weights = resp,
    mean = cf_mean,
    query = query,
    evidence = as.numeric(evidence),
    do = as.numeric(do),
    name = sprintf("counterfactual(%s)", g@name)
  )
}

## gmm_counterfactual_law class --------------------------------------------

#' A per-unit counterfactual law
#'
#' The return type of [gmm_counterfactual()]: a discrete law on `K` atoms (the
#' per-component counterfactual outcomes) with abduction weights, plus the one
#' identified summary -- the counterfactual `mean`. The atom spread is regime
#' (epistemic) uncertainty about which structural component generated the unit;
#' it is **not** the counterfactual outcome's dispersion, which is unidentified.
#' Accordingly [gmm_cf_variance()] and [gmm_cf_tail_prob()] refuse to answer.
#'
#' @param atoms Numeric vector of length `K` -- the per-component
#'   counterfactual outcomes.
#' @param weights Numeric vector of length `K` -- the abduction weights
#'   \eqn{\pi_k(\text{evidence})}.
#' @param mean Numeric scalar -- the identified counterfactual mean.
#' @param query Integer scalar -- the queried coordinate index.
#' @param evidence Numeric vector -- the observed unit.
#' @param do Numeric vector -- the intervention.
#' @param name Human-readable name.
#' @param metadata Optional list of descriptors.
#'
#' @returns An S7 object of class `gmm_counterfactual_law`.
#' @family classes
#' @export
gmm_counterfactual_law <- S7::new_class(
  name = "gmm_counterfactual_law",
  package = "proxymix",
  properties = list(
    atoms = S7::class_double,
    weights = S7::class_double,
    mean = S7::class_double,
    query = S7::class_integer,
    evidence = S7::class_double,
    do = S7::class_double,
    name = S7::new_property(
      class = S7::class_character,
      default = "gmm_counterfactual_law"
    ),
    metadata = S7::new_property(
      class = S7::class_list,
      default = list()
    )
  ),
  validator = function(self) {
    if (length(self@atoms) != length(self@weights)) {
      return("`atoms` and `weights` must have the same length")
    }
    if (length(self@mean) != 1L) {
      return("`mean` must be a length-1 numeric")
    }
    NULL
  }
)

#' @export
S7::method(print, gmm_counterfactual_law) <- function(x, ...) {
  cat(sprintf("<gmm_counterfactual_law>: %d atoms\n", length(x@atoms)))
  cat(sprintf("  counterfactual mean : %.6g  (identified)\n", x@mean))
  ord <- order(x@weights, decreasing = TRUE)
  show_max <- min(length(x@atoms), 4L)
  for (j in seq_len(show_max)) {
    k <- ord[j]
    cat(sprintf("  atom[%d] value = %.4g, weight = %.4f\n",
                k, x@atoms[k], x@weights[k]))
  }
  if (length(x@atoms) > show_max) {
    cat(sprintf("  ... %d more atoms\n", length(x@atoms) - show_max))
  }
  cat("  variance / tail   : not identified (cross-world coupling)\n")
  invisible(x)
}

## Identified accessor ------------------------------------------------------

#' The identified counterfactual mean
#'
#' Returns the one identified summary of a [gmm_counterfactual_law] -- its mean.
#'
#' @param x A `gmm_counterfactual_law`.
#'
#' @returns Numeric scalar.
#' @family decision
#' @export
#' @examples
#' g <- gmm(weights = 1, means = list(c(0, 0, 0)),
#'          covariances = list(diag(3)))
#' cf <- gmm_counterfactual(g, evidence = c(1, 0, 0.2),
#'                          do = c(NA, 1, NA), query = 1L)
#' gmm_cf_mean(cf)
gmm_cf_mean <- function(x) {
  if (!S7::S7_inherits(x, gmm_counterfactual_law)) {
    cli::cli_abort("`x` must be a {.cls gmm_counterfactual_law}.")
  }
  x@mean
}

## Refused accessors (identification limits) --------------------------------

#' Refused: the variance of an individual counterfactual law
#'
#' The per-unit counterfactual *variance* is not identified from the joint
#' density: it depends on the cross-world coupling of the structural residuals
#' under the factual and counterfactual treatments, which no goodness-of-fit
#' can certify. This accessor therefore raises an error rather than return the
#' (misleading) spread of the abduction atoms.
#'
#' @param x A [gmm_counterfactual_law].
#'
#' @returns Never returns; always raises a `proxymix_not_identified` error.
#' @family decision
#' @export
#' @examples
#' g <- gmm(weights = 1, means = list(c(0, 0, 0)),
#'          covariances = list(diag(3)))
#' cf <- gmm_counterfactual(g, evidence = c(1, 0, 0.2),
#'                          do = c(NA, 1, NA), query = 1L)
#' try(gmm_cf_variance(cf))
gmm_cf_variance <- function(x) {
  if (!S7::S7_inherits(x, gmm_counterfactual_law)) {
    cli::cli_abort("`x` must be a {.cls gmm_counterfactual_law}.")
  }
  cli::cli_abort(
    c("The individual counterfactual variance is not identified.",
      "i" = "It depends on the unidentified cross-world coupling of the structural residuals.",
      "i" = "Only the counterfactual mean is identified; see {.fn gmm_cf_mean}."),
    class = "proxymix_not_identified"
  )
}

#' Refused: a tail probability of an individual counterfactual law
#'
#' Like [gmm_cf_variance()], a per-unit counterfactual tail probability
#' \eqn{P(Y_{t'} > c \mid y, t, x)} is not identified -- it is a functional of
#' the unidentified counterfactual law, not of its identified mean. This
#' accessor refuses rather than mislead.
#'
#' @param x A [gmm_counterfactual_law].
#' @param threshold Numeric scalar `c`.
#'
#' @returns Never returns; always raises a `proxymix_not_identified` error.
#' @family decision
#' @export
#' @examples
#' g <- gmm(weights = 1, means = list(c(0, 0, 0)),
#'          covariances = list(diag(3)))
#' cf <- gmm_counterfactual(g, evidence = c(1, 0, 0.2),
#'                          do = c(NA, 1, NA), query = 1L)
#' try(gmm_cf_tail_prob(cf, threshold = 2))
gmm_cf_tail_prob <- function(x, threshold) {
  if (!S7::S7_inherits(x, gmm_counterfactual_law)) {
    cli::cli_abort("`x` must be a {.cls gmm_counterfactual_law}.")
  }
  cli::cli_abort(
    c("Individual counterfactual tail probabilities are not identified.",
      "i" = "They are functionals of the unidentified counterfactual law, not its mean.",
      "i" = "Only the counterfactual mean is identified; see {.fn gmm_cf_mean}."),
    class = "proxymix_not_identified"
  )
}
