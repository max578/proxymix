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
## The counterfactual assumes the unit keeps its structural residual under the
## intervention; the data cannot confirm this, and a redrawn residual gives a
## different per-unit mean from the same joint law. `gmm_counterfactual()`
## returns the mean under that assumption and refuses the spread and tails.

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
  .check_quality(g, "gmm_intervene")
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
#' where \eqn{\beta_k^T} is component `k`'s within-class treatment slope.
#'
#' This law assumes the unit keeps its structural residual when the treatment
#' changes. The joint density cannot confirm that assumption: if the residual
#' were redrawn instead, the same mixture would give a different counterfactual
#' law for the unit, mean included. [gmm_cf_mean()], [gmm_cf_variance()] and
#' [gmm_cf_tail_prob()] summarise the law under the same assumption.
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
#'   abduction weights, and the counterfactual `mean`.
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
  .check_quality(g, "gmm_counterfactual")
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
#' per-component counterfactual outcomes) with abduction weights, plus the
#' counterfactual `mean`. All three assume the unit keeps its structural
#' residual under the intervention, which the data cannot confirm. Under that
#' assumption the atom spread is uncertainty about which component generated
#' the unit, reported by [gmm_cf_variance()] and [gmm_cf_tail_prob()].
#'
#' @param atoms Numeric vector of length `K` -- the per-component
#'   counterfactual outcomes.
#' @param weights Numeric vector of length `K` -- the abduction weights
#'   \eqn{\pi_k(\text{evidence})}.
#' @param mean Numeric scalar -- the counterfactual mean.
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
  cat(sprintf("  counterfactual mean : %.6g  (residual kept)\n", x@mean))
  cat(sprintf("  variance            : %.6g  (residual kept)\n",
              gmm_cf_variance(x)))
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
  invisible(x)
}

## Mean accessor ------------------------------------------------------------

#' The counterfactual mean
#'
#' Returns the mean of a [gmm_counterfactual_law]. Like the law itself, it
#' assumes the unit keeps its structural residual under the intervention.
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

## Variance and tail accessors ----------------------------------------------

#' The counterfactual variance
#'
#' Returns the variance of a [gmm_counterfactual_law], the weighted spread of
#' its atoms around the counterfactual mean. Like the law itself, it assumes the
#' unit keeps its structural residual under the intervention. Under that
#' assumption the only uncertainty left is which component generated the unit;
#' if the residual were redrawn, the same mixture would give a different
#' variance.
#'
#' @param x A `gmm_counterfactual_law`.
#'
#' @returns Numeric scalar. Zero when the mixture has one component.
#' @family decision
#' @export
#' @examples
#' g <- gmm(weights = c(0.6, 0.4),
#'          means = list(c(0, 0, 0), c(2, 1, 1)),
#'          covariances = list(diag(3), diag(3)))
#' cf <- gmm_counterfactual(g, evidence = c(1.2, 0, 0.5),
#'                          do = c(NA, 1, NA), query = 1L)
#' gmm_cf_variance(cf)
gmm_cf_variance <- function(x) {
  if (!S7::S7_inherits(x, gmm_counterfactual_law)) {
    cli::cli_abort("`x` must be a {.cls gmm_counterfactual_law}.")
  }
  sum(x@weights * (x@atoms - x@mean)^2)
}

#' A counterfactual tail probability
#'
#' Returns \eqn{P(Y_{t'} > c \mid y, t, x)} for a [gmm_counterfactual_law]:
#' the total weight of the atoms above the threshold `c`. The law is discrete,
#' so the probability changes only where `c` crosses an atom. It rests on the
#' same kept-residual assumption as [gmm_cf_mean()].
#'
#' @param x A `gmm_counterfactual_law`.
#' @param threshold Numeric scalar `c`.
#'
#' @returns Numeric scalar in `[0, 1]`.
#' @family decision
#' @export
#' @examples
#' g <- gmm(weights = c(0.6, 0.4),
#'          means = list(c(0, 0, 0), c(2, 1, 1)),
#'          covariances = list(diag(3), diag(3)))
#' cf <- gmm_counterfactual(g, evidence = c(1.2, 0, 0.5),
#'                          do = c(NA, 1, NA), query = 1L)
#' gmm_cf_tail_prob(cf, threshold = 2)
gmm_cf_tail_prob <- function(x, threshold) {
  if (!S7::S7_inherits(x, gmm_counterfactual_law)) {
    cli::cli_abort("`x` must be a {.cls gmm_counterfactual_law}.")
  }
  if (!is.numeric(threshold) || length(threshold) != 1L || is.na(threshold)) {
    cli::cli_abort("`threshold` must be a single non-missing number.")
  }
  sum(x@weights[x@atoms > threshold])
}
