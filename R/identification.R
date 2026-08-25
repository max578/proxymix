## identification.R -- The audit layer: what is identified, assumed, and not.
##
## The decision module's differentiator is that every effect carries its
## identification status. These verbs price the unobserved-confounding risk
## (`proxy_confounding_gap()`), value a deployable policy off-line
## (`proxy_policy_value()`), read per-unit retrospective uplift off the
## counterfactual mean (`proxy_retrospective_uplift()`), expose the fitted
## regimes (`proxy_regime_segments()`), and assemble the executive one-pager
## (`proxy_identification_report()`).

## proxy_confounding_gap -----------------------------------------------------

#' Confounding gap: the sensitivity of the effect to the latent regime
#'
#' Per-unit difference between the ignorability-mode and do-mode effects,
#' \eqn{\Delta(x) = \tau_{\mathrm{obs}}(x) - \tau_{\mathrm{do}}(x)}. Under
#' ignorability the two coincide and \eqn{\Delta \equiv 0}; a non-zero gap is a
#' **sensitivity signal** -- how much the estimated effect would move if a
#' fitted regime confounded treatment and outcome beyond `X` -- not a
#' correction the data licenses.
#'
#' @param model An [uplift_model].
#' @param newdata A data frame carrying the covariate columns.
#' @param t1,t0 The treated and control treatment values. Default the treatment levels observed at fit time (`model@treatment_levels`); a value matching neither observed level aborts.
#'
#' @returns A [data.table::data.table] with columns `id`, `tau_obs`, `tau_do`,
#'   `gap`, `overlap_flag`.
#' @family decision
#' @seealso [proxy_cate()], [proxy_identification_report()]
#' @export
#' @examples
#' set.seed(1)
#' n <- 600L
#' x <- stats::rnorm(n)
#' t <- stats::rbinom(n, 1L, 0.5)
#' y <- 0.5 * t + x + stats::rnorm(n, sd = 0.5)
#' dat <- data.frame(y = y, t = t, x = x)
#' m <- fit_uplift(dat, "y", "t", "x", N = 2L, regime = "sample",
#'                 max_iter = 80L, seed = 1L)
#' proxy_confounding_gap(m, data.frame(x = c(-1, 0, 1)))
proxy_confounding_gap <- function(model, newdata, t1 = NULL, t0 = NULL) {
  if (!S7::S7_inherits(model, uplift_model)) {
    cli::cli_abort("`model` must be an {.cls uplift_model}.")
  }
  resolved <- .resolve_treatment_arms(model, t1, t0)
  t1 <- resolved$t1
  t0 <- resolved$t0
  g <- model@fit
  z_idx <- c(model@roles$treatment, model@roles$covariate)
  x_idx <- model@roles$covariate
  X <- .newdata_x(model, newdata)
  cache <- .uplift_cache(model)
  arms <- c(t0, t1)

  n_units <- nrow(X)
  tau_obs <- numeric(n_units)
  tau_do <- numeric(n_units)
  flag <- logical(n_units)
  for (u in seq_len(n_units)) {
    obs <- .cate_unit(cache, g, z_idx, x_idx, X[u, ], t1, t0, "ignorability")
    do_ <- .cate_unit(cache, g, z_idx, x_idx, X[u, ], t1, t0, "latent_confounder")
    tau_obs[u] <- obs$tau
    tau_do[u] <- do_$tau
    flag[u] <- (!obs$covered) ||
      .unit_coverage(g, z_idx, arms, X[u, ]) < 0.01
  }
  dt <- data.table::data.table(
    id = seq_len(n_units),
    tau_obs = tau_obs,
    tau_do = tau_do,
    gap = tau_obs - tau_do,
    overlap_flag = flag
  )
  return(dt[])
}

## proxy_policy_value -------------------------------------------------------

## Resolve a policy specification into a 0/1 action per unit.
.resolve_policy <- function(policy, model, newdata, value, cost, t1, t0) {
  n <- nrow(.newdata_x(model, newdata))
  if (is.character(policy) && length(policy) == 1L) {
    return(switch(policy,
      all = rep(1L, n),
      none = rep(0L, n),
      optimal = proxy_decide(model, newdata, value = value, cost = cost,
                             t1 = t1, t0 = t0)$action,
      cli::cli_abort("`policy` string must be \"all\", \"none\" or \"optimal\".")
    ))
  }
  if (is.function(policy)) {
    ce <- proxy_cate(model, newdata, t1 = t1, t0 = t0, se = FALSE)
    act <- as.integer(policy(ce))
  } else {
    act <- as.integer(policy)
  }
  if (length(act) != n || anyNA(act) || !all(act %in% c(0L, 1L))) {
    cli::cli_abort("`policy` must yield one 0/1 action per row of `newdata`.")
  }
  act
}

#' Off-line value of a targeting policy
#'
#' Estimates the expected value of deploying a per-unit targeting policy,
#' \eqn{V(d) = E_X[\,\text{value}\cdot E[Y \mid do(T = d(X)), X] - \text{cost}
#' \cdot d(X)\,]}, from the fitted model alone -- no live A/B test. Units that
#' fail the overlap diagnostic are excluded by default and their count is
#' reported, never silently dropped.
#'
#' @param model An [uplift_model].
#' @param newdata A data frame carrying the covariate columns -- the population
#'   the policy would be deployed on.
#' @param policy A per-unit action specification: a 0/1 vector of length
#'   `nrow(newdata)`, a function of the [proxy_cate()] table returning actions,
#'   or one of the strings `"all"`, `"none"`, `"optimal"`.
#' @param value Numeric scalar -- the value of one unit of outcome.
#' @param cost Numeric scalar -- the cost of treating one unit. Default `0`.
#' @param t1,t0 The treated and control treatment values. Default the treatment levels observed at fit time (`model@treatment_levels`); a value matching neither observed level aborts.
#' @param exclude_low_overlap Logical -- drop overlap-flagged units from the
#'   average (and report the count). Default `TRUE`.
#'
#' @returns A one-row [data.table::data.table] with columns `policy_value`,
#'   `n_used`, `n_excluded`, `n_treated`.
#' @family decision
#' @seealso [proxy_decide()], [proxy_overlap()]
#' @export
#' @examples
#' set.seed(1)
#' n <- 600L
#' x <- stats::rnorm(n)
#' t <- stats::rbinom(n, 1L, 0.5)
#' y <- 1 + (0.4 + x) * t + stats::rnorm(n, sd = 0.5)
#' dat <- data.frame(y = y, t = t, x = x)
#' m <- fit_uplift(dat, "y", "t", "x", N = 2L, regime = "sample",
#'                 max_iter = 80L, seed = 1L)
#' nd <- data.frame(x = stats::rnorm(200))
#' proxy_policy_value(m, nd, policy = "optimal", value = 1, cost = 0.3)
proxy_policy_value <- function(model,
                               newdata,
                               policy,
                               value,
                               cost = 0,
                               t1 = NULL,
                               t0 = NULL,
                               exclude_low_overlap = TRUE) {
  if (!S7::S7_inherits(model, uplift_model)) {
    cli::cli_abort("`model` must be an {.cls uplift_model}.")
  }
  if (length(value) != 1L || !is.numeric(value)) {
    cli::cli_abort("`value` must be a numeric scalar.")
  }
  resolved <- .resolve_treatment_arms(model, t1, t0)
  t1 <- resolved$t1
  t0 <- resolved$t0
  g <- model@fit
  z_idx <- c(model@roles$treatment, model@roles$covariate)
  x_idx <- model@roles$covariate
  X <- .newdata_x(model, newdata)
  cache <- .uplift_cache(model)
  arms <- c(t0, t1)
  action <- .resolve_policy(policy, model, newdata, value, cost, t1, t0)

  n_units <- nrow(X)
  contrib <- numeric(n_units)
  keep <- rep(TRUE, n_units)
  for (u in seq_len(n_units)) {
    t_u <- if (action[u] == 1L) t1 else t0
    resp <- .response_mean(cache, g, z_idx, x_idx, X[u, ], t_u, model@assume)
    contrib[u] <- value * resp - cost * action[u]
    if (exclude_low_overlap &&
        .unit_coverage(g, z_idx, arms, X[u, ]) < 0.01) {
      keep[u] <- FALSE
    }
  }
  n_excl <- sum(!keep)
  if (n_excl > 0L) {
    cli::cli_inform("Excluded {n_excl} low-overlap unit{?s} from the policy value.")
  }
  dt <- data.table::data.table(
    policy_value = mean(contrib[keep]),
    n_used = sum(keep),
    n_excluded = n_excl,
    n_treated = sum(action[keep] == 1L)
  )
  return(dt[])
}

## proxy_retrospective_uplift -----------------------------------------------

#' Retrospective (counterfactual-mean) uplift for observed units
#'
#' For each observed unit `(y, t, x)`, the counterfactual-mean uplift of moving
#' from `t0` to `t1`, \eqn{E[Y_{t_1} \mid y, t, x] - E[Y_{t_0} \mid y, t, x]},
#' computed by [gmm_counterfactual()]. Unlike [proxy_cate()], the abduction
#' gate uses the observed outcome `y` as well, sharpening the per-unit estimate.
#' Only the counterfactual *mean* is identified; the spread is not (see
#' [gmm_cf_variance()]).
#'
#' @param model An [uplift_model].
#' @param observed A data frame carrying the outcome, treatment, and covariate
#'   columns of the observed units.
#' @param t1,t0 The treated and control treatment values. Default the treatment levels observed at fit time (`model@treatment_levels`); a value matching neither observed level aborts.
#'
#' @returns A [data.table::data.table] with columns `id`, `y_obs`, `t_obs`,
#'   `cf_mean_t1`, `retro_uplift`.
#' @family decision
#' @seealso [gmm_counterfactual()], [proxy_cate()]
#' @export
#' @examples
#' set.seed(1)
#' n <- 400L
#' x <- stats::rnorm(n)
#' t <- stats::rbinom(n, 1L, 0.5)
#' y <- 1 + (0.5 + x) * t + stats::rnorm(n, sd = 0.5)
#' dat <- data.frame(y = y, t = t, x = x)
#' m <- fit_uplift(dat, "y", "t", "x", N = 2L, regime = "sample",
#'                 max_iter = 80L, seed = 1L)
#' proxy_retrospective_uplift(m, observed = dat[1:5, ])
proxy_retrospective_uplift <- function(model, observed, t1 = NULL, t0 = NULL) {
  if (!S7::S7_inherits(model, uplift_model)) {
    cli::cli_abort("`model` must be an {.cls uplift_model}.")
  }
  resolved <- .resolve_treatment_arms(model, t1, t0)
  t1 <- resolved$t1
  t0 <- resolved$t0
  g <- model@fit
  r <- model@roles
  p <- gmm_dim(g)
  if (!is.data.frame(observed)) {
    cli::cli_abort("`observed` must be a data frame.")
  }
  need <- c(r$outcome_name, r$treatment_name, r$covariate_names)
  if (!all(need %in% names(observed))) {
    miss <- setdiff(need, names(observed))
    cli::cli_abort("`observed` is missing column{?s}: {.val {miss}}.")
  }

  y_obs <- as.numeric(observed[[r$outcome_name]])
  t_obs <- as.numeric(observed[[r$treatment_name]])
  X <- as.matrix(observed[r$covariate_names])
  storage.mode(X) <- "double"

  n_units <- nrow(observed)
  cf1 <- numeric(n_units)
  retro <- numeric(n_units)
  for (u in seq_len(n_units)) {
    evidence <- numeric(p)
    evidence[r$outcome] <- y_obs[u]
    evidence[r$treatment] <- t_obs[u]
    evidence[r$covariate] <- X[u, ]
    do1 <- rep(NA_real_, p)
    do1[r$treatment] <- t1
    do0 <- rep(NA_real_, p)
    do0[r$treatment] <- t0
    m1 <- gmm_counterfactual(g, evidence = evidence, do = do1,
                             query = r$outcome)@mean
    m0 <- gmm_counterfactual(g, evidence = evidence, do = do0,
                             query = r$outcome)@mean
    cf1[u] <- m1
    retro[u] <- m1 - m0
  }
  dt <- data.table::data.table(
    id = seq_len(n_units),
    y_obs = y_obs,
    t_obs = t_obs,
    cf_mean_t1 = cf1,
    retro_uplift = retro
  )
  return(dt[])
}

## proxy_regime_segments ----------------------------------------------------

#' The fitted regimes as an interpretable segment table
#'
#' Exposes the `K` mixture components as decision segments: each regime's
#' prevalence (`weight`), its within-segment treatment effect (the within-class
#' treatment slope), its residual standard deviation, and its covariate centre.
#' This is the interpretable by-product the closed-form reading gives for free.
#'
#' @param model An [uplift_model].
#' @param t1,t0 The treated and control treatment values used to scale the
#'   within-segment effect. Default the treatment levels observed at fit time (`model@treatment_levels`).
#'
#' @returns A [data.table::data.table] with columns `regime`, `weight`,
#'   `effect`, `sigma`, and one column per covariate centre.
#' @family decision
#' @export
#' @examples
#' set.seed(1)
#' n <- 600L
#' x <- stats::rnorm(n)
#' t <- stats::rbinom(n, 1L, 0.5)
#' y <- 1 + (0.5 + x) * t + stats::rnorm(n, sd = 0.5)
#' dat <- data.frame(y = y, t = t, x = x)
#' m <- fit_uplift(dat, "y", "t", "x", N = 2L, regime = "sample",
#'                 max_iter = 80L, seed = 1L)
#' proxy_regime_segments(m)
proxy_regime_segments <- function(model, t1 = NULL, t0 = NULL) {
  if (!S7::S7_inherits(model, uplift_model)) {
    cli::cli_abort("`model` must be an {.cls uplift_model}.")
  }
  resolved <- .resolve_treatment_arms(model, t1, t0)
  t1 <- resolved$t1
  t0 <- resolved$t0
  g <- model@fit
  r <- model@roles
  cache <- .uplift_cache(model)
  K <- length(cache)
  dt <- t1 - t0

  eff <- vapply(cache, function(ck) ck$beta_t * dt, numeric(1L))
  sig <- vapply(cache, function(ck) sqrt(ck$sigma2), numeric(1L))
  centres <- matrix(0, nrow = K, ncol = length(r$covariate))
  for (k in seq_len(K)) {
    centres[k, ] <- g@means[[k]][r$covariate]
  }
  colnames(centres) <- r$covariate_names

  out <- data.table::data.table(
    regime = seq_len(K),
    weight = g@weights,
    effect = eff,
    sigma = sig
  )
  out <- cbind(out, data.table::as.data.table(centres))
  return(out[])
}

## proxy_identification_report ----------------------------------------------

#' The identification report (an executive one-pager)
#'
#' The differentiator: a structured audit of what the decision model identifies,
#' what it assumes, and what it cannot answer. Carries the estimand, the
#' identification regime and its requirement, the overlap rate on the supplied
#' population, the confounding-gap magnitude (the value at risk from unobserved
#' confounding), and the explicit non-identification of the individual
#' counterfactual law.
#'
#' @param model An [uplift_model].
#' @param newdata A data frame carrying the covariate columns -- the population
#'   the report is computed over.
#' @param t1,t0 The treated and control treatment values. Default the treatment levels observed at fit time (`model@treatment_levels`); a value matching neither observed level aborts.
#'
#' @returns An S7 object of class `uplift_identification` with a `print` method.
#' @family decision
#' @seealso [proxy_confounding_gap()], [proxy_overlap()]
#' @export
#' @examples
#' set.seed(1)
#' n <- 600L
#' x <- stats::rnorm(n)
#' t <- stats::rbinom(n, 1L, 0.5)
#' y <- 0.5 * t + x + stats::rnorm(n, sd = 0.5)
#' dat <- data.frame(y = y, t = t, x = x)
#' m <- fit_uplift(dat, "y", "t", "x", N = 2L, regime = "sample",
#'                 max_iter = 80L, seed = 1L)
#' proxy_identification_report(m, data.frame(x = stats::rnorm(100)))
proxy_identification_report <- function(model, newdata, t1 = NULL, t0 = NULL) {
  if (!S7::S7_inherits(model, uplift_model)) {
    cli::cli_abort("`model` must be an {.cls uplift_model}.")
  }
  gap <- proxy_confounding_gap(model, newdata, t1 = t1, t0 = t0)
  ov <- proxy_overlap(model, newdata, t1 = t1, t0 = t0)

  uplift_identification(
    estimand = "CATE / uplift: E[Y | do(T=t1), X] - E[Y | do(T=t0), X]",
    assume = model@assume,
    n_units = nrow(gap),
    overlap_pct = 100 * mean(!ov$overlap_flag),
    confounding_gap_mean = mean(abs(gap$gap)),
    confounding_gap_max = max(abs(gap$gap)),
    K = gmm_n_components(model@fit),
    outcome_type = model@outcome_type
  )
}

#' Identification-report object
#'
#' The structured return type of [proxy_identification_report()]. Print it for
#' the executive one-pager.
#'
#' @param estimand Character -- the target estimand.
#' @param assume Character -- the identification regime.
#' @param n_units Integer -- units the report covers.
#' @param overlap_pct Numeric -- percentage of units with adequate overlap.
#' @param confounding_gap_mean,confounding_gap_max Numeric -- mean and maximum
#'   absolute confounding gap over the population.
#' @param K Integer -- the number of fitted regimes.
#' @param outcome_type Character -- the outcome scale.
#'
#' @returns An S7 object of class `uplift_identification`.
#' @family decision
#' @export
uplift_identification <- S7::new_class(
  name = "uplift_identification",
  package = "proxymix",
  properties = list(
    estimand = S7::class_character,
    assume = S7::class_character,
    n_units = S7::class_integer,
    overlap_pct = S7::class_double,
    confounding_gap_mean = S7::class_double,
    confounding_gap_max = S7::class_double,
    K = S7::class_integer,
    outcome_type = S7::new_property(
      class = S7::class_character,
      default = "continuous"
    )
  )
)

#' @export
S7::method(print, uplift_identification) <- function(x, ...) {
  cat("== Identification report ==================================\n")
  cat(sprintf("  Estimand   : %s\n", x@estimand))
  cat(sprintf("  Assumption : %s\n", x@assume))
  if (identical(x@assume, "ignorability")) {
    cat("               requires (Y(0), Y(1)) independent of T given X.\n")
  } else {
    cat("               requires the fitted regime to be the only confounder\n")
    cat("               of T and Y -- NOT certified by the fit.\n")
  }
  cat(sprintf("  Regimes    : K = %d   Outcome scale: %s\n",
              x@K, x@outcome_type))
  cat(sprintf("  Units      : %d\n", x@n_units))
  cat(sprintf("  Overlap    : %.1f%% of units adequately supported\n",
              x@overlap_pct))
  cat(sprintf("  Confounding gap (value at risk if a latent regime confounds):\n"))
  cat(sprintf("               mean |Delta| = %.4g, max |Delta| = %.4g\n",
              x@confounding_gap_mean, x@confounding_gap_max))
  cat("  NOT identified : the individual counterfactual law\n")
  cat("                   (its variance and tail probabilities).\n")
  cat("===========================================================\n")
  invisible(x)
}
