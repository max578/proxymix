## uplift.R -- The decision model: one joint fit read as an uplift engine.
##
## `fit_uplift()` assembles a joint Gaussian-mixture proxy over
## (outcome, treatment, covariates) and wraps it, with the column roles and
## the chosen identification assumption, in an `uplift_model`. Every decision
## verb (`proxy_cate()`, `proxy_decide()`, ...) reads that single fit in closed
## form. Nothing here re-fits per query -- the fit is the model.

## uplift_model class -------------------------------------------------------

#' A fitted uplift / next-best-action model
#'
#' The object returned by [fit_uplift()]: a joint [gmm_fit] over
#' `(outcome, treatment, covariates)` together with the column roles, the
#' identification assumption it will be read under, the outcome type, and the
#' training sample (retained so that resampling standard errors and overlap
#' diagnostics are available). The decision verbs dispatch on this class.
#'
#' @param fit The joint [gmm_fit] over the stacked `(outcome, treatment,
#'   covariates)` coordinates, in that column order.
#' @param roles A list with integer indices `outcome`, `treatment`,
#'   `covariate` and the matching `outcome_name`, `treatment_name`,
#'   `covariate_names`.
#' @param assume One of `"ignorability"` (the default uplift assumption) or
#'   `"latent_confounder"` (the do-operator reading; always flagged).
#' @param outcome_type One of `"continuous"`, `"binary"`, `"count"`.
#' @param data The `n` by `p` numeric training matrix `cbind(outcome,
#'   treatment, covariates)`.
#' @param n_train Integer scalar -- the training sample size.
#' @param treatment_levels Numeric length-2 vector `c(t0, t1)` -- the control
#'   and treated values used at fit time (the observed treatment arms).
#' @param name Human-readable name.
#' @param metadata Optional list of descriptors (e.g. the K-selection trace).
#'
#' @returns An S7 object of class `uplift_model`.
#' @family decision
#' @export
uplift_model <- S7::new_class(
  name = "uplift_model",
  package = "proxymix",
  properties = list(
    fit = S7::class_any,
    roles = S7::class_list,
    assume = S7::new_property(
      class = S7::class_character,
      default = "ignorability"
    ),
    outcome_type = S7::new_property(
      class = S7::class_character,
      default = "continuous"
    ),
    data = S7::class_any,
    n_train = S7::class_integer,
    treatment_levels = S7::new_property(
      class = S7::class_double,
      default = c(0, 1)
    ),
    name = S7::new_property(
      class = S7::class_character,
      default = "uplift_model"
    ),
    metadata = S7::new_property(
      class = S7::class_list,
      default = list()
    )
  ),
  validator = function(self) {
    if (!self@assume %in% c("ignorability", "latent_confounder")) {
      return("`assume` must be \"ignorability\" or \"latent_confounder\"")
    }
    if (!self@outcome_type %in% c("continuous", "binary", "count")) {
      return("`outcome_type` must be \"continuous\", \"binary\" or \"count\"")
    }
    NULL
  }
)

#' @export
S7::method(print, uplift_model) <- function(x, ...) {
  r <- x@roles
  cat(sprintf("<uplift_model>: K = %d regimes, assume = \"%s\"\n",
              gmm_n_components(x@fit), x@assume))
  cat(sprintf("  outcome   : %s (%s)\n", r$outcome_name, x@outcome_type))
  cat(sprintf("  treatment : %s  (arms %g vs %g)\n",
              r$treatment_name, x@treatment_levels[1L], x@treatment_levels[2L]))
  cat(sprintf("  covariates: %s\n",
              paste(r$covariate_names, collapse = ", ")))
  cat(sprintf("  trained on: %d units\n", x@n_train))
  if (identical(x@assume, "latent_confounder")) {
    cat("  NOTE: do-operator reading -- effects assume the fitted regime\n")
    cat("        is the only confounder; not certified by the fit.\n")
  }
  invisible(x)
}

## Validation ---------------------------------------------------------------

.check_uplift_roles <- function(data, outcome, treatment, covariates) {
  if (!is.data.frame(data)) {
    cli::cli_abort("`data` must be a data frame.")
  }
  nm <- names(data)
  if (length(outcome) != 1L || !is.character(outcome) || !outcome %in% nm) {
    cli::cli_abort("`outcome` must name a single column of `data`.")
  }
  if (length(treatment) != 1L || !is.character(treatment) || !treatment %in% nm) {
    cli::cli_abort("`treatment` must name a single column of `data`.")
  }
  if (!is.character(covariates) || length(covariates) < 1L ||
      !all(covariates %in% nm)) {
    cli::cli_abort("`covariates` must name one or more columns of `data`.")
  }
  if (treatment %in% covariates || outcome %in% covariates ||
      outcome == treatment) {
    cli::cli_abort("`outcome`, `treatment` and `covariates` must be distinct columns.")
  }
  invisible(TRUE)
}

## Bayesian-information-criterion sweep for the component count -------------

## A component holding a discrete column constant has a density limited only by
## the ridge, so scores floor such a column's variance at spacing^2 / 12.
.floored_loglik <- function(g, M, floor_var) {
  scored <- gmm(
    weights = g@weights,
    means = g@means,
    covariances = lapply(g@covariances, function(S) S + diag(floor_var, ncol(M)))
  )
  sum(dgmm(M, scored, log = TRUE))
}

## Half k-means, half random: on Project STAR only k-means reaches the best optimum,
## and on a mixed continuous and binary joint only random starts do.
.uplift_starts <- function(samples, N, n_starts, ridge_eps) {
  n_km <- ceiling(n_starts / 2)
  km <- lapply(seq_len(n_km), function(i) {
    tryCatch(init_kmeans(samples, N = N, ridge_eps = ridge_eps, nstart = 1L),
             error = function(e) NULL)
  })
  mu <- colMeans(samples)
  sds <- sqrt(diag(stats::cov(samples)))
  rnd <- lapply(seq_len(n_starts - n_km), function(i) {
    init_random(N = N, p = ncol(samples), centre = mu, scale = mean(sds),
                sigma_diag = mean(sds^2))
  })
  c(km[!vapply(km, is.null, logical(1L))], rnd)
}

## EM from each start, keeping the fit with the best floored score.
.fit_em_floored <- function(target, N, floor_var, max_iter = 500L, tol = 1e-8,
                            ridge_eps = 1e-6, n_starts = NULL, seed = NULL,
                            canonicalise = TRUE) {
  samples <- target@samples
  if (is.null(n_starts)) n_starts <- if (any(floor_var[-2L] > 0)) 100L else 20L
  if (N == 1L) n_starts <- 1L
  inits <- if (is.null(seed)) {
    .uplift_starts(samples, N, n_starts, ridge_eps)
  } else {
    withr::with_seed(as.integer(seed), .uplift_starts(samples, N, n_starts, ridge_eps))
  }
  call <- match.call()
  fits <- lapply(inits, function(init) {
    tryCatch(em_samples_one_run(samples, init, target, max_iter, tol, ridge_eps,
                                call = call),
             error = function(e) NULL)
  })
  fits <- fits[!vapply(fits, is.null, logical(1L))]
  if (!length(fits)) cli::cli_abort("all EM starts failed.")
  scores <- vapply(fits, .floored_loglik, numeric(1L), M = samples,
                   floor_var = floor_var)
  .finish_em_samples(fits[[which.max(scores)]], samples, ridge_eps, canonicalise)
}

.fit_uplift_k <- function(target, k, regime, floor_var, ...) {
  args <- list(...)
  uses_em <- (regime == "sample" || (regime == "auto" && k > 1L)) &&
    is.null(args$init) && !isTRUE(args$anneal)
  if (!uses_em) {
    return(do.call(fit_proxymix, c(list(target = target, N = k, regime = regime), args)))
  }
  args$anneal <- NULL
  args$temp_schedule <- NULL
  do.call(.fit_em_floored, c(list(target = target, N = k, floor_var = floor_var), args))
}

## Every K starts from the same random state, so the optimum reached at K does not
## depend on the earlier K.
.select_k_bic <- function(target, M, n_grid, regime, floor_var, ...) {
  n <- nrow(M)
  p <- ncol(M)
  fit_seed <- sample.int(.Machine$integer.max, 1L)
  fits <- vector("list", length(n_grid))
  bic <- numeric(length(n_grid))
  for (i in seq_along(n_grid)) {
    k <- n_grid[i]
    fits[[i]] <- withr::with_seed(
      fit_seed, .fit_uplift_k(target, k, regime, floor_var, ...)
    )
    n_par <- k * (p + p * (p + 1L) / 2L) + (k - 1L)
    bic[i] <- -2 * .floored_loglik(fits[[i]], M, floor_var) + n_par * log(n)
  }
  best <- which.min(bic)
  list(fit = fits[[best]], K = n_grid[best],
       bic = stats::setNames(bic, n_grid))
}

## fit_uplift ---------------------------------------------------------------

#' Fit an uplift / next-best-action model from a data frame
#'
#' Assembles a joint Gaussian-mixture proxy over the outcome, the treatment,
#' and the covariates, and returns an [uplift_model] that the decision verbs
#' read in closed form. One fit yields prediction, heterogeneous treatment
#' effects, optimal actions, off-line policy value, and an identification
#' audit -- see [proxy_cate()], [proxy_decide()], [proxy_policy_value()] and
#' [proxy_identification_report()].
#'
#' The component count `N` may be fixed or chosen automatically. With
#' `N = "auto"` the function sweeps `n_grid` and keeps the K that minimises the
#' joint BIC; the full BIC trace is stored in the model's metadata. A binary
#' or other discrete column can be held constant inside a component, where its
#' density is limited only by the ridge, so the BIC floors the variance of a
#' column with at most ten distinct values at `spacing^2 / 12`. The same score
#' chooses among the EM starts, half from k-means and half random: 100 starts
#' per fit when the outcome or a covariate is discrete and 20 otherwise. Fits at
#' different `N` start from the same random state, so `N = "auto"` and the
#' same fixed `N` agree at one `seed`.
#'
#' A discrete covariate raises a warning of class `proxymix_discrete_columns`.
#' In a small sample, or with many tied values, different seeds can still reach
#' optima with different effects, even in sign. A mixture
#' can also smooth over a difference in effect between the values of such a column.
#' Compare several seeds, and a regression with the interaction, before reading an
#' effect. The treatment is binary at this version (`{t0, t1}`); a continuous dose
#' is a future extension.
#'
#' @param data A data frame holding the outcome, treatment and covariate
#'   columns.
#' @param outcome A single column name -- the outcome `Y`.
#' @param treatment A single column name -- the binary treatment `T`.
#' @param covariates A character vector of one or more column names -- the
#'   pre-treatment covariates `X`.
#' @param N Number of mixture components, or `"auto"` (the default) to select
#'   by BIC over `n_grid`.
#' @param regime One of `"auto"`, `"moment"`, `"sample"`, `"kld"`, forwarded to
#'   [fit_proxymix()]. The default `"auto"` uses classical EM on the supplied
#'   rows.
#' @param assume One of `"ignorability"` (the default) or
#'   `"latent_confounder"` -- the identification regime the effects are read
#'   under. See [proxy_cate()] and [proxy_confounding_gap()].
#' @param outcome_type One of `"continuous"` (the default), `"binary"` or
#'   `"count"`. Effects are reported on the response scale via a discretised
#'   predictive for the non-continuous types; see [proxy_cate()].
#' @param n_grid Integer vector of candidate component counts used when
#'   `N = "auto"`. Default `1:4`.
#' @param seed Optional integer. When supplied, the fitting (including any
#'   random EM starts) runs under a fixed seed and the global RNG state is
#'   restored on exit, so the fit is reproducible without disturbing the
#'   caller's stream.
#' @param ... Additional arguments forwarded to [fit_proxymix()] (e.g.
#'   `max_iter`, `n_starts`, `ridge_eps`).
#'
#' @returns An [uplift_model].
#' @family decision
#' @seealso [proxy_cate()], [proxy_decide()], [proxy_identification_report()]
#' @export
#' @examples
#' set.seed(1)
#' n <- 400L
#' x <- stats::rnorm(n)
#' t <- stats::rbinom(n, 1L, 0.5)
#' y <- 1 + 0.5 * x + (1 + x) * t + stats::rnorm(n, sd = 0.5)
#' dat <- data.frame(y = y, t = t, x = x)
#' m <- fit_uplift(dat, outcome = "y", treatment = "t", covariates = "x",
#'                 N = 2L, regime = "sample", max_iter = 50L, seed = 1L)
#' m
fit_uplift <- function(data,
                       outcome,
                       treatment,
                       covariates,
                       N = "auto",
                       regime = "auto",
                       assume = c("ignorability", "latent_confounder"),
                       outcome_type = c("continuous", "binary", "count"),
                       n_grid = 1:4,
                       seed = NULL,
                       ...) {
  .check_uplift_roles(data, outcome, treatment, covariates)
  assume <- rlang::arg_match(assume)
  outcome_type <- rlang::arg_match(outcome_type)

  if (!is.null(seed)) {
    if (exists(".Random.seed", envir = .GlobalEnv)) {
      old_seed <- get(".Random.seed", envir = .GlobalEnv)
      on.exit(assign(".Random.seed", old_seed, envir = .GlobalEnv), add = TRUE)
    }
    set.seed(as.integer(seed))
  }

  # Assemble the joint matrix in (outcome, treatment, covariates) order --------
  y <- as.numeric(data[[outcome]])
  t <- as.numeric(data[[treatment]])
  X <- as.matrix(data[covariates])
  storage.mode(X) <- "double"
  M <- cbind(y, t, X)
  colnames(M) <- c(outcome, treatment, covariates)
  if (anyNA(M)) {
    cli::cli_abort("`data` has missing values in the model columns; drop or impute first.")
  }

  t_levels <- sort(unique(t))
  if (length(t_levels) != 2L) {
    cli::cli_warn(c(
      "Treatment is not binary (found {length(t_levels)} level{?s}).",
      "i" = "This version targets a binary treatment; effects use t0 / t1 = the two extreme observed values."
    ))
    t_levels <- range(t)
  }

  roles <- list(
    outcome = 1L,
    treatment = 2L,
    covariate = seq.int(3L, 2L + length(covariates)),
    outcome_name = outcome,
    treatment_name = treatment,
    covariate_names = covariates
  )

  target <- gmm_target_from_samples(M, name = "uplift_joint")

  .warn_discrete_columns(
    M[, -(1:2), drop = FALSE], covariates,
    advice = paste(
      "Different seeds can still reach optima with different effects, and a",
      "mixture can smooth over a difference in effect between these values.",
      "Compare several seeds and a regression with the interaction."
    )
  )
  floor_var <- .discrete_variance_floor(M)

  # Fit (fixed K or BIC-selected) ---------------------------------------------
  if (identical(N, "auto")) {
    sel <- .select_k_bic(target, M, n_grid = as.integer(n_grid),
                         regime = regime, floor_var = floor_var, ...)
    fit <- sel$fit
    K <- sel$K
    bic_trace <- sel$bic
    cli::cli_inform("Selected K = {K} by BIC.")
  } else {
    K <- as.integer(N)
    fit <- withr::with_seed(
      sample.int(.Machine$integer.max, 1L),
      .fit_uplift_k(target, K, regime, floor_var, ...)
    )
    bic_trace <- NULL
  }

  uplift_model(
    fit = fit,
    roles = roles,
    assume = assume,
    outcome_type = outcome_type,
    data = M,
    n_train = nrow(M),
    treatment_levels = c(t_levels[1L], t_levels[2L]),
    name = "uplift_model",
    metadata = list(bic_trace = bic_trace)
  )
}
