## decide.R -- The scoring verbs: CATE, uplift, optimal action, overlap.
##
## All four read the fitted joint mixture in closed form, O(K) per unit. The
## point estimates use two identification modes (ignorability vs the
## latent-confounder do-operator); the default delta-method standard error is
## the within-component prediction variance, which reduces to the ordinary
## least-squares standard error at K = 1. A resampling alternative
## (`se_method = "mc"`) captures the regime-gate uncertainty the delta method
## holds fixed.

## Per-component cache -------------------------------------------------------

## Pre-computes, once per model, every per-component quantity the scoring loop
## needs: the (treatment, covariate) mean / precision, the within-class
## regression coefficients of the outcome, the residual variance, and the
## effective sample size. The treatment sits at position 1 of the Z block.
.uplift_cache <- function(model) {
  g <- model@fit
  r <- model@roles
  y_idx <- r$outcome
  z_idx <- c(r$treatment, r$covariate)
  x_local <- seq.int(2L, length(z_idx))
  K <- gmm_n_components(g)
  n <- model@n_train

  lapply(seq_len(K), function(k) {
    mu <- g@means[[k]]
    S <- g@covariances[[k]]
    mZ <- mu[z_idx]
    Szz <- S[z_idx, z_idx, drop = FALSE]
    Pz <- chol2inv(chol(Szz))
    Szy <- S[z_idx, y_idx, drop = FALSE]
    beta <- as.numeric(Pz %*% Szy)
    sigma2 <- as.numeric(S[y_idx, y_idx] - crossprod(Szy, beta))
    list(
      mY = mu[y_idx],
      mZ = mZ,
      Pz = Pz,
      beta = beta,
      beta_t = beta[1L],
      var_t = as.numeric(Pz[1L, 1L]),
      sigma2 = max(sigma2, .Machine$double.eps),
      n_eff = max(g@weights[k] * n, 1),
      x_local = x_local
    )
  })
}

## Per-unit CATE + variance under either identification mode.
.cate_unit <- function(cache, g, z_idx, x_idx, x, t1, t0, mode) {
  z1 <- c(t1, x)
  z0 <- c(t0, x)
  K <- length(cache)
  dt <- t1 - t0

  if (identical(mode, "latent_confounder")) {
    # do-operator: gate on X only, contrast on the within-class slope ----------
    r_x <- .responsibilities_safe(g, x_idx, x)
    tau <- 0
    v <- 0
    for (k in seq_len(K)) {
      ck <- cache[[k]]
      tau <- tau + r_x[k] * ck$beta_t * dt
      v <- v + (r_x[k] * dt)^2 * ck$sigma2 * ck$var_t / ck$n_eff
    }
    return(list(tau = tau, var = v, covered = attr(r_x, "covered")))
  }

  # ignorability: gate on (T, X), contrast two conditional means ---------------
  r1 <- .responsibilities_safe(g, z_idx, z1)
  r0 <- .responsibilities_safe(g, z_idx, z0)
  m1 <- 0
  m0 <- 0
  v <- 0
  for (k in seq_len(K)) {
    ck <- cache[[k]]
    d1 <- z1 - ck$mZ
    d0 <- z0 - ck$mZ
    mu1 <- ck$mY + sum(ck$beta * d1)
    mu0 <- ck$mY + sum(ck$beta * d0)
    m1 <- m1 + r1[k] * mu1
    m0 <- m0 + r0[k] * mu0
    h11 <- 1 / ck$n_eff + as.numeric(crossprod(d1, ck$Pz %*% d1)) / ck$n_eff
    h00 <- 1 / ck$n_eff + as.numeric(crossprod(d0, ck$Pz %*% d0)) / ck$n_eff
    h10 <- 1 / ck$n_eff + as.numeric(crossprod(d1, ck$Pz %*% d0)) / ck$n_eff
    v <- v + ck$sigma2 *
      (r1[k]^2 * h11 - 2 * r1[k] * r0[k] * h10 + r0[k]^2 * h00)
  }
  list(tau = m1 - m0, var = max(v, 0),
       covered = attr(r1, "covered") && attr(r0, "covered"))
}

## Per-unit interventional response level (not a contrast).
.response_mean <- function(cache, g, z_idx, x_idx, x, t, mode) {
  z <- c(t, x)
  K <- length(cache)
  r <- if (identical(mode, "latent_confounder")) {
    .responsibilities_safe(g, x_idx, x)
  } else {
    .responsibilities_safe(g, z_idx, z)
  }
  mu <- 0
  for (k in seq_len(K)) {
    ck <- cache[[k]]
    mu <- mu + r[k] * (ck$mY + sum(ck$beta * (z - ck$mZ)))
  }
  mu
}

## Per-unit response-scale class probability P(Y > threshold | do/see T=t, X=x)
## under the fitted conditional Gaussian mixture -- the discretised predictive
## used for binary outcomes (decision 4: latent scale + discretised scoring).
.response_prob <- function(cache, g, z_idx, x_idx, x, t, mode, threshold) {
  z <- c(t, x)
  K <- length(cache)
  r <- if (identical(mode, "latent_confounder")) {
    .responsibilities_safe(g, x_idx, x)
  } else {
    .responsibilities_safe(g, z_idx, z)
  }
  p <- 0
  for (k in seq_len(K)) {
    ck <- cache[[k]]
    mu <- ck$mY + sum(ck$beta * (z - ck$mZ))
    p <- p + r[k] * stats::pnorm((mu - threshold) / sqrt(ck$sigma2))
  }
  p
}

## Per-unit CATE point estimate on the requested scale (link or response). The
## response scale differs from the link only for a binary outcome, where the
## contrast is on the discretised predictive probability.
.tau_point <- function(cache, g, z_idx, x_idx, x, t1, t0, mode,
                       scale, threshold, outcome_type) {
  if (identical(scale, "response") && identical(outcome_type, "binary")) {
    p1 <- .response_prob(cache, g, z_idx, x_idx, x, t1, mode, threshold)
    p0 <- .response_prob(cache, g, z_idx, x_idx, x, t0, mode, threshold)
    return(p1 - p0)
  }
  .cate_unit(cache, g, z_idx, x_idx, x, t1, t0, mode)$tau
}

## Vectorised CATE + variance over a matrix of units. This is the serving hot
## path: responsibilities, conditional means and the delta variance are all
## computed in matrix form (O(K) passes, no per-unit R loop), realising the
## closed-form scoring edge.
.cate_batch <- function(cache, g, z_idx, x_idx, X, t1, t0, mode) {
  n <- nrow(X)
  K <- length(cache)
  dt <- t1 - t0

  if (identical(mode, "latent_confounder")) {
    rx <- .responsibilities_batch(g, x_idx, X)
    beta_t <- vapply(cache, function(ck) ck$beta_t, numeric(1L))
    vk <- vapply(cache, function(ck) ck$sigma2 * ck$var_t / ck$n_eff,
                 numeric(1L))
    tau <- as.numeric(rx %*% beta_t) * dt
    v <- as.numeric((rx^2) %*% vk) * dt^2
    return(list(tau = tau, var = v, covered = attr(rx, "covered")))
  }

  Z1 <- cbind(t1, X)
  Z0 <- cbind(t0, X)
  r1 <- .responsibilities_batch(g, z_idx, Z1)
  r0 <- .responsibilities_batch(g, z_idx, Z0)
  m1 <- numeric(n)
  m0 <- numeric(n)
  v <- numeric(n)
  for (k in seq_len(K)) {
    ck <- cache[[k]]
    D1 <- sweep(Z1, 2L, ck$mZ)
    D0 <- sweep(Z0, 2L, ck$mZ)
    m1 <- m1 + r1[, k] * (ck$mY + as.numeric(D1 %*% ck$beta))
    m0 <- m0 + r0[, k] * (ck$mY + as.numeric(D0 %*% ck$beta))
    D1P <- D1 %*% ck$Pz
    D0P <- D0 %*% ck$Pz
    h11 <- (1 + rowSums(D1P * D1)) / ck$n_eff
    h00 <- (1 + rowSums(D0P * D0)) / ck$n_eff
    h10 <- (1 + rowSums(D1P * D0)) / ck$n_eff
    v <- v + ck$sigma2 *
      (r1[, k]^2 * h11 - 2 * r1[, k] * r0[, k] * h10 + r0[, k]^2 * h00)
  }
  list(tau = m1 - m0, var = pmax(v, 0),
       covered = attr(r1, "covered") & attr(r0, "covered"))
}

## Vectorised response-scale class probability over a matrix of units.
.response_prob_batch <- function(cache, g, z_idx, x_idx, X, t, mode, threshold) {
  n <- nrow(X)
  K <- length(cache)
  Z <- cbind(t, X)
  r <- if (identical(mode, "latent_confounder")) {
    .responsibilities_batch(g, x_idx, X)
  } else {
    .responsibilities_batch(g, z_idx, Z)
  }
  p <- numeric(n)
  for (k in seq_len(K)) {
    ck <- cache[[k]]
    mu <- ck$mY + as.numeric(sweep(Z, 2L, ck$mZ) %*% ck$beta)
    p <- p + r[, k] * stats::pnorm((mu - threshold) / sqrt(ck$sigma2))
  }
  p
}

## Resampling (mc) standard errors ------------------------------------------

## Non-parametric bootstrap of the training matrix: refit B times and take the
## empirical SD of tau(x) per unit. Captures the regime-gate uncertainty the
## delta method holds fixed; the documented "sampling-based draws" option.
.cate_mc_se <- function(model, x_mat, t1, t0, B, scale, threshold, ...) {
  M <- model@data
  n <- nrow(M)
  K <- gmm_n_components(model@fit)
  regime <- model@fit@regime
  z_idx <- c(model@roles$treatment, model@roles$covariate)
  x_idx <- model@roles$covariate
  n_units <- nrow(x_mat)
  draws <- matrix(NA_real_, nrow = B, ncol = n_units)
  for (b in seq_len(B)) {
    idx <- sample.int(n, n, replace = TRUE)
    tgt_b <- gmm_target_from_samples(M[idx, , drop = FALSE], name = "boot")
    fit_b <- tryCatch(
      fit_proxymix(tgt_b, N = K, regime = regime, ...),
      error = function(e) NULL
    )
    if (is.null(fit_b)) next
    cache_b <- .uplift_cache(
      uplift_model(fit = fit_b, roles = model@roles, assume = model@assume,
                   outcome_type = model@outcome_type, data = M[idx, , drop = FALSE],
                   n_train = n, treatment_levels = model@treatment_levels)
    )
    for (u in seq_len(n_units)) {
      draws[b, u] <- .tau_point(cache_b, fit_b, z_idx, x_idx, x_mat[u, ],
                                t1, t0, model@assume, scale, threshold,
                                model@outcome_type)
    }
  }
  apply(draws, 2L, stats::sd, na.rm = TRUE)
}

## Covariate-matrix extraction ----------------------------------------------

.newdata_x <- function(model, newdata) {
  cov_names <- model@roles$covariate_names
  if (!is.data.frame(newdata)) {
    cli::cli_abort("`newdata` must be a data frame.")
  }
  if (!all(cov_names %in% names(newdata))) {
    miss <- setdiff(cov_names, names(newdata))
    cli::cli_abort("`newdata` is missing covariate column{?s}: {.val {miss}}.")
  }
  X <- as.matrix(newdata[cov_names])
  storage.mode(X) <- "double"
  X
}

## proxy_cate ---------------------------------------------------------------

#' Heterogeneous treatment effects (CATE / uplift)
#'
#' Per-unit conditional average treatment effect \eqn{\tau(x) = E[Y \mid
#' do(T = t_1), X = x] - E[Y \mid do(T = t_0), X = x]}, read in closed form off
#' the fitted joint mixture. Under the model's default `"ignorability"`
#' assumption this is the contrast of two component-gated conditional means;
#' under `"latent_confounder"` it is the regime-gated within-class slope (the
#' do-operator). The two coincide when treatment carries no information about
#' the regime beyond `X`; their difference is [proxy_confounding_gap()].
#'
#' The default delta-method standard error is the within-component prediction
#' variance, holding the regime gate fixed; it reduces to the ordinary
#' least-squares standard error of the treatment effect at K = 1. Set
#' `se_method = "mc"` for a resampling standard error that also reflects gate
#' uncertainty.
#'
#' @param model An [uplift_model].
#' @param newdata A data frame carrying the covariate columns.
#' @param t1,t0 The treated and control treatment values. Default `1` and `0`.
#' @param se Logical -- compute standard errors and confidence intervals.
#' @param se_method One of `"delta"` (closed form, the default) or `"mc"`
#'   (resampling).
#' @param level Confidence level for the interval. Default `0.95`.
#' @param B Number of bootstrap refits when `se_method = "mc"`. Default `200`.
#' @param scale One of `"link"` (the latent / continuous scale, the default) or
#'   `"response"`. For a binary outcome the response scale reports the effect on
#'   the discretised predictive probability `P(Y > threshold)`; for continuous
#'   and count outcomes the two scales coincide.
#' @param threshold Decision threshold for the binary discretised predictive.
#'   Default `0.5`.
#' @param ... Forwarded to [fit_proxymix()] inside the `"mc"` refits.
#'
#' @returns A [data.table::data.table] with columns `id`, `tau`, `se`,
#'   `ci_lo`, `ci_hi`, `overlap_flag`.
#' @family decision
#' @seealso [proxy_decide()], [proxy_confounding_gap()], [proxy_overlap()]
#' @export
#' @examples
#' set.seed(1)
#' n <- 400L
#' x <- stats::rnorm(n)
#' t <- stats::rbinom(n, 1L, 0.5)
#' y <- 1 + 0.5 * x + (1 + x) * t + stats::rnorm(n, sd = 0.5)
#' dat <- data.frame(y = y, t = t, x = x)
#' m <- fit_uplift(dat, "y", "t", "x", N = 2L, regime = "sample",
#'                 max_iter = 50L, seed = 1L)
#' proxy_cate(m, newdata = data.frame(x = c(-1, 0, 1)))
proxy_cate <- function(model,
                       newdata,
                       t1 = 1,
                       t0 = 0,
                       se = TRUE,
                       se_method = c("delta", "mc"),
                       level = 0.95,
                       B = 200L,
                       scale = c("link", "response"),
                       threshold = 0.5,
                       ...) {
  if (!S7::S7_inherits(model, uplift_model)) {
    cli::cli_abort("`model` must be an {.cls uplift_model}.")
  }
  se_method <- rlang::arg_match(se_method)
  scale <- rlang::arg_match(scale)
  response_binary <- identical(scale, "response") &&
    identical(model@outcome_type, "binary")
  g <- model@fit
  z_idx <- c(model@roles$treatment, model@roles$covariate)
  x_idx <- model@roles$covariate
  X <- .newdata_x(model, newdata)
  cache <- .uplift_cache(model)
  arms <- c(t0, t1)

  n_units <- nrow(X)
  res <- .cate_batch(cache, g, z_idx, x_idx, X, t1, t0, model@assume)
  cov_vec <- .coverage_batch(g, z_idx, arms, X)
  flag <- (!res$covered) | (cov_vec < 0.01)
  if (response_binary) {
    tau <- .response_prob_batch(cache, g, z_idx, x_idx, X, t1,
                                model@assume, threshold) -
      .response_prob_batch(cache, g, z_idx, x_idx, X, t0,
                           model@assume, threshold)
    v <- rep(NA_real_, n_units)
  } else {
    tau <- res$tau
    v <- res$var
  }

  se_vec <- rep(NA_real_, n_units)
  ci_lo <- rep(NA_real_, n_units)
  ci_hi <- rep(NA_real_, n_units)
  if (isTRUE(se)) {
    if (identical(se_method, "delta") && response_binary) {
      cli::cli_inform(c(
        "Response-scale standard errors are not available from the delta method.",
        "i" = "Re-run with `se_method = \"mc\"` for response-scale intervals."
      ))
    }
    se_vec <- if (identical(se_method, "mc")) {
      .cate_mc_se(model, X, t1, t0, B = as.integer(B),
                  scale = scale, threshold = threshold, ...)
    } else if (response_binary) {
      rep(NA_real_, n_units)
    } else {
      sqrt(v)
    }
    zq <- stats::qnorm(1 - (1 - level) / 2)
    ci_lo <- tau - zq * se_vec
    ci_hi <- tau + zq * se_vec
  }

  dt <- data.table::data.table(
    id = seq_len(n_units),
    tau = tau,
    se = se_vec,
    ci_lo = ci_lo,
    ci_hi = ci_hi,
    overlap_flag = flag
  )
  return(dt[])
}

#' Uplift (alias of [proxy_cate()] for a binary treatment)
#'
#' For a binary treatment, the uplift is exactly the conditional average
#' treatment effect. This is a thin alias of [proxy_cate()] kept for the
#' next-best-action vocabulary.
#'
#' @inheritParams proxy_cate
#'
#' @returns A [data.table::data.table] -- see [proxy_cate()].
#' @family decision
#' @export
#' @examples
#' set.seed(1)
#' dat <- data.frame(y = stats::rnorm(200), t = stats::rbinom(200, 1L, 0.5),
#'                   x = stats::rnorm(200))
#' m <- fit_uplift(dat, "y", "t", "x", N = 1L, regime = "moment")
#' proxy_uplift(m, newdata = data.frame(x = 0))
proxy_uplift <- function(model, newdata, ...) {
  proxy_cate(model, newdata, ...)
}

## proxy_predict ------------------------------------------------------------

#' Predicted outcome under a treatment (the seeing rung)
#'
#' Per-unit predicted outcome \eqn{E[Y \mid do(T = t), X = x]} -- the first rung
#' of the ladder, risk / response scoring. Under `"ignorability"` this is the
#' component-gated conditional mean; under `"latent_confounder"` it is the
#' regime-gated interventional mean. For a binary outcome with
#' `scale = "response"` the prediction is the discretised predictive
#' probability `P(Y > threshold)`.
#'
#' @param model An [uplift_model].
#' @param newdata A data frame carrying the covariate columns.
#' @param t The treatment value to predict the outcome under.
#' @param scale One of `"link"` (default) or `"response"`.
#' @param threshold Decision threshold for the binary discretised predictive.
#'   Default `0.5`.
#'
#' @returns A [data.table::data.table] with columns `id` and `prediction`.
#' @family decision
#' @seealso [proxy_cate()]
#' @export
#' @examples
#' set.seed(1)
#' dat <- data.frame(y = stats::rnorm(200), t = stats::rbinom(200, 1L, 0.5),
#'                   x = stats::rnorm(200))
#' m <- fit_uplift(dat, "y", "t", "x", N = 1L, regime = "moment")
#' proxy_predict(m, data.frame(x = c(-1, 0, 1)), t = 1)
proxy_predict <- function(model, newdata, t,
                          scale = c("link", "response"), threshold = 0.5) {
  if (!S7::S7_inherits(model, uplift_model)) {
    cli::cli_abort("`model` must be an {.cls uplift_model}.")
  }
  if (length(t) != 1L || !is.numeric(t)) {
    cli::cli_abort("`t` must be a single treatment value.")
  }
  scale <- rlang::arg_match(scale)
  g <- model@fit
  z_idx <- c(model@roles$treatment, model@roles$covariate)
  x_idx <- model@roles$covariate
  X <- .newdata_x(model, newdata)
  cache <- .uplift_cache(model)
  response_binary <- identical(scale, "response") &&
    identical(model@outcome_type, "binary")

  pred <- vapply(seq_len(nrow(X)), function(u) {
    if (response_binary) {
      .response_prob(cache, g, z_idx, x_idx, X[u, ], t, model@assume, threshold)
    } else {
      .response_mean(cache, g, z_idx, x_idx, X[u, ], t, model@assume)
    }
  }, numeric(1L))
  dt <- data.table::data.table(id = seq_len(nrow(X)), prediction = pred)
  return(dt[])
}

## proxy_overlap ------------------------------------------------------------

#' Per-unit overlap / positivity diagnostic
#'
#' Flags units whose `(treatment, covariate)` configuration is poorly covered
#' by the fitted joint -- the proxy's mass coverage is the positivity
#' diagnostic. For each treatment arm the squared Mahalanobis distance to the
#' nearest regime centre is converted to an upper-tail chi-square coverage
#' probability; the reported `coverage` is the minimum across arms, since the
#' treatment effect needs both arms supported. Units below `floor` are flagged
#' and excluded from [proxy_policy_value()] by default.
#'
#' @param model An [uplift_model].
#' @param newdata A data frame carrying the covariate columns.
#' @param t1,t0 The treated and control treatment values. Default `1` and `0`.
#' @param floor Coverage probability below which a unit is flagged. Default
#'   `0.01`.
#'
#' @returns A [data.table::data.table] with columns `id`, `coverage`,
#'   `overlap_flag`.
#' @family decision
#' @export
#' @examples
#' set.seed(1)
#' dat <- data.frame(y = stats::rnorm(200), t = stats::rbinom(200, 1L, 0.5),
#'                   x = stats::rnorm(200))
#' m <- fit_uplift(dat, "y", "t", "x", N = 1L, regime = "moment")
#' proxy_overlap(m, newdata = data.frame(x = c(0, 8)))
proxy_overlap <- function(model, newdata, t1 = 1, t0 = 0, floor = 0.01) {
  if (!S7::S7_inherits(model, uplift_model)) {
    cli::cli_abort("`model` must be an {.cls uplift_model}.")
  }
  g <- model@fit
  z_idx <- c(model@roles$treatment, model@roles$covariate)
  X <- .newdata_x(model, newdata)
  arms <- c(t0, t1)
  cov_vec <- .coverage_batch(g, z_idx, arms, X)
  dt <- data.table::data.table(
    id = seq_len(nrow(X)),
    coverage = cov_vec,
    overlap_flag = cov_vec < floor
  )
  return(dt[])
}

## proxy_decide -------------------------------------------------------------

#' Optimal action and expected incremental value per unit
#'
#' Turns the per-unit treatment effect into a next-best action under a linear
#' value model: treat when the value of the effect exceeds the cost, i.e.
#' \eqn{d^*(x) = 1\{\,\text{value} \cdot \tau(x) > \text{cost}\,\}}, with
#' expected incremental value \eqn{\text{value} \cdot \tau(x) - \text{cost}}.
#' The standard error of `tau` is propagated to an action-flip probability --
#' the chance the recommended action would reverse under sampling noise.
#'
#' @param model An [uplift_model].
#' @param newdata A data frame carrying the covariate columns.
#' @param value Numeric scalar -- the value of one unit of outcome.
#' @param cost Numeric scalar -- the cost of treating one unit. Default `0`.
#' @param t1,t0 The treated and control treatment values. Default `1` and `0`.
#' @param se_method One of `"delta"` (default) or `"mc"`.
#' @param ... Forwarded to [proxy_cate()].
#'
#' @returns A [data.table::data.table] with columns `id`, `action`,
#'   `expected_value`, `tau`, `se`, `flip_prob`, `overlap_flag`.
#' @family decision
#' @seealso [proxy_cate()], [proxy_policy_value()]
#' @export
#' @examples
#' set.seed(1)
#' n <- 400L
#' x <- stats::rnorm(n)
#' t <- stats::rbinom(n, 1L, 0.5)
#' y <- 1 + (x > 0) * t + stats::rnorm(n, sd = 0.5)
#' dat <- data.frame(y = y, t = t, x = x)
#' m <- fit_uplift(dat, "y", "t", "x", N = 2L, regime = "sample",
#'                 max_iter = 50L, seed = 1L)
#' proxy_decide(m, data.frame(x = c(-1, 1)), value = 1, cost = 0.2)
proxy_decide <- function(model,
                         newdata,
                         value,
                         cost = 0,
                         t1 = 1,
                         t0 = 0,
                         se_method = c("delta", "mc"),
                         ...) {
  if (length(value) != 1L || !is.numeric(value)) {
    cli::cli_abort("`value` must be a numeric scalar.")
  }
  if (length(cost) != 1L || !is.numeric(cost)) {
    cli::cli_abort("`cost` must be a numeric scalar.")
  }
  se_method <- rlang::arg_match(se_method)
  ce <- proxy_cate(model, newdata, t1 = t1, t0 = t0, se = TRUE,
                   se_method = se_method, ...)
  gain <- value * ce$tau - cost
  action <- as.integer(gain > 0)
  ## Action-flip probability: the gain threshold is at value * tau = cost.
  flip_prob <- stats::pnorm(-abs(gain) / (abs(value) * ce$se + .Machine$double.eps))
  dt <- data.table::data.table(
    id = ce$id,
    action = action,
    expected_value = gain,
    tau = ce$tau,
    se = ce$se,
    flip_prob = flip_prob,
    overlap_flag = ce$overlap_flag
  )
  return(dt[])
}
