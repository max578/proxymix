## Multiple imputation by Gaussian-mixture conditioning.
##
## A Gaussian mixture fitted to a data matrix with holes turns imputation
## into conditioning: for a row with observed coordinates x_o, the
## distribution of the missing coordinates x_m is the closed-form mixture
## conditional p(x_m | x_o) -- the same Schur-complement algebra as
## `gmm_conditionalise()`. Because the mixture can be multimodal and
## heteroscedastic, this fills cases a single-Gaussian model (`Amelia`) or a
## linear-Gaussian conditional (`mice`, `norm = "norm"`) cannot represent.
##
## Fitting to incomplete data is expectation-maximisation where the E-step
## responsibilities use each row's observed margin and the M-step fills the
## missing entries with their per-component conditional moments, adding back
## the conditional covariance so the component variances are not
## under-estimated (Ghahramani and Jordan 1994). Proper multiple imputation
## then draws the fitting parameters from a bootstrap of the rows, so the
## pooled inference reflects both imputation and parameter uncertainty.

# ---------------------------------------------------------------------------
# internal numerics
# ---------------------------------------------------------------------------

## Row-wise log-sum-exp of a numeric matrix.
.logsumexp_rows <- function(M) {
  mx <- apply(M, 1L, max)
  mx + log(rowSums(exp(M - mx)))
}

## Log multivariate-normal density at the rows of `X` (n by d), for mean
## `mu` (length d) and SPD covariance `Sigma` (d by d). Cholesky-based, so
## it is robust for d = 1 as well.
.dmvn_rows <- function(X, mu, Sigma) {
  d <- ncol(X)
  R <- chol(Sigma)
  cen <- t(X) - mu
  y <- backsolve(R, cen, transpose = TRUE)
  quad <- colSums(y * y)
  -0.5 * (d * log(2 * pi) + 2 * sum(log(diag(R))) + quad)
}

## Distinct missingness patterns of an observed-indicator matrix, as a list
## mapping a pattern key to the row indices that carry it.
.pattern_rows <- function(obs) {
  key <- apply(obs, 1L, function(r) paste0(as.integer(r), collapse = ""))
  split(seq_len(nrow(obs)), key)
}
.pattern_observed <- function(key) which(strsplit(key, "", fixed = TRUE)[[1]] == "1")

# ---------------------------------------------------------------------------
# EM from incomplete data
# ---------------------------------------------------------------------------

## Fit a K-component Gaussian mixture to `X` (n by p, with NA), returning a
## list(gmm, loglik, iterations, converged). Touches the ambient RNG only
## through kmeans initialisation when K > 1; the caller owns RNG state.
.fit_em_missing <- function(X, K, max_iter = 100L, tol = 1e-6,
                            ridge_eps = 1e-6, fallback_mean = NULL) {
  n <- nrow(X)
  p <- ncol(X)
  obs <- !is.na(X)
  pr <- .pattern_rows(obs)
  patterns <- names(pr)

  ## mean-impute for initialisation (use a supplied fallback for any column
  ## that is entirely missing in this (possibly bootstrapped) matrix)
  col_means <- colMeans(X, na.rm = TRUE)
  if (!is.null(fallback_mean)) {
    bad <- !is.finite(col_means)
    col_means[bad] <- fallback_mean[bad]
  }
  Xf <- X
  for (j in seq_len(p)) Xf[is.na(Xf[, j]), j] <- col_means[j]

  if (K == 1L) {
    cl <- rep(1L, n)
  } else {
    cl <- stats::kmeans(Xf, centers = K, nstart = 5L, iter.max = 50L)$cluster
  }
  w <- as.numeric(table(factor(cl, levels = seq_len(K)))) / n
  mu <- lapply(seq_len(K), function(k) colMeans(Xf[cl == k, , drop = FALSE]))
  S <- lapply(seq_len(K), function(k) {
    Sk <- if (sum(cl == k) > 1L) stats::cov(Xf[cl == k, , drop = FALSE]) else diag(p)
    Sk + ridge_eps * diag(p)
  })

  prev <- -Inf
  converged <- FALSE
  iter_done <- 0L
  for (iter in seq_len(max_iter)) {
    iter_done <- iter
    ## E-step: log responsibilities from the observed margin of each row.
    logr <- matrix(0, n, K)
    for (key in patterns) {
      rows <- pr[[key]]
      o <- .pattern_observed(key)
      if (length(o) == 0L) {
        for (k in seq_len(K)) logr[rows, k] <- log(w[k])
        next
      }
      Xo <- X[rows, o, drop = FALSE]
      for (k in seq_len(K)) {
        logr[rows, k] <- log(w[k]) +
          .dmvn_rows(Xo, mu[[k]][o], S[[k]][o, o, drop = FALSE])
      }
    }
    row_ll <- .logsumexp_rows(logr)
    Rsp <- exp(logr - row_ll)
    ll <- sum(row_ll)

    ## M-step: weighted sufficient statistics with the conditional-covariance
    ## correction on the missing block.
    Nk <- numeric(K)
    T1 <- matrix(0, K, p)
    T2 <- lapply(seq_len(K), function(k) matrix(0, p, p))
    for (key in patterns) {
      rows <- pr[[key]]
      o <- .pattern_observed(key)
      mv <- setdiff(seq_len(p), o)
      Xo <- if (length(o)) X[rows, o, drop = FALSE] else NULL
      for (k in seq_len(K)) {
        rk <- Rsp[rows, k]
        Xc <- matrix(0, length(rows), p)
        Cc <- matrix(0, p, p)
        if (length(o)) Xc[, o] <- Xo
        if (length(mv)) {
          if (length(o)) {
            Sooi <- chol2inv(chol(S[[k]][o, o, drop = FALSE]))
            reg <- S[[k]][mv, o, drop = FALSE] %*% Sooi
            cmean <- sweep(Xo, 2L, mu[[k]][o]) %*% t(reg)
            Xc[, mv] <- sweep(cmean, 2L, mu[[k]][mv], "+")
            Cc[mv, mv] <- S[[k]][mv, mv, drop = FALSE] -
              reg %*% S[[k]][o, mv, drop = FALSE]
          } else {
            Xc[, mv] <- matrix(mu[[k]][mv], length(rows), length(mv), byrow = TRUE)
            Cc[mv, mv] <- S[[k]][mv, mv, drop = FALSE]
          }
        }
        Nk[k] <- Nk[k] + sum(rk)
        T1[k, ] <- T1[k, ] + colSums(rk * Xc)
        T2[[k]] <- T2[[k]] + crossprod(Xc * sqrt(rk)) + sum(rk) * Cc
      }
    }
    for (k in seq_len(K)) {
      if (Nk[k] < .Machine$double.eps) next
      w[k] <- Nk[k] / n
      mu[[k]] <- T1[k, ] / Nk[k]
      Sk <- T2[[k]] / Nk[k] - tcrossprod(mu[[k]]) + ridge_eps * diag(p)
      S[[k]] <- (Sk + t(Sk)) / 2
    }
    w <- w / sum(w)

    if (is.finite(prev) && abs(ll - prev) < tol * (abs(prev) + tol)) {
      converged <- TRUE
      break
    }
    prev <- ll
  }
  list(
    gmm = gmm(weights = w, means = mu, covariances = S, name = "em_missing"),
    loglik = prev, iterations = iter_done, converged = converged
  )
}

## Bayesian information criterion for a fitted incomplete-data mixture.
.mixture_bic <- function(loglik, K, p, n) {
  npar <- (K - 1L) + K * p + K * p * (p + 1L) / 2L
  -2 * loglik + npar * log(n)
}

## Choose K over `candidates` by BIC on the full (incomplete) data.
.select_k <- function(X, candidates, max_iter, tol, ridge_eps) {
  p <- ncol(X)
  n <- nrow(X)
  bics <- vapply(candidates, function(K) {
    fit <- .fit_em_missing(X, K, max_iter = max_iter, tol = tol, ridge_eps = ridge_eps)
    .mixture_bic(fit$loglik, K, p, n)
  }, numeric(1))
  candidates[which.min(bics)]
}

# ---------------------------------------------------------------------------
# conditional draw (the imputation step)
# ---------------------------------------------------------------------------

## Draw one set of missing entries for every row of `X` from the mixture
## conditional p(x_missing | x_observed) under `g`. Vectorised by pattern;
## equivalent to `gmm_conditionalise()` per row (asserted in the tests).
.draw_missing <- function(g, X, obs) {
  n <- nrow(X)
  p <- ncol(X)
  K <- length(g@weights)
  out <- X
  pr <- .pattern_rows(obs)
  for (key in names(pr)) {
    rows <- pr[[key]]
    o <- .pattern_observed(key)
    mv <- setdiff(seq_len(p), o)
    if (length(mv) == 0L) next
    nr <- length(rows)
    ## per-component conditional mean (nr by |mv|), conditional cov, evidence
    logw <- matrix(0, nr, K)
    cmeans <- vector("list", K)
    cchol <- vector("list", K)
    for (k in seq_len(K)) {
      if (length(o)) {
        Xo <- X[rows, o, drop = FALSE]
        Sooi <- chol2inv(chol(g@covariances[[k]][o, o, drop = FALSE]))
        reg <- g@covariances[[k]][mv, o, drop = FALSE] %*% Sooi
        cm <- sweep(sweep(Xo, 2L, g@means[[k]][o]) %*% t(reg), 2L, g@means[[k]][mv], "+")
        Cc <- g@covariances[[k]][mv, mv, drop = FALSE] -
          reg %*% g@covariances[[k]][o, mv, drop = FALSE]
        logw[, k] <- log(g@weights[k]) +
          .dmvn_rows(Xo, g@means[[k]][o], g@covariances[[k]][o, o, drop = FALSE])
      } else {
        cm <- matrix(g@means[[k]][mv], nr, length(mv), byrow = TRUE)
        Cc <- g@covariances[[k]][mv, mv, drop = FALSE]
        logw[, k] <- log(g@weights[k])
      }
      cmeans[[k]] <- cm
      cchol[[k]] <- chol((Cc + t(Cc)) / 2 + 1e-10 * diag(length(mv)))
    }
    gamma <- exp(logw - .logsumexp_rows(logw))
    comp <- apply(gamma, 1L, function(pr_row) sample.int(K, 1L, prob = pr_row))
    draw <- matrix(0, nr, length(mv))
    for (k in seq_len(K)) {
      idx <- which(comp == k)
      if (!length(idx)) next
      z <- matrix(stats::rnorm(length(idx) * length(mv)), length(idx), length(mv))
      draw[idx, ] <- cmeans[[k]][idx, , drop = FALSE] + z %*% cchol[[k]]
    }
    out[rows, mv] <- draw
  }
  out
}

# ---------------------------------------------------------------------------
# gmm_imputation
# ---------------------------------------------------------------------------

#' A Gaussian-mixture multiple-imputation result
#'
#' The object returned by [gmm_impute()]. It carries the `m` completed data
#' matrices, the bootstrap-fitted mixtures behind them (used by the analytic
#' pooling in [proxy_pool()]), the mixture fitted to the full data, and a
#' record of the missingness. Pass it to [gmm_complete()] to extract the
#' completed datasets and to [proxy_pool()] / [proxy_fmi()] for inference.
#'
#' @param data The numeric data matrix supplied to [gmm_impute()], with `NA`
#'   for missing entries.
#' @param completions List of `m` completed data matrices.
#' @param fits List of `m` bootstrap-fitted [gmm] objects behind the
#'   completions.
#' @param point_fit The [gmm] fitted to the full data.
#' @param n_components Integer number of mixture components.
#' @param m Integer number of completions.
#' @param mechanism Missingness mechanism (currently `"mar"`).
#' @param observed Logical matrix marking the observed entries.
#' @param var_names Character vector of column names.
#' @param is_data_frame Logical; whether the input was a data frame.
#' @param diagnostics List of fit diagnostics (per-column missing rates,
#'   convergence, iterations).
#' @param call The matched call.
#'
#' @returns An S7 object of class `gmm_imputation`.
#' @family imputation
#' @export
gmm_imputation <- S7::new_class(
  name = "gmm_imputation",
  package = "proxymix",
  properties = list(
    data = S7::class_any,
    completions = S7::class_list,
    fits = S7::class_list,
    point_fit = S7::class_any,
    n_components = S7::class_integer,
    m = S7::class_integer,
    mechanism = S7::new_property(class = S7::class_any, default = "mar"),
    observed = S7::class_any,
    var_names = S7::class_character,
    is_data_frame = S7::new_property(class = S7::class_logical, default = FALSE),
    diagnostics = S7::new_property(class = S7::class_list, default = list()),
    call = S7::class_any
  )
)

#' Multiple imputation by Gaussian-mixture conditioning
#'
#' Fits a Gaussian mixture to a numeric dataset that contains missing values
#' and draws `m` completed datasets from the mixture conditional
#' \eqn{p(x_{\mathrm{missing}} \mid x_{\mathrm{observed}})}. Because the
#' mixture can be multimodal and heteroscedastic, the imputations follow the
#' shape of the joint distribution rather than a single Gaussian, which keeps
#' downstream inference valid on data that a single-Gaussian or
#' linear-Gaussian imputer mis-specifies.
#'
#' Imputation is conditioning. For a row with observed coordinates the
#' missing coordinates follow the closed-form mixture conditional (the same
#' Schur-complement algebra as [gmm_conditionalise()]). The mixture is fitted
#' to the incomplete data by expectation-maximisation whose E-step uses each
#' row's observed margin and whose M-step restores the conditional covariance
#' of the filled entries, so component variances are not under-estimated.
#'
#' Proper multiple imputation requires the fitting parameters themselves to
#' carry uncertainty, otherwise the pooled intervals are too narrow. Each of
#' the `m` imputations is therefore drawn under a mixture fitted to an
#' independent bootstrap resample of the rows, so [proxy_pool()] reflects both
#' imputation and parameter uncertainty.
#'
#' The `mechanism` says how an entry came to be missing, which sets the
#' conditional the missing value is drawn from: [mar()] (the default) for missing
#' at random, [censored()] for a known interval such as a detection limit, or
#' [mnar()] for a value-dependent selection model. The interval and
#' value-dependent gates act on a single coordinate, and a row missing that
#' coordinate must have its other coordinates observed. Numeric data only;
#' categorical variables are out of scope.
#'
#' @param data A numeric matrix or data frame with `NA` for missing entries.
#' @param N Number of mixture components. `NULL` (the default) selects it by
#'   the Bayesian information criterion over `1:6`.
#' @param m Number of completed datasets to draw. Default `20L`.
#' @param mechanism A missingness mechanism: [mar()], [censored()], or [mnar()].
#'   The string `"mar"` is also accepted. Default [mar()].
#' @param seed Optional integer seed. When supplied the result is
#'   reproducible and the ambient random-number state is restored on exit.
#' @param max_iter Maximum EM iterations per fit. Default `100L`.
#' @param tol Relative log-likelihood tolerance for EM convergence. Default
#'   `1e-6`.
#' @param ridge_eps Ridge added to each component covariance at every M-step.
#'   Default `1e-6`.
#'
#' @returns A [gmm_imputation] object.
#' @family imputation
#' @seealso [gmm_complete()] to extract completions, [proxy_pool()] to pool an
#'   estimand across them, [gmm_conditionalise()] for the conditioning algebra.
#' @export
#' @examples
#' set.seed(1)
#' x1 <- rnorm(200)
#' x2 <- x1 + rnorm(200)
#' x2[runif(200) < plogis(x1)] <- NA          # missing at random on x1
#' imp <- gmm_impute(cbind(x1, x2), N = 1L, m = 10L, seed = 1L)
#' proxy_pool(imp, "x2")$estimate             # pooled mean of x2
gmm_impute <- function(data, N = NULL, m = 20L, mechanism = mar(),
                       seed = NULL, max_iter = 100L, tol = 1e-6,
                       ridge_eps = 1e-6) {
  gate <- .as_gate(mechanism)
  is_df <- is.data.frame(data)
  if (is_df) {
    if (!all(vapply(data, is.numeric, logical(1)))) {
      cli::cli_abort(c(
        "Every column of `data` must be numeric.",
        "i" = "Categorical variables are out of scope for {.fn gmm_impute}."
      ))
    }
    var_names <- names(data)
    X <- as.matrix(data)
  } else {
    if (!is.matrix(data) || !is.numeric(data)) {
      cli::cli_abort("`data` must be a numeric matrix or data frame.")
    }
    var_names <- colnames(data)
    X <- data
  }
  if (is.null(var_names)) var_names <- paste0("V", seq_len(ncol(X)))
  storage.mode(X) <- "double"
  if (!anyNA(X)) {
    cli::cli_warn("`data` has no missing values; returning the data unchanged.")
  }
  if (all(is.na(X))) cli::cli_abort("`data` is entirely missing.")
  n <- nrow(X)
  p <- ncol(X)
  obs <- !is.na(X)
  m <- as.integer(m)
  if (length(m) != 1L || is.na(m) || m < 2L) {
    cli::cli_abort("`m` must be an integer >= 2.")
  }

  ## A value-dependent (MNAR) or interval (censored) gate acts on one coordinate,
  ## and a row missing that coordinate must have its other coordinates observed.
  gated <- gate$type != "mar"
  cj <- NA_integer_
  if (gated) {
    cj <- .gate_coord_index(gate, var_names)
    if (ncol(X) > 1L && anyNA(X[, -cj, drop = FALSE])) {
      cli::cli_abort(c(
        "A {.val {gate$type}} mechanism needs {.field coord} to be the only \\
         column with missing values.",
        "i" = "Columns other than {.val {var_names[cj]}} must be fully observed."))
    }
    if (!anyNA(X[, cj])) cli::cli_warn("`coord` has no missing values.")
  }

  ## RNG hygiene: a supplied seed is reproducible and side-effect-free.
  if (!is.null(seed)) {
    if (exists(".Random.seed", envir = globalenv())) {
      old_seed <- get(".Random.seed", envir = globalenv())
      on.exit(assign(".Random.seed", old_seed, envir = globalenv()), add = TRUE)
    } else {
      on.exit(suppressWarnings(rm(".Random.seed", envir = globalenv())), add = TRUE)
    }
    set.seed(seed)
  }

  col_means <- colMeans(X, na.rm = TRUE)
  if (is.null(N)) {
    kmax <- min(6L, max(1L, floor(n / (5L * p))))
    N <- .select_k(X, seq_len(kmax), max_iter, tol, ridge_eps)
  }
  N <- as.integer(N)
  if (length(N) != 1L || is.na(N) || N < 1L) {
    cli::cli_abort("`N` must be a positive integer scalar.")
  }

  fits <- vector("list", m)
  completions <- vector("list", m)
  accept <- rep(NA_real_, m)
  if (gated) {
    point <- .fit_em_gated(X, N, cj, gate, max_iter = max_iter, tol = tol,
                           ridge_eps = ridge_eps, fallback_mean = col_means)
    for (b in seq_len(m)) {
      bi <- sample.int(n, n, replace = TRUE)
      fb <- .fit_em_gated(X[bi, , drop = FALSE], N, cj, gate, max_iter = max_iter,
                          tol = tol, ridge_eps = ridge_eps, fallback_mean = col_means)
      fits[[b]] <- fb$gmm
      dr <- .draw_gated(fb$gmm, X, cj, gate, fb$alpha)
      cb <- dr$data
      accept[b] <- dr$accept
      colnames(cb) <- var_names
      completions[[b]] <- cb
    }
  } else {
    point <- .fit_em_missing(X, N, max_iter = max_iter, tol = tol,
                             ridge_eps = ridge_eps, fallback_mean = col_means)
    for (b in seq_len(m)) {
      bi <- sample.int(n, n, replace = TRUE)
      fit_b <- .fit_em_missing(X[bi, , drop = FALSE], N, max_iter = max_iter,
                               tol = tol, ridge_eps = ridge_eps,
                               fallback_mean = col_means)$gmm
      fits[[b]] <- fit_b
      cb <- .draw_missing(fit_b, X, obs)
      colnames(cb) <- var_names
      completions[[b]] <- cb
    }
  }

  gmm_imputation(
    data = X, completions = completions, fits = fits, point_fit = point$gmm,
    n_components = N, m = m, mechanism = gate, observed = obs,
    var_names = var_names, is_data_frame = is_df,
    diagnostics = list(
      missing_rate = colMeans(!obs),
      converged = point$converged, iterations = point$iterations,
      mnar_alpha = if (gated) point$alpha else NA_real_,
      min_accept = if (gated) min(accept) else NA_real_
    ),
    call = match.call()
  )
}

#' Extract completed datasets from a `gmm_imputation`
#'
#' @param object A [gmm_imputation].
#' @param which Either an integer vector of imputation indices, or `"all"`
#'   for every completion. Default `1L`.
#'
#' @returns When `which` selects one completion, a single completed dataset
#'   (matrix, or data frame if the input was one); otherwise a list of them.
#' @family imputation
#' @export
gmm_complete <- function(object, which = 1L) {
  if (!S7::S7_inherits(object, gmm_imputation)) {
    cli::cli_abort("`object` must be a {.cls gmm_imputation}.")
  }
  m <- object@m
  if (identical(which, "all")) which <- seq_len(m)
  if (!is.numeric(which) || any(which < 1L) || any(which > m)) {
    cli::cli_abort("`which` must be in {.val 1}:{.val {m}} or \"all\".")
  }
  out <- lapply(which, function(i) {
    cm <- object@completions[[i]]
    if (object@is_data_frame) as.data.frame(cm) else cm
  })
  if (length(out) == 1L) out[[1L]] else out
}

S7::method(print, gmm_imputation) <- function(x, ...) {
  cat(sprintf("<gmm_imputation>: m = %d completions, K = %d components, p = %d\n",
              x@m, x@n_components, ncol(x@data)))
  mech <- x@mechanism
  mech_str <- if (inherits(mech, "proxymix_gate")) {
    switch(mech$type,
      mar = "missing at random",
      censored = sprintf("censored on %s to [%g, %g]", mech$coord, mech$lower, mech$upper),
      mnar = sprintf("MNAR on %s (%s, beta = %g)", mech$coord, mech$link, mech$beta))
  } else {
    as.character(mech)
  }
  cat(sprintf("  mechanism  : %s\n", mech_str))
  mr <- x@diagnostics$missing_rate
  cat(sprintf("  missing    : %s\n",
              paste(sprintf("%s %.0f%%", x@var_names, 100 * mr), collapse = ", ")))
  invisible(x)
}
