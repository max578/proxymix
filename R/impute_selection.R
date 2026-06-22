## The gated imputation engine: fitting and drawing when a single coordinate is
## missing through a value-dependent or interval mechanism.
##
## Both the missing-not-at-random and the censored paths are the same algebra: a
## mixture fitted to the data, with the designated coordinate's conditional
## multiplied by a gate delta(y) before it is integrated, filled, or drawn from.
## The gate is a smooth selection weight 1 - g(alpha + beta y) for MNAR, computed
## by Gauss-Hermite quadrature, or the indicator of an interval for censoring,
## computed in closed form as truncated-Gaussian moments. The expectation-
## maximisation that fits the mixture therefore accounts for the mechanism rather
## than ignoring it, which is what lets the missing-not-at-random fit recover an
## estimand that an ignorable fit is biased on.
##
## Scope: the gate acts on one coordinate, and a row missing that coordinate has
## its other coordinates observed (the detection-limit / single-outcome case).

# ---------------------------------------------------------------------------
# Gauss-Hermite nodes (cached), and gated 1-D moments
# ---------------------------------------------------------------------------

.gh_env <- new.env(parent = emptyenv())

## 64-point Gauss-Hermite nodes/weights for weight exp(-x^2), via Golub-Welsch.
## Cached on first use so the eigen-decomposition runs once per session.
.gh_nodes <- function(n = 64L) {
  key <- as.character(n)
  if (is.null(.gh_env[[key]])) {
    i <- seq_len(n - 1L)
    b <- sqrt(i / 2)
    J <- matrix(0, n, n)
    J[cbind(i, i + 1L)] <- b
    J[cbind(i + 1L, i)] <- b
    e <- eigen(J, symmetric = TRUE)
    o <- order(e$values)
    .gh_env[[key]] <- list(x = e$values[o], w = (sqrt(pi) * e$vectors[1L, ]^2)[o])
  }
  .gh_env[[key]]
}

## Moments of the smooth-gated conditional, vectorised over the per-row means
## `mu` at a single component variance `s2`. Gate is g(a + b y), the probability
## of being missing at value y, so the gated conditional is the imputation law
## p(y | missing). Returns the gate normaliser `I` (the per-row probability the
## coordinate is missing under component k) and the mean/variance of that law.
.gated_smooth <- function(mu, s2, a, b, link) {
  gh <- .gh_nodes()
  Y <- outer(mu, sqrt(2 * s2) * gh$x, `+`)            # length(mu) x Q
  d <- if (link == "logit") stats::plogis(a + b * Y) else stats::pnorm(a + b * Y)
  W <- matrix(gh$w, nrow = length(mu), ncol = length(gh$x), byrow = TRUE)
  I  <- rowSums(W * d) / sqrt(pi)
  m1 <- rowSums(W * Y * d) / sqrt(pi) / I
  m2 <- rowSums(W * Y * Y * d) / sqrt(pi) / I
  list(I = I, m1 = m1, v = pmax(m2 - m1 * m1, 0))
}

## Moments of the interval-gated (truncated) conditional, closed form,
## vectorised over the per-row means `mu`. Gate is 1[L <= y <= U].
.gated_trunc <- function(mu, s2, L, U) {
  sdv <- sqrt(s2)
  aL <- (L - mu) / sdv
  aU <- (U - mu) / sdv
  Z <- pmax(stats::pnorm(aU) - stats::pnorm(aL), 1e-300)
  lam <- (stats::dnorm(aL) - stats::dnorm(aU)) / Z
  tL <- ifelse(is.finite(aL), aL * stats::dnorm(aL), 0)
  tU <- ifelse(is.finite(aU), aU * stats::dnorm(aU), 0)
  v <- s2 * (1 + (tL - tU) / Z - lam * lam)
  ## the truncated mean lies in [L, U]; clamp to guard the numerical overshoot
  ## when a component barely overlaps the interval (Z near zero).
  loC <- if (is.finite(L)) L else mu - 10 * sdv
  hiC <- if (is.finite(U)) U else mu + 10 * sdv
  m1 <- pmin(pmax(mu + sdv * lam, loC), hiC)
  list(I = Z, m1 = m1, v = pmax(pmin(v, (hiC - loC)^2), 0))
}

## Expected gate at marginal level: E_{N(mu, s2)}[g(a + b y)] for one component,
## used to calibrate the MNAR intercept to the observed missingness rate.
.expected_g <- function(mu, s2, a, b, link) {
  gh <- .gh_nodes()
  y <- mu + sqrt(2 * s2) * gh$x
  g <- if (link == "logit") stats::plogis(a + b * y) else stats::pnorm(a + b * y)
  sum(gh$w * g) / sqrt(pi)
}

# ---------------------------------------------------------------------------
# per-component conditional of the gated coordinate given the others
# ---------------------------------------------------------------------------

## For component k with mean `mu` and covariance `S`, the conditional of
## coordinate `cj` given the other coordinates `Xo` (n by p-1, all observed):
## per-row mean (length n), the component conditional variance (scalar), and the
## log density of the observed margin at each row.
.coord_cond_k <- function(mu, S, cj, Xo, rest) {
  if (length(rest)) {
    Soo <- S[rest, rest, drop = FALSE]
    Sjo <- S[cj, rest, drop = FALSE]
    Sooi <- chol2inv(chol(Soo))
    reg <- Sjo %*% Sooi
    cmean <- as.numeric(mu[cj] + sweep(Xo, 2L, mu[rest]) %*% t(reg))
    cvar <- as.numeric(S[cj, cj] - reg %*% S[rest, cj, drop = FALSE])
    logo <- .dmvn_rows(Xo, mu[rest], Soo)
  } else {
    cmean <- rep(mu[cj], nrow(Xo))
    cvar <- S[cj, cj]
    logo <- rep(0, nrow(Xo))
  }
  list(cmean = cmean, cvar = max(cvar, 1e-12), logo = logo)
}

# ---------------------------------------------------------------------------
# the gated EM
# ---------------------------------------------------------------------------

## Fit a K-component mixture to `X` (only column `cj` carries NA) under a gate.
## `gate` is the resolved mechanism: list(type, link, beta) for MNAR or
## list(type, lower, upper) for censored. Returns the fitted gmm, the calibrated
## MNAR intercept (NA for censored), and convergence diagnostics. The caller owns
## the RNG (kmeans init only).
.fit_em_gated <- function(X, K, cj, gate, max_iter = 200L, tol = 1e-6,
                          ridge_eps = 1e-6, fallback_mean = NULL) {
  n <- nrow(X)
  p <- ncol(X)
  rest <- setdiff(seq_len(p), cj)
  miss <- is.na(X[, cj])
  p_miss <- mean(miss)
  is_mnar <- identical(gate$type, "mnar")

  ## initialise the missing coordinate. For censoring the missing values lie in
  ## the censored interval, so the observed mean is the wrong side of the bound; a
  ## representative point inside the interval is a stable start. For the smooth
  ## mechanisms the observed mean is the neutral start.
  obs_mean <- if (any(!miss)) mean(X[!miss, cj]) else
    (if (!is.null(fallback_mean)) fallback_mean[cj] else 0)
  obs_sd <- if (sum(!miss) > 1L) stats::sd(X[!miss, cj]) else 1
  fillv <- if (identical(gate$type, "censored")) {
    if (is.finite(gate$lower) && is.finite(gate$upper)) (gate$lower + gate$upper) / 2
    else if (is.finite(gate$upper)) gate$upper - obs_sd
    else gate$lower + obs_sd
  } else {
    obs_mean
  }
  Xf <- X
  Xf[miss, cj] <- fillv
  cl <- if (K == 1L) rep(1L, n) else
    stats::kmeans(Xf, centers = K, nstart = 5L, iter.max = 50L)$cluster
  w <- as.numeric(table(factor(cl, levels = seq_len(K)))) / n
  mu <- lapply(seq_len(K), function(k) colMeans(Xf[cl == k, , drop = FALSE]))
  S <- lapply(seq_len(K), function(k) {
    Sk <- if (sum(cl == k) > 1L) stats::cov(Xf[cl == k, , drop = FALSE]) else diag(p)
    Sk + ridge_eps * diag(p)
  })
  beta <- if (is_mnar) gate$beta else NA_real_
  link <- if (is_mnar) gate$link else NA_character_
  ## intercept calibrated to the marginal missingness rate at the initial fit
  alpha <- if (is_mnar) .calibrate_alpha(w, mu, S, cj, beta, link, p_miss) else NA_real_

  Xo <- X[, rest, drop = FALSE]
  prev <- -Inf
  converged <- FALSE
  iter_done <- 0L
  for (iter in seq_len(max_iter)) {
    iter_done <- iter
    logr <- matrix(0, n, K)
    cm_list <- vector("list", K)        # gated conditional means (missing rows)
    cv_list <- numeric(K)               # gated conditional variances
    for (k in seq_len(K)) {
      cc <- .coord_cond_k(mu[[k]], S[[k]], cj, Xo, rest)
      lk <- numeric(n)
      ## coordinate-observed rows: full joint density
      if (any(!miss)) {
        lk[!miss] <- log(w[k]) + cc$logo[!miss] +
          stats::dnorm(X[!miss, cj], cc$cmean[!miss], sqrt(cc$cvar), log = TRUE)
      }
      ## coordinate-missing rows: integrate the gated conditional
      if (any(miss)) {
        gm <- if (is_mnar) {
          .gated_smooth(cc$cmean[miss], cc$cvar, alpha, beta, link)
        } else {
          .gated_trunc(cc$cmean[miss], cc$cvar, gate$lower, gate$upper)
        }
        lk[miss] <- log(w[k]) + cc$logo[miss] + log(pmax(gm$I, 1e-300))
        cm_list[[k]] <- gm$m1
        cv_list[k] <- NA       # carried per-row below via attribute
        attr(cm_list[[k]], "v") <- gm$v
      } else {
        cm_list[[k]] <- numeric(0)
      }
      logr[, k] <- lk
    }
    logr[is.nan(logr)] <- -Inf          # a vanished component must not poison the row
    row_ll <- .logsumexp_rows(logr)
    Rsp <- exp(logr - row_ll)
    ll <- sum(row_ll)

    ## M-step: weighted sufficient statistics, filling the missing coordinate
    ## with its gated conditional mean and restoring its gated variance.
    Nk <- colSums(Rsp)
    for (k in seq_len(K)) {
      if (!is.finite(Nk[k]) || Nk[k] < .Machine$double.eps) next
      Xk <- X
      vk <- numeric(n)
      if (any(miss)) {
        Xk[miss, cj] <- cm_list[[k]]
        vk[miss] <- attr(cm_list[[k]], "v")
      }
      rk <- Rsp[, k]
      muk <- colSums(rk * Xk) / Nk[k]
      d <- sweep(Xk, 2L, muk)
      Sk <- crossprod(d * sqrt(rk)) / Nk[k]
      Sk[cj, cj] <- Sk[cj, cj] + sum(rk * vk) / Nk[k]      # gated variance correction
      Sk <- Sk + ridge_eps * diag(p)
      w[k] <- Nk[k] / n
      mu[[k]] <- muk
      S[[k]] <- (Sk + t(Sk)) / 2
    }
    w <- w / sum(w)
    if (is_mnar) alpha <- .calibrate_alpha(w, mu, S, cj, beta, link, p_miss)
    if (is.finite(prev) && abs(ll - prev) < tol * (abs(prev) + tol)) {
      converged <- TRUE
      break
    }
    prev <- ll
  }
  list(gmm = gmm(weights = w, means = mu, covariances = S, name = "em_gated"),
       alpha = alpha, beta = beta, link = link, loglik = prev,
       iterations = iter_done, converged = converged, p_miss = p_miss)
}

## Calibrate the MNAR intercept so the mixture-marginal missingness rate of the
## gated coordinate equals the observed rate `p_miss`. Monotone in alpha.
.calibrate_alpha <- function(w, mu, S, cj, beta, link, p_miss) {
  rate <- function(a) {
    sum(vapply(seq_along(w), function(k)
      w[k] * .expected_g(mu[[k]][cj], S[[k]][cj, cj], a, beta, link), numeric(1)))
  }
  if (p_miss <= 0) return(-Inf)
  if (p_miss >= 1) return(Inf)
  lo <- -50; hi <- 50
  if (rate(lo) > p_miss) return(lo)
  if (rate(hi) < p_miss) return(hi)
  stats::uniroot(function(a) rate(a) - p_miss, c(lo, hi), tol = 1e-8)$root
}

# ---------------------------------------------------------------------------
# drawing imputations from the gated conditional
# ---------------------------------------------------------------------------

## Draw the missing coordinate for every missing row from the gated mixture
## conditional under the fitted `g`. MNAR uses rejection from the untilted
## conditional with acceptance 1 - g (the gate); censoring inverts the truncated
## normal. Returns the completed matrix and the minimum per-draw acceptance
## (an effective-sample diagnostic for the MNAR gate).
.draw_gated <- function(g, X, cj, gate, alpha) {
  p <- ncol(X)
  rest <- setdiff(seq_len(p), cj)
  miss <- which(is.na(X[, cj]))
  out <- X
  if (!length(miss)) return(list(data = out, accept = 1))
  K <- length(g@weights)
  Xo <- X[miss, rest, drop = FALSE]
  is_mnar <- identical(gate$type, "mnar")
  logw <- matrix(0, length(miss), K)
  cmean <- matrix(0, length(miss), K)
  cvar <- numeric(K)
  for (k in seq_len(K)) {
    cc <- .coord_cond_k(g@means[[k]], g@covariances[[k]], cj, Xo, rest)
    cmean[, k] <- cc$cmean
    cvar[k] <- cc$cvar
    gI <- if (is_mnar) {
      .gated_smooth(cc$cmean, cc$cvar, alpha, gate$beta, gate$link)$I
    } else {
      .gated_trunc(cc$cmean, cc$cvar, gate$lower, gate$upper)$I
    }
    logw[, k] <- log(g@weights[k]) + cc$logo + log(pmax(gI, 1e-300))
  }
  gamma <- exp(logw - .logsumexp_rows(logw))
  comp <- apply(gamma, 1L, function(pr) sample.int(K, 1L, prob = pr))
  draw <- numeric(length(miss))
  accept <- 1
  for (k in seq_len(K)) {
    idx <- which(comp == k)
    if (!length(idx)) next
    mu_i <- cmean[idx, k]
    sdv <- sqrt(cvar[k])
    if (is_mnar) {
      acc <- .draw_smooth_reject(mu_i, sdv, alpha, gate$beta, gate$link)
      draw[idx] <- acc$y
      accept <- min(accept, acc$rate)
    } else {
      draw[idx] <- .draw_trunc(mu_i, sdv, gate$lower, gate$upper)
    }
  }
  out[miss, cj] <- draw
  list(data = out, accept = accept)
}

## Rejection draw from N(mu, sd^2) * g(alpha + beta y) / I, the imputation law
## for a value known to be missing (gate g = probability of missing).
.draw_smooth_reject <- function(mu, sdv, alpha, beta, link, max_round = 200L) {
  n <- length(mu)
  y <- numeric(n)
  todo <- seq_len(n)
  rounds <- 0L
  while (length(todo) && rounds < max_round) {
    rounds <- rounds + 1L
    prop <- stats::rnorm(length(todo), mu[todo], sdv)
    gp <- if (link == "logit") stats::plogis(alpha + beta * prop) else
      stats::pnorm(alpha + beta * prop)
    keep <- stats::runif(length(todo)) < gp
    y[todo[keep]] <- prop[keep]
    todo <- todo[!keep]
  }
  if (length(todo)) {            # guardrail: fall back to the gated mean
    gm <- .gated_smooth(mu[todo], sdv^2, alpha, beta, link)
    y[todo] <- gm$m1
  }
  list(y = y, rate = (n - length(todo)) / n)
}

## Inverse-CDF draw from N(mu, sd^2) truncated to [L, U].
.draw_trunc <- function(mu, sdv, L, U) {
  pL <- stats::pnorm((L - mu) / sdv)
  pU <- stats::pnorm((U - mu) / sdv)
  u <- pL + stats::runif(length(mu)) * (pU - pL)
  mu + sdv * stats::qnorm(pmin(pmax(u, 1e-12), 1 - 1e-12))
}
