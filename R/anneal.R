## Deterministic-annealing EM and phase-transition component discovery.
##
## The EM objective is a variational free energy F = <E> - T H. Softening the
## E-step responsibilities by a temperature T,
##   gamma_ik proportional to (pi_k N(x_i; mu_k, Sigma_k))^(1/T),
## and lowering T from a high value toward 1 on a geometric schedule, turns the
## fit into a homotopy: at high T the annealed objective has a single smooth
## basin (every component collapses toward the global centroid), and as T falls
## the centroids bifurcate at critical temperatures -- genuine phase transitions
## (Rose, Gurewitz and Fox 1990).
##
## Two facilities are built on this. First, a robust warm-start: components are
## located by shared-isotropic deterministic annealing (the regime in which the
## phase transitions have a clean closed-form critical temperature), and the
## annealed centroids then seed the existing, unchanged cold EM loop for the
## final free-covariance polish. This attacks the package's EM-local-optima
## hard-risk without altering the cold path. Second, `gmm_anneal_path()` tracks
## the number of distinct centroids as T falls -- a physics-derived
## component-count diagnostic.
##
## For a single cluster with data covariance C and shared reference covariance
## Sigma, the first bifurcation occurs at T_c = lambda_max(Sigma^{-1} C). With
## the (pi_k N_k)^(1/T) convention used here this is the Rose-Gurewitz-Fox
## critical temperature 2 lambda_max(C) once the factor-of-two from the Gaussian
## exponent (Sigma = sigma^2 I gives exp(-||x - mu||^2 / (2 sigma^2 T))) and the
## squared-Euclidean distortion are reconciled.

# ---------------------------------------------------------------------------
# Schedule, reference statistics, and the shared-covariance annealed sweep
# ---------------------------------------------------------------------------

## Geometric (log-linear) temperature schedule, strictly descending from
## `t_high` to `t_low` in `n_steps` steps, ending exactly at `t_low`.
.anneal_schedule <- function(t_high, t_low = 1, n_steps = 20L) {
  n_steps <- as.integer(n_steps)
  if (n_steps < 1L) {
    cli::cli_abort("`n_steps` must be a positive integer.")
  }
  if (!is.finite(t_high) || !is.finite(t_low) ||
        t_high <= 0 || t_low <= 0) {
    cli::cli_abort("temperatures must be finite and positive.")
  }
  if (t_high < t_low) {
    cli::cli_abort("`t_high` ({t_high}) must be at least `t_low` ({t_low}).")
  }
  if (n_steps == 1L) {
    return(t_low)
  }
  exp(seq(log(t_high), log(t_low), length.out = n_steps))
}

## Row-weighted centroid, covariance, isotropic scale and leading eigenvalue of
## a data matrix. `rw` is a length-n vector of non-negative weights summing to
## one (1/n for i.i.d. samples; importance weights for regime (iii)).
.anneal_reference <- function(x, rw) {
  n <- nrow(x)
  p <- ncol(x)
  mu0 <- as.numeric(colSums(rw * x))
  diff0 <- x - matrix(mu0, nrow = n, ncol = p, byrow = TRUE)
  C <- symmetrise(crossprod(diff0 * sqrt(rw)))
  eig <- eigen(C, symmetric = TRUE)
  list(mu0 = mu0, C = C, scale0 = sqrt(mean(diag(C))),
       lambda_max = eig$values[1L], pc1 = eig$vectors[, 1L])
}

## One annealed E + M step with a fixed shared covariance Sigma: the
## responsibilities are softened by the temperature, the means and weights
## update, and the covariance is held as the annealing reference scale.
##
## This is the Rose-Gurewitz-Fox *mass-constrained* association,
##   gamma_ik proportional to pi_k * N(x_i; mu_k, Sigma)^(1/T),
## in which the prior mass pi_k enters linearly and only the Gaussian distortion
## is tempered. Raising the prior itself to the 1/T power (the literal
## (pi_k N_k)^(1/T) form) makes a small mass imbalance explode into winner-take-
## all collapse as T falls; the linear form is stable and, because the priors
## cancel at the symmetric bifurcation point, leaves the closed-form critical
## temperature T_c = lambda_max(Sigma^{-1} C) unchanged. At T = 1 it coincides
## with the cold E-step. Returns the updated weights and means and the
## deterministic-annealing free energy
##   F(T) = -T sum_i w_i log sum_k pi_k N(x_i; mu_k, Sigma)^(1/T).
.anneal_centroid_sweep <- function(x, w, weights, means, Sigma, temperature) {
  N <- length(weights)
  covs <- rep(list(Sigma), N)
  gauss_log <- gmm_log_unnorm(x, rep(1, N), means, covs)
  la <- sweep(gauss_log / temperature, 2L, log(weights), `+`)
  lse <- logsumexp_rows(la)
  resp <- exp(la - lse)
  wr <- resp * w
  Nk <- colSums(wr)
  new_weights <- pmax(as.numeric(Nk), 1e-300)
  new_weights <- new_weights / sum(new_weights)
  new_means <- means
  for (k in seq_len(N)) {
    if (Nk[k] >= 1e-12) {
      new_means[[k]] <- as.numeric(colSums(wr[, k] * x) / Nk[k])
    }
  }
  list(weights = new_weights, means = new_means,
       free_energy = -temperature * sum(w * lse))
}

## Drive the shared-covariance annealing dynamics with `n_codevectors`
## codevectors down `schedule`, applying a symmetry-breaking perturbation of
## standard deviation `perturb_abs` to the centroids at each temperature. When
## `track` is `TRUE`, also record the distinct-centroid count (merging at
## `merge_abs`) and the free energy at each temperature. The caller is
## responsible for seeding (the perturbations use `stats::rnorm`).
.anneal_da_run <- function(x, w, n_codevectors, ref, Sigma, schedule,
                           perturb_abs, n_inner, merge_abs = NULL,
                           track = FALSE, init_spread = 0.1) {
  p <- ncol(x)
  N <- n_codevectors
  ## Symmetry-agnostic initialisation: a small isotropic cloud around the
  ## centroid. Stringing the codevectors along a principal axis would bias the
  ## first bifurcation and collapse onto one direction when the data covariance
  ## is (near-)isotropic; an isotropic cloud lets the per-temperature
  ## perturbation break symmetry in every direction.
  means <- lapply(seq_len(N), function(k) {
    if (N == 1L) {
      as.numeric(ref$mu0)
    } else {
      as.numeric(ref$mu0 + stats::rnorm(p, sd = init_spread * ref$scale0))
    }
  })
  weights <- rep(1 / N, N)
  n_eff <- if (track) integer(length(schedule)) else NULL
  free_energy <- if (track) numeric(length(schedule)) else NULL
  for (s in seq_along(schedule)) {
    temperature <- schedule[s]
    if (N > 1L && perturb_abs > 0) {
      for (k in seq_len(N)) {
        means[[k]] <- means[[k]] + stats::rnorm(p, sd = perturb_abs)
      }
    }
    last_fe <- NA_real_
    for (inner in seq_len(n_inner)) {
      step <- .anneal_centroid_sweep(x, w, weights, means, Sigma, temperature)
      weights <- step$weights
      means <- step$means
      last_fe <- step$free_energy
    }
    if (track) {
      n_eff[s] <- .count_distinct_centroids(means, merge_abs)
      free_energy[s] <- last_fe
    }
  }
  list(weights = weights, means = means, n_eff = n_eff,
       free_energy = free_energy)
}

## Default data-adaptive warm-start schedule, or validation of a user override.
## The temperatures live in shared-isotropic units (Sigma = mean-variance * I),
## so the defaults bracket the bifurcation cascade: from 3 * T_c down to a
## fraction of T_c that resolves the requested components without deep
## over-fragmentation.
.anneal_warmstart_schedule <- function(temp_schedule, t_c, n_steps = 40L) {
  if (is.null(temp_schedule)) {
    return(.anneal_schedule(t_high = 3 * t_c, t_low = 0.02 * t_c,
                            n_steps = n_steps))
  }
  if (!is.numeric(temp_schedule) || length(temp_schedule) < 1L ||
        any(!is.finite(temp_schedule)) || any(temp_schedule <= 0)) {
    cli::cli_abort("`temp_schedule` must be a vector of positive finite temperatures.")
  }
  as.numeric(temp_schedule)
}

## Deterministic-annealing warm-start. Locate `N` component centroids by
## shared-isotropic annealing, then estimate free per-component covariances and
## masses for the cold-EM handoff. The temperature schedule actually used is
## stored in the returned `gmm`'s metadata so the caller can record it.
.anneal_em_warmstart <- function(x, rw, N, temp_schedule = NULL,
                                 ridge_eps = 1e-6, n_inner = 10L,
                                 perturb = 0.02, seed = NULL) {
  n <- nrow(x)
  p <- ncol(x)
  N <- as.integer(N)
  ref <- .anneal_reference(x, rw)
  sigma_sq <- mean(diag(ref$C))
  Sigma <- ridge(diag(sigma_sq, p), ridge_eps)
  t_c <- ref$lambda_max / sigma_sq
  schedule <- .anneal_warmstart_schedule(temp_schedule, t_c)
  perturb_abs <- perturb * ref$scale0

  run <- function() {
    da <- .anneal_da_run(x, rw, N, ref, Sigma, schedule, perturb_abs, n_inner)
    means <- da$means
    ## Free-covariance handoff: soft responsibilities under the shared Sigma at
    ## T = 1, then a responsibility-weighted covariance per component.
    covs_ref <- rep(list(Sigma), N)
    log_unnorm <- gmm_log_unnorm(x, da$weights, means, covs_ref)
    resp <- exp(log_unnorm - logsumexp_rows(log_unnorm))
    wr <- resp * rw
    Nk <- colSums(wr)
    new_covs <- vector("list", N)
    for (k in seq_len(N)) {
      if (Nk[k] < 1e-9) {
        new_covs[[k]] <- symmetrise(ridge(ref$C, ridge_eps))
        next
      }
      diff <- x - matrix(means[[k]], nrow = n, ncol = p, byrow = TRUE)
      S_k <- crossprod(diff * sqrt(wr[, k])) / Nk[k]
      new_covs[[k]] <- symmetrise(ridge(S_k, ridge_eps))
    }
    weights <- pmax(as.numeric(Nk), 1e-300)
    weights <- weights / sum(weights)
    gmm(weights = weights, means = means, covariances = new_covs,
        name = "anneal_warmstart",
        metadata = list(temp_schedule = schedule))
  }
  if (is.null(seed)) run() else withr::with_seed(seed, run())
}

# ---------------------------------------------------------------------------
# Phase-transition critical temperature and component discovery
# ---------------------------------------------------------------------------

## Closed-form first critical temperature of a single cluster: the largest
## eigenvalue of Sigma^{-1} C, where C is the (weighted) data covariance and
## Sigma the shared reference covariance. This is the temperature at which the
## degenerate single-centroid solution loses stability and bifurcates.
.critical_temperature <- function(C, Sigma) {
  M <- solve(Sigma, C)
  vals <- eigen(M, symmetric = FALSE, only.values = TRUE)$values
  max(Re(vals))
}

## Count distinct centroids by connected components of the "closer than `tol`"
## graph (transitive single-linkage merging).
.count_distinct_centroids <- function(means, tol) {
  N <- length(means)
  if (N == 0L) {
    return(0L)
  }
  M <- do.call(rbind, means)
  D <- as.matrix(stats::dist(M))
  adj <- D < tol
  comp <- integer(N)
  cid <- 0L
  for (i in seq_len(N)) {
    if (comp[i] == 0L) {
      cid <- cid + 1L
      stack <- i
      while (length(stack) > 0L) {
        v <- stack[length(stack)]
        stack <- stack[-length(stack)]
        if (comp[v] == 0L) {
          comp[v] <- cid
          nb <- which(adj[v, ] & comp == 0L)
          stack <- c(stack, nb)
        }
      }
    }
  }
  cid
}

## The component count occupying the widest cumulative log-temperature span on
## the n_eff(T) staircase. This is the robust readout of the discovered
## component number: well-separated clusters resist fragmentation over a broad
## temperature range, whereas the deep-cold over-fragmentation toward `k_max`
## occupies only narrow ranges.
.widest_plateau <- function(temps, n_eff) {
  log_t <- log(temps)
  m <- length(temps)
  width <- numeric(m)
  for (s in seq_len(m)) {
    lo <- if (s > 1L) (log_t[s - 1L] + log_t[s]) / 2 else log_t[s]
    hi <- if (s < m) (log_t[s] + log_t[s + 1L]) / 2 else log_t[s]
    width[s] <- abs(lo - hi)
  }
  agg <- tapply(width, n_eff, sum)
  as.integer(names(agg)[which.max(agg)])
}

#' Phase-transition component discovery by deterministic annealing
#'
#' Tracks the number of distinct mixture centroids as a function of temperature
#' under mass-constrained deterministic annealing (Rose, Gurewitz and Fox 1990),
#' a physics-derived alternative to information-criterion model selection. The
#' system starts at a high temperature where all `k_max` centroids collapse to
#' the data centroid (a single effective component) and is cooled along a
#' geometric schedule; at each critical temperature a centroid bifurcates, so the
#' number of distinct centroids grows in steps. The temperatures at which it
#' grows are the phase transitions, and the count occupying the widest
#' temperature range is the discovered component number.
#'
#' The first bifurcation has a closed-form critical temperature
#' \eqn{T_c = \lambda_{\max}(\Sigma^{-1} C)}, where \eqn{C} is the (weighted)
#' data covariance and \eqn{\Sigma = \sigma^2 I} the shared reference covariance.
#' This value is returned as `t_critical_analytic` and serves as an independent
#' analytic check on the empirically detected first transition. Subsequent
#' transitions have no comparably simple closed form, and the count is a
#' diagnostic rather than a guarantee.
#'
#' Annealing fixes the component covariance to the reference \eqn{\Sigma} so the
#' temperature is the only scale; this is the clean isotropic regime in which the
#' critical temperature is exact. For robust *fitting* under free covariances,
#' use `anneal = TRUE` on [fit_em_samples()] or [fit_kld_em()] instead.
#'
#' @param x A numeric `n` by `p` matrix of samples, or a [gmm_target] carrying a
#'   `samples` matrix. For regime (iii) targets, pass an importance-resampled
#'   draw.
#' @param k_max Maximum number of centroids tracked (the discovered count is at
#'   most `k_max`).
#' @param sigma Reference scale: the shared covariance is `sigma^2 * I`. When
#'   `NULL` (the default) `sigma` is `1`, so the first critical temperature is
#'   the largest eigenvalue of the data covariance.
#' @param t_high,t_low Top and bottom of the temperature schedule. When `NULL`
#'   they default to `3 * t_critical_analytic` and `0.05 * t_critical_analytic`,
#'   bracketing the bifurcation cascade.
#' @param n_steps Number of temperatures on the geometric schedule.
#' @param n_inner Fixed-point iterations run at each temperature.
#' @param w Optional length-`n` vector of non-negative observation weights
#'   (e.g. importance weights). Defaults to uniform.
#' @param perturb Symmetry-breaking perturbation, as a fraction of the data
#'   scale, applied to the centroids at each temperature.
#' @param merge_tol Two centroids count as distinct when their distance exceeds
#'   `merge_tol` times the data scale.
#' @param ridge_eps Ridge added to the reference covariance for stability.
#' @param seed Optional integer seed for the perturbations (the result is
#'   deterministic given a seed).
#'
#' @returns A list with elements `path` (a data frame of `temperature`,
#'   `n_effective` and `free_energy`), `critical_temperatures` (the temperatures
#'   at which the count increased), `first_critical_temperature` (the first such,
#'   or `NA` if none was detected), `t_critical_analytic`
#'   (\eqn{\lambda_{\max}(\Sigma^{-1} C)}), `k_selected` (the widest-plateau
#'   component count), `lambda_max` and `sigma`.
#' @family diagnostics
#' @references Rose, K., Gurewitz, E. and Fox, G. C. (1990) Statistical
#'   mechanics and phase transitions in clustering. *Physical Review Letters*
#'   65(8), 945--948. \doi{10.1103/PhysRevLett.65.945}
#' @export
#' @examples
#' set.seed(1)
#' x <- rbind(
#'   matrix(stats::rnorm(120, mean = -4), ncol = 2),
#'   matrix(stats::rnorm(120, mean =  4), ncol = 2)
#' )
#' path <- gmm_anneal_path(x, k_max = 4L, n_steps = 40L)
#' path$k_selected
#' path$first_critical_temperature
gmm_anneal_path <- function(x, k_max = 8L,
                            sigma = NULL,
                            t_high = NULL, t_low = NULL,
                            n_steps = 80L,
                            n_inner = 30L,
                            w = NULL,
                            perturb = 0.02,
                            merge_tol = 0.1,
                            ridge_eps = 1e-6,
                            seed = 1L) {
  if (S7::S7_inherits(x, gmm_target)) {
    if (is.null(x@samples)) {
      cli::cli_abort(c(
        "`x` is a {.cls gmm_target} without `samples`.",
        "i" = "Pass a numeric matrix or an importance-resampled draw."
      ))
    }
    x <- x@samples
  }
  if (!is.matrix(x) || !is.numeric(x)) {
    cli::cli_abort("`x` must be a numeric matrix or a {.cls gmm_target} with samples.")
  }
  n <- nrow(x)
  p <- ncol(x)
  k_max <- as.integer(k_max)
  if (k_max < 1L) {
    cli::cli_abort("`k_max` must be a positive integer.")
  }
  if (is.null(w)) {
    w <- rep(1 / n, n)
  } else {
    if (length(w) != n || any(!is.finite(w)) || any(w < 0)) {
      cli::cli_abort("`w` must be a length-{n} vector of non-negative finite weights.")
    }
    w <- w / sum(w)
  }

  ref <- .anneal_reference(x, w)
  sigma_sq <- if (is.null(sigma)) 1 else as.numeric(sigma)^2
  Sigma <- ridge(diag(sigma_sq, p), ridge_eps)
  t_c_analytic <- .critical_temperature(ref$C, Sigma)

  if (is.null(t_high)) t_high <- 3 * t_c_analytic
  if (is.null(t_low)) t_low <- 0.05 * t_c_analytic
  schedule <- .anneal_schedule(t_high, t_low, n_steps)
  merge_abs <- merge_tol * ref$scale0
  perturb_abs <- perturb * ref$scale0

  da <- withr::with_seed(seed, {
    .anneal_da_run(x, w, k_max, ref, Sigma, schedule, perturb_abs,
                   n_inner, merge_abs = merge_abs, track = TRUE)
  })

  n_eff <- pmin(da$n_eff, k_max)
  increases <- which(diff(n_eff) > 0L) + 1L
  critical_temperatures <- schedule[increases]
  first_critical <- if (length(increases) > 0L) {
    schedule[increases[1L]]
  } else {
    NA_real_
  }

  list(
    path = data.frame(
      temperature = schedule,
      n_effective = n_eff,
      free_energy = da$free_energy
    ),
    critical_temperatures = critical_temperatures,
    first_critical_temperature = first_critical,
    t_critical_analytic = t_c_analytic,
    k_selected = .widest_plateau(schedule, n_eff),
    lambda_max = ref$lambda_max,
    sigma = sqrt(sigma_sq)
  )
}
