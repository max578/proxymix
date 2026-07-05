## Regime (iii) applied to an objective: a mixture proxy for its optima.
##
## An objective f(x) to be minimised over a box defines an unnormalised
## density, the Gibbs (Boltzmann) measure exp(-f(x) / T), whose mass
## concentrates on the low regions of f as the temperature T falls. That
## density can be evaluated point-wise but not directly sampled -- exactly
## regime (iii) of Hoek and Elliott (2024). `from_objective()` therefore
## fits a Gaussian-mixture proxy to exp(-f / T) by cooling a short
## temperature ladder through importance-sampled KLD-EM (`fit_kld_em()`),
## warm-starting each step from the previous fit. The result is a closed-
## form mixture whose components settle over the basins of f; `gmm_modes()`
## then resolves it into the distinct optima. Multimodal objectives are
## mapped together rather than one optimum at a time.
##
## The engine is the package's own `fit_kld_em()` M-step throughout: the
## rank-weighted Gaussian update that minimises KL from a mixture to a
## peaked target is the same update an evolution strategy performs, so this
## is the existing regime-(iii) fit driven against a Gibbs target, not a
## separate optimiser.

# ---------------------------------------------------------------------------
# internal helpers
# ---------------------------------------------------------------------------

## Row-wise objective evaluation: `objective` takes a length-p numeric
## vector and returns a scalar; this returns one value per row of `X`.
## Errors and non-scalar / non-numeric returns are mapped to `NA_real_`
## (the caller replaces those with a finite penalty), so objectives defined
## only on part of the box never abort the fit.
.objective_rows <- function(objective, X, p) {
  if (is.null(dim(X))) X <- matrix(X, ncol = p)
  apply(X, 1L, function(z) {
    val <- tryCatch(objective(z), error = function(e) NA_real_)
    if (length(val) != 1L || !is.numeric(val)) NA_real_ else val
  })
}

## Are the rows of `X` inside the box [lower, upper] (closed)?
.in_box <- function(X, lower, upper) {
  lo <- matrix(lower, nrow = nrow(X), ncol = ncol(X), byrow = TRUE)
  hi <- matrix(upper, nrow = nrow(X), ncol = ncol(X), byrow = TRUE)
  rowSums(X >= lo & X <= hi) == ncol(X)
}

## Numerically stable log(exp(a) + exp(b)) element-wise.
.log_add_exp <- function(a, b) pmax(a, b) + log1p(exp(-abs(a - b)))

## A defensive adaptive importance proposal: with probability
## `exploration` draw uniformly over the box (so no basin is ever starved),
## otherwise draw from the (inflated) current mixture `g` to exploit the
## basins found so far.
.defensive_proposal <- function(g, lower, upper, p, log_vol, exploration) {
  sampler <- function(n) {
    n <- as.integer(n)
    n_unif <- stats::rbinom(1L, n, exploration)
    ru <- matrix(0, nrow = n_unif, ncol = p)
    if (n_unif > 0L) {
      for (j in seq_len(p)) ru[, j] <- stats::runif(n_unif, lower[j], upper[j])
    }
    if (n_unif >= n) return(ru)
    rbind(ru, rgmm(n - n_unif, g))
  }
  log_dens <- function(X) {
    if (is.null(dim(X))) X <- matrix(X, ncol = p)
    inside <- .in_box(X, lower, upper)
    log_unif <- ifelse(inside, log(exploration) - log_vol, -Inf)
    .log_add_exp(log_unif, log1p(-exploration) + dgmm(X, g, log = TRUE))
  }
  is_proposal(
    n_dim = p, sample = sampler, log_density = log_dens,
    name = "from_objective_defensive"
  )
}

# ---------------------------------------------------------------------------
# from_objective
# ---------------------------------------------------------------------------

#' Map the optima of an objective with a Gaussian-mixture proxy
#'
#' Fits a Gaussian-mixture proxy to the Gibbs measure
#' \eqn{\exp(-f(x) / T)} of a user-supplied objective `f` over a bounded
#' box, by cooling a short temperature ladder through regime-(iii)
#' importance-sampled KLD-EM ([fit_kld_em()]). As the temperature falls the
#' mixture mass concentrates on the low regions of `f`, so the fitted
#' mixture is a closed-form *map* over the optima rather than a single
#' point estimate. Pair it with [gmm_modes()] to read off the distinct
#' optima.
#'
#' The Gibbs measure can be evaluated point-wise but not directly sampled,
#' which is precisely the setting of regime (iii): minimising
#' the Kullback-Leibler divergence from a Gaussian mixture to a peaked
#' target is the rank-weighted Gaussian update at the heart of [fit_kld_em()].
#' This function is that fit, driven against a sequence of cooling Gibbs
#' targets and warm-started from the previous fit at each step. Because a
#' multimodal `f` produces a multimodal target, the components spread across
#' the basins and recover the optima together.
#'
#' Recovery is most reliable with component headroom -- a number of
#' components `N` comfortably larger than the number of optima you expect,
#' so the defensive proposal can keep a component on each basin. Symmetric
#' landscapes (where several optima are exchangeable) need the most headroom.
#'
#' Dimensional scope. As with the rest of regime (iii), the importance-
#' sampling effective sample size falls sharply with dimension; the guard is
#' `p <= 5` (recommended), `p <= 10` (allowed with a warning), `p > 10`
#' (rejected).
#'
#' @param objective Function taking a length-`p` numeric vector and
#'   returning a finite numeric scalar -- the objective to minimise (or
#'   maximise; see `minimise`).
#' @param lower,upper Numeric vectors of equal length `p` giving the box
#'   over which the optima are sought. Every `upper` must exceed its
#'   `lower`.
#' @param N Number of mixture components. Default `max(10L, 5L * p)`. Raise
#'   it for objectives with many optima or strong symmetry.
#' @param minimise Logical. If `TRUE` (the default) lower `objective`
#'   values are better (the proxy concentrates on the minima); if `FALSE`
#'   the proxy concentrates on the maxima.
#' @param temperature Optional control of the cooling ladder. `NULL` (the
#'   default) derives a ladder automatically from a uniform probe of the
#'   landscape. A positive scalar sets the final (lowest) temperature; a
#'   length-2 numeric `c(high, low)` sets both ends explicitly.
#' @param n_steps Number of temperatures in the cooling ladder. Default
#'   `6L`.
#' @param exploration Probability mass the importance proposal places on
#'   uniform exploration of the box at each step, in `[0, 1]`. Default
#'   `0.5`. Larger values explore more (and starve no basin); smaller values
#'   exploit the basins found so far.
#' @param inflate Factor by which the current mixture covariances are
#'   inflated when used as the exploitation part of the proposal. Default
#'   `1.8`.
#' @param is_size Importance-sample size per cooling step. Default `1e4L`.
#' @param max_iter Maximum EM iterations per cooling step. Default `70L`.
#' @param ridge_eps Ridge added to each component covariance at every
#'   M-step. Forwarded to [fit_kld_em()].
#' @param seed Optional integer seed for reproducibility.
#'
#' @returns A [gmm_fit] (the fitted proxy) carrying a `from_objective`
#'   metadata record with the temperature ladder and box. Pass it to
#'   [gmm_modes()] to extract the distinct optima.
#' @family fitting
#' @seealso [gmm_modes()] to resolve the fitted map into distinct optima.
#' @export
#' @examples
#' ## A bimodal 1-D objective with minima at +/- 2.
#' f <- function(v) (v[1]^2 - 4)^2
#' fit <- from_objective(f, lower = -5, upper = 5, N = 6L,
#'                       is_size = 2000L, n_steps = 5L, seed = 1L)
#' gmm_modes(fit)$modes
from_objective <- function(objective, lower, upper,
                           N = NULL, minimise = TRUE,
                           temperature = NULL, n_steps = 6L,
                           exploration = 0.5, inflate = 1.8,
                           is_size = 1e4L, max_iter = 70L,
                           ridge_eps = 1e-4, seed = NULL) {
  if (!is.function(objective)) {
    cli::cli_abort("`objective` must be a function of a length-p numeric vector.")
  }
  lower <- as.numeric(lower)
  upper <- as.numeric(upper)
  if (length(lower) != length(upper)) {
    cli::cli_abort("`lower` and `upper` must have the same length.")
  }
  if (anyNA(lower) || anyNA(upper) || any(!is.finite(c(lower, upper)))) {
    cli::cli_abort("`lower` and `upper` must be finite.")
  }
  if (any(upper <= lower)) {
    cli::cli_abort("`upper` must exceed `lower` in every coordinate.")
  }
  p <- length(lower)
  if (p > 10L) {
    cli::cli_abort(c(
      "{.fn from_objective} currently supports {.code p <= 10} (got p = {p}).",
      "i" = "Regime-(iii) importance sampling degrades sharply in high dimensions."
    ))
  }
  if (p > 5L) {
    cli::cli_warn(c(
      "{.fn from_objective} is well-tested for {.code p <= 5} (got p = {p}).",
      "i" = "Expect importance-sampling ESS to fall sharply beyond p = 5."
    ))
  }
  if (is.null(N)) N <- max(10L, 5L * p)
  N <- as.integer(N)
  if (length(N) != 1L || is.na(N) || N < 1L) {
    cli::cli_abort("`N` must be a positive integer scalar.")
  }
  n_steps <- as.integer(n_steps)
  if (length(n_steps) != 1L || is.na(n_steps) || n_steps < 1L) {
    cli::cli_abort("`n_steps` must be a positive integer scalar.")
  }
  if (length(exploration) != 1L || is.na(exploration) ||
      exploration < 0 || exploration > 1) {
    cli::cli_abort("`exploration` must be a single number in [0, 1].")
  }
  if (length(inflate) != 1L || is.na(inflate) || inflate < 1) {
    cli::cli_abort("`inflate` must be a single number >= 1.")
  }
  sgn <- if (isTRUE(minimise)) -1 else 1
  vol <- prod(upper - lower)
  log_vol <- sum(log(upper - lower))

  ## Restore the RNG state on exit so a supplied `seed` is side-effect-free.
  if (!is.null(seed)) {
    if (exists(".Random.seed", envir = globalenv())) {
      old_seed <- get(".Random.seed", envir = globalenv())
      on.exit(assign(".Random.seed", old_seed, envir = globalenv()), add = TRUE)
    } else {
      on.exit(
        suppressWarnings(rm(".Random.seed", envir = globalenv())), add = TRUE
      )
    }
    set.seed(seed)
  }

  ## Probe the landscape to scale the temperature ladder.
  probe <- matrix(0, nrow = 2000L, ncol = p)
  for (j in seq_len(p)) probe[, j] <- stats::runif(2000L, lower[j], upper[j])
  o_probe <- .objective_rows(objective, probe, p)
  o_probe <- o_probe[is.finite(o_probe)]
  if (length(o_probe) == 0L) {
    cli::cli_abort(c(
      "`objective` returned no finite values over the box.",
      "i" = "Check that `objective` is finite on `[lower, upper]`."
    ))
  }
  spread <- diff(range(o_probe))
  ## A finite penalty for points where the objective is non-finite or
  ## undefined: strongly unattractive at every temperature, but never -Inf
  ## (a -Inf log-density would give the fit a 0 * -Inf NaN in its KLD trace).
  penalty_spread <- if (is.finite(spread) && spread > 0) spread else 1
  penalty <- if (isTRUE(minimise)) {
    max(o_probe) + 10 * penalty_spread
  } else {
    min(o_probe) - 10 * penalty_spread
  }
  if (is.null(temperature)) {
    t_high <- max(spread / 2, 1e-6)
    t_low <- t_high / 40
  } else if (length(temperature) == 1L) {
    t_low <- as.numeric(temperature)
    t_high <- max(spread / 2, t_low)
  } else if (length(temperature) == 2L) {
    t_high <- max(temperature)
    t_low <- min(temperature)
  } else {
    cli::cli_abort("`temperature` must be NULL, a scalar, or a length-2 vector.")
  }
  if (!is.finite(t_low) || t_low <= 0) {
    cli::cli_abort("`temperature` must be positive and finite.")
  }
  ladder <- exp(seq(log(t_high), log(t_low), length.out = n_steps))

  ## The Gibbs target exp(sgn * objective / T). The box is the search
  ## region (it scales the ladder, the uniform exploration, and the init);
  ## the target itself is evaluated everywhere so the importance weights
  ## stay finite, with non-finite objective values mapped to `penalty`.
  make_target <- function(temp) {
    gmm_target(
      n_dim = p,
      log_density = function(X) {
        if (is.null(dim(X))) X <- matrix(X, ncol = p)
        v <- .objective_rows(objective, X, p)
        v[!is.finite(v)] <- penalty
        sgn * v / temp
      },
      normalised = FALSE,
      name = "from_objective"
    )
  }

  ## Diverse space-filling initialisation: one component near each of N
  ## random points, with covariances at a fifth of the box span.
  centres <- matrix(0, nrow = N, ncol = p)
  for (j in seq_len(p)) centres[, j] <- stats::runif(N, lower[j], upper[j])
  span <- (upper - lower) / 5
  fit <- gmm(
    weights = rep(1 / N, N),
    means = lapply(seq_len(N), function(i) centres[i, ]),
    covariances = rep(list(diag(span^2, nrow = p)), N)
  )

  ## Cool through the ladder, warm-starting each step from the last fit.
  for (temp in ladder) {
    g <- gmm(
      weights = fit@weights, means = fit@means,
      covariances = lapply(fit@covariances,
                           function(S) inflate * S + 1e-3 * diag(p))
    )
    proposal <- .defensive_proposal(g, lower, upper, p, log_vol, exploration)
    fit <- fit_kld_em(
      make_target(temp), N = N, proposal = proposal,
      is_size = is_size, init = gmm(weights = fit@weights, means = fit@means,
                                    covariances = fit@covariances),
      max_iter = max_iter, ridge_eps = ridge_eps, seed = seed,
      support_warn = FALSE, canonicalise = FALSE
    )
  }

  fit@metadata <- modifyList(fit@metadata, list(
    from_objective = list(
      lower = lower, upper = upper, minimise = isTRUE(minimise),
      temperatures = ladder, N = N
    )
  ))
  fit
}

# ---------------------------------------------------------------------------
# gmm_modes
# ---------------------------------------------------------------------------

#' Modes of a Gaussian mixture
#'
#' Returns the distinct local modes of a Gaussian-mixture density by
#' Gaussian mean-shift (the fixed-point hill-climb of Carreira-Perpinan
#' 2000) started from each component mean, together with the mixture density
#' at each mode. Nearby converged points are merged so that each genuine
#' mode is reported once.
#'
#' A mixture of `K` components has at most `K` modes, and fewer when
#' components overlap; mean-shift from every component mean finds them
#' robustly without a grid. The companion of [from_objective()]: applied to
#' the fitted map it returns the recovered optima (ordered by density, which
#' for a Gibbs proxy ranks the deepest optima first).
#'
#' @param object A [gmm] or [gmm_fit].
#' @param starts Optional `m`-by-`p` matrix of starting points for the
#'   mean-shift. Default `NULL` uses the component means.
#' @param tol Convergence tolerance on the mean-shift step length. Default
#'   `1e-5`.
#' @param dedup Distance below which two converged points are treated as the
#'   same mode. Default `NULL` derives a small fraction of the spread of the
#'   component means.
#' @param max_iter Maximum mean-shift iterations per start. Default `200L`.
#'
#' @returns A list with `modes` (an `n`-by-`p` matrix of distinct modes,
#'   ordered by descending mixture density), `density` (the mixture density
#'   at each mode), and `n` (the number of modes).
#' @family ops
#' @seealso [from_objective()], whose fitted map this resolves into optima.
#' @export
#' @examples
#' g <- gmm(
#'   weights = rep(1 / 3, 3),
#'   means = list(c(-3, 0), c(3, 0), c(0, 4)),
#'   covariances = rep(list(0.3 * diag(2)), 3)
#' )
#' gmm_modes(g)$modes
gmm_modes <- function(object, starts = NULL, tol = 1e-5,
                      dedup = NULL, max_iter = 200L) {
  if (!S7::S7_inherits(object, gmm)) {
    cli::cli_abort("`object` must be a {.cls gmm} object.")
  }
  weights <- object@weights
  means <- do.call(rbind, object@means)
  K <- length(weights)
  p <- ncol(means)
  chols <- lapply(object@covariances, chol)
  precisions <- lapply(chols, chol2inv)
  log_dets <- vapply(chols, function(R) 2 * sum(log(diag(R))), numeric(1))
  log_weights <- log(weights)
  const <- -0.5 * p * log(2 * pi)

  ## Component log-densities at a single point `x` (length p).
  comp_log_dens <- function(x) {
    vapply(seq_len(K), function(k) {
      d <- x - means[k, ]
      const - 0.5 * log_dets[k] - 0.5 * sum(d * (precisions[[k]] %*% d))
    }, numeric(1))
  }

  ## One Gaussian mean-shift step (responsibility-weighted precision pool).
  shift_once <- function(x) {
    lr <- log_weights + comp_log_dens(x)
    r <- exp(lr - max(lr))
    r <- r / sum(r)
    a_mat <- matrix(0, nrow = p, ncol = p)
    b_vec <- numeric(p)
    for (k in seq_len(K)) {
      rp <- r[k] * precisions[[k]]
      a_mat <- a_mat + rp
      b_vec <- b_vec + as.numeric(rp %*% means[k, ])
    }
    as.numeric(solve(a_mat, b_vec))
  }

  if (is.null(starts)) starts <- means
  if (is.null(dim(starts))) starts <- matrix(starts, ncol = p)
  n_start <- nrow(starts)

  converged <- matrix(0, nrow = n_start, ncol = p)
  for (i in seq_len(n_start)) {
    x <- starts[i, ]
    for (it in seq_len(max_iter)) {
      x_new <- shift_once(x)
      step <- sqrt(sum((x_new - x)^2))
      x <- x_new
      if (step < tol) break
    }
    converged[i, ] <- x
  }

  if (is.null(dedup)) {
    scale <- if (K > 1L) sqrt(mean(apply(means, 2L, stats::var))) else 1
    dedup <- max(1e-4, 0.02 * scale)
  }
  dens <- dgmm(converged, object)
  ord <- order(dens, decreasing = TRUE)
  modes <- matrix(numeric(0), ncol = p)
  keep <- integer(0)
  for (i in ord) {
    if (nrow(modes) == 0L) {
      modes <- converged[i, , drop = FALSE]
      keep <- i
    } else {
      gaps <- sqrt(rowSums(
        (modes - matrix(converged[i, ], nrow(modes), p, byrow = TRUE))^2
      ))
      if (all(gaps > dedup)) {
        modes <- rbind(modes, converged[i, , drop = FALSE])
        keep <- c(keep, i)
      }
    }
  }
  list(modes = modes, density = dens[keep], n = nrow(modes))
}
