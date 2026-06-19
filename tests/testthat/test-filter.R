## G3: gmm_filter -- bounded Gaussian-sum filtering over an observation series.
##
## Independent oracles (none of them call gmm_affine / gmm_observe / gmm_reduce):
##   * a hand-coded Kalman / Gaussian-sum bank `.hand_filter()` that expands
##     process- and measurement-noise Gaussian sums explicitly and runs the
##     textbook predict / update recursion, returning the filtered mixtures and
##     the per-step log marginal evidence;
##   * an analytic random-walk simulation for the bounded-horizon tracking check.

## ---------------------------------------------------------------------------
## Independent oracle: a textbook predict/update bank.
## ---------------------------------------------------------------------------

## `Qspec` / `Rspec` are each either a covariance matrix (a single Gaussian) or
## a `gmm` (a Gaussian-sum noise). The bank expands every prior component
## through every process-noise component at predict, and every predicted
## component through every measurement-noise component at update -- the exact
## Gaussian-sum filter, with no reduction.
.noise_parts <- function(spec, dim) {
  if (S7::S7_inherits(spec, gmm)) {
    list(w = spec@weights, mu = spec@means, S = spec@covariances)
  } else {
    list(w = 1, mu = list(rep(0, dim)), S = list(spec))
  }
}

.hand_filter <- function(prior, A, b, Qspec, C, d, Rspec, Y) {
  p <- length(prior@means[[1L]])
  m <- nrow(C)
  w <- prior@weights
  mu <- prior@means
  S <- prior@covariances
  qp <- .noise_parts(Qspec, p)
  rp <- .noise_parts(Rspec, m)
  n <- nrow(Y)
  filtered <- vector("list", n)
  log_ev <- numeric(n)
  for (t in seq_len(n)) {
    ## Predict: expand each component through each process-noise component.
    pw <- numeric(0); pmu <- list(); pS <- list()
    for (i in seq_along(w)) {
      for (jq in seq_along(qp$w)) {
        pw <- c(pw, w[i] * qp$w[jq])
        pmu <- c(pmu, list(as.numeric(A %*% mu[[i]] + b + qp$mu[[jq]])))
        pS <- c(pS, list(A %*% S[[i]] %*% t(A) + qp$S[[jq]]))
      }
    }
    ## Update: expand each predicted component through each measurement-noise
    ## component, Kalman-update, and reweight by the marginal evidence.
    yt <- Y[t, ]
    uw_log <- numeric(0); umu <- list(); uS <- list()
    for (k in seq_along(pw)) {
      for (jr in seq_along(rp$w)) {
        mean_obs <- as.numeric(C %*% pmu[[k]] + d + rp$mu[[jr]])
        Sx <- C %*% pS[[k]] %*% t(C) + rp$S[[jr]]
        Sx <- 0.5 * (Sx + t(Sx))
        gain <- pS[[k]] %*% t(C) %*% solve(Sx)
        uw_log <- c(uw_log, log(pw[k]) + log(rp$w[jr]) +
                      mvnfast::dmvn(matrix(yt, nrow = 1L), mu = mean_obs,
                                    sigma = Sx, log = TRUE))
        umu <- c(umu, list(as.numeric(pmu[[k]] + gain %*% (yt - mean_obs))))
        Su <- (diag(p) - gain %*% C) %*% pS[[k]]
        uS <- c(uS, list(0.5 * (Su + t(Su))))
      }
    }
    log_ev[t] <- log(sum(exp(uw_log - max(uw_log)))) + max(uw_log)
    uw <- exp(uw_log - max(uw_log))
    uw <- uw / sum(uw)
    filtered[[t]] <- gmm(weights = uw, means = umu, covariances = uS)
    w <- uw; mu <- umu; S <- uS
  }
  list(filtered = filtered, log_ev = log_ev)
}

## Order-invariant parameter distance between two mixtures (sort by the first
## mean coordinate, then weight).
.mix_dist <- function(a, b) {
  oa <- order(vapply(a@means, function(x) x[1L], numeric(1L)), a@weights)
  ob <- order(vapply(b@means, function(x) x[1L], numeric(1L)), b@weights)
  max(c(max(abs(a@weights[oa] - b@weights[ob])),
        max(abs(unlist(a@means[oa]) - unlist(b@means[ob]))),
        max(abs(unlist(a@covariances[oa]) - unlist(b@covariances[ob])))))
}

## ---------------------------------------------------------------------------
## GS5: filter parity against the independent bank.
## ---------------------------------------------------------------------------

test_that("GS5a: gmm_filter at K = 1 matches a textbook Kalman filter", {
  ## A two-dimensional constant-velocity track, one-dimensional position obs.
  prior <- gmm(weights = 1, means = list(c(0, 1)),
               covariances = list(diag(c(1, 1))))
  A <- matrix(c(1, 1, 0, 1), nrow = 2L, byrow = TRUE)
  Q <- 0.01 * diag(2)
  C <- matrix(c(1, 0), nrow = 1L)
  R <- matrix(0.2, 1L, 1L)
  set.seed(11)
  Y <- matrix(cumsum(stats::rnorm(15, sd = 0.4)) + stats::rnorm(15, sd = 0.4),
              ncol = 1L)

  ref <- .hand_filter(prior, A, b = 0, Q, C, d = 0, R, Y)
  out <- gmm_filter(prior,
                    dynamics = list(A = A, Q = Q),
                    measurement = list(C = C, R = R),
                    y = Y, ridge_eps = 0)

  for (t in seq_len(nrow(Y))) {
    expect_lt(.mix_dist(out$filtered[[t]], ref$filtered[[t]]), 1e-7)
  }
  expect_lt(max(abs(out$summary$log_evidence - ref$log_ev)), 1e-7)
})

test_that("GS5b: gmm_filter matches a Gaussian-sum bank under mixture process noise", {
  prior <- gmm(weights = 1, means = list(0), covariances = list(matrix(1)))
  A <- matrix(1, 1L, 1L)
  ## A well-separated two-component process noise, so the growing bank stays
  ## distinct and an order-invariant comparison is unambiguous.
  Q <- gmm(weights = c(0.5, 0.5), means = list(-3, 3),
           covariances = list(matrix(0.1), matrix(0.1)))
  C <- matrix(1, 1L, 1L)
  R <- matrix(0.5, 1L, 1L)
  Y <- matrix(c(0.2, -0.5, 1.1), ncol = 1L)

  ref <- .hand_filter(prior, A, b = 0, Q, C, d = 0, R, Y)
  out <- gmm_filter(prior,
                    dynamics = list(A = A, Q = Q),
                    measurement = list(C = C, R = R),
                    y = Y, ridge_eps = 0)

  ## No reduction: the component count is the exact Gaussian-sum bank (2, 4, 8).
  expect_equal(out$summary$n_components, c(2L, 4L, 8L))
  for (t in seq_len(nrow(Y))) {
    expect_lt(.mix_dist(out$filtered[[t]], ref$filtered[[t]]), 1e-7)
  }
  expect_lt(max(abs(out$summary$log_evidence - ref$log_ev)), 1e-7)
})

test_that("GS5c: gmm_filter matches the bank under mixture measurement noise", {
  ## A glint-style two-component measurement noise (tight core + heavy tail).
  prior <- gmm(weights = 1, means = list(0), covariances = list(matrix(1)))
  A <- matrix(1, 1L, 1L)
  Q <- matrix(0.05, 1L, 1L)
  C <- matrix(1, 1L, 1L)
  R <- gmm(weights = c(0.8, 0.2), means = list(0, 0),
           covariances = list(matrix(0.2), matrix(4)))
  Y <- matrix(c(0.5, 2.5), ncol = 1L)

  ref <- .hand_filter(prior, A, b = 0, Q, C, d = 0, R, Y)
  out <- gmm_filter(prior,
                    dynamics = list(A = A, Q = Q),
                    measurement = list(C = C, R = R),
                    y = Y, ridge_eps = 0)

  expect_equal(out$summary$n_components, c(2L, 4L))
  for (t in seq_len(nrow(Y))) {
    expect_lt(.mix_dist(out$filtered[[t]], ref$filtered[[t]]), 1e-7)
  }
  expect_lt(max(abs(out$summary$log_evidence - ref$log_ev)), 1e-7)
})

test_that("GS5d: a non-zero dynamics/measurement offset matches the bank", {
  prior <- gmm(weights = 1, means = list(0), covariances = list(matrix(1)))
  A <- matrix(0.9, 1L, 1L)
  Q <- matrix(0.1, 1L, 1L)
  C <- matrix(2, 1L, 1L)
  R <- matrix(0.3, 1L, 1L)
  Y <- matrix(c(1.0, 0.4, -0.3), ncol = 1L)

  ref <- .hand_filter(prior, A, b = 0.5, Q, C, d = -1, R, Y)
  out <- gmm_filter(prior,
                    dynamics = list(A = A, b = 0.5, Q = Q),
                    measurement = list(C = C, R = R, d = -1),
                    y = Y, ridge_eps = 0)
  for (t in seq_len(nrow(Y))) {
    expect_lt(.mix_dist(out$filtered[[t]], ref$filtered[[t]]), 1e-7)
  }
})

## ---------------------------------------------------------------------------
## GS6: bounded horizon.
## ---------------------------------------------------------------------------

test_that("GS6: a 200-step GSF with mixture noise stays bounded and tracks", {
  set.seed(42)
  p <- 1L
  prior <- gmm(weights = 1, means = list(0), covariances = list(matrix(1)))
  ## Heavy-tailed process noise: mostly small steps, occasional large jumps.
  Q <- gmm(weights = c(0.85, 0.15), means = list(0, 0),
           covariances = list(matrix(0.02), matrix(0.5)))
  A <- matrix(1, 1L, 1L)
  C <- matrix(1, 1L, 1L)
  R <- matrix(0.3, 1L, 1L)

  n <- 200L
  truth <- numeric(n)
  x <- 0
  for (t in seq_len(n)) {
    x <- x + as.numeric(rgmm(1L, Q))
    truth[t] <- x
  }
  y <- truth + stats::rnorm(n, sd = sqrt(0.3))

  out <- gmm_filter(prior,
                    dynamics = list(A = A, Q = Q),
                    measurement = list(C = C, R = R),
                    y = y, k_max = 5L, reduce = "merge")

  ## Bounded: the component count never exceeds the cap.
  expect_true(all(out$summary$n_components <= 5L))
  expect_length(out$filtered, n)
  ## Finite throughout: no blow-up, no NaN.
  expect_true(all(is.finite(out$mean)))
  expect_true(all(is.finite(out$summary$log_evidence)))
  expect_true(all(vapply(out$cov, function(S) all(is.finite(S)), logical(1L))))
  ## Tracks: the filtered mean is closer to the truth than the raw observation.
  rmse_filter <- sqrt(mean((out$mean[, 1L] - truth)^2))
  rmse_obs <- sqrt(mean((y - truth)^2))
  expect_lt(rmse_filter, rmse_obs)
})

test_that("GS6: reduction keeps the moment-preserving mean within tolerance", {
  ## With a single Gaussian process noise the count is constant, so a generous
  ## cap is a no-op and the filter equals the unreduced Kalman filter.
  prior <- gmm(weights = 1, means = list(0), covariances = list(matrix(1)))
  set.seed(7)
  y <- cumsum(stats::rnorm(30, sd = 0.3)) + stats::rnorm(30, sd = 0.5)
  base <- gmm_filter(prior,
                     dynamics = list(A = matrix(1), Q = matrix(0.09)),
                     measurement = list(C = matrix(1), R = matrix(0.25)),
                     y = y, ridge_eps = 0)
  capped <- gmm_filter(prior,
                       dynamics = list(A = matrix(1), Q = matrix(0.09)),
                       measurement = list(C = matrix(1), R = matrix(0.25)),
                       y = y, k_max = 10L, ridge_eps = 0)
  expect_equal(base$mean, capped$mean)
  expect_equal(base$summary$n_components, rep(1L, 30L))
})

## ---------------------------------------------------------------------------
## Composition, return shape, and convenience input forms.
## ---------------------------------------------------------------------------

test_that("the return value has the documented shape", {
  prior <- gmm(weights = c(0.5, 0.5), means = list(c(-1, 0), c(1, 0)),
               covariances = list(diag(2), diag(2)))
  A <- diag(2)
  Q <- 0.1 * diag(2)
  C <- diag(2)
  R <- 0.2 * diag(2)
  Y <- matrix(stats::rnorm(8), ncol = 2L)

  out <- gmm_filter(prior,
                    dynamics = list(A = A, Q = Q),
                    measurement = list(C = C, R = R),
                    y = Y)
  expect_named(out, c("filtered", "mean", "cov", "summary"))
  expect_length(out$filtered, 4L)
  expect_true(all(vapply(out$filtered, function(g) S7::S7_inherits(g, gmm),
                         logical(1L))))
  expect_equal(dim(out$mean), c(4L, 2L))
  expect_length(out$cov, 4L)
  expect_equal(dim(out$cov[[1L]]), c(2L, 2L))
  expect_equal(out$summary$step, seq_len(4L))
  expect_true(all(c("mean_1", "mean_2", "sd_1", "sd_2",
                    "n_components", "log_evidence") %in% names(out$summary)))
  ## The summary mean equals the matrix mean.
  expect_equal(as.matrix(out$summary[, c("mean_1", "mean_2")]),
               out$mean, ignore_attr = TRUE)
})

test_that("a scalar observation series is accepted as a plain vector", {
  prior <- gmm(weights = 1, means = list(0), covariances = list(matrix(1)))
  y <- c(0.5, -0.2, 0.1)
  out_vec <- gmm_filter(prior,
                        dynamics = list(A = matrix(1), Q = matrix(0.1)),
                        measurement = list(C = matrix(1), R = matrix(0.25)),
                        y = y, ridge_eps = 0)
  out_mat <- gmm_filter(prior,
                        dynamics = list(A = matrix(1), Q = matrix(0.1)),
                        measurement = list(C = matrix(1), R = matrix(0.25)),
                        y = matrix(y, ncol = 1L), ridge_eps = 0)
  out_list <- gmm_filter(prior,
                         dynamics = list(A = matrix(1), Q = matrix(0.1)),
                         measurement = list(C = matrix(1), R = matrix(0.25)),
                         y = as.list(y), ridge_eps = 0)
  expect_equal(out_vec$mean, out_mat$mean)
  expect_equal(out_vec$mean, out_list$mean)
})

test_that("time-varying dynamics/measurement functions equal the constant form", {
  prior <- gmm(weights = 1, means = list(c(0, 0)),
               covariances = list(diag(2)))
  A <- matrix(c(1, 0.5, 0, 1), nrow = 2L, byrow = TRUE)
  Q <- 0.05 * diag(2)
  C <- matrix(c(1, 0), nrow = 1L)
  R <- matrix(0.3, 1L, 1L)
  Y <- matrix(stats::rnorm(6), ncol = 1L)

  const <- gmm_filter(prior,
                      dynamics = list(A = A, Q = Q),
                      measurement = list(C = C, R = R),
                      y = Y, ridge_eps = 0)
  varying <- gmm_filter(prior,
                        dynamics = function(t) list(A = A, Q = Q),
                        measurement = function(t) list(C = C, R = R),
                        y = Y, ridge_eps = 0)
  expect_equal(const$mean, varying$mean)
  for (t in seq_len(nrow(Y))) {
    expect_lt(.mix_dist(const$filtered[[t]], varying$filtered[[t]]), 1e-10)
  }
})

## ---------------------------------------------------------------------------
## Vanishing-evidence handling.
## ---------------------------------------------------------------------------

test_that("a far observation under Gaussian noise yields a no-update step", {
  prior <- gmm(weights = 1, means = list(0), covariances = list(matrix(1)))
  expect_warning(
    out <- gmm_filter(prior,
                      dynamics = list(A = matrix(1), Q = matrix(1e-6)),
                      measurement = list(C = matrix(1), R = matrix(1e-10)),
                      y = 1e6),
    regexp = "Marginal evidence"
  )
  expect_equal(out$summary$log_evidence, -Inf)
})

test_that("a far observation under mixture noise yields a no-update step", {
  prior <- gmm(weights = 1, means = list(0), covariances = list(matrix(1)))
  R <- gmm(weights = c(0.5, 0.5), means = list(0, 0),
           covariances = list(matrix(1e-10), matrix(1e-10)))
  suppressWarnings(
    out <- gmm_filter(prior,
                      dynamics = list(A = matrix(1), Q = matrix(1e-6)),
                      measurement = list(C = matrix(1), R = R),
                      y = 1e6)
  )
  expect_equal(out$summary$log_evidence, -Inf)
  ## The belief is the predicted prior, unchanged by the (refused) update.
  expect_equal(gmm_n_components(out$filtered[[1L]]), 1L)
})

## ---------------------------------------------------------------------------
## Input validation.
## ---------------------------------------------------------------------------

test_that("input validation rejects malformed arguments", {
  prior <- gmm(weights = 1, means = list(c(0, 0)), covariances = list(diag(2)))
  good_dyn <- list(A = diag(2), Q = 0.1 * diag(2))
  good_meas <- list(C = matrix(c(1, 0), nrow = 1L), R = matrix(0.2, 1L, 1L))
  Y <- matrix(stats::rnorm(4), ncol = 1L)

  expect_error(
    gmm_filter("not a gmm", good_dyn, good_meas, Y),
    regexp = "must be a"
  )
  expect_error(
    gmm_filter(prior, list(A = diag(3), Q = 0.1 * diag(3)), good_meas, Y),
    regexp = "state-to-state"
  )
  expect_error(
    gmm_filter(prior, list(Q = 0.1 * diag(2)), good_meas, Y),
    regexp = "at least an `A`"
  )
  expect_error(
    gmm_filter(prior, list(A = diag(2), Q = list(1)), good_meas, Y),
    regexp = "dynamics\\$Q"
  )
  expect_error(
    gmm_filter(prior, good_dyn, list(C = matrix(1, 1L, 1L), R = matrix(0.2, 1L, 1L)), Y),
    regexp = "columns"
  )
  expect_error(
    gmm_filter(prior, good_dyn, list(C = matrix(c(1, 0), nrow = 1L)), Y),
    regexp = "`C` matrix and an `R`"
  )
  expect_error(
    gmm_filter(prior, good_dyn, good_meas, matrix(stats::rnorm(4), ncol = 2L)),
    regexp = "to match"
  )
  expect_error(
    gmm_filter(prior, good_dyn,
               list(C = matrix(c(1, 0), nrow = 1L), R = 0.2 * diag(2)), Y),
    regexp = "numeric matrix or a"
  )
  expect_error(
    gmm_filter(prior, good_dyn,
               list(C = matrix(c(1, 0), nrow = 1L),
                    R = gmm(weights = 1, means = list(c(0, 0)),
                            covariances = list(diag(2)))), Y),
    regexp = "must live in"
  )
  expect_error(
    gmm_filter(prior, good_dyn, good_meas, "not numeric"),
    regexp = "numeric matrix or a length-n list"
  )
  expect_error(
    gmm_filter(prior, good_dyn, good_meas, list()),
    regexp = "at least one observation"
  )
  expect_error(
    gmm_filter(prior, good_dyn, good_meas, Y, k_max = 0L),
    regexp = "positive integer"
  )
  expect_error(
    gmm_filter(prior, good_dyn, good_meas, Y, k_max = -3L),
    regexp = "positive integer"
  )
  expect_error(
    gmm_filter(prior, good_dyn, good_meas, list(c(1), NA_real_)),
    regexp = "finite length"
  )
})

test_that("a measurement dimension that changes mid-series is rejected", {
  prior <- gmm(weights = 1, means = list(0), covariances = list(matrix(1)))
  ## One-dimensional observation at step 1, two-dimensional at step 2.
  meas_fun <- function(t) {
    if (t == 1L) {
      list(C = matrix(1, 1L, 1L), R = matrix(0.2, 1L, 1L))
    } else {
      list(C = matrix(c(1, 1), nrow = 2L), R = 0.2 * diag(2))
    }
  }
  expect_error(
    gmm_filter(prior, list(A = matrix(1), Q = matrix(0.1)),
               meas_fun, y = c(0.3, 0.4)),
    regexp = "stay constant"
  )
})

test_that("a mixture process noise on the wrong dimension is rejected", {
  prior <- gmm(weights = 1, means = list(c(0, 0)), covariances = list(diag(2)))
  bad_Q <- gmm(weights = 1, means = list(0), covariances = list(matrix(1)))
  expect_error(
    gmm_filter(prior, list(A = diag(2), Q = bad_Q),
               list(C = matrix(c(1, 0), nrow = 1L), R = matrix(0.2, 1L, 1L)),
               matrix(stats::rnorm(3), ncol = 1L)),
    regexp = "must live in"
  )
})
