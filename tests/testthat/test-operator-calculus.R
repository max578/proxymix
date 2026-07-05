## Tests for the affine-Gaussian operator calculus — v0.3.0 graduation.
##
## Test obligations per `docs/design/operator_calculus_v0.3.md`:
##   A0 identity channel; A1 affine of moments; A2 sum-of-Gaussians;
##   O0 Kalman parity; O1 vanishing-evidence guard; O2 Bayes consistency;
##   G0 aggregate alias; M0 missing vs conditionalise;
##   C0 composition with marginalise.

mk_prior <- function() {
  gmm(
    weights     = c(0.4, 0.6),
    means       = list(c(-1, 0), c(1.5, 0.5)),
    covariances = list(diag(c(0.6, 0.8)), diag(c(0.7, 0.5)))
  )
}

test_that("A0: identity channel preserves the mixture", {
  g <- mk_prior()
  g_id <- gmm_affine(g, diag(2), b = c(0, 0), noise_cov = NULL)
  expect_equal(g_id@weights, g@weights)
  for (k in seq_len(gmm_n_components(g))) {
    expect_equal(g_id@means[[k]], g@means[[k]])
    ## Allow for the small numerical-hygiene ridge.
    expect_lt(max(abs(g_id@covariances[[k]] - g@covariances[[k]])),
              1e-5)
  }
})

test_that("A1: affine of moments matches by-hand", {
  g <- mk_prior()
  A <- matrix(c(1, 0,
                0, 1,
                1, -1), nrow = 3L, byrow = TRUE)
  b <- c(0.5, 0, -1)
  R <- 0.1 * diag(3)
  ga <- gmm_affine(g, A, b, noise_cov = R)
  for (k in seq_len(gmm_n_components(g))) {
    expect_equal(ga@means[[k]],
                 as.numeric(A %*% g@means[[k]] + b))
    expect_lt(max(abs(ga@covariances[[k]] -
                      (A %*% g@covariances[[k]] %*% t(A) + R))),
              1e-5)
  }
  expect_equal(ga@weights, g@weights)
})

test_that("A2: sum-of-Gaussians under identity-with-noise", {
  g <- gmm(weights = 1, means = list(c(0, 0)),
           covariances = list(diag(c(1, 2))))
  R <- diag(c(0.3, 0.7))
  ga <- gmm_affine(g, A = diag(2), noise_cov = R)
  expect_lt(max(abs(ga@covariances[[1L]] - (diag(c(1, 2)) + R))), 1e-5)
})

test_that("O0: Kalman parity vs single-component update", {
  g <- gmm(weights = 1, means = list(c(0, 0)),
           covariances = list(diag(c(1, 2))))
  H <- matrix(c(1, 0), nrow = 1L)
  R <- matrix(0.5, 1L, 1L)
  y <- 0.8
  g_post <- gmm_observe(g, A = H, y = y, noise_cov = R)

  ## Hand Kalman:
  S    <- diag(c(1, 2))
  Sx   <- H %*% S %*% t(H) + R
  gain <- S %*% t(H) %*% solve(Sx)
  mu_hat <- as.numeric(gain * y)
  S_hat  <- S - gain %*% H %*% S

  expect_lt(max(abs(mu_hat - g_post@means[[1L]])), 1e-5)
  expect_lt(max(abs(S_hat  - g_post@covariances[[1L]])), 1e-5)
})

test_that("O1: vanishing evidence triggers warn and no-update", {
  g <- mk_prior()
  H <- matrix(c(1, 0), nrow = 1L)
  ## Place y enormously far from any component on a vanishing noise.
  expect_warning(
    g_no <- gmm_observe(g, A = H, y = 1e6, noise_cov = matrix(1e-12, 1L, 1L)),
    regexp = "Marginal evidence"
  )
  expect_true(isTRUE(g_no@metadata$gmm_observe_no_update))
  expect_equal(g_no@weights, g@weights)
})

test_that("O2: sequential observations equal stacked observation", {
  g <- mk_prior()
  g_a <- gmm_observe(g,   A = matrix(c(1, 0), nrow = 1L),
                     y = 0.5, noise_cov = matrix(0.25, 1L, 1L))
  g_ab <- gmm_observe(g_a, A = matrix(c(0, 1), nrow = 1L),
                      y = 0.2, noise_cov = matrix(0.25, 1L, 1L))
  g_stack <- gmm_observe(
    g,
    A         = diag(2),
    y         = c(0.5, 0.2),
    noise_cov = 0.25 * diag(2)
  )
  expect_lt(max(abs(g_ab@weights - g_stack@weights)), 1e-5)
  for (k in seq_len(gmm_n_components(g))) {
    expect_lt(max(abs(g_ab@means[[k]] - g_stack@means[[k]])), 1e-5)
    expect_lt(max(abs(g_ab@covariances[[k]] -
                      g_stack@covariances[[k]])), 1e-4)
  }
})

test_that("O3: log marginal evidence metadata records the full log-sum", {
  g <- mk_prior()
  H <- matrix(c(1, 0), nrow = 1L)
  R <- matrix(0.25, 1L, 1L)
  y <- 0.5
  g_post <- gmm_observe(g, A = H, y = y, noise_cov = R)

  log_terms <- vapply(seq_len(gmm_n_components(g)), function(k) {
    Sx <- H %*% g@covariances[[k]] %*% t(H) + R
    log(g@weights[k]) + mvnfast::dmvn(
      matrix(y, nrow = 1L),
      mu = as.numeric(H %*% g@means[[k]]),
      sigma = Sx,
      log = TRUE
    )
  }, numeric(1L))
  expected <- max(log_terms) + log(sum(exp(log_terms - max(log_terms))))
  expect_equal(g_post@metadata$log_marginal_evidence, expected,
               tolerance = 1e-12)
})

test_that("G0: gmm_aggregate equals gmm_affine for an aggregation matrix", {
  g <- gmm(weights = c(0.5, 0.5),
           means = list(c(-1, 0, 1), c(1, 0, -1)),
           covariances = list(diag(3), diag(3)))
  A <- matrix(c(1, 1, 0,
                0, 0, 1), nrow = 2L, byrow = TRUE)
  g_agg  <- gmm_aggregate(g, A)
  g_aff  <- gmm_affine(g, A)
  expect_equal(g_agg@weights, g_aff@weights)
  for (k in seq_len(gmm_n_components(g))) {
    expect_lt(max(abs(g_agg@means[[k]] - g_aff@means[[k]])), 1e-12)
    expect_lt(max(abs(g_agg@covariances[[k]] - g_aff@covariances[[k]])),
              1e-12)
  }
  expect_identical(utils::tail(g_agg@metadata$provenance, 1L), "aggregate")
})

test_that("M0: gmm_missing equals gmm_conditionalise", {
  g <- mk_prior()
  g_m <- gmm_missing(g, observed = 2L, values = 0.5)
  g_c <- gmm_conditionalise(g, given = c(NA, 0.5))
  expect_lt(max(abs(g_m@weights - g_c@weights)), 1e-12)
  for (k in seq_len(gmm_n_components(g))) {
    expect_lt(max(abs(g_m@means[[k]] - g_c@means[[k]])), 1e-12)
    expect_lt(max(abs(g_m@covariances[[k]] - g_c@covariances[[k]])),
              1e-12)
  }
})

test_that("C0: gmm_marginalise commutes with gmm_affine on coordinate sub-channel", {
  ## A channel that only uses coordinates {1, 3} of a 3D mixture maps
  ## to (a, c) -> (a + c). Pushforward, then marginalise over the
  ## (trivially singleton) output, should equal marginalising the
  ## input over {2} and then pushing through the same A.
  g <- gmm(weights = c(0.4, 0.6),
           means = list(c(0, 0, 0), c(2, 1, -1)),
           covariances = list(diag(3), diag(3)))
  A_full <- matrix(c(1, 0, 1), nrow = 1L)   # uses coords 1 and 3
  A_sub  <- matrix(c(1, 1),    nrow = 1L)   # acts on coords (1, 3)
  push_then_keep <- gmm_affine(g, A_full)
  keep_then_push <- gmm_affine(
    gmm_marginalise(g, keep = c(1L, 3L)), A_sub
  )
  expect_lt(max(abs(push_then_keep@weights - keep_then_push@weights)),
            1e-12)
  for (k in seq_len(gmm_n_components(g))) {
    expect_lt(max(abs(push_then_keep@means[[k]] -
                      keep_then_push@means[[k]])), 1e-12)
    expect_lt(max(abs(push_then_keep@covariances[[k]] -
                      keep_then_push@covariances[[k]])), 1e-5)
  }
})

test_that("Input validation: A shape, b length, noise_cov shape", {
  g <- mk_prior()
  expect_error(gmm_affine(g, A = matrix(1:6, nrow = 2L, ncol = 3L)),
               regexp = "2 columns")
  expect_error(gmm_affine(g, A = diag(2), b = c(1, 2, 3)),
               regexp = "scalar or length-2")
  expect_error(gmm_affine(g, A = diag(2), noise_cov = diag(3)),
               regexp = "2 x 2 matrix")
  expect_error(gmm_observe(g, A = diag(2), y = 0.5,
                           noise_cov = diag(2)),
               regexp = "length-2 numeric vector")
  expect_error(gmm_observe(g, A = diag(2), y = c(0.5, 0.2),
                           noise_cov = NULL),
               regexp = "noise_cov.*required")
})

test_that("gmm_missing input validation: indices, lengths, uniqueness", {
  g <- mk_prior()
  expect_error(gmm_missing(g, observed = 3L, values = 0.5),
               regexp = "indices in 1..2")
  expect_error(gmm_missing(g, observed = c(1L, 2L), values = c(0.5, 0.2)),
               regexp = "leave at least one")
  expect_error(gmm_missing(g, observed = c(1L, 1L), values = c(0.5, 0.2)),
               regexp = "unique")
  expect_error(gmm_missing(g, observed = 1L, values = c(0.5, 0.2)),
               regexp = "length 1")
})
