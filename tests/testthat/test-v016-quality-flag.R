## The fit-quality flag fires on non-convergence, degeneracy, or an
## approximate KL above 0.3 nats, and not on low relative ESS.

.chain_target <- function() {
  energy <- function(X) {
    X <- matrix(X, ncol = 3)
    rowSums((X^2 - 1)^2) - 0.7 * (X[, 1] * X[, 2] + X[, 2] * X[, 3])
  }
  gmm_target(n_dim = 3L, log_density = function(X) -energy(X))
}

.chain_fit <- function(N, is_size, seed) {
  fit_kld_em(.chain_target(), N = N, proposal = is_uniform(3L, -3, 3),
             is_size = is_size, anneal = N > 1L, seed = seed,
             validation_size = 0L, support_warn = FALSE)
}

## The advisory is shown once per verb and session, so each call uses a new
## verb name.
.flag_calls <- new.env()
.flag_calls$n <- 0L
.is_flagged <- function(g) {
  .flag_calls$n <- .flag_calls$n + 1L
  flagged <- FALSE
  withCallingHandlers(
    .check_quality(g, paste0("flag_test_", .flag_calls$n)),
    proxymix_low_quality = function(m) {
      flagged <<- TRUE
      invokeRestart("muffleMessage")
    }
  )
  flagged
}

test_that("kld_approx tracks the KL divergence on a normalised target", {
  fit <- fit_proxymix(banana_target(), N = 3L, regime = "kld",
                      proposal = proposal_mvt(n_dim = 2L, mean = c(0, 0),
                                              sigma = 4 * diag(2), df = 5),
                      is_size = 2000L, max_iter = 60L, seed = 1L)
  q <- gmm_fit_quality(fit)
  ## Exact KL of this fit by quadrature is 0.016.
  expect_lt(q$kld_approx, 0.03)
  expect_gt(q$kld_approx, 0.005)
  expect_false(.is_flagged(fit))
})

test_that("a one-component fit to the chain target is flagged", {
  fit <- .chain_fit(N = 1L, is_size = 8000L, seed = 1L)
  q <- gmm_fit_quality(fit)
  expect_true(q$converged)
  expect_false(q$degenerate)
  ## Exact KL of this fit is about 0.86.
  expect_gt(q$kld_approx, 0.3)
  expect_true(.is_flagged(fit))
})

test_that("good chain fits are not flagged for low relative ESS", {
  fit <- .chain_fit(N = 8L, is_size = 8000L, seed = 1L)
  q <- gmm_fit_quality(fit)
  ## The uniform proposal keeps relative ESS near 0.05 however good the fit.
  expect_lt(q$ess_relative, 0.1)
  expect_lt(q$kld_approx, 0.3)
  expect_false(.is_flagged(fit))
  expect_true(all(gmm_independence_graph(fit) ==
                    matrix(c(0, 1, 0, 1, 0, 1, 0, 1, 0), 3L)))
})

test_that("flagging does not increase with is_size on the chain target", {
  seeds <- 1:5
  small <- vapply(seeds, function(s) .is_flagged(.chain_fit(8L, 6000L, s)),
                  logical(1L))
  large <- vapply(seeds, function(s) .is_flagged(.chain_fit(8L, 12000L, s)),
                  logical(1L))
  bad <- .chain_fit(N = 1L, is_size = 8000L, seed = 1L)
  expect_true(.is_flagged(bad) && .is_flagged(bad))
  expect_lte(sum(large), sum(small))
  expect_lte(sum(large), 1L)
})
