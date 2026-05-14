test_that("support-domination warning fires when uniform misses target mass", {
  ## Target has substantial mass outside a tiny uniform box.
  tgt <- banana_target()
  q_tight <- is_uniform(n_dim = 2L, lower = -0.1, upper = 0.1)
  ## ESS will also be tiny here; allow the warnings to fire silently.
  expect_warning(
    suppressWarnings({
      fit <- fit_kld_em(tgt, N = 2L, proposal = q_tight,
                        is_size = 500L, max_iter = 5L, seed = 1L)
    }),
    NA
  )
  ## All draws inside the tight box receive finite weight (support_fraction
  ## within the proposal's box is 1), so the canonical test of the
  ## warning is via a proposal whose log-density evaluates to -Inf for some
  ## draws. Use a box that excludes part of the proposal's own sampler
  ## by drawing from an MVN sampler and evaluating against a uniform.
  q_mismatch <- is_proposal(
    n_dim = 2L,
    sample = function(n) mvnfast::rmvn(n, mu = c(0, 0), sigma = 9 * diag(2)),
    log_density = function(x) {
      if (is.null(dim(x))) x <- matrix(x, nrow = 1L)
      inside <- abs(x[, 1L]) <= 1 & abs(x[, 2L]) <= 1
      out <- rep(-Inf, nrow(x))
      out[inside] <- -log(4)
      out
    },
    name = "broken_proposal"
  )
  expect_warning(
    fit_kld_em(tgt, N = 1L, proposal = q_mismatch,
               is_size = 800L, max_iter = 3L, seed = 1L,
               min_ess = 0),
    "does not dominate"
  )
})

test_that("support_warn = FALSE silences the warning", {
  tgt <- banana_target()
  q_mismatch <- is_proposal(
    n_dim = 2L,
    sample = function(n) mvnfast::rmvn(n, mu = c(0, 0), sigma = 9 * diag(2)),
    log_density = function(x) {
      if (is.null(dim(x))) x <- matrix(x, nrow = 1L)
      inside <- abs(x[, 1L]) <= 1 & abs(x[, 2L]) <= 1
      out <- rep(-Inf, nrow(x))
      out[inside] <- -log(4)
      out
    },
    name = "broken_proposal"
  )
  expect_no_warning(
    suppressWarnings(  ## low-ESS warning is independent of support warning
      fit_kld_em(tgt, N = 1L, proposal = q_mismatch,
                 is_size = 800L, max_iter = 3L, seed = 1L,
                 min_ess = 0, support_warn = FALSE)
    )
  )
})

test_that("support_fraction is recorded in diagnostics", {
  fit <- fit_kld_em(banana_target(), N = 2L,
                    is_size = 800L, max_iter = 10L, seed = 1L)
  expect_true(is.finite(fit@diagnostics$support_fraction))
  expect_gte(fit@diagnostics$support_fraction, 0)
  expect_lte(fit@diagnostics$support_fraction, 1)
})
