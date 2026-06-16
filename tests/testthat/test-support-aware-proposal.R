## Support-aware importance-proposal selection (regime-(iii) robustness layer).
##
## Locks out the NaN-weight failure mode on bounded and one-sided targets
## (regress spec criterion S10): the default multivariate-t proposal places
## importance mass where a compact or one-sided target's log-density is -Inf,
## which produced non-finite weights. A declared `support` now selects a
## matched uniform proposal automatically.

# ---- epanechnikov_target fixture -------------------------------------------

test_that("epanechnikov_target is a normalised compact-support target", {
  e <- epanechnikov_target()
  expect_true(S7::S7_inherits(e, gmm_target))
  expect_equal(e@n_dim, 1L)
  expect_equal(e@support$lower, -1)
  expect_equal(e@support$upper, 1)

  ## Integrates to one over its support.
  val <- stats::integrate(
    function(z) exp(e@log_density(matrix(z, ncol = 1L))), -1, 1
  )$value
  expect_equal(val, 1, tolerance = 1e-6)

  ## -Inf outside, finite strictly inside, no NaN warning.
  expect_equal(e@log_density(matrix(c(-1.5, 1.5), ncol = 1L)), c(-Inf, -Inf))
  expect_true(all(is.finite(e@log_density(matrix(c(-0.9, 0, 0.9), ncol = 1L)))))
  expect_no_warning(e@log_density(matrix(c(-2, 0, 2), ncol = 1L)))
})

test_that("epanechnikov_target Devroye samples match the Epanechnikov law", {
  es <- epanechnikov_target(with_samples = TRUE, n = 5000L, seed = 1L)
  expect_equal(dim(es@samples), c(5000L, 1L))
  expect_true(all(abs(es@samples) <= 1))
  ## Variance of the Epanechnikov kernel on [-1, 1] is 1/5.
  expect_equal(stats::var(as.numeric(es@samples)), 0.2, tolerance = 0.03)
})

# ---- support validation ----------------------------------------------------

test_that("gmm_target validates the support declaration", {
  lf <- function(x) rep(-1, nrow(x))
  expect_no_error(
    gmm_target(n_dim = 1L, log_density = lf, support = list(lower = 0, upper = 1))
  )
  expect_error(
    gmm_target(n_dim = 1L, log_density = lf, support = list(0, 1)),
    "lower.*upper"
  )
  expect_error(
    gmm_target(n_dim = 1L, log_density = lf, support = list(lower = 1, upper = 0)),
    "exceed"
  )
  expect_error(
    gmm_target(n_dim = 2L, log_density = lf,
               support = list(lower = c(0, 0, 0), upper = c(1, 1, 1))),
    "length 1 or 2"
  )
  expect_error(
    gmm_target(n_dim = 1L, log_density = lf,
               support = list(lower = NA_real_, upper = 1)),
    "NA"
  )
})

# ---- selector --------------------------------------------------------------

test_that("support-matched proposal is NULL for unbounded, uniform for compact", {
  ## Unbounded target keeps the heavy-tailed default (selector returns NULL).
  expect_null(.support_matched_proposal(banana_target()))

  ## Compact target selects an inset uniform.
  q <- .support_matched_proposal(epanechnikov_target())
  expect_true(S7::S7_inherits(q, is_proposal))
  expect_match(q@name, "support-matched")
  expect_gt(q@metadata$lower, -1)
  expect_lt(q@metadata$upper, 1)
})

test_that("one-sided support needs samples and yields finite working bounds", {
  lf <- function(x) stats::dgamma(x[, 1L], shape = 2, rate = 1, log = TRUE)

  ## No samples -> informative abort, never a silent crashing fallback.
  t_nodata <- gmm_target(n_dim = 1L, log_density = lf,
                         support = list(lower = 0, upper = Inf))
  expect_error(.support_matched_proposal(t_nodata), "without samples")

  ## With samples -> a finite uniform on roughly [0, max + margin].
  smp <- withr::with_seed(1L,
    matrix(stats::rgamma(2000L, shape = 2, rate = 1), ncol = 1L))
  t_data <- gmm_target(n_dim = 1L, log_density = lf, samples = smp,
                       support = list(lower = 0, upper = Inf),
                       normalised = TRUE, log_normalizer = 0)
  q <- .support_matched_proposal(t_data)
  expect_true(S7::S7_inherits(q, is_proposal))
  expect_true(is.finite(q@metadata$lower) && is.finite(q@metadata$upper))
  expect_gte(q@metadata$lower, 0)
})

# ---- end-to-end: no NaN weights on the bounded battery (S10) ---------------

test_that("fit_kld_em auto-selects a uniform on a compact target, no NaN weights", {
  e <- epanechnikov_target()
  expect_message(
    fit <- fit_kld_em(e, N = 2L, is_size = 3000L, max_iter = 20L, seed = 1L),
    "support-matched"
  )
  expect_equal(fit@diagnostics$proposal_name, "is_uniform[support-matched]")
  expect_true(is.finite(fit@diagnostics$kld_final))
  expect_gt(fit@diagnostics$ess, 50)
  expect_equal(fit@diagnostics$support_fraction, 1, tolerance = 1e-9)
  expect_false(any(vapply(fit@means, anyNA, logical(1L))))
})

test_that("fit_kld_em handles a one-sided Gamma target via the auto uniform", {
  lf <- function(x) stats::dgamma(x[, 1L], shape = 3, rate = 1, log = TRUE)
  smp <- withr::with_seed(2L,
    matrix(stats::rgamma(3000L, shape = 3, rate = 1), ncol = 1L))
  g <- gmm_target(n_dim = 1L, log_density = lf, samples = smp,
                  support = list(lower = 0, upper = Inf),
                  normalised = TRUE, log_normalizer = 0, name = "gamma")
  expect_message(
    fit <- fit_kld_em(g, N = 2L, is_size = 4000L, max_iter = 30L, seed = 1L),
    "support-matched"
  )
  expect_true(is.finite(fit@diagnostics$kld_final))
  expect_gt(fit@diagnostics$ess, 50)
  expect_false(any(vapply(fit@means, anyNA, logical(1L))))
})

# ---- backward compatibility ------------------------------------------------

test_that("a NULL-support target keeps the heavy-tailed default (no message)", {
  expect_no_message(
    fit <- suppressWarnings(
      fit_kld_em(banana_target(), N = 2L, is_size = 1500L, max_iter = 10L, seed = 1L)
    )
  )
  expect_equal(fit@diagnostics$proposal_name, "is_mvt")
})
