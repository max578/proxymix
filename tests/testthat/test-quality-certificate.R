## The fit-quality certificate: stamped by every fitter, carried unchanged
## through the operator calculus, and read by downstream verbs.

.certificate_fields <- c("regime", "converged", "degenerate", "ess",
                         "ess_relative", "min_component_ess", "max_weight",
                         "support_fraction", "kld_final", "validation_gap")

test_that("all three regimes stamp a certificate with the shared schema", {
  x <- withr::with_seed(11, matrix(stats::rnorm(240), ncol = 2))
  tgt <- gmm_target_from_samples(x)
  fits <- list(
    moment = fit_moment_match(tgt, N = 1L),
    sample = fit_em_samples(tgt, N = 2L, max_iter = 25L, n_starts = 2L,
                            seed = 3L),
    kld = fit_kld_em(banana_target(), N = 2L, is_size = 1200L,
                     max_iter = 15L, seed = 1L)
  )
  for (nm in names(fits)) {
    q <- gmm_fit_quality(fits[[nm]])
    expect_named(q, .certificate_fields, info = nm)
    expect_equal(q$regime, unname(c(moment = "moment", sample = "sample",
                                    kld = "kld")[nm]))
    expect_false(isTRUE(q$degenerate), info = nm)
  }
  ## The kld certificate carries the live IS diagnostics and the held-out
  ## validation gap (validation is on by default now).
  qk <- gmm_fit_quality(fits$kld)
  expect_true(is.finite(qk$ess))
  expect_true(is.finite(qk$validation_gap))
  expect_true(is.finite(qk$min_component_ess))
})

test_that("the certificate and provenance survive an operator chain", {
  fit <- fit_kld_em(banana_target(), N = 2L, is_size = 1200L,
                    max_iter = 15L, seed = 1L)
  chained <- gmm_marginalise(
    gmm_affine(fit, matrix(c(1, 0.3, -0.2, 1), 2), b = c(1, -1)),
    keep = 1L
  )
  expect_identical(gmm_fit_quality(chained), gmm_fit_quality(fit))
  expect_identical(chained@metadata$provenance, c("affine", "marginalise"))
  ## Names no longer nest without bound.
  expect_lt(nchar(chained@name), 60L)
})

test_that("a plain mixture has no certificate and passes checks silently", {
  g <- gmm(weights = c(0.5, 0.5), means = list(c(-1, 0), c(1, 0)),
           covariances = list(diag(2), diag(2)))
  expect_null(gmm_fit_quality(g))
  expect_silent(hush <- gmm_entropy(g, order = "renyi2"))
})

test_that("ESS collapse flags the fit as degenerate with a classed warning", {
  far <- gmm_target(
    n_dim = 2L,
    log_density = function(x) {
      if (is.null(dim(x))) x <- matrix(x, ncol = 2L)
      -0.5 * rowSums((x - 50)^2) / 0.01 - log(2 * pi * 0.01)
    },
    normalised = TRUE, name = "far_spike"
  )
  q <- is_mvn(2L, cov = diag(2))
  expect_warning(
    fit <- suppressWarnings(
      fit_kld_em(far, N = 2L, proposal = q, is_size = 400L,
                 max_iter = 5L, seed = 1L, validation_size = 0L),
      classes = c("proxymix_support", "proxymix_nonmonotone")
    ),
    class = "proxymix_low_ess"
  )
  expect_false(fit@converged)
  expect_true(fit@diagnostics$degenerate)
  expect_true(gmm_fit_quality(fit)$degenerate)

  ## The abort path refuses to return the degenerate fit at all.
  expect_error(
    suppressWarnings(
      fit_kld_em(far, N = 2L, proposal = q, is_size = 400L,
                 max_iter = 5L, seed = 1L, validation_size = 0L,
                 on_low_ess = "abort")
    ),
    class = "proxymix_degenerate_fit"
  )

  ## A downstream verb reading the flagged certificate raises the one-shot
  ## advisory.
  expect_message(
    gmm_entropy(fit, order = "renyi2"),
    class = "proxymix_low_quality"
  )
})

test_that("fit_kld_em is reproducible end-to-end under a seed", {
  tgt <- banana_target()
  q <- is_mvt(n_dim = 2L, mean = c(0, 0), sigma = 4 * diag(2), df = 5)
  ## Different ambient RNG states; the seed alone must pin the fit.
  f1 <- withr::with_seed(1, fit_kld_em(tgt, N = 2L, proposal = q,
                                       is_size = 1200L, max_iter = 15L,
                                       seed = 7L))
  f2 <- withr::with_seed(999, fit_kld_em(tgt, N = 2L, proposal = q,
                                         is_size = 1200L, max_iter = 15L,
                                         seed = 7L))
  expect_equal(f1@weights, f2@weights, tolerance = 1e-14)
  expect_equal(f1@means, f2@means, tolerance = 1e-14)
  expect_equal(f1@covariances, f2@covariances, tolerance = 1e-14)
  expect_equal(f1@diagnostics$kld_final, f2@diagnostics$kld_final,
               tolerance = 1e-14)
})

test_that("the dimension guard fires at the core fitter", {
  p <- 12L
  tgt <- gmm_target(
    n_dim = p,
    log_density = function(x) {
      if (is.null(dim(x))) x <- matrix(x, ncol = p)
      -0.5 * rowSums(x^2) - 0.5 * p * log(2 * pi)
    },
    normalised = TRUE, name = "iso12"
  )
  q <- is_mvn(p, cov = 4 * diag(p))
  expect_warning(
    fit <- suppressWarnings(
      fit_kld_em(tgt, N = 1L, proposal = q, is_size = 2000L,
                 max_iter = 10L, seed = 1L, validation_size = 0L),
      classes = c("proxymix_low_ess", "proxymix_nonmonotone")
    ),
    class = "proxymix_high_dimension"
  )
  ## Whatever happened, the certificate reports it rather than hiding it.
  q12 <- gmm_fit_quality(fit)
  expect_true(is.finite(q12$ess_relative))
})
