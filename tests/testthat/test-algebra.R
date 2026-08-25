## Oracles for the completed mixture algebra: moments, product,
## convolution, mixing, and the one-dimensional p/q functions, plus the
## evidence estimator. Reference quantities are hand-derived or computed
## by quadrature / sampling cross-checks independent of the implementations.

test_that("gmm_mean / gmm_cov match the exact mixture-moment formulas", {
  w <- c(0.3, 0.7)
  m1 <- c(-1, 0.5); m2 <- c(2, -1)
  S1 <- matrix(c(1, 0.3, 0.3, 0.8), 2)
  S2 <- matrix(c(0.6, -0.1, -0.1, 1.2), 2)
  g <- gmm(weights = w, means = list(m1, m2), covariances = list(S1, S2))
  mu <- w[1] * m1 + w[2] * m2
  Sg <- w[1] * (S1 + tcrossprod(m1)) + w[2] * (S2 + tcrossprod(m2)) -
    tcrossprod(mu)
  expect_equal(gmm_mean(g), mu, tolerance = 1e-12)
  expect_equal(gmm_cov(g), Sg, tolerance = 1e-12)
  ## And against large-sample moments.
  x <- withr::with_seed(7, rgmm(200000L, g))
  expect_lt(max(abs(colMeans(x) - mu)), 0.02)
  expect_lt(max(abs(stats::cov(x) - Sg)), 0.05)
})

test_that("gmm_product satisfies the exact density identity p*q = Z * product-density", {
  g1 <- gmm(weights = c(0.4, 0.6), means = list(c(-1, 0), c(1.5, 1)),
            covariances = list(diag(2), matrix(c(0.8, 0.2, 0.2, 1), 2)))
  g2 <- gmm(weights = c(0.5, 0.5), means = list(c(0, 0.5), c(-0.5, -1)),
            covariances = list(0.7 * diag(2), 1.3 * diag(2)))
  prod12 <- gmm_product(g1, g2)
  Z <- exp(prod12@metadata$log_integral)
  xs <- withr::with_seed(3, matrix(stats::rnorm(20), ncol = 2))
  ## Pointwise: dgmm(x, g1) * dgmm(x, g2) == Z * dgmm(x, product).
  expect_equal(dgmm(xs, g1) * dgmm(xs, g2), Z * dgmm(xs, prod12),
               tolerance = 1e-10)
  ## Commutativity at the density level.
  prod21 <- gmm_product(g2, g1)
  expect_equal(dgmm(xs, prod12), dgmm(xs, prod21), tolerance = 1e-10)
  expect_equal(prod12@metadata$log_integral, prod21@metadata$log_integral,
               tolerance = 1e-10)
})

test_that("gmm_product at K = 1 matches the hand conjugate-update formulas", {
  m0 <- c(1, -1); S0 <- matrix(c(2, 0.5, 0.5, 1), 2)
  m1 <- c(0, 0.5); S1 <- matrix(c(1, 0, 0, 0.5), 2)
  g0 <- gmm(weights = 1, means = list(m0), covariances = list(S0))
  g1 <- gmm(weights = 1, means = list(m1), covariances = list(S1))
  pr <- gmm_product(g0, g1)
  P <- solve(solve(S0) + solve(S1))
  mu <- as.numeric(P %*% (solve(S0) %*% m0 + solve(S1) %*% m1))
  expect_equal(pr@covariances[[1L]], P, tolerance = 1e-10)
  expect_equal(pr@means[[1L]], mu, tolerance = 1e-10)
})

test_that("gmm_product agrees with gmm_observe on the identity channel", {
  ## Conditioning g on y = x + eps, eps ~ N(0, R), equals the product of g
  ## with the Gaussian likelihood N(y; x, R) read as a mixture in x.
  g <- gmm(weights = c(0.4, 0.6), means = list(c(-1, 0), c(1, 1)),
           covariances = list(diag(2), 0.5 * diag(2)))
  y <- c(0.3, -0.2)
  R <- 0.4 * diag(2)
  obs <- gmm_observe(g, A = diag(2), y = y, noise_cov = R, ridge_eps = 0)
  lik <- gmm(weights = 1, means = list(y), covariances = list(R))
  pr <- gmm_product(g, lik)
  xs <- withr::with_seed(5, matrix(stats::rnorm(16), ncol = 2))
  expect_equal(dgmm(xs, obs), dgmm(xs, pr), tolerance = 1e-8)
})

test_that("gmm_convolve is the distribution of the sum of independent draws", {
  g1 <- gmm(weights = c(0.5, 0.5), means = list(-1, 1),
            covariances = list(matrix(0.5), matrix(0.5)))
  g2 <- gmm(weights = c(0.3, 0.7), means = list(0, 2),
            covariances = list(matrix(1), matrix(0.3)))
  cv <- gmm_convolve(g1, g2)
  ## Exact moment additivity.
  expect_equal(gmm_mean(cv), gmm_mean(g1) + gmm_mean(g2), tolerance = 1e-12)
  expect_equal(gmm_cov(cv), gmm_cov(g1) + gmm_cov(g2), tolerance = 1e-12)
  ## KS of simulated sums against the convolution's own CDF.
  s <- withr::with_seed(11, as.numeric(rgmm(4000L, g1)) +
                          as.numeric(rgmm(4000L, g2)))
  ks <- suppressWarnings(stats::ks.test(s, function(q) pgmm(q, cv)))
  expect_gt(ks$p.value, 0.01)
})

test_that("gmm_mix averages densities with the supplied weights", {
  g1 <- gmm(weights = 1, means = list(-1), covariances = list(matrix(1)))
  g2 <- gmm(weights = 1, means = list(2), covariances = list(matrix(0.5)))
  mx <- gmm_mix(list(g1, g2), weights = c(0.7, 0.3))
  xs <- matrix(seq(-3, 4, by = 0.5), ncol = 1L)
  expect_equal(dgmm(xs, mx), 0.7 * dgmm(xs, g1) + 0.3 * dgmm(xs, g2),
               tolerance = 1e-12)
})

test_that("pgmm matches quadrature of dgmm and qgmm inverts it", {
  g <- gmm(weights = c(0.4, 0.6), means = list(-2, 1),
           covariances = list(matrix(0.5), matrix(1)))
  for (q in c(-2.5, 0, 1.5)) {
    num <- stats::integrate(function(x) dgmm(matrix(x), g), -20, q,
                            rel.tol = 1e-10)$value
    expect_equal(pgmm(q, g), num, tolerance = 1e-8)
  }
  ps <- c(0.05, 0.25, 0.5, 0.9, 0.99)
  expect_equal(pgmm(qgmm(ps, g), g), ps, tolerance = 1e-8)
  ## Upper tail and errors.
  expect_equal(pgmm(0, g, lower.tail = FALSE), 1 - pgmm(0, g))
  g2d <- gmm(weights = 1, means = list(c(0, 0)), covariances = list(diag(2)))
  expect_error(pgmm(0, g2d), "one-dimensional")
  expect_error(qgmm(1.2, g), "inside")
})

test_that("gmm_evidence recovers a known normalising constant", {
  offsets <- c(0, 3)
  for (off in offsets) {
    tgt <- gmm_target(
      n_dim = 2L,
      log_density = local({
        o <- off
        function(x) {
          if (is.null(dim(x))) x <- matrix(x, ncol = 2L)
          -0.5 * rowSums(x^2) - log(2 * pi) + o
        }
      }),
      normalised = off == 0, name = sprintf("gauss_offset_%g", off)
    )
    fit <- fit_kld_em(tgt, N = 1L, is_size = 3000L, max_iter = 50L,
                      seed = 1L, validation_size = 0L,
                      proposal = is_mvt(2L, sigma = 4 * diag(2), df = 5))
    ev <- gmm_evidence(fit, n = 4000L, seed = 2L)
    expect_lt(abs(ev$log_z - off), 5 * ev$se_log_z + 0.01)
    expect_false(ev$flagged)
  }
})

test_that("gmm_evidence recovers log Z for a genuinely multimodal target", {
  gtrue <- gmm(weights = c(0.35, 0.65), means = list(c(-2, 0), c(2, 1)),
               covariances = list(0.6 * diag(2), diag(2)))
  const <- -1.7
  tgt <- gmm_target(
    n_dim = 2L,
    log_density = function(x) {
      if (is.null(dim(x))) x <- matrix(x, ncol = 2L)
      dgmm(x, gtrue, log = TRUE) + const
    },
    normalised = FALSE, name = "mm_offset"
  )
  fit <- fit_kld_em(tgt, N = 2L, is_size = 4000L, max_iter = 60L, seed = 3L,
                    proposal = is_mvt(2L, sigma = 9 * diag(2), df = 5))
  ev <- gmm_evidence(fit, n = 6000L, seed = 4L)
  expect_lt(abs(ev$log_z - const), 5 * ev$se_log_z + 0.02)
})

test_that("gmm_evidence flags a proposal with too-light tails", {
  ## Target: heavy-ish mixture; proxy deliberately too narrow.
  tgt <- gmm_target(
    n_dim = 1L,
    log_density = function(x) {
      xx <- if (is.matrix(x)) x[, 1L] else x
      stats::dnorm(xx, 0, 3, log = TRUE)
    },
    normalised = TRUE, name = "wide"
  )
  narrow <- gmm_fit(
    weights = 1, means = list(0), covariances = list(matrix(0.25)),
    target = tgt, regime = "kld", diagnostics = list(),
    converged = TRUE, iterations = 1L, call = quote(f()),
    name = "narrow"
  )
  expect_warning(ev <- gmm_evidence(narrow, n = 3000L, seed = 5L),
                 class = "proxymix_heavy_tail")
  expect_true(ev$flagged)
})

test_that("the accessors return the component parameters verbatim", {
  g <- gmm(weights = c(0.3, 0.7), means = list(-1, 2),
           covariances = list(matrix(1), matrix(0.5)))
  expect_identical(gmm_weights(g), g@weights)
  expect_identical(gmm_means(g), g@means)
  expect_identical(gmm_covariances(g), g@covariances)
})

test_that("tidiers produce the component table and the fit summary", {
  skip_if_not_installed("generics")
  fit <- fit_kld_em(banana_target(), N = 2L, is_size = 1000L,
                    max_iter = 10L, seed = 1L, validation_size = 0L)
  td <- generics::tidy(fit)
  expect_equal(nrow(td), 2L)
  expect_true(all(c("component", "weight", "mean_1", "mean_2",
                    "var_1", "var_2") %in% names(td)))
  expect_equal(sum(td$weight), 1, tolerance = 1e-12)
  gl <- generics::glance(fit)
  expect_equal(nrow(gl), 1L)
  expect_equal(gl$regime, "kld")
  expect_true(is.finite(gl$ess))
})

test_that("a sealed mechanism rejects unknown gate types", {
  fake <- structure(list(type = "shadow"), class = "proxymix_gate")
  x <- withr::with_seed(9, {
    m <- matrix(stats::rnorm(60), ncol = 2, dimnames = list(NULL, c("a", "b")))
    m[c(3, 7), 2] <- NA
    m
  })
  expect_error(gmm_impute(x, mechanism = fake, m = 2L, seed = 1L),
               "sealed")
})

# ---- quality-certificate composition across binary/n-ary operators (F2) ----

.mk_cert_gmm <- function(nm, degenerate, converged, ess_rel) {
  gmm(weights = 1, means = list(0), covariances = list(matrix(1)), name = nm,
      metadata = list(
        quality = list(regime = "kld", converged = converged,
                       degenerate = degenerate, ess = NA_real_,
                       ess_relative = ess_rel, min_component_ess = NA_real_,
                       max_weight = 0.1, support_fraction = 1,
                       kld_final = 0.01, validation_gap = NA_real_),
        provenance = nm))
}

test_that("binary/n-ary operators combine operand certificates conservatively", {
  good <- .mk_cert_gmm("good", degenerate = FALSE, converged = TRUE, ess_rel = 0.6)
  bad  <- .mk_cert_gmm("bad",  degenerate = TRUE,  converged = FALSE, ess_rel = 0.02)

  ops <- list(
    product_gb = function() gmm_product(good, bad),
    product_bg = function() gmm_product(bad, good),
    convolve   = function() gmm_convolve(good, bad),
    mix        = function() gmm_mix(list(good, bad))
  )
  for (nm in names(ops)) {
    q <- gmm_fit_quality(ops[[nm]]())
    expect_false(is.null(q), info = nm)
    expect_identical(q$regime, "composite", info = nm)
    expect_true(q$degenerate, info = nm)        # any operand degenerate
    expect_false(q$converged, info = nm)        # not all operands converged
    expect_equal(q$ess_relative, 0.02, info = nm)  # worst (min) relative ESS
    expect_length(q$quality_sources, 2L)
  }
})

test_that("mixing / multiplying plain (uncertified) mixtures makes no false certificate", {
  a <- gmm(weights = 1, means = list(-1), covariances = list(matrix(1)))
  b <- gmm(weights = 1, means = list(2),  covariances = list(matrix(0.5)))
  expect_null(gmm_fit_quality(gmm_mix(list(a, b))))
  expect_null(gmm_fit_quality(gmm_product(a, b)))
})
