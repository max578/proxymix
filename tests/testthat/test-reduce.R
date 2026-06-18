## G1: gmm_reduce -- moment-preserving greedy mixture reduction.
##
## Independent oracles: the analytic mixture moments (computed directly, not via
## the package), an independently hand-coded Runnalls (2007) reduction, the
## planted redundant structure, and the closed-form Cauchy-Schwarz divergence
## (validated in the entropy layer).

## Analytic global mean and covariance of a mixture (independent of gmm_reduce).
.mix_moments <- function(g) {
  w <- g@weights
  m <- g@means
  S <- g@covariances
  mu <- Reduce(`+`, Map(function(wk, mk) wk * mk, w, m))
  m2 <- Reduce(`+`, Map(function(wk, mk, Sk) wk * (Sk + tcrossprod(mk)), w, m, S))
  list(mean = as.numeric(mu), cov = m2 - tcrossprod(mu))
}

## Independent hand-coded Runnalls greedy reduction (determinant(), no ridge).
.hand_runnalls <- function(g, k_max) {
  w <- g@weights; m <- g@means; S <- g@covariances; K <- length(w)
  ld <- function(M) as.numeric(determinant(M, logarithm = TRUE)$modulus)
  mg2 <- function(wi, mi, Si, wj, mj, Sj) {
    ww <- wi + wj; mu <- (wi * mi + wj * mj) / ww; d <- mi - mj
    list(w = ww, mu = as.numeric(mu),
         S = (wi * Si + wj * Sj) / ww + (wi * wj / ww^2) * tcrossprod(d))
  }
  while (K > k_max) {
    best <- Inf; bi <- NA; bj <- NA; bm <- NULL
    for (i in seq_len(K - 1L)) for (j in (i + 1L):K) {
      mg <- mg2(w[i], m[[i]], S[[i]], w[j], m[[j]], S[[j]])
      b <- 0.5 * ((w[i] + w[j]) * ld(mg$S) - w[i] * ld(S[[i]]) - w[j] * ld(S[[j]]))
      if (b < best) { best <- b; bi <- i; bj <- j; bm <- mg }
    }
    w[bi] <- bm$w; m[[bi]] <- bm$mu; S[[bi]] <- bm$S
    w <- w[-bj]; m <- m[-bj]; S <- S[-bj]; K <- K - 1L
  }
  gmm(weights = w / sum(w), means = m, covariances = S)
}

## Canonicalised parameter distance between two mixtures (order-invariant).
.mix_dist <- function(a, b) {
  ca <- order(vapply(a@means, function(x) x[1L], numeric(1L)), a@weights)
  cb <- order(vapply(b@means, function(x) x[1L], numeric(1L)), b@weights)
  max(c(max(abs(a@weights[ca] - b@weights[cb])),
        max(abs(unlist(a@means[ca]) - unlist(b@means[cb]))),
        max(abs(unlist(a@covariances[ca]) - unlist(b@covariances[cb])))))
}

.random_mixture <- function(K, p = 2L, seed = 7L) {
  withr::with_seed(seed, {
    u <- stats::runif(K)
    gmm(
      weights = u / sum(u),
      means = lapply(seq_len(K), function(.) stats::rnorm(p, sd = 3)),
      covariances = lapply(seq_len(K), function(.) {
        A <- matrix(stats::rnorm(p * p), p, p)
        crossprod(A) / p + 0.4 * diag(p)
      })
    )
  })
}

# ---- GS1: moment preservation ---------------------------------------------

test_that("reduction preserves the global mean and covariance exactly", {
  g <- .random_mixture(7L)
  o <- .mix_moments(g)
  for (km in c(1L, 2L, 3L, 4L)) {
    r <- gmm_reduce(g, km)
    rm_ <- .mix_moments(r)
    expect_equal(rm_$mean, o$mean, tolerance = 1e-10)
    expect_equal(rm_$cov, o$cov, tolerance = 1e-10)
    expect_lte(gmm_n_components(r), km)
  }
  ## the cs cost preserves moments too
  rc <- gmm_reduce(g, 3L, cost = "cs")
  expect_equal(.mix_moments(rc)$mean, o$mean, tolerance = 1e-10)
  expect_equal(.mix_moments(rc)$cov, o$cov, tolerance = 1e-10)
})

test_that("reducing to one component returns the moment-matched Gaussian", {
  g <- .random_mixture(6L, seed = 3L)
  o <- .mix_moments(g)
  r1 <- gmm_reduce(g, 1L)
  expect_equal(gmm_n_components(r1), 1L)
  expect_equal(r1@means[[1L]], o$mean, tolerance = 1e-10)
  expect_equal(r1@covariances[[1L]], o$cov, tolerance = 1e-10)
})

# ---- GS2: identity when the budget is not binding -------------------------

test_that("k_max at or above the component count leaves the mixture unchanged", {
  g <- .random_mixture(5L)
  expect_equal(.mix_dist(gmm_reduce(g, 5L), g), 0)
  expect_equal(gmm_n_components(gmm_reduce(g, 9L)), 5L)
  expect_equal(.mix_dist(gmm_reduce(g, 9L), g), 0)
})

# ---- Runnalls parity against an independent reference ---------------------

test_that("the kl-cost reduction matches an independent Runnalls reduction", {
  g <- .random_mixture(8L, seed = 11L)
  for (km in c(2L, 3L, 4L, 5L)) {
    pm <- gmm_reduce(g, km, cost = "kl", ridge_eps = 0)
    hr <- .hand_runnalls(g, km)
    expect_equal(.mix_dist(pm, hr), 0, tolerance = 1e-10)
  }
})

# ---- GS3: fidelity on a redundant mixture ---------------------------------

test_that("a redundant mixture reduces near-losslessly and recovers its clusters", {
  ## three well-separated clusters, each split into two near-identical components
  gr <- gmm(
    weights = rep(1 / 6, 6),
    means = list(c(-5, 0), c(-5, 0.15), c(5, 0), c(5.1, -0.1), c(0, 6), c(0.1, 6.1)),
    covariances = rep(list(0.5 * diag(2)), 6)
  )
  rk <- gmm_reduce(gr, 3L, cost = "kl")
  rc <- gmm_reduce(gr, 3L, cost = "cs")
  expect_equal(gmm_n_components(rk), 3L)
  ## merging near-duplicates is near-lossless: tiny divergence from the original
  dk <- gmm_divergence(gr, rk, type = "cs")
  dc <- gmm_divergence(gr, rc, type = "cs")
  expect_lt(dk, 1e-3)
  expect_lt(dc, 1e-3)
  ## the three recovered centres sit on the three planted clusters
  cx <- sort(vapply(rk@means, function(m) m[1L], numeric(1L)))
  expect_equal(cx, c(-5, 0.05, 5.05), tolerance = 0.2)
})

test_that("reduction divergence is monotone in the component budget", {
  g <- .random_mixture(7L, seed = 21L)
  divs <- vapply(5:1, function(km) {
    gmm_divergence(g, gmm_reduce(g, as.integer(km)), type = "cs")
  }, numeric(1L))
  ## smaller budget -> larger (or equal) divergence from the original
  expect_true(all(diff(divs) >= -1e-9))
  expect_gt(divs[length(divs)], divs[1L])
})

# ---- reduced object is a valid mixture ------------------------------------

test_that("the reduction returns a valid gmm with normalised weights", {
  g <- .random_mixture(6L)
  r <- gmm_reduce(g, 3L)
  expect_true(S7::S7_inherits(r, gmm))
  expect_equal(sum(r@weights), 1, tolerance = 1e-12)
  expect_true(all(r@weights >= 0))
})

test_that("gmm_reduce accepts a gmm_fit and preserves its moments", {
  x <- withr::with_seed(1L, rbind(
    matrix(stats::rnorm(200, -4), ncol = 2L),
    matrix(stats::rnorm(200, 4), ncol = 2L)
  ))
  fit <- fit_em_samples(gmm_target_from_samples(x), N = 5L, max_iter = 60L)
  o <- .mix_moments(fit)
  r <- gmm_reduce(fit, 2L)
  expect_true(S7::S7_inherits(r, gmm))
  expect_equal(gmm_n_components(r), 2L)
  expect_equal(.mix_moments(r)$mean, o$mean, tolerance = 1e-9)
  expect_equal(.mix_moments(r)$cov, o$cov, tolerance = 1e-9)
})

# ---- input validation -----------------------------------------------------

test_that("gmm_reduce validates its inputs", {
  g <- .random_mixture(4L)
  expect_error(gmm_reduce(g, 0L), "positive integer")
  expect_error(gmm_reduce(g, -1L), "positive integer")
  expect_error(gmm_reduce("not a gmm", 2L), "gmm")
  expect_error(gmm_reduce(g, 2L, cost = "nope"))
})

test_that("a single-component mixture is returned unchanged", {
  g1 <- gmm(weights = 1, means = list(c(0, 0)), covariances = list(diag(2)))
  expect_equal(gmm_n_components(gmm_reduce(g1, 1L)), 1L)
})
