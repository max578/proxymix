## E3: deterministic-annealing EM and phase-transition component discovery.
##
## Independent oracles: the closed-form first critical temperature
## lambda_max(Sigma^{-1} C) (pure linear algebra, not produced by the fitter);
## the analytic Gaussian/k-means fixed points; and a planted DGP whose global
## optimum is known by construction.

# ---- temperature schedule ------------------------------------------------

test_that(".anneal_schedule is geometric, descending, and lands on t_low", {
  s <- proxymix:::.anneal_schedule(10, 1, 20L)
  expect_length(s, 20L)
  expect_equal(s[1L], 10)
  expect_equal(s[length(s)], 1)
  expect_true(all(diff(s) < 0))
  ## geometric: constant ratio between successive temperatures.
  ratios <- s[-1L] / s[-length(s)]
  expect_lt(stats::sd(ratios), 1e-9)
  expect_equal(proxymix:::.anneal_schedule(5, 2, 1L), 2)
})

test_that(".anneal_schedule rejects malformed temperatures", {
  expect_error(proxymix:::.anneal_schedule(1, 10, 5L), "at least")
  expect_error(proxymix:::.anneal_schedule(-1, 1, 5L), "positive")
  expect_error(proxymix:::.anneal_schedule(10, 1, 0L), "positive integer")
})

# ---- annealed E-step reduces to the cold step at T = 1 -------------------

test_that("annealed sweep at T = 1 equals the cold shared-covariance E/M step", {
  set.seed(1)
  x <- rbind(matrix(stats::rnorm(80, -3), ncol = 2),
             matrix(stats::rnorm(80, 3), ncol = 2))
  n <- nrow(x)
  rw <- rep(1 / n, n)
  weights <- c(0.5, 0.5)
  means <- list(c(-2, -2), c(2, 2))
  Sigma <- diag(2)

  ## Cold reference computed by hand from log(pi_k) + log N_k.
  gauss_log <- proxymix:::gmm_log_unnorm(x, c(1, 1), means, rep(list(Sigma), 2L))
  la <- sweep(gauss_log, 2L, log(weights), `+`)
  resp <- exp(la - proxymix:::logsumexp_rows(la))
  wr <- resp * rw
  Nk <- colSums(wr)
  mu_cold <- lapply(seq_len(2L), function(k) as.numeric(colSums(wr[, k] * x) / Nk[k]))

  step <- proxymix:::.anneal_centroid_sweep(x, rw, weights, means, Sigma, 1)
  expect_equal(step$means, mu_cold, tolerance = 1e-12)
  expect_equal(step$weights, as.numeric(Nk / sum(Nk)), tolerance = 1e-12)
})

test_that("at high temperature the annealed responsibilities are uniform", {
  set.seed(2)
  x <- matrix(stats::rnorm(200), ncol = 2)
  n <- nrow(x)
  rw <- rep(1 / n, n)
  weights <- c(0.5, 0.5)
  means <- list(c(-3, 0), c(3, 0))
  ## At very high T every observation is shared equally, so both means collapse
  ## to the global centroid.
  step <- proxymix:::.anneal_centroid_sweep(x, rw, weights, means, diag(2), 1e8)
  centroid <- colMeans(x)
  expect_equal(step$means[[1L]], as.numeric(centroid), tolerance = 1e-3)
  expect_equal(step$means[[2L]], as.numeric(centroid), tolerance = 1e-3)
})

# ---- closed-form critical temperature (independent linear-algebra oracle) -

test_that(".critical_temperature equals lambda_max(Sigma^{-1} C)", {
  C <- matrix(c(6.25, 1.2, 1.2, 0.49), 2, 2)
  Sigma <- diag(c(1, 1))
  expect_equal(proxymix:::.critical_temperature(C, Sigma),
               max(eigen(solve(Sigma, C), only.values = TRUE)$values),
               tolerance = 1e-12)
  ## Scaling Sigma scales the critical temperature inversely.
  expect_equal(proxymix:::.critical_temperature(C, 2 * Sigma),
               proxymix:::.critical_temperature(C, Sigma) / 2,
               tolerance = 1e-12)
})

test_that("the first phase transition matches the analytic critical temperature", {
  ## A single anisotropic Gaussian blob: the first bifurcation should occur at
  ## T_c = lambda_max(C) (with Sigma = I), checked against the empirically
  ## detected first transition.
  x <- withr::with_seed(42L,
    cbind(stats::rnorm(800, sd = 2.5), stats::rnorm(800, sd = 0.7)))
  path <- gmm_anneal_path(x, k_max = 4L, n_steps = 70L, n_inner = 25L, seed = 1L)
  ## analytic critical temperature == leading eigenvalue of the data covariance
  ## (up to the stabilising ridge on the reference covariance).
  expect_equal(path$t_critical_analytic, path$lambda_max, tolerance = 1e-4)
  ratio <- path$first_critical_temperature / path$t_critical_analytic
  expect_gt(ratio, 0.8)
  expect_lt(ratio, 1.1)
})

# ---- widest-plateau readout and K discovery -------------------------------

test_that(".widest_plateau picks the longest log-temperature plateau", {
  temps <- exp(seq(log(8), log(0.1), length.out = 9L))
  ## counts: 1 over a long high-T stretch, then a brief 2, 3.
  n_eff <- c(1L, 1L, 1L, 1L, 1L, 1L, 2L, 2L, 3L)
  expect_equal(proxymix:::.widest_plateau(temps, n_eff), 1L)
  n_eff2 <- c(1L, 2L, 3L, 3L, 3L, 3L, 3L, 4L, 5L)
  expect_equal(proxymix:::.widest_plateau(temps, n_eff2), 3L)
})

test_that("phase-transition tracking recovers the planted component count", {
  x <- withr::with_seed(7L, rbind(
    matrix(stats::rnorm(200), ncol = 2) +
      matrix(rep(c(-6, -6), each = 100L), ncol = 2),
    matrix(stats::rnorm(200), ncol = 2) +
      matrix(rep(c(6, -6), each = 100L), ncol = 2),
    matrix(stats::rnorm(200), ncol = 2) +
      matrix(rep(c(0, 7), each = 100L), ncol = 2)
  ))
  path <- gmm_anneal_path(x, k_max = 6L, n_steps = 90L, n_inner = 25L, seed = 1L)
  expect_equal(path$k_selected, 3L)
  expect_true(is.data.frame(path$path))
  expect_named(path$path, c("temperature", "n_effective", "free_energy"))
  expect_true(all(path$path$n_effective >= 1L))
})

# ---- robustness: annealing beats cold single-start (SE6) ------------------

test_that("annealing recovers the global optimum where cold EM lands in local optima", {
  skip_on_cran()
  ## Five closely-spaced unequal clusters: a poor central random init regularly
  ## merges neighbours. The planted global optimum is recovered by a heavy
  ## multi-start reference; annealing should match it from a single start.
  plant <- function(seed) withr::with_seed(seed, {
    centres <- list(c(0, 0), c(2.6, 0), c(1.3, 2.3), c(-1.3, 2.3), c(1.3, -2.3))
    nk <- c(160L, 120L, 100L, 80L, 120L)
    do.call(rbind, lapply(seq_along(centres), function(k) {
      matrix(stats::rnorm(nk[k] * 2L, 0, 0.55), ncol = 2) +
        matrix(rep(centres[[k]], each = nk[k]), ncol = 2)
    }))
  })
  x <- plant(404L)
  tgt <- gmm_target_from_samples(x)
  best <- fit_em_samples(tgt, N = 5L, n_starts = 30L, max_iter = 300L)
  ll_star <- best@diagnostics$loglik_final

  n_seeds <- 12L
  cold_fail <- 0L
  ann_fail <- 0L
  for (s in seq_len(n_seeds)) {
    ini <- init_random(N = 5L, p = 2L, centre = colMeans(x),
                       scale = 0.4, sigma_diag = 1, seed = s)
    cold <- fit_em_samples(tgt, N = 5L, init = ini, max_iter = 300L)
    if (ll_star - cold@diagnostics$loglik_final > 1) cold_fail <- cold_fail + 1L
    ann <- fit_em_samples(tgt, N = 5L, anneal = TRUE, max_iter = 300L, seed = s)
    if (ll_star - ann@diagnostics$loglik_final > 1) ann_fail <- ann_fail + 1L
  }
  ## Cold single-start genuinely struggles; annealing essentially never does.
  expect_gt(cold_fail, 3L)
  expect_lt(ann_fail, cold_fail)
  expect_lte(ann_fail, 1L)
})

# ---- determinism and diagnostics -----------------------------------------

test_that("annealed fits are deterministic given a seed and carry the flag", {
  x <- withr::with_seed(1L, rbind(
    matrix(stats::rnorm(120, -4), ncol = 2),
    matrix(stats::rnorm(120, 4), ncol = 2)
  ))
  tgt <- gmm_target_from_samples(x)
  a1 <- fit_em_samples(tgt, N = 2L, anneal = TRUE, seed = 7L, max_iter = 100L)
  a2 <- fit_em_samples(tgt, N = 2L, anneal = TRUE, seed = 7L, max_iter = 100L)
  expect_equal(a1@diagnostics$loglik_final, a2@diagnostics$loglik_final)
  expect_equal(a1@means, a2@means)
  expect_true(isTRUE(a1@diagnostics$annealed))
  expect_true(is.numeric(a1@diagnostics$temp_schedule))
  expect_gt(length(a1@diagnostics$temp_schedule), 1L)
})

test_that("fit_proxymix forwards anneal to the sample regime", {
  x <- withr::with_seed(3L, rbind(
    matrix(stats::rnorm(120, -4), ncol = 2),
    matrix(stats::rnorm(120, 4), ncol = 2)
  ))
  tgt <- gmm_target_from_samples(x)
  fit <- fit_proxymix(tgt, N = 2L, regime = "sample", anneal = TRUE, seed = 3L)
  expect_true(isTRUE(fit@diagnostics$annealed))
  expect_equal(fit@regime, "sample")
})

test_that("KLD-EM accepts the annealing warm-start on a bimodal target", {
  logf <- function(x) {
    if (is.null(dim(x))) x <- matrix(x, nrow = 1L)
    m1 <- c(-3, 0)
    m2 <- c(3, 0)
    l1 <- -0.5 * rowSums((x - matrix(m1, nrow(x), 2L, byrow = TRUE))^2)
    l2 <- -0.5 * rowSums((x - matrix(m2, nrow(x), 2L, byrow = TRUE))^2)
    mx <- pmax(l1, l2)
    log(0.5) + mx + log(exp(l1 - mx) + exp(l2 - mx)) - log(2 * pi)
  }
  tgt <- gmm_target(n_dim = 2L, log_density = logf, normalised = TRUE,
                    name = "two_modes")
  q <- is_mvt(n_dim = 2L, mean = c(0, 0), sigma = 16 * diag(2), df = 5)
  fit <- fit_kld_em(tgt, N = 2L, proposal = q, is_size = 2000L,
                    seed = 1L, max_iter = 40L, anneal = TRUE)
  expect_true(isTRUE(fit@diagnostics$annealed))
  expect_true(is.numeric(fit@diagnostics$temp_schedule))
  ## Both modes recovered (component means near +/- 3 on the first coordinate).
  m1 <- sort(vapply(fit@means, function(m) m[1L], numeric(1L)))
  expect_lt(m1[1L], -2)
  expect_gt(m1[2L], 2)
})

# ---- input validation -----------------------------------------------------

test_that("gmm_anneal_path validates its inputs", {
  expect_error(gmm_anneal_path("not a matrix"), "numeric matrix")
  tgt <- gmm_target(n_dim = 2L, log_density = function(x) rep(0, nrow(x)))
  expect_error(gmm_anneal_path(tgt), "without `samples`")
  x <- matrix(stats::rnorm(40), ncol = 2)
  expect_error(gmm_anneal_path(x, w = rep(1, 5L)), "length")
  expect_error(gmm_anneal_path(x, k_max = 0L), "positive integer")
})

test_that("gmm_anneal_path accepts a gmm_target carrying samples", {
  x <- withr::with_seed(5L, rbind(
    matrix(stats::rnorm(120, -5), ncol = 2),
    matrix(stats::rnorm(120, 5), ncol = 2)
  ))
  tgt <- gmm_target_from_samples(x)
  path <- gmm_anneal_path(tgt, k_max = 4L, n_steps = 60L, n_inner = 20L)
  expect_equal(path$k_selected, 2L)
})

test_that("an explicit temp_schedule overrides the default warm-start path", {
  x <- withr::with_seed(8L, rbind(
    matrix(stats::rnorm(120, -4), ncol = 2),
    matrix(stats::rnorm(120, 4), ncol = 2)
  ))
  tgt <- gmm_target_from_samples(x)
  sched <- c(8, 4, 2, 1, 0.5)
  fit <- fit_em_samples(tgt, N = 2L, anneal = TRUE, temp_schedule = sched,
                        seed = 1L, max_iter = 100L)
  expect_equal(fit@diagnostics$temp_schedule, sched)
  expect_error(
    fit_em_samples(tgt, N = 2L, anneal = TRUE, temp_schedule = c(-1, 1)),
    "positive finite"
  )
})
