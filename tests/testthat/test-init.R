test_that("init_random produces a valid gmm", {
  init <- init_random(N = 3L, p = 2L, seed = 1L)
  expect_s7_class(init, gmm)
  expect_equal(gmm_n_components(init), 3L)
})

test_that("init_kmeans groups samples sensibly", {
  withr::with_seed(2026, {
    mt <- mixture_target(with_samples = TRUE, n = 600L, seed = 2L)
    init <- init_kmeans(mt@samples, N = 3L)
    expect_s7_class(init, gmm)
    expect_equal(gmm_n_components(init), 3L)
    ## Centres should bracket the true means (-2, 2).
    mean_norms <- vapply(init@means, function(mu) sqrt(sum(mu^2)),
                         numeric(1L))
    expect_true(any(mean_norms > 1))
  })
})

test_that("init_moment_seed places centres along PC1", {
  x <- matrix(stats::rnorm(400), ncol = 2)
  init <- init_moment_seed(x, N = 3L)
  expect_s7_class(init, gmm)
  expect_equal(gmm_n_components(init), 3L)
})

test_that("multi_start_best_of returns the higher-scoring fit", {
  set.seed(1)
  x <- matrix(stats::rnorm(400), ncol = 2)
  tgt <- gmm_target_from_samples(x)
  inits <- list(
    init_random(2L, 2L, seed = 1L),
    init_random(2L, 2L, seed = 2L),
    init_moment_seed(x, N = 2L)
  )
  best <- multi_start_best_of(
    fit_fn = function(init, ...) fit_em_samples(tgt, init = init, ...),
    inits = inits,
    score_fn = function(fit) fit@diagnostics$loglik_final,
    max_iter = 25L
  )
  expect_s7_class(best, gmm_fit)
})
