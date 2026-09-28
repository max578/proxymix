## Tail-probability warning in proxy_functional_ci().

t3_ensemble <- function() {
  tgt <- gmm_target(
    n_dim = 1L,
    log_density = function(x) stats::dt(as.numeric(x), df = 3, log = TRUE),
    normalised = TRUE, name = "student_t3"
  )
  fit <- fit_kld_em(tgt, N = 2L, is_size = 1500L, max_iter = 30L, seed = 1L,
                    validation_size = 0L,
                    proposal = is_mvt(1L, sigma = matrix(9), df = 3))
  gmm_fit_ensemble(fit, B = 20L, seed = 2L)
}

test_that("a tail probability below tail_warn raises a classed warning", {
  ens <- t3_ensemble()
  ## Exact P(X > 15) for t_3 is 3.2e-4; the fitted mixture's value is smaller.
  expect_warning(
    proxy_functional_ci(ens, function(g) pgmm(15, g, lower.tail = FALSE)),
    class = "proxymix_tail_functional"
  )
  ## The complement of a small tail triggers it too.
  expect_warning(
    proxy_functional_ci(ens, function(g) pgmm(15, g)),
    class = "proxymix_tail_functional"
  )
})

test_that("central probabilities and non-probability functionals stay silent", {
  ens <- t3_ensemble()
  expect_no_warning(
    proxy_functional_ci(ens, function(g) pgmm(0.5, g, lower.tail = FALSE))
  )
  expect_no_warning(proxy_functional_ci(ens, gmm_mean))
  expect_no_warning(proxy_functional_ci(ens, function(g) qgmm(0.999, g)))
})

test_that("tail_warn = NULL or 0 silences the warning", {
  ens <- t3_ensemble()
  tail_fn <- function(g) pgmm(15, g, lower.tail = FALSE)
  expect_no_warning(proxy_functional_ci(ens, tail_fn, tail_warn = NULL))
  expect_no_warning(proxy_functional_ci(ens, tail_fn, tail_warn = 0))
  expect_error(proxy_functional_ci(ens, tail_fn, tail_warn = -1), "tail_warn")
})

test_that("gmm_eos_test runs with no measurement noise and matches the AR(1) innovations", {
  y <- withr::with_seed(3, as.numeric(stats::arima.sim(list(ar = 0.6), n = 80L))) + 2
  n <- length(y)
  prior <- gmm(weights = 1, means = list(2),
               covariances = list(matrix(1 / (1 - 0.36))))
  dynamics <- list(A = matrix(0.6), b = 2 * 0.4, Q = matrix(1))
  measurement <- list(C = matrix(1), R = matrix(0))
  res <- gmm_eos_test(prior, dynamics, measurement, y, m = 3L)
  z2 <- (y[-1L] - 2 - 0.6 * (y[-n] - 2))^2
  expect_equal(res$statistic, sum(z2[(n - 3L):(n - 1L)]), tolerance = 1e-10)
})
