test_that("banana_target log-density is normalised by Monte Carlo", {
  withr::with_seed(2026, {
    b <- banana_target()
    q <- is_mvt(n_dim = 2L, mean = c(0, 0),
                sigma = 4 * diag(2), df = 5)
    x <- q@sample(20000L)
    log_w <- b@log_density(x) - q@log_density(x)
    z_est <- mean(exp(log_w))
    ## Should be approximately 1 since the banana density is normalised.
    expect_equal(z_est, 1, tolerance = 0.1)
  })
})

test_that("donut_target log-density is normalised by Monte Carlo", {
  withr::with_seed(2026, {
    d <- donut_target()
    q <- is_uniform(n_dim = 2L, lower = -5, upper = 5)
    x <- q@sample(20000L)
    log_w <- d@log_density(x) - q@log_density(x)
    z_est <- mean(exp(log_w))
    expect_equal(z_est, 1, tolerance = 0.15)
  })
})

test_that("mixture_target log-density agrees with dgmm on the same weights", {
  m <- mixture_target()
  g <- gmm(weights = m@metadata$weights,
           means = m@metadata$means,
           covariances = m@metadata$covariances)
  x <- matrix(c(0, 0, -2, -2, 2, 2, 5, 5), ncol = 2, byrow = TRUE)
  expect_equal(m@log_density(x), dgmm(x, g, log = TRUE),
               tolerance = 1e-10)
})

test_that("mixture_target with_samples generates correct moments", {
  withr::with_seed(2026, {
    m <- mixture_target(with_samples = TRUE, n = 5000L, seed = 1L)
    ## Mixture mean = 0.3 * (-2,-2) + 0.4 * (0,0) + 0.3 * (2,2) = (0, 0).
    expect_equal(colMeans(m@samples), c(0, 0), tolerance = 0.1)
  })
})
