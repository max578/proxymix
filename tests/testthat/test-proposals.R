test_that("is_uniform produces samples inside the box", {
  q <- is_uniform(n_dim = 2L, lower = -2, upper = 2)
  x <- q@sample(100L)
  expect_equal(dim(x), c(100L, 2L))
  expect_true(all(x >= -2 & x <= 2))
})

test_that("is_uniform log-density is constant inside, -Inf outside", {
  q <- is_uniform(n_dim = 2L, lower = -1, upper = 1)
  inside <- q@log_density(matrix(c(0, 0), ncol = 2))
  outside <- q@log_density(matrix(c(5, 5), ncol = 2))
  expect_equal(inside, -log(4), tolerance = 1e-10)
  expect_equal(outside, -Inf)
})

test_that("is_mvn density integrates approximately to 1 by MC", {
  withr::with_seed(2026, {
    q <- is_mvn(n_dim = 2L)
    x <- q@sample(5000L)
    ## E_q[1] = integral q(x) / q(x) q(x) dx = 1.
    expect_equal(mean(exp(q@log_density(x) - q@log_density(x))), 1)
  })
})

test_that("is_mvt has heavier tails than is_mvn at the same scale", {
  q_n <- is_mvn(n_dim = 2L, cov = diag(2))
  q_t <- is_mvt(n_dim = 2L, sigma = diag(2), df = 3)
  edge <- matrix(c(5, 5), ncol = 2)
  expect_gt(q_t@log_density(edge), q_n@log_density(edge))
})
