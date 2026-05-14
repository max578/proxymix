test_that("dgmm matches a hand-computed density at the origin", {
  g <- gmm(weights = c(0.5, 0.5),
           means = list(c(-1, 0), c(1, 0)),
           covariances = list(diag(2), diag(2)))
  d_origin <- dgmm(c(0, 0), g)
  ## Symmetric: each component contributes equally.
  expected <- 0.5 * mvnfast::dmvn(matrix(c(0, 0), nrow = 1),
                                  mu = c(-1, 0), sigma = diag(2)) +
    0.5 * mvnfast::dmvn(matrix(c(0, 0), nrow = 1),
                        mu = c(1, 0), sigma = diag(2))
  expect_equal(d_origin, as.numeric(expected), tolerance = 1e-10)
})

test_that("dgmm works for matrix inputs", {
  g <- gmm(weights = c(0.4, 0.6),
           means = list(c(-1, 0), c(1, 0)),
           covariances = list(diag(2), diag(2)))
  x <- matrix(c(0, 0, 1, 0, -1, 0, 5, 5), ncol = 2, byrow = TRUE)
  d <- dgmm(x, g)
  expect_length(d, 4L)
  expect_true(all(d > 0))
  expect_true(d[4L] < d[1L])
})

test_that("dgmm log = TRUE returns log-densities", {
  g <- gmm(weights = 1, means = list(c(0, 0)),
           covariances = list(diag(2)))
  expect_equal(dgmm(c(0, 0), g, log = TRUE),
               -log(2 * pi), tolerance = 1e-10)
})

test_that("rgmm samples have correct shape and approximate moments", {
  withr::with_seed(2026, {
    g <- gmm(weights = c(0.4, 0.6),
             means = list(c(-2, 0), c(2, 0)),
             covariances = list(diag(2), diag(2)))
    x <- rgmm(5000L, g)
    expect_equal(dim(x), c(5000L, 2L))
    ## Mixture mean = 0.4 * (-2) + 0.6 * 2 = 0.4 in x1; 0 in x2.
    expect_equal(colMeans(x)[1L], 0.4, tolerance = 0.1)
    expect_equal(colMeans(x)[2L], 0, tolerance = 0.1)
  })
})

test_that("gmm_marginalise reduces dimensionality and preserves weights", {
  g <- gmm(weights = c(0.5, 0.5),
           means = list(c(-1, 0, 2), c(1, 0, -2)),
           covariances = list(diag(3), diag(3)))
  m <- gmm_marginalise(g, keep = c(1L, 3L))
  expect_equal(gmm_dim(m), 2L)
  expect_equal(m@weights, g@weights)
  expect_equal(m@means[[1L]], c(-1, 2))
  expect_equal(m@means[[2L]], c(1, -2))
})

test_that("gmm_conditionalise gives a sane Schur complement", {
  ## Single Gaussian: conditional mean shifts correctly.
  g <- gmm(weights = 1,
           means = list(c(0, 0)),
           covariances = list(matrix(c(1, 0.5, 0.5, 1), 2, 2)))
  c1 <- gmm_conditionalise(g, given = c(NA, 1))
  ## Conditional mean: 0 + 0.5/1 * (1 - 0) = 0.5.
  ## Conditional variance: 1 - 0.5^2 / 1 = 0.75.
  expect_equal(c1@means[[1L]], 0.5, tolerance = 1e-10)
  expect_equal(c1@covariances[[1L]][1, 1], 0.75, tolerance = 1e-10)
})

test_that("gmm_kld is approximately zero when p == q", {
  withr::with_seed(2026, {
    g <- gmm(weights = c(0.4, 0.6),
             means = list(c(-1, 0), c(1, 0)),
             covariances = list(diag(2), diag(2)))
    kld <- gmm_kld(g, g, n_mc = 5000L, variational = TRUE)
    expect_lt(abs(kld$mc), 0.1)
    expect_lt(abs(kld$variational), 0.1)
  })
})

test_that("gmm_kld is positive when q is too narrow", {
  withr::with_seed(2026, {
    p <- gmm(weights = 1,
             means = list(c(0, 0)),
             covariances = list(diag(2) * 2))
    q <- gmm(weights = 1,
             means = list(c(0, 0)),
             covariances = list(diag(2) * 0.5))
    kld <- gmm_kld(p, q, n_mc = 5000L, variational = TRUE)
    expect_gt(kld$mc, 0)
    expect_gt(kld$variational, 0)
  })
})
