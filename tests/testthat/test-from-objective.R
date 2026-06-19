## from_objective() maps the optima of an objective via its Gibbs measure
## (regime iii); gmm_modes() resolves the fitted map into distinct optima.
## Ground truth is known analytically for every scenario below.

# ---------------------------------------------------------------------------
# gmm_modes: pure, exact cases (no fitting)
# ---------------------------------------------------------------------------

test_that("gmm_modes recovers the means of a well-separated mixture", {
  g <- gmm(
    weights = rep(1 / 3, 3),
    means = list(c(-3, 0), c(3, 0), c(0, 4)),
    covariances = rep(list(0.25 * diag(2)), 3)
  )
  res <- gmm_modes(g)
  expect_equal(res$n, 3L)
  ## each true mean is matched by a recovered mode
  truth <- rbind(c(-3, 0), c(3, 0), c(0, 4))
  nearest <- apply(truth, 1L, function(tm) {
    min(sqrt(rowSums((res$modes - matrix(tm, nrow(res$modes), 2L,
                                         byrow = TRUE))^2)))
  })
  expect_true(all(nearest < 1e-2))
})

test_that("gmm_modes on a single component returns its mean", {
  g <- gmm(weights = 1, means = list(c(1.5, -2)), covariances = list(diag(2)))
  res <- gmm_modes(g)
  expect_equal(res$n, 1L)
  expect_equal(as.numeric(res$modes), c(1.5, -2), tolerance = 1e-6)
})

test_that("gmm_modes works in one dimension (no transpose bug)", {
  g <- gmm(
    weights = c(0.5, 0.5),
    means = list(-2, 2),
    covariances = list(matrix(0.2), matrix(0.2))
  )
  res <- gmm_modes(g)
  expect_equal(res$n, 2L)
  expect_equal(ncol(res$modes), 1L)
  expect_equal(sort(as.numeric(res$modes)), c(-2, 2), tolerance = 1e-3)
})

test_that("gmm_modes rejects non-gmm input", {
  expect_error(gmm_modes(list(a = 1)), class = "rlang_error")
})

# ---------------------------------------------------------------------------
# from_objective: ground-truth recovery
# ---------------------------------------------------------------------------

test_that("from_objective recovers both minima of a 1-D bimodal objective", {
  f <- function(v) (v[1]^2 - 4)^2          # minima at +/- 2
  fit <- from_objective(f, lower = -5, upper = 5, N = 6L,
                        is_size = 2500L, n_steps = 5L, seed = 1L)
  expect_s7_class(fit, gmm_fit)
  modes <- gmm_modes(fit)$modes
  ## both true minima recovered within tolerance
  for (m in c(-2, 2)) {
    expect_true(min(abs(modes[, 1] - m)) < 0.25)
  }
})

test_that("from_objective is reproducible under a fixed seed", {
  f <- function(v) (v[1]^2 - 4)^2
  a <- from_objective(f, lower = -5, upper = 5, N = 6L,
                      is_size = 2000L, n_steps = 4L, seed = 7L)
  b <- from_objective(f, lower = -5, upper = 5, N = 6L,
                      is_size = 2000L, n_steps = 4L, seed = 7L)
  expect_equal(do.call(rbind, a@means), do.call(rbind, b@means))
  expect_equal(a@weights, b@weights)
})

test_that("from_objective leaves the global RNG state untouched", {
  set.seed(99L)
  before <- .Random.seed
  f <- function(v) (v[1]^2 - 4)^2
  invisible(from_objective(f, lower = -5, upper = 5, N = 6L,
                           is_size = 1500L, n_steps = 3L, seed = 1L))
  expect_identical(.Random.seed, before)
})

test_that("from_objective maximises when minimise = FALSE", {
  ## peaks at +/- 2 (negated bimodal well)
  f <- function(v) -((v[1]^2 - 4)^2)
  fit <- from_objective(f, lower = -5, upper = 5, N = 6L, minimise = FALSE,
                        is_size = 2500L, n_steps = 5L, seed = 3L)
  modes <- gmm_modes(fit)$modes
  for (m in c(-2, 2)) {
    expect_true(min(abs(modes[, 1] - m)) < 0.3)
  }
})

test_that("from_objective recovers the four Himmelblau minima", {
  skip_on_cran()
  himmelblau <- function(v) {
    x <- v[1]; y <- v[2]
    (x * x + y - 11)^2 + (x + y * y - 7)^2
  }
  truth <- rbind(c(3, 2), c(-2.805118, 3.131312),
                 c(-3.779310, -3.283186), c(3.584428, -1.848126))
  fit <- from_objective(himmelblau, lower = c(-5, -5), upper = c(5, 5),
                        N = 10L, is_size = 4000L, n_steps = 6L, seed = 1L)
  modes <- gmm_modes(fit)$modes
  nearest <- apply(truth, 1L, function(tm) {
    min(sqrt(rowSums((modes - matrix(tm, nrow(modes), 2L, byrow = TRUE))^2)))
  })
  expect_true(all(nearest < 0.4))
})

# ---------------------------------------------------------------------------
# from_objective: input validation
# ---------------------------------------------------------------------------

test_that("from_objective validates its inputs", {
  f <- function(v) sum(v^2)
  expect_error(from_objective("nope", -1, 1), "must be a function")
  expect_error(from_objective(f, c(-1, -1), 1), "same length")
  expect_error(from_objective(f, 1, -1), "must exceed")
  expect_error(from_objective(f, -Inf, 1), "finite")
  expect_error(
    from_objective(f, rep(-1, 11L), rep(1, 11L)), "p <= 10"
  )
})

test_that("from_objective abstains on a non-finite objective", {
  f <- function(v) NaN
  expect_error(
    from_objective(f, lower = -1, upper = 1, N = 4L, is_size = 500L, seed = 1L),
    "no finite values"
  )
})
