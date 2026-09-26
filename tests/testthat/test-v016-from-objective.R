## from_objective(): the box bounds the map, and one fit raises at most one
## low effective-sample-size warning.

test_that("from_objective keeps the map inside the box when the optimum lies outside it", {
  f <- function(v) (v[1] + 2)^2            # unconstrained minimum at -2
  fit <- from_objective(f, lower = 0, upper = 4, N = 6L,
                        is_size = 2000L, n_steps = 5L, seed = 1L)
  set.seed(11L)
  x <- rgmm(20000L, fit)
  expect_gt(mean(x >= 0 & x <= 4), 0.9)
  modes <- gmm_modes(fit)$modes
  expect_lt(abs(modes[1L, 1L]), 0.3)
})

test_that("from_objective does not evaluate the objective outside the box", {
  outside <- 0L
  f <- function(v) {
    if (v[1] < 0 || v[1] > 4) outside <<- outside + 1L
    (v[1] + 2)^2
  }
  invisible(from_objective(f, lower = 0, upper = 4, N = 4L,
                           is_size = 1000L, n_steps = 3L, seed = 1L))
  expect_identical(outside, 0L)
})

test_that("from_objective raises one low-ESS warning per call, not one per cooling step", {
  f <- function(v) sum((v - 0.3)^2) * 1e3
  warns <- list()
  fit <- withCallingHandlers(
    from_objective(f, rep(-1, 5L), rep(1, 5L), N = 4L, is_size = 400L,
                   n_steps = 6L, max_iter = 20L, seed = 3L),
    proxymix_low_ess = function(cnd) {
      warns[[length(warns) + 1L]] <<- cnd
      invokeRestart("muffleWarning")
    }
  )
  expect_length(warns, 1L)
  expect_match(conditionMessage(warns[[1L]]), "cooling step")
  ess <- fit@metadata$from_objective$ess
  expect_length(ess, 6L)
  expect_equal(ess[6L], fit@diagnostics$ess)
})

test_that("from_objective stays silent when only an earlier cooling step had low ESS", {
  f <- function(v) v[1]^2 - 10 * cos(2 * pi * v[1]) + 10
  fit <- expect_no_warning(
    from_objective(f, lower = -5, upper = 5, N = 4L, is_size = 500L,
                   n_steps = 6L, max_iter = 20L, seed = 2L),
    class = "proxymix_low_ess"
  )
  ess <- fit@metadata$from_objective$ess
  expect_true(any(ess[-6L] < 50))
  expect_gte(ess[6L], 50)
})

## a logistic-regression slope with 150 observations; its likelihood is the
## reference, by grid quadrature on the box
logit_nll <- local({
  set.seed(20260926L)
  x <- stats::rnorm(150L)
  y <- stats::rbinom(150L, 1L, stats::plogis(1.2 * x))
  function(b) sum(log1p(exp(b[1] * x))) - sum(y * b[1] * x)
})

hellinger_to_likelihood <- function(fit, lower, upper) {
  grid <- seq(lower, upper, length.out = 2001L)
  log_lik <- -vapply(grid, logit_nll, numeric(1))
  r <- exp(log_lik - max(log_lik))
  r <- r / sum(r)
  q <- dgmm(matrix(grid, ncol = 1L), fit)
  q <- q / sum(q)
  sqrt(max(0, 1 - sum(sqrt(r * q))))
}

test_that("from_objective with scale = 'loglik' returns the likelihood surface", {
  fit <- from_objective(logit_nll, lower = -5, upper = 5, N = 4L,
                        is_size = 2000L, seed = 1L, scale = "loglik")
  temps <- fit@metadata$from_objective$temperatures
  expect_equal(temps[length(temps)], 1)
  expect_identical(fit@metadata$from_objective$scale, "loglik")
  expect_lt(hellinger_to_likelihood(fit, -5, 5), 0.05)
  ## the default ladder ends well above 1 here and gives a wider map
  wide <- from_objective(logit_nll, lower = -5, upper = 5, N = 4L,
                         is_size = 2000L, seed = 1L)
  expect_gt(min(wide@metadata$from_objective$temperatures), 1)
  expect_gt(hellinger_to_likelihood(wide, -5, 5), 0.2)
})

test_that("from_objective keeps the default ladder and a supplied temperature unchanged", {
  f <- function(v) (v[1]^2 - 4)^2
  fit <- from_objective(f, lower = -5, upper = 5, N = 4L,
                        is_size = 1000L, n_steps = 4L, seed = 1L)
  temps <- fit@metadata$from_objective$temperatures
  expect_length(temps, 4L)
  expect_equal(temps[1L] / temps[4L], 40)
  a <- from_objective(f, lower = -5, upper = 5, N = 4L, temperature = c(20, 0.5),
                      is_size = 1000L, n_steps = 4L, seed = 1L)
  b <- from_objective(f, lower = -5, upper = 5, N = 4L, temperature = c(20, 0.5),
                      is_size = 1000L, n_steps = 4L, seed = 1L, scale = "loglik")
  expect_equal(a@metadata$from_objective$temperatures,
               exp(seq(log(20), log(0.5), length.out = 4L)))
  expect_equal(a@means, b@means)
  expect_equal(a@weights, b@weights)
})

test_that("from_objective with scale = 'loglik' lengthens a steep default ladder", {
  f <- function(v) 1e4 * v[1]^2
  fit <- suppressWarnings(
    from_objective(f, lower = -1, upper = 1, N = 2L, is_size = 300L,
                   max_iter = 5L, seed = 1L, scale = "loglik")
  )
  temps <- fit@metadata$from_objective$temperatures
  expect_gt(length(temps), 6L)
  expect_lte(max(temps[-length(temps)] / temps[-1L]), 40^(1 / 5) + 1e-8)
  expect_error(from_objective(f, -1, 1, scale = "nope"), class = "rlang_error")
})
