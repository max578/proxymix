## Tests for gmm_target_from_posterior() — v0.2.0 Contract A constructor.

test_that("function method accepts a vectorised log-posterior", {
  log_post <- function(theta) -0.5 * rowSums(theta^2)
  tgt <- gmm_target_from_posterior(
    log_post, parameter_names = c("mu", "sigma")
  )
  expect_s7_class(tgt, gmm_target)
  expect_equal(tgt@n_dim, 2L)
  expect_false(isTRUE(tgt@normalised))
  expect_true(is.na(tgt@log_normalizer))
  expect_equal(tgt@metadata$parameter_names, c("mu", "sigma"))
  ## Attribute preserved on the stored callable.
  expect_equal(attr(tgt@log_density, "parameter_names"),
               c("mu", "sigma"))
})

test_that("function method passes log_normalizer when known", {
  log_post <- function(theta) -0.5 * rowSums(theta^2)
  tgt <- gmm_target_from_posterior(
    log_post, parameter_names = c("a", "b"),
    log_normalizer = -log(2 * pi)
  )
  expect_false(isTRUE(tgt@normalised))
  expect_equal(tgt@log_normalizer, -log(2 * pi))
})

test_that("vectorisation contract is enforced by probe", {
  scalar_only <- function(theta) {
    sum(theta^2) ## returns a single number, not a vector
  }
  expect_error(
    gmm_target_from_posterior(scalar_only, parameter_names = c("x", "y")),
    regexp = "length 1"
  )
  raises <- function(theta) stop("boom")
  expect_error(
    gmm_target_from_posterior(raises, parameter_names = c("x", "y")),
    regexp = "boom"
  )
})

test_that("default method points users at correct alternatives", {
  expect_error(
    gmm_target_from_posterior(list(opaque = TRUE)),
    regexp = "No .* method registered"
  )
})

test_that("constructor rejects empty / duplicate / NA parameter names", {
  log_post <- function(theta) rep(0, nrow(theta))
  expect_error(
    gmm_target_from_posterior(log_post, parameter_names = character()),
    regexp = "non-empty character vector"
  )
  expect_error(
    gmm_target_from_posterior(log_post, parameter_names = c("a", "a")),
    regexp = "unique"
  )
  expect_error(
    gmm_target_from_posterior(log_post, parameter_names = c("a", NA_character_)),
    regexp = "non-empty character"
  )
})

test_that("compiled posterior round-trips through fit_proxymix(regime = 'kld')", {
  ## Standard normal in 2D, supplied as an unnormalised log-posterior.
  log_post <- function(theta) -0.5 * rowSums(theta^2)
  tgt <- gmm_target_from_posterior(
    log_post, parameter_names = c("a", "b"),
    log_normalizer = -log(2 * pi)
  )
  fit <- fit_proxymix(
    tgt, N = 1L, regime = "kld",
    proposal = is_mvn(2L, cov = 4 * diag(2)),
    is_size = 1500L, max_iter = 25L, seed = 1L
  )
  expect_equal(fit@regime, "kld")
  expect_true(fit@diagnostics$kld_is_shifted)
  expect_false(is.na(fit@diagnostics$kld_final_absolute))
  ## Absolute KL of the standard normal against itself should be near 0.
  expect_lt(abs(fit@diagnostics$kld_final_absolute), 0.05)
  ## Recovered single-component proxy is near-standard-normal.
  mu_hat <- fit@means[[1L]]
  expect_lt(max(abs(mu_hat)), 0.2)
})

test_that("parameter_names can be carried as an attribute on the callable", {
  log_post <- function(theta) rep(-1, nrow(theta))
  attr(log_post, "parameter_names") <- c("p", "q")
  tgt <- gmm_target_from_posterior(log_post)
  expect_equal(tgt@n_dim, 2L)
  expect_equal(tgt@metadata$parameter_names, c("p", "q"))
})
