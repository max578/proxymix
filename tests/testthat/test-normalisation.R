test_that("built-in targets are declared normalised = TRUE", {
  expect_true(banana_target()@normalised)
  expect_true(donut_target()@normalised)
  expect_true(mixture_target()@normalised)
})

test_that("gmm_target_from_samples is NA-normalised (unknown)", {
  x <- matrix(stats::rnorm(20), ncol = 2)
  expect_true(is.na(gmm_target_from_samples(x)@normalised))
})

test_that("explicit normalised flag round-trips", {
  f <- function(x) -0.5 * rowSums(x^2)
  un <- gmm_target(n_dim = 2L, log_density = f, normalised = FALSE)
  expect_false(un@normalised)
  expect_true(is.na(un@log_normalizer))

  norm <- gmm_target(n_dim = 2L, log_density = f, normalised = TRUE,
                     log_normalizer = 0)
  expect_true(norm@normalised)
  expect_equal(norm@log_normalizer, 0)
})

test_that("validator rejects non-scalar normalised / log_normalizer", {
  expect_error(
    gmm_target(n_dim = 2L,
               log_density = function(x) 0,
               normalised = c(TRUE, FALSE)),
    "length-1 logical"
  )
  expect_error(
    gmm_target(n_dim = 2L,
               log_density = function(x) 0,
               log_normalizer = c(0, 1)),
    "length-1 numeric"
  )
})

test_that("KLD on a normalised target reports kld_is_shifted = FALSE", {
  fit <- fit_kld_em(banana_target(), N = 2L,
                    is_size = 1000L, max_iter = 10L, seed = 1L)
  expect_false(fit@diagnostics$kld_is_shifted)
  expect_true(is.na(fit@diagnostics$kld_shift_explanation))
  expect_equal(fit@diagnostics$kld_final_absolute, fit@diagnostics$kld_final)
})

test_that("KLD on an unnormalised target reports kld_is_shifted = TRUE", {
  ## Standard 2-D Gaussian log-density but missing the 2*pi normaliser.
  unnormalised_log_f <- function(x) -0.5 * rowSums(x^2)
  ut <- gmm_target(n_dim = 2L,
                   log_density = unnormalised_log_f,
                   normalised = FALSE)
  fit <- fit_kld_em(ut, N = 1L, is_size = 1000L,
                    max_iter = 10L, seed = 1L)
  expect_true(fit@diagnostics$kld_is_shifted)
  expect_true(is.na(fit@diagnostics$kld_final_absolute))
})

test_that("supplying log_normalizer corrects the absolute KLD", {
  ## The same log-density as above, but now declared unnormalised with
  ## the correct log_normalizer = log(2 * pi) for a 2-D N(0, I).
  unnormalised_log_f <- function(x) -0.5 * rowSums(x^2)
  ut <- gmm_target(n_dim = 2L,
                   log_density = unnormalised_log_f,
                   normalised = FALSE,
                   log_normalizer = log(2 * pi))
  fit <- fit_kld_em(ut, N = 1L, is_size = 1500L,
                    max_iter = 15L, seed = 1L)
  expect_true(is.finite(fit@diagnostics$kld_final_absolute))
  ## The shift uses target@log_normalizer additively.
  expect_equal(
    fit@diagnostics$kld_final_absolute,
    fit@diagnostics$kld_final + log(2 * pi),
    tolerance = 1e-12
  )
})

test_that("hellinger_mc warns when target is not normalised", {
  unnormalised_log_f <- function(x) -0.5 * rowSums(x^2)
  ut <- gmm_target(n_dim = 2L,
                   log_density = unnormalised_log_f,
                   normalised = FALSE)
  fit <- fit_kld_em(ut, N = 2L, is_size = 800L, max_iter = 10L, seed = 1L)
  expect_warning(hellinger_mc(fit, n_mc = 500L, seed = 1L),
                 "not declared")
})
