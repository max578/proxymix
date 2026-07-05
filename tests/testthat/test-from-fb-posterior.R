## Internal as of v0.13.0; tested through `:::`.
fb_log_posterior_spec <- proxymix:::fb_log_posterior_spec
fb_producer_available <- proxymix:::fb_producer_available
mock_fb_posterior <- proxymix:::mock_fb_posterior
from_fb_posterior <- proxymix:::from_fb_posterior

## Tests for the posterior-producer seam: fb_log_posterior_spec(),
## fb_producer_available(), mock_fb_posterior(), and from_fb_posterior().
## No producer package is assumed to be installed, so every test runs
## against a bare callable or the synthetic mock producer.

# ---- fb_log_posterior_spec() --------------------------------------------------

test_that("spec is built from a vectorised bare callable", {
  log_post <- function(theta) -0.5 * rowSums(theta^2)
  spec <- fb_log_posterior_spec(
    log_post,
    parameter_names = c("mu", "log_sigma"),
    log_normalizer = -log(2 * pi)
  )
  expect_s3_class(spec, "fb_log_posterior_spec")
  expect_equal(spec$n_dim, 2L)
  expect_equal(spec$parameter_names, c("mu", "log_sigma"))
  expect_equal(spec$log_normalizer, -log(2 * pi))
  ## Attribute preserved on the stored callable.
  expect_equal(attr(spec$log_density, "parameter_names"),
               c("mu", "log_sigma"))
})

test_that("spec accepts parameter_names attached to the callable", {
  log_post <- function(theta) rep(-1, nrow(theta))
  attr(log_post, "parameter_names") <- c("p", "q")
  spec <- fb_log_posterior_spec(log_post)
  expect_equal(spec$n_dim, 2L)
  expect_equal(spec$parameter_names, c("p", "q"))
})

test_that("an existing spec is returned unchanged (idempotent)", {
  spec <- mock_fb_posterior(shape = "gaussian", n_dim = 2L)
  expect_identical(fb_log_posterior_spec(spec), spec)
})

test_that("the vectorisation contract is enforced by a probe", {
  scalar_only <- function(theta) sum(theta^2)
  expect_error(
    fb_log_posterior_spec(scalar_only, parameter_names = c("x", "y")),
    regexp = "expected 2"
  )
  raises <- function(theta) stop("boom")
  expect_error(
    fb_log_posterior_spec(raises, parameter_names = c("x", "y")),
    regexp = "boom"
  )
  non_numeric <- function(theta) rep("a", nrow(theta))
  expect_error(
    fb_log_posterior_spec(non_numeric, parameter_names = c("x", "y")),
    regexp = "expected numeric"
  )
})

test_that("invalid parameter_names are rejected", {
  log_post <- function(theta) rep(0, nrow(theta))
  expect_error(
    fb_log_posterior_spec(log_post, parameter_names = character()),
    regexp = "non-empty character"
  )
  expect_error(
    fb_log_posterior_spec(log_post, parameter_names = c("a", "a")),
    regexp = "unique"
  )
  expect_error(
    fb_log_posterior_spec(log_post, parameter_names = c("a", NA_character_)),
    regexp = "non-empty character"
  )
})

test_that("support bounds and draws are validated", {
  log_post <- function(theta) -0.5 * rowSums(theta^2)
  expect_error(
    fb_log_posterior_spec(log_post, parameter_names = c("a", "b"),
                          support_lower = 0),
    regexp = "length-2"
  )
  expect_error(
    fb_log_posterior_spec(log_post, parameter_names = c("a", "b"),
                          support_lower = c(0, 0),
                          support_upper = c(-1, 1)),
    regexp = "must exceed"
  )
  expect_error(
    fb_log_posterior_spec(log_post, parameter_names = c("a", "b"),
                          draws = matrix(0, nrow = 3L, ncol = 3L)),
    regexp = "2 columns"
  )
})

# ---- fb_producer_available() --------------------------------------------------

test_that("the capability probe is a non-erroring logical scalar", {
  res <- fb_producer_available()
  expect_type(res, "logical")
  expect_length(res, 1L)
  expect_false(is.na(res))
})

# ---- mock_fb_posterior() ------------------------------------------------------

test_that("the gaussian mock is a normalised conforming spec", {
  spec <- mock_fb_posterior(shape = "gaussian", n_dim = 3L, sd = 2)
  expect_s3_class(spec, "fb_log_posterior_spec")
  expect_equal(spec$n_dim, 3L)
  expect_equal(spec$log_normalizer, 0)
  ## A normalised density integrates (the log-density at the mode equals the
  ## log normalising constant of an isotropic Gaussian with sd = 2).
  mode_logdens <- spec$log_density(matrix(0, nrow = 1L, ncol = 3L))
  expect_equal(mode_logdens, -3 * (0.5 * log(2 * pi) + log(2)))
})

test_that("the unnormalised gaussian mock records its offset as log Z", {
  spec <- mock_fb_posterior(shape = "gaussian", n_dim = 2L,
                            unnormalised = TRUE)
  expect_false(spec$log_normalizer == 0)
  expect_true(is.finite(spec$log_normalizer))
})

test_that("the banana mock is a 2-D conforming spec", {
  spec <- mock_fb_posterior(shape = "banana")
  expect_equal(spec$n_dim, 2L)
  expect_equal(spec$parameter_names, c("z1", "z2"))
})

# ---- from_fb_posterior(): the consumer verb -----------------------------------

test_that("a fitted model object with no producer raises a clear seam error", {
  ## An opaque 'fitted object' has no registered producer.
  fake_fit <- structure(list(), class = "some_bayes_fit")
  expect_error(
    from_fb_posterior(fake_fit, N = 1L),
    regexp = "posterior producer"
  )
})

test_that("a Gaussian posterior compresses to a near-exact single component", {
  ## Unnormalised isotropic 2-D Gaussian: a single component should recover
  ## it, and the absolute KLD against the truth should be near zero.
  spec <- mock_fb_posterior(shape = "gaussian", n_dim = 2L,
                            unnormalised = TRUE)
  fit <- from_fb_posterior(
    spec, N = 1L,
    proposal = is_mvn(2L, cov = 4 * diag(2)),
    is_size = 2000L, max_iter = 30L, seed = 1L
  )
  expect_s7_class(fit, gmm_fit)
  expect_equal(fit@regime, "kld")
  ## The producer supplied a finite log Z, so the absolute KLD is recoverable.
  expect_true(fit@diagnostics$kld_is_shifted)
  expect_false(is.na(fit@diagnostics$kld_final_absolute))
  expect_lt(abs(fit@diagnostics$kld_final_absolute), 0.05)
  ## Recovered single component is near standard normal.
  expect_lt(max(abs(fit@means[[1L]])), 0.2)
  expect_equal(diag(fit@covariances[[1L]]), c(1, 1), tolerance = 0.15)
  ## Provenance stamped on the fit metadata.
  expect_equal(fit@metadata$from_fb_posterior$producer_name,
               spec$name)
})

test_that("a banana posterior compresses with good moment agreement", {
  ## Quality bar: the proxy's first two moments should match the target's,
  ## estimated from a large independent sample drawn from the true banana.
  spec <- mock_fb_posterior(shape = "banana")
  fit <- from_fb_posterior(
    spec, N = 3L,
    proposal = is_mvt(2L, sigma = 6 * diag(2), df = 5),
    is_size = 6000L, max_iter = 60L, seed = 2L,
    validation_size = 6000L
  )
  expect_s7_class(fit, gmm_fit)

  ## Ground-truth moments of the banana (z1 ~ N(0,1); z2 = w + 0.5(z1^2 - 1),
  ## w ~ N(0,1) independent). E[z1] = 0, Var[z1] = 1, E[z2] = 0,
  ## Var[z2] = Var[w] + 0.25 Var[z1^2] = 1 + 0.25 * 2 = 1.5.
  truth_mean <- c(0, 0)
  truth_var <- c(1, 1.5)

  ## Proxy moments via a large draw from the fitted mixture.
  set.seed(99L)
  proxy_draws <- rgmm(40000L, fit)
  proxy_mean <- colMeans(proxy_draws)
  proxy_var <- apply(proxy_draws, 2L, stats::var)

  expect_lt(max(abs(proxy_mean - truth_mean)), 0.15)
  expect_lt(max(abs(proxy_var - truth_var)), 0.40)

  ## Held-out validation KLD is finite and small -- the fit generalises
  ## beyond the single fitting IS draw.
  expect_true(is.finite(fit@diagnostics$validation_kld))
  expect_lt(fit@diagnostics$validation_kld, 0.5)
})

test_that("the bare-callable path matches the spec path", {
  log_post <- function(theta) -0.5 * rowSums(theta^2) - log(2 * pi)
  fit_callable <- from_fb_posterior(
    log_post, N = 1L, parameter_names = c("a", "b"),
    log_normalizer = 0,
    proposal = is_mvn(2L, cov = 4 * diag(2)),
    is_size = 1500L, max_iter = 25L, seed = 7L
  )
  spec <- fb_log_posterior_spec(log_post, parameter_names = c("a", "b"),
                                log_normalizer = 0)
  fit_spec <- from_fb_posterior(
    spec, N = 1L,
    proposal = is_mvn(2L, cov = 4 * diag(2)),
    is_size = 1500L, max_iter = 25L, seed = 7L
  )
  expect_equal(fit_callable@means[[1L]], fit_spec@means[[1L]])
  expect_equal(fit_callable@covariances[[1L]], fit_spec@covariances[[1L]])
})

test_that("the default proposal is seeded from supplied draws", {
  ## A posterior offset far from the origin: the default proposal must follow
  ## the supplied draws, not sit on the origin, or the IS ESS collapses.
  log_post <- function(theta) {
    -0.5 * rowSums((theta - matrix(c(10, 10), nrow = nrow(theta),
                                   ncol = 2L, byrow = TRUE))^2) - log(2 * pi)
  }
  draws <- mvnfast::rmvn(500L, mu = c(10, 10), sigma = diag(2))
  spec <- fb_log_posterior_spec(log_post, parameter_names = c("a", "b"),
                                log_normalizer = 0, draws = draws)
  fit <- from_fb_posterior(spec, N = 1L, is_size = 2000L,
                           max_iter = 30L, seed = 3L)
  ## Recovered mean near (10, 10) -- only possible if the proposal tracked
  ## the draws.
  expect_lt(max(abs(fit@means[[1L]] - c(10, 10))), 0.3)
  ## ESS not collapsed.
  expect_gt(fit@diagnostics$ess_relative, 0.05)
})

test_that("the dimensional guard rejects p > 10 and warns for p > 5", {
  spec_big <- mock_fb_posterior(shape = "gaussian", n_dim = 11L)
  expect_error(
    from_fb_posterior(spec_big, N = 1L),
    regexp = "n_dim <= 10"
  )
  spec_mid <- mock_fb_posterior(shape = "gaussian", n_dim = 6L)
  ## A wide, heavy-tailed proposal keeps the IS ESS healthy, so the only
  ## warning is the dimensional "well-tested" caution.
  wide <- is_mvt(6L, sigma = 9 * diag(6L), df = 5)
  expect_warning(
    from_fb_posterior(spec_mid, N = 1L, proposal = wide,
                      is_size = 1500L, max_iter = 5L,
                      seed = 1L, support_warn = FALSE, min_ess = 1L),
    regexp = "well-tested"
  )
})
