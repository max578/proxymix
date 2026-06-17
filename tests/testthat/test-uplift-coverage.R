## Branch coverage for the decision layer: constructors, prints, validation
## errors, automatic component selection, prediction, response-scale paths and
## the out-of-support fallbacks.

.toy_uplift <- function(seed = 1L, n = 400L, assume = "ignorability") {
  dat <- withr::with_seed(seed, {
    x <- stats::rnorm(n)
    t <- stats::rbinom(n, 1L, 0.5)
    y <- 1 + 0.5 * x + (0.5 + x) * t + stats::rnorm(n, sd = 0.5)
    data.frame(y = y, t = t, x = x)
  })
  fit_uplift(dat, "y", "t", "x", N = 2L, regime = "sample",
             assume = assume, max_iter = 80L, seed = seed)
}

test_that("fit_uplift selects K automatically by BIC", {
  dat <- withr::with_seed(1L, {
    n <- 600L
    x <- stats::rnorm(n)
    t <- stats::rbinom(n, 1L, 0.5)
    y <- 0.5 * t + x + stats::rnorm(n, sd = 0.5)
    data.frame(y = y, t = t, x = x)
  })
  expect_message(
    m <- fit_uplift(dat, "y", "t", "x", N = "auto", n_grid = 1:2,
                    regime = "sample", max_iter = 60L, seed = 1L),
    "Selected K"
  )
  expect_true(gmm_n_components(m@fit) %in% 1:2)
  expect_false(is.null(m@metadata$bic_trace))
})

test_that("uplift_model prints both identification regimes", {
  m_ign <- .toy_uplift(assume = "ignorability")
  m_do <- .toy_uplift(assume = "latent_confounder")
  expect_output(print(m_ign), "uplift_model")
  expect_output(print(m_do), "do-operator")
})

test_that("fit_uplift validates its role arguments", {
  dat <- data.frame(y = stats::rnorm(50), t = stats::rbinom(50, 1L, 0.5),
                    x = stats::rnorm(50))
  expect_error(fit_uplift(as.list(dat), "y", "t", "x"), "data frame")
  expect_error(fit_uplift(dat, "nope", "t", "x"), "outcome")
  expect_error(fit_uplift(dat, "y", "nope", "x"), "treatment")
  expect_error(fit_uplift(dat, "y", "t", "nope"), "covariates")
  expect_error(fit_uplift(dat, "y", "t", "y"), "distinct")
})

test_that("fit_uplift errors on missing values and warns on non-binary T", {
  dat <- data.frame(y = c(NA, stats::rnorm(49)),
                    t = stats::rbinom(50, 1L, 0.5), x = stats::rnorm(50))
  expect_error(fit_uplift(dat, "y", "t", "x", N = 1L, regime = "moment"),
               "missing")

  dat2 <- withr::with_seed(2L, {
    n <- 300L
    x <- stats::rnorm(n)
    tt <- sample(0:2, n, replace = TRUE)            # three levels
    y <- 0.5 * tt + x + stats::rnorm(n, sd = 0.5)
    data.frame(y = y, t = tt, x = x)
  })
  expect_warning(
    fit_uplift(dat2, "y", "t", "x", N = 1L, regime = "moment"),
    "not binary"
  )
})

test_that("proxy_predict returns predictions on link and response scales", {
  m <- .toy_uplift()
  pr <- proxy_predict(m, data.frame(x = c(-1, 0, 1)), t = 1)
  expect_equal(nrow(pr), 3L)
  expect_true(all(is.finite(pr$prediction)))
  expect_error(proxy_predict(m, data.frame(x = 0), t = c(0, 1)), "single")

  ## Binary response scale yields probabilities in [0, 1].
  db <- withr::with_seed(3L, {
    n <- 600L
    x <- stats::rnorm(n)
    t <- stats::rbinom(n, 1L, 0.5)
    y <- stats::rbinom(n, 1L, stats::plogis(-0.3 + x + t))
    data.frame(y = y, t = t, x = x)
  })
  mb <- fit_uplift(db, "y", "t", "x", N = 2L, regime = "sample",
                   outcome_type = "binary", max_iter = 80L, seed = 3L)
  pb <- proxy_predict(mb, data.frame(x = c(-1, 1)), t = 1, scale = "response")
  expect_true(all(pb$prediction >= 0 & pb$prediction <= 1))
})

test_that("response-scale CATE: delta declines, mc delivers", {
  db <- withr::with_seed(4L, {
    n <- 800L
    x <- stats::rnorm(n)
    t <- stats::rbinom(n, 1L, 0.5)
    y <- stats::rbinom(n, 1L, stats::plogis(-0.2 + 0.5 * x + t))
    data.frame(y = y, t = t, x = x)
  })
  mb <- fit_uplift(db, "y", "t", "x", N = 2L, regime = "sample",
                   outcome_type = "binary", max_iter = 80L, seed = 4L)
  expect_message(
    ce <- proxy_cate(mb, data.frame(x = c(0, 1)), scale = "response"),
    "Response-scale"
  )
  expect_true(all(is.na(ce$se)))

  ce_mc <- proxy_cate(mb, data.frame(x = c(0, 1)), scale = "response",
                      se_method = "mc", B = 25L)
  expect_true(all(is.finite(ce_mc$se)))
})

test_that("decide validates value and cost; policy_value resolves all forms", {
  m <- .toy_uplift()
  expect_error(proxy_decide(m, data.frame(x = 0), value = c(1, 2)), "scalar")
  expect_error(proxy_decide(m, data.frame(x = 0), value = 1, cost = c(1, 2)),
               "scalar")

  nd <- data.frame(x = stats::rnorm(100))
  expect_error(proxy_policy_value(m, nd, "bogus", value = 1), "all")
  v_fun <- proxy_policy_value(m, nd, function(ce) as.integer(ce$tau > 0),
                              value = 1, cost = 0.2)$policy_value
  v_vec <- proxy_policy_value(m, nd, rep(c(0L, 1L), 50L),
                              value = 1, cost = 0.2)$policy_value
  expect_true(is.finite(v_fun) && is.finite(v_vec))
  expect_error(proxy_policy_value(m, nd, rep(0L, 3L), value = 1), "0/1 action")
})

test_that("intervene rejects an empty do and accepts an all-NA given", {
  g <- gmm(weights = 1, means = list(c(0, 0, 0)), covariances = list(diag(3)))
  expect_error(gmm_intervene(g, do = c(NA, NA, NA)), "at least one")
  ## All-NA `given` (no conditioning) is valid: Y is left free.
  expect_no_error(gmm_intervene(g, do = c(NA, 1, NA), given = c(NA, NA, NA)))
  ## counterfactual rejects a query that is also intervened.
  expect_error(
    gmm_counterfactual(g, evidence = c(1, 0, 0.2), do = c(1, NA, NA),
                       query = 1L),
    "intervened"
  )
})

test_that("the counterfactual aborts when evidence has zero density", {
  g <- gmm(weights = 1, means = list(c(0, 0, 0)), covariances = list(diag(3)))
  ## A finite but astronomically far observation underflows every component's
  ## density to -Inf -> abduction has no support -> informative abort.
  expect_error(
    gmm_counterfactual(g, evidence = c(1e160, 1, 0.2),
                       do = c(NA, 1, NA), query = 1L)
  )
})

test_that("out-of-support units fall back and are flagged, not errored", {
  m <- .toy_uplift()
  ## A covariate far outside the fitted support underflows the gate -> the
  ## safe responsibilities fall back to the prior and the unit is flagged.
  ce <- proxy_cate(m, data.frame(x = 60), se = FALSE)
  expect_true(ce$overlap_flag[1L])
  expect_true(is.finite(ce$tau[1L]))
})

test_that("the identification report prints the latent-confounder branch", {
  m_do <- .toy_uplift(assume = "latent_confounder")
  rep <- proxy_identification_report(m_do, data.frame(x = stats::rnorm(80)))
  expect_output(print(rep), "only confounder")
})

test_that("a counterfactual law prints with its identified mean", {
  g <- gmm(weights = c(0.5, 0.5), means = list(c(0, 0, 0), c(2, 1, 1)),
           covariances = list(diag(3), diag(3)))
  cf <- gmm_counterfactual(g, evidence = c(1, 0, 0.3),
                           do = c(NA, 1, NA), query = 1L)
  expect_output(print(cf), "counterfactual mean")
  expect_output(print(cf), "not identified")
})
