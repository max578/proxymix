## --- helpers ---------------------------------------------------------------

## A correlated bivariate dataset with x2 missing-at-random on x1.
.mar_data <- function(n = 400L, beta = 1.5, seed = 1L) {
  withr::with_seed(seed, {
    x1 <- stats::rnorm(n)
    x2 <- beta * x1 + stats::rnorm(n)
    X <- cbind(x1 = x1, x2 = x2)
    X[stats::runif(n) < stats::plogis(0.4 + 0.6 * x1), 2L] <- NA
    X
  })
}

## --- EM from incomplete data -----------------------------------------------

test_that("EM with holes recovers the generating parameters (K = 1)", {
  X <- withr::with_seed(2L, {
    z <- mvnfast::rmvn(800L, mu = c(0.5, -1), sigma = matrix(c(1, 0.7, 0.7, 1.5), 2))
    z[stats::runif(800L) < stats::plogis(z[, 1L]), 2L] <- NA
    z
  })
  fit <- proxymix:::.fit_em_missing(X, K = 1L, max_iter = 200L)$gmm
  expect_equal(fit@means[[1L]], c(0.5, -1), tolerance = 0.15)
  expect_equal(fit@covariances[[1L]][1L, 2L], 0.7, tolerance = 0.2)
  expect_true(fit@covariances[[1L]][2L, 2L] > 1)        # variance not collapsed
})

test_that("EM with holes fits a two-component mixture", {
  X <- withr::with_seed(3L, {
    lab <- sample(1:2, 600L, replace = TRUE)
    mu <- list(c(-2, -2), c(2, 2))
    z <- t(vapply(lab, function(l) mvnfast::rmvn(1L, mu[[l]], 0.6 * diag(2)), numeric(2)))
    z[stats::runif(600L) < 0.35, 2L] <- NA
    z
  })
  fit <- proxymix:::.fit_em_missing(X, K = 2L, max_iter = 200L)$gmm
  centres <- sort(vapply(fit@means, `[`, numeric(1), 1L))
  expect_equal(centres, c(-2, 2), tolerance = 0.4)
})

## --- gmm_impute ------------------------------------------------------------

test_that("gmm_impute returns m complete datasets with no missing values", {
  X <- .mar_data()
  imp <- gmm_impute(X, N = 1L, m = 12L, seed = 1L)
  expect_s3_class(gmm_complete(imp, 1L), NA)            # a matrix, not data frame
  expect_equal(imp@m, 12L)
  comps <- gmm_complete(imp, "all")
  expect_length(comps, 12L)
  expect_false(any(vapply(comps, anyNA, logical(1))))
  ## observed entries are left untouched
  obs <- !is.na(X[, 2L])
  expect_equal(comps[[1L]][obs, 2L], X[obs, 2L])
})

test_that("gmm_impute is reproducible and does not disturb the global RNG", {
  X <- .mar_data()
  set.seed(99L)
  before <- .Random.seed
  a <- gmm_impute(X, N = 1L, m = 6L, seed = 5L)
  b <- gmm_impute(X, N = 1L, m = 6L, seed = 5L)
  expect_identical(a@completions, b@completions)        # reproducible
  expect_identical(before, .Random.seed)                # side-effect-free
})

test_that("a data frame round-trips through gmm_impute / gmm_complete", {
  X <- .mar_data()
  df <- as.data.frame(X)
  imp <- gmm_impute(df, N = 1L, m = 5L, seed = 1L)
  out <- gmm_complete(imp, 1L)
  expect_s3_class(out, "data.frame")
  expect_named(out, c("x1", "x2"))
})

test_that("gmm_impute validates its inputs", {
  X <- .mar_data()
  expect_error(gmm_impute(data.frame(a = letters[1:5]), m = 5L), "numeric")
  expect_error(gmm_impute(X, m = 1L), "m")
  expect_error(gmm_impute(X, mechanism = "mnar"), "arg")
})

## --- the conditioning algebra is the public operator -----------------------

test_that("the internal conditional matches gmm_conditionalise (oracle)", {
  X <- .mar_data()
  imp <- gmm_impute(X, N = 2L, m = 4L, seed = 1L)
  g <- imp@point_fit
  obs <- !is.na(X)
  cm <- proxymix:::.col_cond_moments(g, X, obs, 2L)
  row <- which(!obs[, 2L])[1L]
  gc <- gmm_conditionalise(g, given = c(X[row, 1L], NA))
  analytic_mean <- sum(gc@weights * vapply(gc@means, `[`, numeric(1), 1L))
  expect_equal(cm$mean[row], analytic_mean, tolerance = 1e-10)
})

## --- pooling ---------------------------------------------------------------

test_that("analytic column-mean pooling is sensible and splits the variance", {
  X <- .mar_data(beta = 1.5)
  imp <- gmm_impute(X, N = 1L, m = 20L, seed = 1L)
  pooled <- proxy_pool(imp, "x2")
  expect_equal(attr(pooled, "method"), "analytic")
  expect_true(pooled$conf.low < 0 && pooled$conf.high > 0)   # truth E[x2] = 0
  expect_true(pooled$fmi > 0 && pooled$fmi < 1)
  an <- proxymix:::.analytic_mean(imp, 2L)
  expect_named(an$components, c("complete", "imputation", "parameter"))
  expect_equal(sum(an$components), an$var, tolerance = 1e-12)
})

test_that("analytic and Rubin column-mean pooling agree closely", {
  X <- .mar_data()
  imp <- gmm_impute(X, N = 1L, m = 40L, seed = 2L)
  a <- proxy_pool(imp, "x2", method = "analytic")
  r <- proxy_pool(imp, "x2", method = "rubin")
  expect_lt(abs(a$estimate - r$estimate), 0.03)              # differ by m=40 MC noise
  expect_equal(a$std.error, r$std.error, tolerance = 0.25)   # analytic = m->inf limit
})

test_that("as_mids hands the completions to mice for model pooling", {
  skip_if_not_installed("mice")
  X <- .mar_data(beta = 1.5)
  imp <- gmm_impute(X, N = 1L, m = 20L, seed = 3L)
  md <- as_mids(imp)
  expect_s3_class(md, "mids")
  expect_equal(md$m, 20L)
  slope <- summary(mice::pool(with(md, stats::lm(x2 ~ x1))))
  slope <- slope[slope$term == "x1", "estimate"]
  expect_equal(slope, 1.5, tolerance = 0.35)   # robust, not a single-seed CI
})

test_that("proxy_pool refuses a non-column estimand and points to mice", {
  X <- .mar_data()
  imp <- gmm_impute(X, N = 1L, m = 5L, seed = 1L)
  expect_error(proxy_pool(imp, x2 ~ x1), "column")
})

test_that("proxy_fmi returns a column-mean fraction in [0, 1]", {
  X <- .mar_data()
  imp <- gmm_impute(X, N = 1L, m = 10L, seed = 1L)
  f <- proxy_fmi(imp, "x2")
  expect_length(f, 1L)
  expect_true(f >= 0 && f <= 1)
})

## --- multimodal: the gap single-Gaussian imputers leave --------------------

test_that("a mixture imputer keeps a bimodal column bimodal", {
  X <- withr::with_seed(7L, {
    lab <- sample(1:2, 500L, replace = TRUE)
    x1 <- ifelse(lab == 1L, -2, 2) + stats::rnorm(500L, 0, 0.5)
    x2 <- ifelse(lab == 1L, -2, 2) + stats::rnorm(500L, 0, 0.5)
    X <- cbind(x1 = x1, x2 = x2)
    X[stats::runif(500L) < 0.4, 2L] <- NA
    X
  })
  imp <- gmm_impute(X, N = 2L, m = 5L, seed = 1L)
  filled <- gmm_complete(imp, 1L)[, 2L]
  ## imputed values should populate both modes, not the gap between them
  expect_true(mean(filled < 0) > 0.2 && mean(filled > 0) > 0.2)
  expect_true(mean(abs(filled) < 0.8) < 0.35)            # few land in the gap
})
