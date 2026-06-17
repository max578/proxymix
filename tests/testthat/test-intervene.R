## The do-operator and counterfactual, graded against hand-built linear SCMs.
##
## Each component of a joint Gaussian over (Y, T, X) is an affine structural
## equation, so the interventional and counterfactual means have closed forms
## we can write down independently of any proxymix output. The fixtures below
## fix the covariance so the within-component regression coefficients are known.

## A single-Gaussian SCM with Y = a T + b X + noise, (T, X) independent --------
.linear_scm_k1 <- function(a = 2, b = -0.7, s2 = 0.09) {
  var_y <- a^2 + b^2 + s2
  S <- matrix(c(var_y, a,   b,
                a,     1,   0,
                b,     0,   1), nrow = 3L, byrow = TRUE)
  gmm(weights = 1,
      means = list(c(0, 0, 0)),
      covariances = list(S),
      name = "linear_scm_k1")
}

## Marginal mean of the Y coordinate (position 1) of a mixture.
.y_mean <- function(g, pos = 1L) {
  sum(g@weights * vapply(g@means, function(m) m[pos], numeric(1L)))
}

test_that("do-CATE at K = 1 equals the within-component treatment slope", {
  a <- 2
  g <- .linear_scm_k1(a = a)

  do1 <- gmm_intervene(g, do = c(NA, 1, NA))
  do0 <- gmm_intervene(g, do = c(NA, 0, NA))

  ## free coordinates are (Y, X); Y is position 1.
  tau_do <- .y_mean(do1) - .y_mean(do0)
  expect_equal(tau_do, a, tolerance = 1e-8)
})

test_that("intervening while conditioning on X gives the X-specific do-response", {
  a <- 2
  b <- -0.7
  g <- .linear_scm_k1(a = a, b = b)
  x0 <- 0.5

  do1 <- gmm_intervene(g, do = c(NA, 1, NA), given = c(NA, NA, x0))
  ## With X conditioned, the only free coordinate is Y.
  expect_equal(gmm_dim(do1), 1L)
  expect_equal(.y_mean(do1), a * 1 + b * x0, tolerance = 1e-8)
})

test_that("the do-gate ignores the intervened coordinate but honours `given`", {
  ## Two regimes with different treatment centres: observing T is informative
  ## about the regime (seeing), but do() must not let it move the gate.
  g <- gmm(
    weights = c(0.5, 0.5),
    means = list(c(0, -1, 0), c(3, 1, 0)),
    covariances = list(diag(3), diag(3))
  )

  ## do(T = 1), no conditioning: gate stays at the prior.
  di <- gmm_intervene(g, do = c(NA, 1, NA))
  expect_equal(di@weights, g@weights, tolerance = 1e-12)

  ## do(T = 1) while conditioning on X = 0.4: gate is the X-only posterior,
  ## which equals conditioning on X alone via gmm_conditionalise().
  di_x <- gmm_intervene(g, do = c(NA, 1, NA), given = c(NA, NA, 0.4))
  cond_x <- gmm_conditionalise(g, given = c(NA, NA, 0.4))
  expect_equal(di_x@weights, cond_x@weights, tolerance = 1e-10)
})

test_that("seeing and doing diverge when treatment is confounded with regime", {
  ## Regime drives both T and the Y baseline -> a back-door path. Seeing T = 1
  ## up-weights the high-Y regime; do(T = 1) does not.
  g <- gmm(
    weights = c(0.5, 0.5),
    means = list(c(0, -2, 0), c(4, 2, 0)),
    covariances = list(diag(3), diag(3))
  )

  see1 <- gmm_conditionalise(g, given = c(NA, 1, NA))   # p(Y, X | T = 1)
  do1 <- gmm_intervene(g, do = c(NA, 1, NA))            # p(Y, X | do(T = 1))

  ## do() keeps the prior 50/50 gate; seeing tilts heavily to regime 2.
  expect_equal(do1@weights, c(0.5, 0.5), tolerance = 1e-12)
  expect_gt(see1@weights[2L], 0.95)
  ## The two interventional/observational Y-means must differ materially.
  expect_gt(abs(.y_mean(see1) - .y_mean(do1)), 1)
})

test_that("counterfactual mean at K = 1 is y + slope * (t' - t)", {
  a <- 2
  g <- .linear_scm_k1(a = a)

  y_obs <- 1.2
  cf <- gmm_counterfactual(g, evidence = c(y_obs, 0, 0.5),
                           do = c(NA, 1, NA), query = 1L)

  ## One component -> one atom -> mean equals the atom.
  expect_length(cf@atoms, 1L)
  expect_equal(cf@mean, y_obs + a * (1 - 0), tolerance = 1e-8)
  expect_equal(cf@atoms[1L], cf@mean, tolerance = 1e-12)
})

test_that("counterfactual atoms preserve the observed residual within a regime", {
  ## Two parallel regimes with the SAME treatment slope a: whichever regime
  ## generated the unit, the counterfactual shift is a * (t' - t), so the mean
  ## is y + a (t' - t) regardless of the abduction weights.
  a <- 1.5
  S <- matrix(c(a^2 + 1, a, 0,
                a,       1, 0,
                0,       0, 1), nrow = 3L, byrow = TRUE)
  g <- gmm(weights = c(0.4, 0.6),
           means = list(c(0, 0, -1), c(2, 0, 1)),
           covariances = list(S, S))

  y_obs <- 0.8
  cf <- gmm_counterfactual(g, evidence = c(y_obs, 0, 0.3),
                           do = c(NA, 1, NA), query = 1L)
  expect_equal(cf@mean, y_obs + a * 1, tolerance = 1e-8)
})

test_that("counterfactual variance and tail accessors refuse (not identified)", {
  g <- .linear_scm_k1()
  cf <- gmm_counterfactual(g, evidence = c(1, 0, 0.2),
                           do = c(NA, 1, NA), query = 1L)

  expect_error(gmm_cf_variance(cf), class = "proxymix_not_identified")
  expect_error(gmm_cf_tail_prob(cf, threshold = 2),
               class = "proxymix_not_identified")
  expect_equal(gmm_cf_mean(cf), cf@mean)
})

test_that("intervene / counterfactual validate their coordinate arguments", {
  g <- .linear_scm_k1()

  ## A coordinate in both do and given.
  expect_error(
    gmm_intervene(g, do = c(NA, 1, NA), given = c(NA, 0, NA)),
    "both"
  )
  ## No free coordinate left.
  expect_error(
    gmm_intervene(g, do = c(NA, 1, NA), given = c(NA, NA, 0)),
    NA  # this one is valid (Y free); ensure it does NOT error
  )
  ## query not observed in evidence.
  expect_error(
    gmm_counterfactual(g, evidence = c(NA, 0, 0.2),
                       do = c(NA, 1, NA), query = 1L),
    "observed"
  )
  ## intervened coordinate not observed.
  expect_error(
    gmm_counterfactual(g, evidence = c(1, NA, 0.2),
                       do = c(NA, 1, NA), query = 1L),
    "observed"
  )
})
