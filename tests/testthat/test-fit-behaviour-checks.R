# Each check compares with a reference the fit did not produce: a closed form,
# a quadrature, the moments supplied, or a constructed data set.

collect_classes <- function(expr) {
  classes <- character(0)
  value <- withCallingHandlers(expr, warning = function(w) {
    classes <<- c(classes, class(w)[1L])
    invokeRestart("muffleWarning")
  }, message = function(m) invokeRestart("muffleMessage"))
  list(value = value, classes = classes)
}

mixture_moments <- function(g) {
  m <- Reduce(`+`, Map(`*`, g@weights, g@means))
  v <- Reduce(`+`, Map(function(w, mu, S) w * (S + tcrossprod(mu - m)),
                       g@weights, g@means, g@covariances))
  list(mean = m, cov = v)
}

# A flat-prior logistic posterior with a constant known by quadrature.
logistic_case <- function() {
  x <- c(-17, -13, -12, -11, -10, -8, -7, -5, -4, -3, -2, -1, 0, 1, 2, 3, 5, 6, 8, 9, 10, 11, 12)
  y <- c(1, 1, 1, 1, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0)
  lp <- function(th) {
    if (is.null(dim(th))) th <- matrix(th, ncol = 2L)
    eta <- outer(rep(1, length(x)), th[, 1L]) + outer(x, th[, 2L])
    colSums(y * eta - (pmax(eta, 0) + log1p(exp(-abs(eta)))))
  }
  op <- stats::optim(c(0, 0), function(b) -lp(b), hessian = TRUE)
  list(lp = lp, mode = op$par, sigma = solve(op$hessian))
}

test_that("moment regime matches the supplied moments at every N", {
  S <- matrix(c(4, 1.5, 1.5, 2), 2L)
  tgt <- maxent_target(moments = list(mean = c(1, 2), cov = S))
  for (N in 2:4) {
    r <- collect_classes(fit_proxymix(tgt, N = N, regime = "moment"))
    mo <- mixture_moments(r$value)
    expect_equal(as.numeric(mo$mean), c(1, 2), tolerance = 1e-10)
    expect_equal(mo$cov, S, tolerance = 1e-10)
  }
  expect_message(fit_proxymix(tgt, N = 2L, regime = "moment"), "adds nothing")
})

test_that("set.seed reaches the random starts of the sample regime", {
  x <- withr::with_seed(8L, matrix(c(stats::rnorm(70, 0), stats::rnorm(50, 3), stats::rnorm(30, 6))))
  tgt <- gmm_target_from_samples(x)
  fit_at <- function(s) {
    set.seed(s)
    collect_classes(fit_em_samples(tgt, N = 5L))$value@diagnostics$loglik_final
  }
  expect_identical(fit_at(1L), fit_at(1L))
  expect_gt(length(unique(round(vapply(1:12, fit_at, numeric(1L)), 3))), 1L)
})

test_that("select_N does not choose a collapsed fit when another is available", {
  set.seed(11)
  x <- matrix(c(stats::rnorm(120), rep(4, 14), rep(-4, 14)))
  sel <- collect_classes(select_N(gmm_target_from_samples(x), candidates = 1:4, seed = 1L))$value
  tab <- sel$table
  expect_false(tab$collapsed[tab$chosen])
  expect_true(any(tab$collapsed))
  expect_lt(min(tab$bic), tab$bic[tab$chosen])
})

test_that("select_N defaults reach the optimum that many starts reach", {
  set.seed(3)
  x <- matrix(c(stats::rnorm(60, 40, 3), stats::rnorm(50, 58, 3), stats::rnorm(40, 75, 2)))
  tgt <- gmm_target_from_samples(x)
  many <- collect_classes(select_N(tgt, candidates = 3L, seed = 1L, n_starts = 60L))$value
  def <- collect_classes(select_N(tgt, candidates = 3L, seed = 1L))$value
  expect_equal(def$table$bic, many$table$bic, tolerance = 1e-3)
})

test_that("heldout_kld_se describes the error of heldout_kld", {
  skip_on_cran()
  case <- logistic_case()
  tgt <- gmm_target(n_dim = 2L, log_density = case$lp, normalised = FALSE)
  q <- is_mvt(2L, mean = case$mode, sigma = 4 * case$sigma, df = 5)
  res <- do.call(rbind, lapply(1:12, function(s) {
    d <- collect_classes(fit_kld_em(tgt, N = 2L, proposal = q, is_size = 3000L, seed = s,
                                    validation_size = 5000L))$value@diagnostics
    c(held = d$heldout_kld, se = d$heldout_kld_se, mc = d$validation_mc_se)
  }))
  ## All twelve estimate nearly the same KL, so their spread is the error.
  ratio <- stats::sd(res[, "held"]) / mean(res[, "se"])
  expect_gt(ratio, 0.5)
  expect_lt(ratio, 2)
  expect_gt(mean(res[, "se"]), 2 * mean(res[, "mc"]))
})

test_that("adapt = pmc stops on its own and raises no spurious objective warning", {
  case <- logistic_case()
  tgt <- gmm_target(n_dim = 2L, log_density = case$lp, normalised = FALSE)
  q <- is_mvt(2L, mean = case$mode, sigma = 4 * case$sigma, df = 5)
  r <- collect_classes(fit_kld_em(tgt, N = 2L, proposal = q, is_size = 2000L, seed = 1L,
                                  adapt = "pmc", max_iter = 1000L, validation_size = 0L))
  expect_true(r$value@converged)
  expect_lt(r$value@iterations, 1000L)
  expect_false("proxymix_nonmonotone" %in% r$classes)
})

test_that("the default iteration limit lets an easy target converge", {
  case <- logistic_case()
  tgt <- gmm_target(n_dim = 2L, log_density = case$lp, normalised = FALSE)
  q <- is_mvt(2L, mean = case$mode, sigma = 4 * case$sigma, df = 5)
  f <- collect_classes(fit_kld_em(tgt, N = 3L, proposal = q, is_size = 3000L, seed = 1L,
                                  validation_size = 0L))$value
  expect_true(f@converged)
})

test_that("the heavy-tail warning separates mixed scales from stable laws", {
  cf_normals <- function(t) 0.99 * exp(-t[, 1L]^2 / 2) + 0.01 * exp(-(100 * t[, 1L])^2 / 2)
  cf_exps <- function(t) 0.5 / (1 - 1i * t[, 1L]) + 0.5 / (1 - 100i * t[, 1L])
  cf_stable <- function(a) function(t) exp(-abs(t[, 1L])^a)
  run <- function(cf) {
    collect_classes(from_cf(cf, n_dim = 1L, N = 2L, seed = 1L, n_freq = 512L, bracket = FALSE))
  }
  for (cf in list(cf_normals, cf_exps)) {
    r <- run(cf)
    expect_false("proxymix_cf_heavy_tail" %in% r$classes)
    expect_gt(min(r$value@diagnostics$tail_index), 1.9)
  }
  for (a in c(0.8, 1.5)) {
    r <- run(cf_stable(a))
    expect_true("proxymix_cf_heavy_tail" %in% r$classes)
    expect_equal(min(r$value@diagnostics$tail_index), a, tolerance = 1e-3)
  }
})

test_that("adapt = pmc stops late enough to match a long run by quadrature", {
  skip_on_cran()
  ## DAAG::orings: temperature minus 70, any damage on that flight.
  temp <- c(53, 57, 58, 63, 66, 67, 67, 67, 68, 69, 70, 70, 70, 70, 72, 73, 75,
            75, 76, 76, 78, 79, 81) - 70
  dam <- c(1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0)
  lp <- function(th) {
    if (is.null(dim(th))) th <- matrix(th, ncol = 2L)
    eta <- outer(rep(1, length(temp)), th[, 1L]) + outer(temp, th[, 2L])
    colSums(dam * eta - (pmax(eta, 0) + log1p(exp(-abs(eta)))))
  }
  op <- stats::optim(c(-1, -0.2), function(b) -lp(b), hessian = TRUE)
  sig <- solve(op$hessian)
  sds <- sqrt(diag(sig))
  g1 <- seq(op$par[1L] - 14 * sds[1L], op$par[1L] + 24 * sds[1L], length.out = 241L)
  g2 <- seq(op$par[2L] - 14 * sds[2L], op$par[2L] + 24 * sds[2L], length.out = 241L)
  gr <- as.matrix(expand.grid(g1, g2))
  lp_gr <- lp(gr)
  da <- (g1[2L] - g1[1L]) * (g2[2L] - g2[1L])
  log_z <- log(sum(exp(lp_gr)) * da)
  true_kl <- function(f) {
    keep <- lp_gr - log_z > -700
    sum(exp(lp_gr[keep] - log_z) * (lp_gr[keep] - log_z - dgmm(gr[keep, ], f, log = TRUE))) * da
  }
  tgt <- gmm_target(n_dim = 2L, log_density = lp, normalised = FALSE)
  q <- is_mvt(2L, mean = op$par, sigma = 4 * sig, df = 5)
  kl <- vapply(1:3, function(s) {
    f <- collect_classes(fit_kld_em(tgt, N = 3L, proposal = q, is_size = 5000L, seed = s,
                                    adapt = "pmc", max_iter = 3000L,
                                    validation_size = 0L))$value
    expect_true(f@converged)
    true_kl(f)
  }, numeric(1L))
  expect_lt(mean(kl), 0.0035)
})

test_that("from_cf_sources accepts a finite-variance source that mixes scales", {
  cf <- function(t) {
    s1 <- t[, 1L] + 0.5 * t[, 2L]
    s2 <- 0.5 * t[, 1L] + t[, 2L]
    (0.5 / (1 - 1i * s1) + 0.5 / (1 - 100i * s1)) * (1 / (1 + s2^2))
  }
  r <- collect_classes(from_cf_sources(cf, n_dim = 2L, N = 2L, seed = 1L))
  expect_type(r$value, "list")
})

test_that("gmm_impute warns on discrete columns and avoids a censored spike", {
  set.seed(5)
  obs <- stats::rnorm(60, 0, 1.5)
  low <- stats::quantile(obs, 0.3)
  y <- matrix(ifelse(obs < low, NA, obs), dimnames = list(NULL, "y"))
  im <- collect_classes(gmm_impute(y, N = NULL, m = 3L, seed = 1L,
                                   mechanism = censored("y", upper = low)))$value
  sds <- sqrt(unlist(im@point_fit@covariances))
  expect_gt(min(sds), 0.05)
  m <- cbind(y = stats::rnorm(200), b = stats::rbinom(200, 1L, 0.4))
  m[1:30, "y"] <- NA
  r <- collect_classes(gmm_impute(m, N = NULL, m = 3L, seed = 1L))
  expect_true("proxymix_discrete_columns" %in% r$classes)
  ridges <- vapply(c(1e-6, 1e-3), function(rd) {
    collect_classes(gmm_impute(m, N = NULL, m = 3L, seed = 1L, ridge_eps = rd))$value@n_components
  }, integer(1L))
  expect_identical(ridges[1L], ridges[2L])
})

test_that("fit_uplift warns on discrete covariates and not on the treatment", {
  set.seed(6)
  n <- 200L
  dat <- data.frame(t = stats::rbinom(n, 1L, 0.5), b = stats::rbinom(n, 1L, 0.4), x = stats::rnorm(n))
  dat$y <- dat$x + dat$t + stats::rnorm(n)
  with_b <- collect_classes(fit_uplift(dat, "y", "t", c("b", "x"), N = 1L, seed = 1L))
  without <- collect_classes(fit_uplift(dat, "y", "t", "x", N = 1L, seed = 1L))
  expect_true("proxymix_discrete_columns" %in% with_b$classes)
  expect_false("proxymix_discrete_columns" %in% without$classes)
})

test_that("uplift selection is reproducible across seeds and agrees with a fixed N", {
  set.seed(9)
  n <- 500L
  dat <- data.frame(t = stats::rbinom(n, 1L, 0.5), b = stats::rbinom(n, 1L, 0.3), x = stats::rnorm(n))
  dat$y <- dat$x + dat$t * (1 + dat$b) + stats::rnorm(n)
  auto <- collect_classes(fit_uplift(dat, "y", "t", c("b", "x"), N = "auto", n_grid = 1:3,
                                     regime = "sample", seed = 4L))$value
  fixed <- collect_classes(fit_uplift(dat, "y", "t", c("b", "x"), N = gmm_n_components(auto@fit),
                                      regime = "sample", seed = 4L))$value
  expect_equal(gmm_n_components(fixed@fit), gmm_n_components(auto@fit))
  ## The score does not run away with the ridge.
  expect_gt(min(auto@metadata$bic_trace), -1e5 / 10)
})

test_that("fit_uplift reaches the true effects from seeds that used to miss them", {
  skip_on_cran()
  set.seed(20261107)
  n <- 3000L
  dat <- data.frame(x = stats::rnorm(n), b1 = stats::rbinom(n, 1L, 0.4),
                    b2 = stats::rbinom(n, 1L, 0.3), t = stats::rbinom(n, 1L, 0.5))
  dat$y <- 1 + dat$x + dat$b1 + dat$t * (0.5 + dat$b1) + stats::rnorm(n, sd = 0.5)
  truth <- 0.5 + dat$b1
  covs <- c("x", "b1", "b2")
  dist <- vapply(c(11L, 18L, 19L, 20L), function(s) {
    m <- collect_classes(fit_uplift(dat, "y", "t", covs, N = 4L, regime = "sample", seed = s))$value
    tau <- collect_classes(proxy_cate(m, newdata = dat[, covs], se = FALSE))$value$tau
    sqrt(mean((tau - truth)^2))
  }, numeric(1L))
  expect_lt(max(dist), 0.06)
})

test_that("proxy_overlap flags a profile whose components hold one arm", {
  ## Treatment is determined by the covariate, so positivity fails by design.
  set.seed(2)
  n <- 400L
  x <- stats::rnorm(n)
  t <- as.numeric(x > 0)
  dat <- data.frame(y = x + t + stats::rnorm(n, sd = 0.3), t = t, x = x)
  m <- collect_classes(fit_uplift(dat, "y", "t", "x", N = 2L, regime = "sample", seed = 1L))$value
  ov <- collect_classes(proxy_overlap(m, data.frame(x = c(-2, 2))))$value
  expect_true(all(ov$overlap_flag))
  ## Randomised treatment: both arms supported, propensity close to the truth.
  t2 <- stats::rbinom(n, 1L, 0.5)
  dat2 <- data.frame(y = x + t2 + stats::rnorm(n, sd = 0.3), t = t2, x = x)
  m2 <- collect_classes(fit_uplift(dat2, "y", "t", "x", N = 1L, regime = "sample", seed = 1L))$value
  ov2 <- collect_classes(proxy_overlap(m2, data.frame(x = c(-1, 0, 1))))$value
  expect_false(any(ov2$overlap_flag))
  expect_true(all(abs(ov2$propensity - 0.5) < 0.1))
})

test_that("mutual information offers the Shannon value, checked in closed form", {
  rho <- 0.7
  g <- gmm(weights = 1, means = list(c(0, 0)),
           covariances = list(matrix(c(1, rho, rho, 1), 2L)))
  r <- gmm_mutual_information(g, 1L, 2L, measure = "shannon", n_mc = 40000L, seed = 1L)
  truth <- -0.5 * log(1 - rho^2)
  expect_lt(abs(r$mc - truth), 4 * r$mc_se)
  expect_lt(gmm_mutual_information(g, 1L, 2L), truth)
})

test_that("gmm_modes filters by height and keeps every mode by default", {
  g <- gmm(weights = c(0.99, 0.01), means = list(0, 8),
           covariances = list(matrix(1), matrix(0.05)))
  expect_equal(gmm_modes(g)$n, 2L)
  kept <- gmm_modes(g, min_height = 0.05)
  expect_equal(kept$n, 1L)
  expect_equal(as.numeric(kept$modes), 0, tolerance = 1e-3)
  expect_error(gmm_modes(g, min_height = 2), "between 0 and 1")
})

test_that("a negative variational divergence raises a classed warning", {
  p <- gmm(weights = c(0.5, 0.5), means = list(163, 176),
           covariances = list(matrix(50), matrix(50)))
  q <- gmm(weights = 1, means = list(169.5), covariances = list(matrix(92)))
  expect_warning(gmm_kld(p, q, n_mc = 200L), class = "proxymix_negative_variational")
})
