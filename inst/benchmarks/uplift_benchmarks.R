# uplift_benchmarks.R -- External benchmarks for the proxymix decision layer.
#
# Grades the closed-form uplift reader against independent implementations
# (grf causal forest, S/T/X meta-learners, DoubleML) on synthetic ground truth
# and on two public uplift datasets (Hillstrom, Criteo). Not part of the
# package check surface and not a package dependency: it requires `grf` and,
# for the DoubleML ATE baseline, `DoubleML` + `mlr3` + `mlr3learners` +
# `ranger` to be installed separately, and reads data downloaded outside the
# package tree. The DoubleML baseline degrades to `NA` if those are absent.
#
# Data directory (no hard-coded paths): set the environment variable
# PROXYMIX_BENCH_DATA to the folder holding `hillstrom.csv` and
# `criteo-uplift-v2.1.csv.gz`. Defaults to "benchmarks_data".
#
# Metrics. Qini coefficient (ranking quality, higher is better) and AUUC for
# the real datasets; PEHE (lower is better) against known ground truth for the
# synthetic set; per-unit scoring throughput for the speed comparison.

suppressMessages({
  library(proxymix)
  library(grf)
})

data_dir <- Sys.getenv("PROXYMIX_BENCH_DATA", unset = "benchmarks_data")

# Metrics -------------------------------------------------------------------

## Qini curve value at every prefix when units are ranked by `uplift`.
qini_curve <- function(uplift, treat, y) {
  o <- order(uplift, decreasing = TRUE)
  treat <- treat[o]
  y <- y[o]
  cum_t <- cumsum(y * treat)
  cum_c <- cumsum(y * (1 - treat))
  n_t <- cumsum(treat)
  n_c <- cumsum(1 - treat)
  ratio <- ifelse(n_c > 0, n_t / n_c, 0)
  cum_t - cum_c * ratio
}

## Qini coefficient: area between the Qini curve and the random-targeting line,
## normalised by the number of units. Higher is better.
qini_coef <- function(uplift, treat, y) {
  q <- qini_curve(uplift, treat, y)
  n <- length(q)
  rand <- seq_len(n) / n * q[n]
  sum(q - rand) / n
}

## Area under the uplift curve (AUUC), normalised likewise.
auuc <- function(uplift, treat, y) {
  q <- qini_curve(uplift, treat, y)
  sum(q) / length(q)
}

# Meta-learners (grf regression forests as base learners) -------------------

.rf <- function(x, y) grf::regression_forest(as.matrix(x), y, num.trees = 500L)

learner_T <- function(xtr, ttr, ytr, xte) {
  m1 <- .rf(xtr[ttr == 1, , drop = FALSE], ytr[ttr == 1])
  m0 <- .rf(xtr[ttr == 0, , drop = FALSE], ytr[ttr == 0])
  as.numeric(predict(m1, as.matrix(xte))$predictions) -
    as.numeric(predict(m0, as.matrix(xte))$predictions)
}

learner_S <- function(xtr, ttr, ytr, xte) {
  m <- .rf(cbind(xtr, t = ttr), ytr)
  p1 <- predict(m, as.matrix(cbind(xte, t = 1)))$predictions
  p0 <- predict(m, as.matrix(cbind(xte, t = 0)))$predictions
  as.numeric(p1) - as.numeric(p0)
}

learner_X <- function(xtr, ttr, ytr, xte, e = 0.5) {
  m1 <- .rf(xtr[ttr == 1, , drop = FALSE], ytr[ttr == 1])
  m0 <- .rf(xtr[ttr == 0, , drop = FALSE], ytr[ttr == 0])
  d1 <- ytr[ttr == 1] -
    as.numeric(predict(m0, as.matrix(xtr[ttr == 1, , drop = FALSE]))$predictions)
  d0 <- as.numeric(predict(m1, as.matrix(xtr[ttr == 0, , drop = FALSE]))$predictions) -
    ytr[ttr == 0]
  t1 <- .rf(xtr[ttr == 1, , drop = FALSE], d1)
  t0 <- .rf(xtr[ttr == 0, , drop = FALSE], d0)
  tau1 <- as.numeric(predict(t1, as.matrix(xte))$predictions)
  tau0 <- as.numeric(predict(t0, as.matrix(xte))$predictions)
  e * tau0 + (1 - e) * tau1
}

# proxymix uplift scoring ---------------------------------------------------

proxymix_uplift <- function(train, xte, outcome_type = "continuous",
                            N = 4L, seed = 1L) {
  cov_names <- setdiff(names(train), c("y", "t"))
  m <- fit_uplift(train, "y", "t", cov_names, N = N, regime = "sample",
                  outcome_type = outcome_type, max_iter = 150L, seed = seed)
  nd <- as.data.frame(xte)
  names(nd) <- cov_names
  list(model = m,
       tau = proxy_cate(m, nd, se = FALSE)$tau)
}

# Real-dataset benchmark (Qini / AUUC vs grf and learners; ATE vs DoubleML) --

bench_real <- function(train, test, label, outcome_type = "binary", N = 4L) {
  cov_names <- setdiff(names(train), c("y", "t"))
  xtr <- train[cov_names]
  xte <- test[cov_names]
  ttr <- train$t
  ytr <- train$y

  cat(sprintf("\n[%s] train n=%d, test n=%d, %d covariates\n",
              label, nrow(train), nrow(test), length(cov_names)))

  scores <- list()
  ## proxymix
  pm <- proxymix_uplift(train, xte, outcome_type = outcome_type, N = N)
  scores[["proxymix"]] <- pm$tau
  ## grf causal forest
  cf <- grf::causal_forest(as.matrix(xtr), ytr, ttr, num.trees = 1000L)
  scores[["grf"]] <- as.numeric(predict(cf, as.matrix(xte))$predictions)
  ## meta-learners
  scores[["T-learner"]] <- learner_T(xtr, ttr, ytr, xte)
  scores[["S-learner"]] <- learner_S(xtr, ttr, ytr, xte)
  scores[["X-learner"]] <- learner_X(xtr, ttr, ytr, xte)

  qini <- vapply(scores, function(s) qini_coef(s, test$t, test$y), numeric(1L))
  au <- vapply(scores, function(s) auuc(s, test$t, test$y), numeric(1L))

  ## DoubleML interactive regression model: a debiased ATE oracle.
  ate <- tryCatch({
    if (requireNamespace("DoubleML", quietly = TRUE) &&
        requireNamespace("mlr3", quietly = TRUE) &&
        requireNamespace("mlr3learners", quietly = TRUE)) {
      df <- data.frame(xtr, y = ytr, t = ttr)
      dml_data <- DoubleML::double_ml_data_from_data_frame(
        df, y_col = "y", d_cols = "t", x_cols = cov_names)
      lrn <- mlr3::lrn("regr.ranger", num.trees = 200L)
      lrn_c <- mlr3::lrn("classif.ranger", num.trees = 200L,
                         predict_type = "prob")
      obj <- DoubleML::DoubleMLIRM$new(dml_data, ml_g = lrn, ml_m = lrn_c)
      suppressMessages(obj$fit())
      obj$coef
    } else NA_real_
  }, error = function(e) NA_real_)

  out <- data.frame(
    method = names(qini),
    qini = round(qini, 5),
    auuc = round(au, 5),
    row.names = NULL
  )
  attr(out, "ate_proxymix") <- mean(pm$tau)
  attr(out, "ate_dml") <- ate
  attr(out, "ate_grf") <- mean(scores[["grf"]])
  print(out[order(-out$qini), ])
  cat(sprintf("[%s] ATE: proxymix=%.5f  grf=%.5f  DoubleML(IRM)=%.5f\n",
              label, mean(pm$tau), mean(scores[["grf"]]), ate))
  invisible(out)
}

# Synthetic PEHE benchmark (proxymix vs grf, ground truth) ------------------

bench_synthetic <- function(seed = 99L) {
  set.seed(seed)
  dgps <- list(
    linear      = function(x) 0.5 + 0 * x,
    crossing    = function(x) 2 * x,
    nonlinear   = function(x) 1 / (1 + exp(-2 * x)),
    interaction = function(x) 0.4 + 1.2 * x
  )
  rows <- lapply(names(dgps), function(nm) {
    f <- dgps[[nm]]
    n <- 4000L
    x <- rnorm(n)
    t <- rbinom(n, 1L, 0.5)
    y <- 1 + 0.3 * x + f(x) * t + rnorm(n, sd = 0.5)
    tr <- seq_len(3000L)
    te <- setdiff(seq_len(n), tr)
    train <- data.frame(y = y[tr], t = t[tr], x = x[tr])
    truth <- f(x[te])

    pm <- proxymix_uplift(train, data.frame(x = x[te]), N = 3L)
    cf <- grf::causal_forest(matrix(x[tr], ncol = 1L), y[tr], t[tr],
                             num.trees = 1000L)
    tau_cf <- as.numeric(predict(cf, matrix(x[te], ncol = 1L))$predictions)

    data.frame(
      dgp = nm,
      pehe_proxymix = sqrt(mean((pm$tau - truth)^2)),
      pehe_grf = sqrt(mean((tau_cf - truth)^2))
    )
  })
  res <- do.call(rbind, rows)
  res$ratio <- round(res$pehe_proxymix / res$pehe_grf, 3)
  res$pehe_proxymix <- round(res$pehe_proxymix, 4)
  res$pehe_grf <- round(res$pehe_grf, 4)
  cat("\n== Synthetic PEHE (proxymix vs grf, known ground truth) ==\n")
  print(res, row.names = FALSE)
  invisible(res)
}

# Speed benchmark -----------------------------------------------------------

bench_speed <- function(seed = 7L, n_score = 1e5L) {
  set.seed(seed)
  n <- 4000L
  x <- rnorm(n)
  t <- rbinom(n, 1L, 0.5)
  y <- 1 + (0.4 + x) * t + rnorm(n, sd = 0.5)
  train <- data.frame(y = y, t = t, x = x)
  m <- fit_uplift(train, "y", "t", "x", N = 3L, regime = "sample",
                  max_iter = 120L, seed = 1L)
  cf <- grf::causal_forest(matrix(x, ncol = 1L), y, t, num.trees = 1000L)

  xs <- data.frame(x = rnorm(n_score))
  t_pm <- system.time(proxy_cate(m, xs, se = FALSE))[["elapsed"]]
  t_cf <- system.time(predict(cf, matrix(xs$x, ncol = 1L)))[["elapsed"]]
  cat(sprintf("\n== Scoring %d units ==\n  proxymix: %.3fs (%.0f units/s)\n  grf:      %.3fs (%.0f units/s)\n  speedup:  %.1fx\n",
              n_score, t_pm, n_score / t_pm, t_cf, n_score / t_cf, t_cf / t_pm))
  invisible(c(proxymix_s = t_pm, grf_s = t_cf, speedup = t_cf / t_pm))
}

# Hillstrom loader ----------------------------------------------------------

load_hillstrom <- function(path) {
  d <- utils::read.csv(path, stringsAsFactors = FALSE)
  ## Binary treatment: Womens E-Mail (1) vs No E-Mail (0); drop Mens arm.
  d <- d[d$segment %in% c("Womens E-Mail", "No E-Mail"), ]
  t <- as.integer(d$segment == "Womens E-Mail")
  zip <- model.matrix(~ factor(zip_code) - 1, d)
  chan <- model.matrix(~ factor(channel) - 1, d)
  X <- data.frame(
    recency = d$recency, history = d$history,
    mens = d$mens, womens = d$womens, newbie = d$newbie,
    zip[, -1, drop = FALSE], chan[, -1, drop = FALSE]
  )
  names(X) <- make.names(names(X))
  data.frame(y = d$visit, t = t, X)
}

# Criteo loader (subsampled -- 25M rows is intractable to fit whole) ---------

load_criteo <- function(path, n_sub = 250000L, seed = 1L) {
  d <- data.table::fread(cmd = paste("gzcat", shQuote(path)))
  set.seed(seed)
  idx <- sample.int(nrow(d), min(n_sub, nrow(d)))
  d <- d[idx, ]
  cat(sprintf("[criteo] subsampled %d of %d rows (logged, not silent)\n",
              nrow(d), 25309483L))
  feats <- grep("^f", names(d), value = TRUE)
  out <- data.frame(y = d$visit, t = d$treatment, as.data.frame(d[, ..feats]))
  out
}

# Runner --------------------------------------------------------------------

run_all <- function(which = c("synthetic", "speed", "hillstrom", "criteo")) {
  res <- list()
  if ("synthetic" %in% which) res$synthetic <- bench_synthetic()
  if ("speed" %in% which) res$speed <- bench_speed()
  if ("hillstrom" %in% which) {
    h <- load_hillstrom(file.path(data_dir, "hillstrom.csv"))
    set.seed(1L)
    tr <- sample.int(nrow(h), floor(0.7 * nrow(h)))
    res$hillstrom <- bench_real(h[tr, ], h[-tr, ], "hillstrom",
                                outcome_type = "binary", N = 4L)
  }
  if ("criteo" %in% which) {
    cr <- load_criteo(file.path(data_dir, "criteo-uplift-v2.1.csv.gz"))
    set.seed(1L)
    tr <- sample.int(nrow(cr), floor(0.6 * nrow(cr)))
    res$criteo <- bench_real(cr[tr, ], cr[-tr, ], "criteo",
                             outcome_type = "binary", N = 4L)
  }
  invisible(res)
}

if (sys.nframe() == 0L) {
  run_all()
}
