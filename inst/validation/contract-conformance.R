## Durable two-sided conformance case set for proxymix (the authored, non-
## regenerable validation IP). Each function carries POSITIVE cases (valid /
## extreme-but-valid input runs, checked against a cheap oracle) and NEGATIVE
## cases (malformed input is REJECTED with a clean error), plus metamorphic
## invariants. This file is the durable home (ships with the package, tracked on
## the package remote); it is executed by the /pkg-validation dossier, which adds
## the verdict engine and writes the api-coverage gate the release relies on.
##
## Contract-conformance sweep -- two-sided "corner to corner" over the documented
## function surface of the sixteen capabilities. For each function:
##   * POSITIVE cases  -- a valid input (including extreme-but-valid: single
##                        component, tiny variance, large N) runs, and where an
##                        oracle is cheap the result is checked against it;
##   * NEGATIVE cases  -- a malformed, out-of-domain, wrong-type or NA/Inf input is
##                        REJECTED with a clean error, never accepted silently.
## Plus metamorphic invariants (identity maps, composition, log/exp, reduction
## mean-preservation, same-seed reproducibility). Writes
## bench/results/contract-conformance.rds for the AI-auditor dossier
## (audit_docs.R reads $verdicts/$summary/$metamorphic/$findings/$n_*).
##
## Run from report/ :  Rscript bench/contract-conformance.R
## The independent oracles here are base-R re-derivations and the verdict engine
## itself (pv_grid/pv_verdict/pv_metamorphic), not the package under study.

## Single-threaded BLAS so the forked sweep is fork-safe: macOS Accelerate (and
## some OpenBLAS builds) are not fork-safe, and a Cholesky/backsolve in a forked
## worker can abort the worker. Set before any matrix op or fork.
Sys.setenv(VECLIB_MAXIMUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1",
           OMP_NUM_THREADS = "1")

if (!"proxymix" %in% loadedNamespaces()) suppressMessages(library(proxymix))
## The two-sided sweep is driven by the /pkg-validation dossier (which sources the
## verdict engine: pv_grid / pv_verdict / pv_assemble_conformance). Sourced
## standalone, this file defines the authored case set for inspection.
.pv_engine <- path.expand("~/.claude/skills/pkg-validation/engine/conformance.R")
if (file.exists(.pv_engine)) source(.pv_engine)
set.seed(20260620)
RESDIR <- "bench/results"; if (!dir.exists(RESDIR)) dir.create(RESDIR, recursive = TRUE)

## ---- valid fixtures --------------------------------------------------------
g1  <- gmm(weights = c(0.6, 0.4), means = list(-1, 2),
           covariances = list(matrix(1), matrix(0.5)))                # 1-D, K = 2
g2  <- gmm(weights = c(0.5, 0.5), means = list(c(0, 0), c(2, 1)),
           covariances = list(diag(2), diag(2)))                      # 2-D, K = 2
g1a <- gmm(weights = 1, means = list(0), covariances = list(matrix(1)))   # single comp
g4  <- gmm(weights = rep(0.25, 4), means = list(-3, -1, 1, 3),
           covariances = rep(list(matrix(0.3)), 4))                   # 1-D, K = 4

tgt1 <- gmm_target(n_dim = 1L,
  log_density = function(x) { xx <- if (is.matrix(x)) x[, 1] else x; dnorm(xx, 0, 1.5, log = TRUE) },
  normalised = TRUE)
Xs  <- matrix(rnorm(200), ncol = 1L); st1 <- gmm_target_from_samples(Xs)
zz  <- sample(1:2, 200L, TRUE); X2 <- cbind(rnorm(200L, c(-2, 2)[zz]), rnorm(200L, c(0, 1)[zz]))
colnames(X2) <- c("x1", "x2"); st2 <- gmm_target_from_samples(X2)

fit_s <- fit_em_samples(st2, N = 2L, n_starts = 3L, seed = 1L)
fit_k <- fit_kld_em(tgt1, N = 2L, proposal = is_mvn(1L, 0, matrix(9)),
                    is_size = 600L, seed = 1L, support_warn = FALSE)

Xi <- X2; Xi[runif(nrow(Xi)) < 0.25, "x2"] <- NA
imp <- gmm_impute(Xi, m = 5L, seed = 1L)

du <- data.frame(x1 = rnorm(200L), t = rbinom(200L, 1, 0.5))
du$y <- 1 + 0.8 * du$x1 + 1.5 * du$t + rnorm(200L, 0, 0.5)
um <- fit_uplift(du, outcome = "y", treatment = "t", covariates = "x1")

prior <- g1a
dyn_f <- list(A = matrix(1), b = 0, Q = matrix(0.1))
meas_f <- list(C = matrix(1), R = matrix(0.5))
yf <- matrix(cumsum(rnorm(30, 0, 0.3)) + rnorm(30, 0, 0.7), ncol = 1L)
dyn_e <- list(A = matrix(1), Q = matrix(0.04)); meas_e <- list(C = matrix(1), R = matrix(1))
ye <- matrix(c(rnorm(60), 5), ncol = 1L)

## -- extended-surface fixtures (for the uncovered functions, v0.11.2 refresh) -
bt  <- banana_target()                                   # a built-in target
g3  <- gmm(weights = 1, means = list(c(0, 0, 0)),
           covariances = list(diag(3)))                  # 3-D single component
gK  <- gmm(weights = c(0.6, 0.4), means = list(c(0, 0, 0), c(2, 1, 1)),
           covariances = list(diag(3), diag(3)))         # 3-D, K = 2 (causal joint)
cf  <- gmm_counterfactual(g3, evidence = c(1, 0, 0.2),
                          do = c(NA, 1, NA), query = 1L)  # a counterfactual law
mfb <- mock_fb_posterior(shape = "gaussian", n_dim = 2L) # a mock FB producer

## MNAR data with a genuine MNAR NA in column "y", taken verbatim from
## proxy_mnar_sensitivity()'s own @examples (its missingness probability
## depends on the unobserved y itself).
set.seed(1L)
mnar_x1 <- rnorm(300L)
mnar_y  <- mnar_x1 + rnorm(300L)
mnar_y[runif(300L) < plogis(-0.4 + 0.8 * mnar_y)] <- NA
mnar_dat <- data.frame(x1 = mnar_x1, y = mnar_y)
set.seed(20260620)

## global mean of a mixture (independent re-derivation, for oracles)
gmean <- function(g) Reduce(`+`, Map(function(w, m) w * m, g@weights, g@means))

## ---- case constructors -----------------------------------------------------
ok_case  <- function(fn, label, call, check = NULL)
  list(fn = fn, label = label, call = call, expect = "ok", check = check)
err_case <- function(fn, label, call)
  list(fn = fn, label = label, call = call, expect = "error", check = NULL)

CASES <- list(

  ## -- regimes: fit_proxymix / fit_moment_match / fit_em_samples / fit_kld_em --
  ok_case("fit_proxymix", "fit_proxymix N=1 moment",
          function() fit_proxymix(st1, N = 1L, regime = "moment")),
  ok_case("fit_proxymix", "fit_proxymix auto on samples",
          function() fit_proxymix(st2, N = 2L, regime = "sample")),
  ok_case("fit_moment_match", "fit_moment_match recovers sample mean",
          function() fit_moment_match(st1, N = 1L),
          check = function(v) abs(v@means[[1]] - mean(Xs)) < 1e-6),
  ok_case("fit_em_samples", "fit_em_samples K=1 extreme",
          function() fit_em_samples(st2, N = 1L, seed = 1L),
          check = function(v) length(v@weights) == 1L),
  ok_case("fit_em_samples", "fit_em_samples K=2",
          function() fit_em_samples(st2, N = 2L, n_starts = 2L, seed = 1L),
          check = function(v) length(v@weights) == 2L),
  ok_case("fit_kld_em", "fit_kld_em evaluate-only target",
          function() fit_kld_em(tgt1, N = 2L, proposal = is_mvn(1L, 0, matrix(9)),
                                is_size = 400L, seed = 1L, support_warn = FALSE),
          check = function(v) length(v@weights) == 2L),
  err_case("fit_proxymix", "fit_proxymix non-target", function() fit_proxymix(42)),
  err_case("fit_proxymix", "fit_proxymix N=0", function() fit_proxymix(tgt1, N = 0L)),
  err_case("fit_proxymix", "fit_proxymix N=-2", function() fit_proxymix(tgt1, N = -2L)),
  err_case("fit_moment_match", "fit_moment_match char target",
           function() fit_moment_match("x", N = 1L)),
  err_case("fit_em_samples", "fit_em_samples N=0", function() fit_em_samples(st2, N = 0L)),
  err_case("fit_em_samples", "fit_em_samples char", function() fit_em_samples("x", N = 2L)),
  err_case("fit_kld_em", "fit_kld_em is_size=0",
           function() fit_kld_em(tgt1, N = 2L, is_size = 0L, seed = 1L)),
  err_case("fit_kld_em", "fit_kld_em non-target", function() fit_kld_em(42, N = 2L)),
  err_case("fit_kld_em", "fit_kld_em N=-1",
           function() fit_kld_em(tgt1, N = -1L, is_size = 200L, seed = 1L)),

  ## -- gmm constructor (out-of-domain) ----------------------------------------
  ok_case("gmm", "gmm valid 2-D", function() gmm(c(0.5, 0.5), list(0, 1),
                                                 list(matrix(1), matrix(1)))),
  err_case("gmm", "gmm negative weight",
           function() gmm(c(-0.5, 1.5), list(0, 1), list(matrix(1), matrix(1)))),
  err_case("gmm", "gmm negative variance",
           function() gmm(c(0.5, 0.5), list(0, 1), list(matrix(-1), matrix(1)))),
  err_case("gmm", "gmm length mismatch",
           function() gmm(c(0.5, 0.3, 0.2), list(0, 1), list(matrix(1), matrix(1)))),

  ## -- affine operator calculus: gmm_affine / gmm_aggregate / gmm_missing / dgmm / rgmm
  ok_case("gmm_affine", "gmm_affine scale+shift (oracle: A mu + b)",
          function() gmm_affine(g2, matrix(c(2, 0, 0, 2), 2), b = c(1, 1)),
          check = function(v) max(abs(unlist(v@means) - (2 * unlist(g2@means) + 1))) < 1e-9),
  ok_case("gmm_affine", "gmm_affine identity",
          function() gmm_affine(g2, diag(2)),
          check = function(v) max(abs(unlist(v@means) - unlist(g2@means))) < 1e-12),
  ok_case("gmm_aggregate", "gmm_aggregate sum channel",
          function() gmm_aggregate(g2, matrix(c(1, 1), 1, 2))),
  ok_case("dgmm", "dgmm = analytic mixture density",
          function() dgmm(matrix(c(0, 0), 1, 2), g2),
          check = function(v) abs(v - (0.5 * dnorm(0)^2 +
                                       0.5 * dnorm(2) * dnorm(1))) < 1e-9),
  ok_case("rgmm", "rgmm draws right shape",
          function() rgmm(10L, g2),
          check = function(v) is.matrix(v) && all(dim(v) == c(10L, 2L))),
  err_case("gmm_affine", "gmm_affine non-gmm", function() gmm_affine(42, diag(2))),
  err_case("gmm_affine", "gmm_affine char A", function() gmm_affine(g2, "A")),
  err_case("gmm_affine", "gmm_affine non-conformable A",
           function() gmm_affine(g2, matrix(1, 3, 3))),
  err_case("gmm_affine", "gmm_affine NA in A",
           function() gmm_affine(g2, matrix(NA_real_, 2, 2))),
  err_case("gmm_aggregate", "gmm_aggregate non-gmm",
           function() gmm_aggregate(42, matrix(c(1, 1), 1, 2))),
  err_case("dgmm", "dgmm non-gmm", function() dgmm(matrix(0, 1, 2), 42)),
  err_case("dgmm", "dgmm char x", function() dgmm("x", g2)),
  err_case("rgmm", "rgmm n=-1", function() rgmm(-1L, g2)),
  err_case("rgmm", "rgmm non-gmm", function() rgmm(10L, 42)),
  err_case("rgmm", "rgmm n=NA", function() rgmm(NA_integer_, g2)),
  ok_case("gmm_missing", "gmm_missing impute one coord",
          function() gmm_missing(g2, observed = 2L, values = 1)),
  err_case("gmm_missing", "gmm_missing non-gmm",
           function() gmm_missing(42, c(TRUE, FALSE), 1)),
  err_case("gmm_missing", "gmm_missing observed length mismatch",
           function() gmm_missing(g2, observed = c(TRUE, FALSE, TRUE), values = 1)),

  ## -- regression / conditioning: gmm_conditionalise / gmm_marginalise --------
  ok_case("gmm_conditionalise", "gmm_conditionalise -> gmm",
          function() gmm_conditionalise(g2, given = c(NA, 0.5)),
          check = function(v) inherits(v, class(g2)[1])),
  ok_case("gmm_marginalise", "gmm_marginalise keep coord 1",
          function() gmm_marginalise(g2, keep = 1L),
          check = function(v) abs(v@means[[1]] - g2@means[[1]][1]) < 1e-12),
  err_case("gmm_conditionalise", "gmm_conditionalise given wrong length",
           function() gmm_conditionalise(g2, given = c(1, 2, 3))),
  err_case("gmm_conditionalise", "gmm_conditionalise non-gmm",
           function() gmm_conditionalise(42, given = c(NA, 1))),
  err_case("gmm_marginalise", "gmm_marginalise keep out of range",
           function() gmm_marginalise(g2, keep = 5L)),
  err_case("gmm_marginalise", "gmm_marginalise keep=0",
           function() gmm_marginalise(g2, keep = 0L)),

  ## -- pca: fit_proxymix single component (covariance eigenbasis) -------------
  ok_case("fit_proxymix", "fit_proxymix K=1 covariance",
          function() fit_proxymix(st2, N = 1L, regime = "moment"),
          check = function(v) length(v@weights) == 1L),

  ## -- bayes observation channel: gmm_observe / from_fb_posterior -------------
  ok_case("gmm_observe", "gmm_observe linear-Gaussian update",
          function() gmm_observe(g2, matrix(c(1, 0), 1, 2), y = 0.5,
                                 noise_cov = matrix(0.2)),
          check = function(v) inherits(v, class(g2)[1])),
  err_case("gmm_observe", "gmm_observe non-gmm",
           function() gmm_observe(42, matrix(c(1, 0), 1, 2), y = 0.5, noise_cov = matrix(0.2))),
  err_case("gmm_observe", "gmm_observe y=NA",
           function() gmm_observe(g2, matrix(c(1, 0), 1, 2), y = NA_real_, noise_cov = matrix(0.2))),
  ok_case("from_fb_posterior", "from_fb_posterior consumes a mock producer",
          function() from_fb_posterior(mock_fb_posterior("gaussian", n_dim = 2L),
                                       N = 2L, is_size = 800L, seed = 1L)),
  err_case("from_fb_posterior", "from_fb_posterior non-producer",
           function() from_fb_posterior(42)),

  ## -- information: entropy / divergence / mutual information / conditional ----
  ok_case("gmm_entropy", "gmm_entropy finite",
          function() gmm_entropy(g1), check = function(v) is.finite(v)),
  ok_case("gmm_divergence", "gmm_divergence self ~ 0",
          function() gmm_divergence(g1, g1, type = "cs"),
          check = function(v) abs(v) < 1e-6),
  ok_case("gmm_mutual_information", "gmm_mutual_information finite",
          function() gmm_mutual_information(g2, block_a = 1L, block_b = 2L),
          check = function(v) is.finite(v)),
  ok_case("gmm_conditional_entropy", "gmm_conditional_entropy finite",
          function() gmm_conditional_entropy(g2, given = c(NA, 0)),
          check = function(v) is.finite(v)),
  ok_case("gmm_independence_graph", "gmm_independence_graph symmetric 0/1",
          function() gmm_independence_graph(g2),
          check = function(v) is.matrix(v) && isSymmetric(unname(v)) &&
            all(v %in% c(0L, 1L))),
  err_case("gmm_independence_graph", "gmm_independence_graph non-gmm",
           function() gmm_independence_graph(42)),
  err_case("gmm_independence_graph", "gmm_independence_graph bad threshold",
           function() gmm_independence_graph(g2, threshold = -1)),
  err_case("gmm_entropy", "gmm_entropy non-gmm", function() gmm_entropy(42)),
  err_case("gmm_entropy", "gmm_entropy bad order",
           function() gmm_entropy(g1, order = "banana")),
  err_case("gmm_divergence", "gmm_divergence non-gmm q",
           function() gmm_divergence(g1, 42, type = "cs")),
  err_case("gmm_divergence", "gmm_divergence bad type",
           function() gmm_divergence(g1, g1, type = "zzz")),
  err_case("gmm_mutual_information", "gmm_mutual_information block out of range",
           function() gmm_mutual_information(g2, block_a = 1L, block_b = 5L)),
  err_case("gmm_conditional_entropy", "gmm_conditional_entropy non-gmm",
           function() gmm_conditional_entropy(42, given = 1L)),

  ## -- filtering: gmm_filter / gmm_reduce -------------------------------------
  ok_case("gmm_filter", "gmm_filter Kalman track",
          function() gmm_filter(prior, dyn_f, meas_f, yf, k_max = 4L),
          check = function(v) length(v$mean) == nrow(yf)),
  ok_case("gmm_reduce", "gmm_reduce 4->2 preserves global mean",
          function() gmm_reduce(g4, k_max = 2L),
          check = function(v) length(v@weights) <= 2L &&
            max(abs(gmean(v) - gmean(g4))) < 1e-8),
  err_case("gmm_filter", "gmm_filter non-gmm prior",
           function() gmm_filter(42, dyn_f, meas_f, yf)),
  err_case("gmm_filter", "gmm_filter char y",
           function() gmm_filter(prior, dyn_f, meas_f, "x")),
  err_case("gmm_reduce", "gmm_reduce k_max=0", function() gmm_reduce(g4, k_max = 0L)),
  err_case("gmm_reduce", "gmm_reduce non-gmm", function() gmm_reduce(42, k_max = 2L)),
  err_case("gmm_reduce", "gmm_reduce bad cost",
           function() gmm_reduce(g4, k_max = 2L, cost = "zzz")),

  ## -- density compression: from_kde ------------------------------------------
  ok_case("from_kde", "from_kde 2-D compress",
          function() from_kde(X2[1:80, ], N = 2L, is_size = 400L, seed = 1L)),
  err_case("from_kde", "from_kde char", function() from_kde("x", N = 3L)),
  err_case("from_kde", "from_kde p>5 scope guard",
           function() from_kde(matrix(rnorm(60), ncol = 6L), N = 2L, seed = 1L)),
  err_case("from_kde", "from_kde N=0", function() from_kde(X2[1:80, ], N = 0L)),

  ## -- optimisation: from_objective / gmm_modes -------------------------------
  ok_case("from_objective", "from_objective bimodal",
          function() from_objective(function(v) (v[1]^2 - 4)^2, lower = -5, upper = 5,
                                    N = 4L, is_size = 400L, seed = 1L)),
  ok_case("gmm_modes", "gmm_modes returns modes",
          function() gmm_modes(fit_s),
          check = function(v) is.list(v) && !is.null(v$modes)),
  err_case("from_objective", "from_objective char objective",
           function() from_objective("f", lower = -1, upper = 1)),
  err_case("from_objective", "from_objective lower>upper",
           function() from_objective(function(v) v[1]^2, lower = 1, upper = -1, seed = 1L)),
  err_case("gmm_modes", "gmm_modes non-object", function() gmm_modes(42)),

  ## -- annealing: gmm_anneal_path ---------------------------------------------
  ok_case("gmm_anneal_path", "gmm_anneal_path discovers K",
          function() gmm_anneal_path(X2[1:120, ], k_max = 3L, n_steps = 20L,
                                     n_inner = 10L, seed = 1L)),
  err_case("gmm_anneal_path", "gmm_anneal_path char", function() gmm_anneal_path("x")),
  err_case("gmm_anneal_path", "gmm_anneal_path k_max=0",
           function() gmm_anneal_path(X2[1:120, ], k_max = 0L)),

  ## -- instability: gmm_eos_test ----------------------------------------------
  ok_case("gmm_eos_test", "gmm_eos_test chisq",
          function() gmm_eos_test(prior, dyn_e, meas_e, ye, m = 1L, method = "chisq"),
          check = function(v) is.finite(v$p_value) && v$p_value >= 0 && v$p_value <= 1),
  ok_case("gmm_eos_test", "gmm_eos_test andrews",
          function() gmm_eos_test(prior, dyn_e, meas_e, ye, m = 1L, method = "andrews"),
          check = function(v) is.finite(v$p_value) && v$p_value >= 0 && v$p_value <= 1),
  err_case("gmm_eos_test", "gmm_eos_test m=0",
           function() gmm_eos_test(prior, dyn_e, meas_e, ye, m = 0L)),
  err_case("gmm_eos_test", "gmm_eos_test bad method",
           function() gmm_eos_test(prior, dyn_e, meas_e, ye, m = 1L, method = "zzz")),

  ## -- imputation (MAR core): gmm_impute / gmm_complete / proxy_pool / proxy_fmi / as_mids
  ok_case("gmm_impute", "gmm_impute MAR",
          function() gmm_impute(Xi, m = 5L, seed = 1L),
          check = function(v) inherits(v, class(imp)[1])),
  ok_case("proxy_pool", "proxy_pool column mean",
          function() proxy_pool(imp, "x2"),
          check = function(v) is.finite(v$estimate)),
  ok_case("proxy_fmi", "proxy_fmi finite",
          function() proxy_fmi(imp, "x2")),
  ok_case("as_mids", "as_mids -> mids",
          function() as_mids(imp), check = function(v) inherits(v, "mids")),
  ok_case("gmm_complete", "gmm_complete first completion",
          function() gmm_complete(imp, which = 1L)),
  err_case("gmm_complete", "gmm_complete which out of range",
           function() gmm_complete(imp, which = 99L)),
  err_case("gmm_impute", "gmm_impute char", function() gmm_impute("x")),
  err_case("gmm_impute", "gmm_impute m=0", function() gmm_impute(Xi, m = 0L)),
  err_case("gmm_impute", "gmm_impute bad mechanism",
           function() gmm_impute(Xi, mechanism = "bad")),
  err_case("proxy_pool", "proxy_pool non-imputation", function() proxy_pool(42, "x2")),
  err_case("proxy_pool", "proxy_pool missing column",
           function() proxy_pool(imp, "nope")),
  err_case("proxy_fmi", "proxy_fmi non-imputation", function() proxy_fmi(42, "x2")),
  err_case("as_mids", "as_mids non-imputation", function() as_mids(42)),

  ## -- missing-not-at-random: mnar / censored / mar mechanisms ----------------
  ok_case("mnar", "mnar mechanism builds",
          function() mnar("x2", beta = 0.5)),
  ok_case("censored", "censored mechanism builds",
          function() censored("x2", lower = -1, upper = 1)),
  ok_case("mar", "mar mechanism builds", function() mar()),
  ok_case("gmm_impute", "gmm_impute MNAR gated",
          function() gmm_impute(Xi, m = 5L, mechanism = mnar("x2", beta = 0.5), seed = 1L),
          check = function(v) inherits(v, class(imp)[1])),
  ok_case("gmm_impute", "gmm_impute censored gated",
          function() gmm_impute(Xi, m = 5L, mechanism = censored("x2", lower = -3), seed = 1L),
          check = function(v) inherits(v, class(imp)[1])),
  err_case("mnar", "mnar beta char", function() mnar("x2", beta = "x")),
  err_case("censored", "censored lower>upper",
           function() censored("x2", lower = 1, upper = -1)),
  ok_case("proxy_mnar_sensitivity", "proxy_mnar_sensitivity over a beta grid",
          function() proxy_mnar_sensitivity(mnar_dat, "y", beta_grid = c(0, 0.5),
                                            m = 5L, seed = 1L)),
  err_case("proxy_mnar_sensitivity", "proxy_mnar_sensitivity non-data",
           function() proxy_mnar_sensitivity(42, "x2")),

  ## -- decision / causal: proxy_cate / proxy_decide / proxy_policy_value ------
  ##    plus the counterfactual surface (negative-side: refusal contract) ------
  ok_case("proxy_cate", "proxy_cate tau vector",
          function() proxy_cate(um, du, se = FALSE),
          check = function(v) length(v$tau) == nrow(du)),
  ok_case("proxy_decide", "proxy_decide recommends",
          function() proxy_decide(um, du, value = 1)),
  ok_case("proxy_policy_value", "proxy_policy_value of a treat-all policy",
          function() proxy_policy_value(um, du, policy = rep(1, nrow(du)), value = 1)),
  ok_case("gmm_intervene", "gmm_intervene do+given on a 3-D joint",
          function() gmm_intervene(gK, do = c(NA, 1, NA), given = c(NA, NA, 0.3))),
  ok_case("gmm_counterfactual", "gmm_counterfactual law on a 3-D joint",
          function() gmm_counterfactual(gK, evidence = c(1.2, 0, 0.5),
                                        do = c(NA, 1, NA), query = 1L)),
  err_case("proxy_cate", "proxy_cate non-model", function() proxy_cate(42, du)),
  err_case("proxy_cate", "proxy_cate char newdata", function() proxy_cate(um, "x")),
  err_case("proxy_decide", "proxy_decide non-model",
           function() proxy_decide(42, du, value = 1)),
  err_case("proxy_policy_value", "proxy_policy_value non-model",
           function() proxy_policy_value(42, du, policy = rep(1, nrow(du)), value = 1)),
  err_case("gmm_intervene", "gmm_intervene non-gmm",
           function() gmm_intervene(42, do = c(NA, 1))),
  err_case("gmm_counterfactual", "gmm_counterfactual non-gmm",
           function() gmm_counterfactual(42, evidence = 1, do = 1, query = 1)),
  err_case("gmm_cf_variance", "gmm_cf_variance non-law", function() gmm_cf_variance(42)),
  err_case("gmm_cf_tail_prob", "gmm_cf_tail_prob non-law",
           function() gmm_cf_tail_prob(42, threshold = 1)),

  ## -- abstention / robustness: ess_summary / bic_aic / maxent_target ---------
  ok_case("ess_summary", "ess_summary on kld fit",
          function() ess_summary(fit_k)),
  ok_case("bic_aic", "bic_aic finite",
          function() bic_aic(fit_s),
          check = function(v) any(is.finite(unlist(v)))),
  ok_case("maxent_target", "maxent_target Gaussian moments",
          function() maxent_target(moments = list(mean = 0, cov = matrix(1)))),
  err_case("ess_summary", "ess_summary non-fit", function() ess_summary(42)),
  err_case("bic_aic", "bic_aic non-fit", function() bic_aic("x")),
  err_case("maxent_target", "maxent_target char moments",
           function() maxent_target(moments = "x")),

  ## ===== extended surface coverage (v0.11.2 refresh) ========================

  ## -- built-in targets -------------------------------------------------------
  ok_case("banana_target", "banana_target builds", function() banana_target()),
  err_case("banana_target", "banana_target negative sample count",
           function() banana_target(with_samples = TRUE, n = -5L)),
  ok_case("donut_target", "donut_target builds", function() donut_target()),
  err_case("donut_target", "donut_target negative sample count",
           function() donut_target(with_samples = TRUE, n = -5L)),
  ok_case("epanechnikov_target", "epanechnikov_target builds",
          function() epanechnikov_target()),
  err_case("epanechnikov_target", "epanechnikov_target negative half_width",
           function() epanechnikov_target(half_width = -1)),
  ok_case("mixture_target", "mixture_target builds", function() mixture_target()),
  err_case("mixture_target", "mixture_target negative sample count",
           function() mixture_target(with_samples = TRUE, n = -5L)),
  ok_case("gmm_target", "gmm_target evaluate-only",
          function() gmm_target(n_dim = 1L,
            log_density = function(x) dnorm(if (is.matrix(x)) x[, 1] else x, log = TRUE),
            normalised = TRUE)),
  err_case("gmm_target", "gmm_target log_density not a function",
           function() gmm_target(n_dim = 1L, log_density = 42)),
  ok_case("gmm_target_from_samples", "gmm_target_from_samples builds",
          function() gmm_target_from_samples(X2)),
  err_case("gmm_target_from_samples", "gmm_target_from_samples char",
           function() gmm_target_from_samples("x")),
  ok_case("gmm_target_from_posterior", "gmm_target_from_posterior from log-posterior",
          function() gmm_target_from_posterior(
            function(theta) -0.5 * rowSums(theta^2), parameter_names = c("a", "b"))),
  err_case("gmm_target_from_posterior", "gmm_target_from_posterior non-producer",
           function() gmm_target_from_posterior(42)),

  ## -- importance proposals ---------------------------------------------------
  ok_case("is_mvn", "is_mvn builds 2-D",
          function() is_mvn(n_dim = 2L, mean = c(0, 0), cov = 4 * diag(2))),
  err_case("is_mvn", "is_mvn char n_dim", function() is_mvn(n_dim = "x")),
  ok_case("is_mvt", "is_mvt builds 2-D", function() is_mvt(n_dim = 2L, df = 5)),
  err_case("is_mvt", "is_mvt df<=0", function() is_mvt(n_dim = 2L, df = -1)),
  ok_case("is_uniform", "is_uniform builds 2-D",
          function() is_uniform(n_dim = 2L, lower = -5, upper = 5)),
  err_case("is_uniform", "is_uniform lower>upper",
           function() is_uniform(n_dim = 2L, lower = 5, upper = -5)),
  ok_case("is_proposal", "is_proposal wraps a sampler",
          function() is_proposal(n_dim = 1L,
            sample = function(n) matrix(rnorm(n), ncol = 1L),
            log_density = function(x) dnorm(x, log = TRUE))),
  err_case("is_proposal", "is_proposal sample not a function",
           function() is_proposal(n_dim = 1L, sample = 42)),

  ## -- initialisers -----------------------------------------------------------
  ok_case("init_kmeans", "init_kmeans seeds N=3", function() init_kmeans(X2, N = 3L)),
  err_case("init_kmeans", "init_kmeans char samples", function() init_kmeans("x", N = 3L)),
  ok_case("init_moment_seed", "init_moment_seed seeds N=3",
          function() init_moment_seed(X2, N = 3L)),
  err_case("init_moment_seed", "init_moment_seed char samples",
           function() init_moment_seed("x", N = 3L)),
  ok_case("init_random", "init_random draws N=3",
          function() init_random(N = 3L, p = 2L, seed = 1L)),
  err_case("init_random", "init_random N=-1", function() init_random(N = -1L, p = 2L)),
  ok_case("init_warm_start", "init_warm_start from gmm", function() init_warm_start(g2)),
  err_case("init_warm_start", "init_warm_start non-gmm", function() init_warm_start(42)),

  ## -- gmm accessors ----------------------------------------------------------
  ok_case("gmm_dim", "gmm_dim of 2-D gmm",
          function() gmm_dim(g2), check = function(v) v == 2L),
  err_case("gmm_dim", "gmm_dim non-gmm", function() gmm_dim(42)),
  ok_case("gmm_n_components", "gmm_n_components of K=2",
          function() gmm_n_components(g2), check = function(v) v == 2L),
  err_case("gmm_n_components", "gmm_n_components non-gmm", function() gmm_n_components(42)),
  ok_case("gmm_canonicalise", "gmm_canonicalise orders K=4",
          function() gmm_canonicalise(g4),
          check = function(v) length(v@weights) == 4L),
  err_case("gmm_canonicalise", "gmm_canonicalise non-gmm", function() gmm_canonicalise(42)),
  ok_case("gmm_kld", "gmm_kld finite",
          function() gmm_kld(g2, gmm(weights = 1, means = list(c(0, 0)),
                                     covariances = list(diag(2))), n_mc = 1000L),
          check = function(v) is.finite(v$mc) && is.finite(v$variational)),
  err_case("gmm_kld", "gmm_kld non-gmm q", function() gmm_kld(g2, 42)),

  ## -- diagnostics traces -----------------------------------------------------
  ok_case("ess_trace", "ess_trace on kld fit", function() ess_trace(fit_k)),
  err_case("ess_trace", "ess_trace non-fit", function() ess_trace(42)),
  ok_case("kld_trace", "kld_trace on kld fit", function() kld_trace(fit_k)),
  err_case("kld_trace", "kld_trace non-fit", function() kld_trace(42)),
  ok_case("hellinger_mc", "hellinger_mc finite",
          function() hellinger_mc(fit_k, n_mc = 1000L, seed = 1L),
          check = function(v) is.finite(unlist(v)[1])),
  err_case("hellinger_mc", "hellinger_mc non-fit", function() hellinger_mc(42)),

  ## -- counterfactual-law accessors -------------------------------------------
  ok_case("gmm_cf_mean", "gmm_cf_mean finite",
          function() gmm_cf_mean(cf), check = function(v) all(is.finite(unlist(v)))),
  err_case("gmm_cf_mean", "gmm_cf_mean non-law", function() gmm_cf_mean(42)),
  ## the individual counterfactual variance / tail probability are not identified,
  ## so the contract is to refuse the query even on a valid law (a principled abstention)
  err_case("gmm_cf_variance", "gmm_cf_variance refuses an unidentified query",
           function() gmm_cf_variance(cf)),
  err_case("gmm_cf_tail_prob", "gmm_cf_tail_prob refuses an unidentified query",
           function() gmm_cf_tail_prob(cf, threshold = 1)),

  ## -- uplift estimation + the uplift_model surface ---------------------------
  ok_case("fit_uplift", "fit_uplift moment regime",
          function() fit_uplift(du, outcome = "y", treatment = "t",
                                covariates = "x1", N = 1L, regime = "moment")),
  err_case("fit_uplift", "fit_uplift missing outcome column",
           function() fit_uplift(du, outcome = "nope", treatment = "t", covariates = "x1")),
  ok_case("proxy_uplift", "proxy_uplift tau on newdata",
          function() proxy_uplift(um, du)),
  err_case("proxy_uplift", "proxy_uplift non-model", function() proxy_uplift(42, du)),
  ok_case("proxy_predict", "proxy_predict potential outcomes",
          function() proxy_predict(um, du, t = 1)),
  err_case("proxy_predict", "proxy_predict non-model", function() proxy_predict(42, du, t = 1)),
  ok_case("proxy_overlap", "proxy_overlap diagnostics", function() proxy_overlap(um, du)),
  err_case("proxy_overlap", "proxy_overlap non-model", function() proxy_overlap(42, du)),
  ok_case("proxy_confounding_gap", "proxy_confounding_gap on model",
          function() proxy_confounding_gap(um, du)),
  err_case("proxy_confounding_gap", "proxy_confounding_gap non-model",
           function() proxy_confounding_gap(42, du)),
  ok_case("proxy_identification_report", "proxy_identification_report builds",
          function() proxy_identification_report(um, du)),
  err_case("proxy_identification_report", "proxy_identification_report non-model",
           function() proxy_identification_report(42, du)),
  ok_case("proxy_regime_segments", "proxy_regime_segments builds",
          function() proxy_regime_segments(um)),
  err_case("proxy_regime_segments", "proxy_regime_segments non-model",
           function() proxy_regime_segments(42)),
  ok_case("proxy_retrospective_uplift", "proxy_retrospective_uplift builds",
          function() proxy_retrospective_uplift(um, observed = du)),
  err_case("proxy_retrospective_uplift", "proxy_retrospective_uplift non-model",
           function() proxy_retrospective_uplift(42, observed = du)),
  ok_case("uplift_identification", "uplift_identification report object",
          function() uplift_identification(estimand = "ATE", assume = "ignorability",
                                           n_units = 100L)),
  err_case("uplift_identification", "uplift_identification char n_units",
           function() uplift_identification(n_units = "x")),
  ok_case("uplift_model", "uplift_model constructs",
          function() uplift_model(roles = list(outcome = "y"), outcome_type = "continuous")),
  err_case("uplift_model", "uplift_model char outcome_type",
           function() uplift_model(outcome_type = 42)),

  ## -- forward / FB producers -------------------------------------------------
  ok_case("fb_log_posterior_spec", "fb_log_posterior_spec builds",
          function() fb_log_posterior_spec(function(theta) -0.5 * rowSums(theta^2),
                                           parameter_names = c("a", "b"))),
  err_case("fb_log_posterior_spec", "fb_log_posterior_spec non-function producer",
           function() fb_log_posterior_spec(42)),
  ok_case("fb_producer_available", "fb_producer_available returns logical",
          function() fb_producer_available(),
          check = function(v) is.logical(v) && length(v) == 1L),
  ok_case("mock_fb_posterior", "mock_fb_posterior gaussian",
          function() mock_fb_posterior(shape = "gaussian", n_dim = 2L)),
  err_case("mock_fb_posterior", "mock_fb_posterior bad shape",
           function() mock_fb_posterior(shape = "zzz")),

  ## -- S7 result constructors -------------------------------------------------
  ok_case("gmm_fit", "gmm_fit constructs a fit object",
          function() gmm_fit(weights = c(0.5, 0.5), means = list(0, 1),
                             covariances = list(matrix(1), matrix(1)))),
  err_case("gmm_fit", "gmm_fit char weights", function() gmm_fit(weights = "x")),
  ok_case("gmm_imputation", "gmm_imputation constructs",
          function() gmm_imputation(m = 5L, mechanism = "mar")),
  err_case("gmm_imputation", "gmm_imputation char m", function() gmm_imputation(m = "x")),
  ok_case("gmm_counterfactual_law", "gmm_counterfactual_law constructs",
          function() gmm_counterfactual_law(atoms = c(0, 1), weights = c(0.5, 0.5),
                                            mean = 0.5)),
  err_case("gmm_counterfactual_law", "gmm_counterfactual_law char weights",
           function() gmm_counterfactual_law(weights = "x")),

  ## -- optimisation / scenario / collider helpers -----------------------------
  ok_case("multi_start_best_of", "multi_start_best_of picks the best start",
          function() multi_start_best_of(
            fit_fn = function(init) fit_em_samples(st2, N = 2L, seed = 1L),
            inits = list(init_random(2L, 2L, seed = 1L), init_moment_seed(X2, N = 2L)),
            score_fn = function(fit) -bic_aic(fit)$bic)),
  err_case("multi_start_best_of", "multi_start_best_of char fit_fn",
           function() multi_start_best_of("x", list(), identity)),
  err_case("to_apsim_scenarios", "to_apsim_scenarios non-fit",
           function() to_apsim_scenarios(42)),
  err_case("fit_kld_em_collider", "fit_kld_em_collider non-target",
           function() fit_kld_em_collider(42, dag = matrix(0, 2, 2)))
)

## ---- deliberately left without TWO-sided coverage (honest documentation) ----
## A function counts as two-sided covered only when it carries BOTH an "ok" and
## an "error" case. The following are exported but intentionally not two-sided,
## and are NOT findings -- there is no rejectable / acceptable counterpart to add
## without fabricating one:
##
##   * Tier-2 stubs in R/stubs.R -- from_aggregate_likelihood, from_simulator,
##     to_apsim_scenarios, fit_kld_em_collider. Every call signals
##     `proxymix_not_yet_implemented`, so there is no valid input that returns
##     ok; only the error side is meaningful (see their err_case rows above).
##   * gmm_cf_variance / gmm_cf_tail_prob -- by-design refusal accessors. The
##     individual counterfactual variance and tail probability are not
##     identified, so the contract is to refuse even a valid law (a principled
##     abstention). Both sides of their surface are therefore "error"; there is
##     no legitimate ok case (see the refusal err_case rows above).
##   * mar / fb_producer_available -- zero-required-argument constructors / a
##     pure environment probe. They take no rejectable input, so there is no
##     malformed call to add; only the ok side exists.

## ---- run the sweep ---------------------------------------------------------
fn_of <- vapply(CASES, function(cs) cs$fn, character(1))
verdicts <- pv_grid(
  CASES,
  fn     = function(cs) cs$call(),
  expect = function(cs) cs$expect,
  label  = function(cs, i) cs$label,
  check  = function(v, cs) if (is.null(cs$check)) TRUE else isTRUE(cs$check(v)),
  cores  = 8L, base_seed = 20260620L)
verdicts$fn <- fn_of[verdicts$idx]
verdicts <- verdicts[, c("idx", "fn", "label", "expect", "status", "seconds",
                         "n_warn", "detail")]

## ---- metamorphic invariants ------------------------------------------------
META <- list(
  pv_metamorphic(unlist(gmm_affine(g2, diag(2))@means), unlist(g2@means),
                 what = "gmm_affine(g, I, 0) is the identity on means"),
  pv_metamorphic(unlist(gmm_affine(g2, diag(2), ridge_eps = 0)@covariances),
                 unlist(g2@covariances),
                 what = "gmm_affine(g, I, 0) is the identity on covariances"),
  pv_metamorphic(unlist(gmm_marginalise(g2, keep = c(1L, 2L))@means), unlist(g2@means),
                 what = "gmm_marginalise(g, all) is the identity"),
  pv_metamorphic(
    unlist(gmm_affine(gmm_affine(g2, matrix(c(1, 0.3, -0.2, 1), 2)),
                      matrix(c(0.5, 0, 0, 2), 2))@means),
    unlist(gmm_affine(g2, matrix(c(0.5, 0, 0, 2), 2) %*%
                          matrix(c(1, 0.3, -0.2, 1), 2))@means),
    what = "affine composition A2(A1 x) equals (A2 A1) x"),
  pv_metamorphic(
    dgmm(matrix(c(0.3, -0.4), 1, 2), g2, log = TRUE),
    log(dgmm(matrix(c(0.3, -0.4), 1, 2), g2, log = FALSE)),
    what = "dgmm log equals log of dgmm density"),
  pv_metamorphic(gmean(gmm_reduce(g4, k_max = 2L)), gmean(g4),
                 what = "gmm_reduce preserves the global mean"),
  local({
    a <- withr::with_seed(7, rgmm(40L, g2)); b <- withr::with_seed(7, rgmm(40L, g2))
    pv_metamorphic(a, b, what = "rgmm is reproducible under a fixed seed")
  })
)
meta_df <- do.call(rbind, lapply(META, function(m) data.frame(
  what = m$what, pass = m$pass, gap = signif(m$gap, 3),
  detail = m$detail, stringsAsFactors = FALSE)))

## ---- assemble + persist ----------------------------------------------------
## Assembled by the engine (engine/conformance.R) so the object carries the
## coverage fields the dossier + pv_api_coverage_gate read; verdicts and the
## metamorphic results are unchanged.
out <- pv_assemble_conformance(verdicts, meta_df,
                               surface = getNamespaceExports("proxymix"),
                               master_seed = 20260620L,
                               provenance = pv_provenance("../proxymix"))
summary_df  <- out$summary
findings_df <- out$findings
saveRDS(out, file.path(RESDIR, "contract-conformance.rds"))

cat(sprintf("contract-conformance: %d cells over %d functions ; %d findings ; metamorphic %d/%d pass\n",
            out$n_cells, out$n_functions, out$n_findings,
            out$n_meta - out$n_meta_fail, out$n_meta))
cat("status tally:\n"); print(summary_df[summary_df$n > 0, ], row.names = FALSE)
if (nrow(findings_df)) {
  cat("\nFINDINGS (status, expect, label):\n")
  print(findings_df[, c("fn", "label", "expect", "status", "detail")], row.names = FALSE)
}
if (out$n_meta_fail) { cat("\nMETAMORPHIC FAILURES:\n"); print(meta_df[!meta_df$pass, ], row.names = FALSE) }
