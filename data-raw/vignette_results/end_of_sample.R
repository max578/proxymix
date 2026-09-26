# =============================================================================
# Summary of the stored end-of-sample benchmark for the short end-of-sample vignette
#
# Author      : Max Moldovan
# Version date: 26 Sep 2026
# Usage       : Rscript data-raw/vignette_results/end_of_sample.R
#               (run from the package root)
#
# Data        : vignettes/articles/extended/end_of_sample_sim.rds, written by
#               data-raw/articles/end_of_sample_sim.R
# Requirements: proxymix, strucchange, changepoint, KFAS, dlm (section 2 only)
#
# Outputs     : vignettes/results/end_of_sample.rds, read by vignettes/end_of_sample.Rmd
#
# Structure   : 1. Read the stored simulation
#               2. Time the tests of each package on one series
#               3. Write the summary
# =============================================================================


rm(list = ls())

ptm1 <- proc.time()	# total execution time


# Configuration: paths
########################################################

file_name_sim <- paste0(getwd(), "/vignettes/articles/extended/end_of_sample_sim.rds")
dir_out <- paste0(getwd(), "/vignettes/results/")
########################################################


# =============================================================================
# 1. Read the stored simulation
# =============================================================================
print("==== 1. Read the stored simulation ====")
ptm_local <- proc.time()

sim <- readRDS(file_name_sim)

# the version stamp is the proxymix version the simulation ran under
stopifnot(all(c("sim_tab", "n", "phi", "delta", "m_vec", "n_rep", "versions",
                "run_date", "elapsed_secs") %in% names(sim)))
stopifnot("proxymix" %in% names(sim$versions))
print(sim$versions)


time_section <- proc.time() - ptm_local
print("Section: 1. Read the stored simulation")
print(time_section)


# =============================================================================
# 2. Time the tests of each package on one series
# =============================================================================
print("==== 2. Time the tests of each package on one series ====")
ptm_local <- proc.time()

require(proxymix)
require(strucchange)
require(changepoint)
require(KFAS)
require(dlm)

# one series, as in the vignette's competitor chunk
set.seed(1L)
n <- 100L
m <- 3L
y <- as.numeric(arima.sim(list(ar = 0.6), n = n))
y[(n - m + 1L):n] <- y[(n - m + 1L):n] + 3
y_fit <- y[seq_len(n - m)]
last <- (n - m + 1L):n

# each expression is the vignette chunk's code for one package
expr_list <- list(
  proxymix = quote({
    ar <- arima(y_fit, order = c(1L, 0L, 0L))
    a1 <- coef(ar)[["ar1"]]
    mu <- coef(ar)[["intercept"]]
    s2 <- ar$sigma2
    prior <- gmm(weights = 1, means = list(mu),
                 covariances = list(matrix(s2 / (1 - a1^2))))
    dynamics <- list(A = matrix(a1), b = mu * (1 - a1), Q = matrix(s2))
    measurement <- list(C = matrix(1), R = matrix(0))
    p_chisq <- gmm_eos_test(prior, dynamics, measurement, y, m = m,
                            method = "chisq")$p_value
    p_andrews <- gmm_eos_test(prior, dynamics, measurement, y, m = m,
                              method = "andrews")$p_value
  }),
  strucchange = quote({
    reg <- data.frame(y = y[-1L], y_lag = y[-n])
    n_reg <- nrow(reg)
    p_supf <- if (m >= 3L) {
      fs <- Fstats(y ~ y_lag, data = reg, from = n_reg - m, to = n_reg - 3L)
      unname(sctest(fs)$p.value)
    } else NA_real_
    p_cusum <- unname(sctest(efp(y ~ y_lag, data = reg,
                                 type = "OLS-CUSUM"))$p.value)
  }),
  changepoint = quote({
    cp <- cpts(cpt.mean(y / sd(y_fit)))
    p_cpt <- if (any(cp >= n - m)) 0 else 1
  }),
  KFAS = quote({
    y_bar <- mean(y_fit)
    kfas_update <- function(pars, model) {
      part <- SSMarima(ar = 0.999 * tanh(pars[1L]), Q = exp(pars[2L]))
      model["T", "arima"] <- part$T
      model["R", "arima"] <- part$R
      model["Q", "arima"] <- part$Q
      model["P1", "arima"] <- part$P1
      model
    }
    kfas_fit <- fitSSM(
      SSModel(I(y_fit - y_bar) ~ -1 + SSMarima(ar = 0.5, Q = 1), H = 0),
      inits = c(0, 0), updatefn = kfas_update, method = "BFGS"
    )
    kfas_par <- kfas_fit$optim.out$par
    kfas_model <- SSModel(
      I(y - y_bar) ~ -1 + SSMarima(ar = 0.999 * tanh(kfas_par[1L]),
                                    Q = exp(kfas_par[2L])),
      H = 0
    )
    kfas_pred <- predict(kfas_model, interval = "prediction", level = 0.95,
                         filtered = TRUE)
    kfas_z <- (y[last] - y_bar - kfas_pred[last, "fit"]) /
      ((kfas_pred[last, "upr"] - kfas_pred[last, "lwr"]) / (2 * qnorm(0.975)))
    p_kfas <- min(1, m * 2 * pnorm(-max(abs(kfas_z))))
  }),
  dlm = quote({
    y_bar <- mean(y_fit)
    dlm_build <- function(pars) {
      dlmModARMA(ar = 0.999 * tanh(pars[1L]), sigma2 = exp(pars[2L]), dV = 0)
    }
    dlm_fit <- dlmMLE(y_fit - y_bar, parm = c(0, 0), build = dlm_build)
    dlm_filt <- dlmFilter(y - y_bar, dlm_build(dlm_fit$par))
    dlm_var <- unlist(dlmSvd2var(dlm_filt$U.R, dlm_filt$D.R))
    dlm_z <- (y[last] - y_bar - dlm_filt$f[last]) / sqrt(dlm_var[last])
    p_dlm <- min(1, m * 2 * pnorm(-max(abs(dlm_z))))
  })
)

# each run repeats a package's code n_call times; the time kept is per call
n_time <- 5L
n_call <- 10L
env_time <- new.env()
time_mat <- matrix(NA_real_, nrow = n_time, ncol = length(expr_list),
                   dimnames = list(NULL, names(expr_list)))
for (i1 in seq_len(n_time)) {
  for (k1 in names(expr_list)) {
    time_mat[i1, k1] <- system.time(
      for (j1 in seq_len(n_call)) eval(expr_list[[k1]], envir = env_time)
    )[["elapsed"]] / n_call
  } # end of for (k1 in names(expr_list))
} # end of for (i1 in seq_len(n_time))

# median seconds for one series at m = 3, per package
time_vec <- apply(time_mat, 2L, median)
print(time_vec)
print(c(env_time$p_chisq, env_time$p_andrews, env_time$p_supf, env_time$p_cusum,
        env_time$p_cpt, env_time$p_kfas, env_time$p_dlm))
time_versions <- vapply(names(expr_list), function(k) {
  as.character(packageVersion(k))
}, character(1L))


time_section <- proc.time() - ptm_local
print("Section: 2. Time the tests of each package on one series")
print(time_section)


# =============================================================================
# 3. Write the summary
# =============================================================================
print("==== 3. Write the summary ====")
ptm_local <- proc.time()

## package versions as each DESCRIPTION writes them (strucchange 1.6-0, not 1.6.0),
## taken from the installed package only when it is the version that was run
desc_version <- function(v) {
  v_out <- v
  pkg_vec <- intersect(setdiff(names(v), "R"), rownames(utils::installed.packages()))
  v_out[pkg_vec] <- vapply(pkg_vec, function(s1) {
    if (identical(as.character(packageVersion(s1)), v[[s1]])) {
      utils::packageDescription(s1)$Version
    } else {
      v[[s1]]
    }
  }, character(1L))
  v_out
}

res_out <- list(sim_tab = sim$sim_tab,
                n = sim$n,
                phi = sim$phi,
                delta = sim$delta,
                m_vec = sim$m_vec,
                n_rep = sim$n_rep,
                proxymix_version = sim$versions[["proxymix"]],
                versions = desc_version(sim$versions),
                run_date = sim$run_date,
                elapsed_secs = sim$elapsed_secs,
                time_secs = time_vec,
                time_m = m,
                time_n_call = n_call,
                time_versions = desc_version(time_versions),
                time_platform = R.version$platform)

dir.create(dir_out, showWarnings = FALSE)
file_name_out <- paste0(dir_out, "end_of_sample.rds")
saveRDS(res_out, file_name_out)

print(paste0("written: ", file_name_out, " (", file.size(file_name_out), " bytes)"))
print(paste0("version stamp: ", res_out$proxymix_version))


time_section <- proc.time() - ptm_local
print("Section: 3. Write the summary")
print(time_section)


# =============================================================================
# Completion
# =============================================================================
exec_time1 <- proc.time() - ptm1
print(exec_time1)
