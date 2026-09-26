# =============================================================================
# Summary of the stored filtering benchmark for the short operator-calculus vignette
#
# Author      : Max Moldovan
# Version date: 26 Sep 2026
# Usage       : Rscript data-raw/vignette_results/operator_calculus.R
#               (run from the package root)
#
# Data        : vignettes/articles/extended/operator_calculus_sim.rds, written by
#               data-raw/articles/operator_calculus_sim.R; datasets::Nile
# Requirements: proxymix, dlm, KFAS (section 3 only)
#
# Outputs     : vignettes/results/operator_calculus.rds, read by vignettes/operator_calculus.Rmd
#
# Structure   : 1. Read the stored simulation
#               2. Paired differences between filters
#               3. The Nile at one and two noise components
#               4. Write the summary
# =============================================================================


rm(list = ls())

ptm1 <- proc.time()	# total execution time


# Configuration: paths
########################################################

file_name_sim <- paste0(getwd(), "/vignettes/articles/extended/operator_calculus_sim.rds")
dir_out <- paste0(getwd(), "/vignettes/results/")
########################################################


# =============================================================================
# 1. Read the stored simulation
# =============================================================================
print("==== 1. Read the stored simulation ====")
ptm_local <- proc.time()

sim <- readRDS(file_name_sim)

# the version stamp is the proxymix version the simulation ran under
stopifnot(all(c("sim_tab", "res", "code", "n_t", "n_rep", "q_state", "c0",
                "r_sd", "r_w", "r_var", "t_scale", "versions", "run_date",
                "elapsed_secs") %in% names(sim)))
stopifnot("proxymix" %in% names(sim$versions))
print(sim$versions)

# settings the simulation code fixes as literals, checked against that code
stopifnot(any(sim$code == "n_particles <- 10000L"))
stopifnot(any(sim$code == "t_df <- 4"))
stopifnot(any(grepl("y = y, k_max = 4L)", sim$code, fixed = TRUE)))
n_particles <- 10000L
t_df <- 4
k_max_sim <- 4L


time_section <- proc.time() - ptm_local
print("Section: 1. Read the stored simulation")
print(time_section)


# =============================================================================
# 2. Paired differences between filters
# =============================================================================
print("==== 2. Paired differences between filters ====")
ptm_local <- proc.time()

# every filter ran on the same series, so RMSE differences are paired by series:
# mean of m1 minus m2, its standard error, and the share of series on which
# m2 has the lower RMSE
paired <- function(m1, m2) {
	a <- sim$res[sim$res$method == m1, ]
	b <- sim$res[sim$res$method == m2, ]
	stopifnot(nrow(a) == sim$n_rep, nrow(b) == sim$n_rep)
	d <- a$rmse[order(a$series)] - b$rmse[order(b$series)]
	c(diff = mean(d), se = sd(d) / sqrt(length(d)), share = mean(d > 0))
}

pm <- "proxymix, two-component noise"
paired_list <- list(dlm = paired("dlm, Gaussian noise", pm),
                    student_t = paired("nimble, Student-t noise", pm),
                    pf = paired(pm, "particle filter, two-component noise"),
                    nimble_mix = paired(pm, "nimble, two-component noise"))
print(paired_list)


time_section <- proc.time() - ptm_local
print("Section: 2. Paired differences between filters")
print(time_section)


# =============================================================================
# 3. The Nile at one and two noise components
# =============================================================================
print("==== 3. The Nile at one and two noise components ====")
ptm_local <- proc.time()

require(proxymix)
require(dlm)
require(KFAS)

# the code of the extended article's chunks nile-fit, nile-three and nile-mixture
y_nile <- as.numeric(datasets::Nile)
n_nile <- length(y_nile)
build_level <- function(p) dlmModPoly(1L, dV = exp(p[1L]), dW = exp(p[2L]))
fit_nile <- dlmMLE(y_nile, parm = c(log(15000), log(1500)),
                   build = build_level)
h_nile <- exp(fit_nile$par[1L])   # observation noise variance
q_nile <- exp(fit_nile$par[2L])   # variance of the level increments
m0_nile <- 1000
c0_nile <- 1e5

mod_dlm <- dlmModPoly(1L, dV = h_nile, dW = q_nile, m0 = m0_nile,
                      C0 = c0_nile)
f_dlm <- dlmFilter(y_nile, mod_dlm)
mean_dlm <- as.numeric(f_dlm$m[-1L])
var_dlm <- unlist(dlmSvd2var(f_dlm$U.C, f_dlm$D.C))[-1L]
ll_dlm <- -dlmLL(y_nile, mod_dlm) - n_nile / 2 * log(2 * pi)

mod_kfas <- SSModel(y_nile ~ SSMtrend(1L, Q = q_nile, a1 = m0_nile,
                                      P1 = c0_nile + q_nile), H = h_nile)
f_kfas <- KFS(mod_kfas, filtering = "state", smoothing = "none")
mean_kfas <- as.numeric(f_kfas$att)
var_kfas <- as.numeric(f_kfas$Ptt)
ll_kfas <- as.numeric(logLik(mod_kfas))

prior_nile <- gmm(weights = 1, means = list(m0_nile),
                  covariances = list(matrix(c0_nile)))
f_pm <- gmm_filter(prior_nile,
                   dynamics = list(A = matrix(1), Q = matrix(q_nile)),
                   measurement = list(C = matrix(1), R = matrix(h_nile)),
                   y = y_nile, ridge_eps = 0)
mean_pm <- f_pm$mean[, 1L]
var_pm <- f_pm$summary$sd_1^2
ll_pm <- sum(f_pm$summary$log_evidence)

r_mix_nile <- gmm(weights = c(0.9, 0.1), means = list(0, 0),
                  covariances = list(matrix(h_nile / 1.9), matrix(10 * h_nile / 1.9)))
f_mix <- gmm_filter(prior_nile,
                    dynamics = list(A = matrix(1), Q = matrix(q_nile)),
                    measurement = list(C = matrix(1), R = r_mix_nile),
                    y = y_nile, k_max = 4L)
ll_mix <- sum(f_mix$summary$log_evidence)

# largest gaps between proxymix at one component and the two Kalman filters
gap <- function(a, b) max(abs(a - b))
nile_list <- list(n = n_nile,
                  years = range(as.numeric(time(datasets::Nile))),
                  h = h_nile,
                  q = q_nile,
                  mean_gap = max(gap(mean_pm, mean_dlm), gap(mean_pm, mean_kfas)),
                  var_gap = max(gap(var_pm, var_dlm), gap(var_pm, var_kfas)),
                  ll_gauss = c(proxymix = ll_pm, dlm = ll_dlm, KFAS = ll_kfas),
                  ll_mix = ll_mix,
                  k_max = 4L)
print(nile_list)

nile_versions <- c(proxymix = as.character(packageVersion("proxymix")),
                   dlm = as.character(packageVersion("dlm")),
                   KFAS = as.character(packageVersion("KFAS")))


time_section <- proc.time() - ptm_local
print("Section: 3. The Nile at one and two noise components")
print(time_section)


# =============================================================================
# 4. Write the summary
# =============================================================================
print("==== 4. Write the summary ====")
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
                paired = paired_list,
                n_t = sim$n_t,
                n_rep = sim$n_rep,
                q_state = sim$q_state,
                c0 = sim$c0,
                r_sd = sim$r_sd,
                r_w = sim$r_w,
                r_var = sim$r_var,
                t_scale = sim$t_scale,
                t_df = t_df,
                n_particles = n_particles,
                k_max = k_max_sim,
                nile = nile_list,
                nile_versions = desc_version(nile_versions),
                proxymix_version = sim$versions[["proxymix"]],
                versions = desc_version(sim$versions),
                run_date = sim$run_date,
                elapsed_secs = sim$elapsed_secs)

dir.create(dir_out, showWarnings = FALSE)
file_name_out <- paste0(dir_out, "operator_calculus.rds")
saveRDS(res_out, file_name_out)

print(paste0("written: ", file_name_out, " (", file.size(file_name_out), " bytes)"))
print(paste0("version stamp: ", res_out$proxymix_version))


time_section <- proc.time() - ptm_local
print("Section: 4. Write the summary")
print(time_section)


# =============================================================================
# Completion
# =============================================================================
exec_time1 <- proc.time() - ptm1
print(exec_time1)
