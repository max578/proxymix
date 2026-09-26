# =============================================================================
# Summary of the stored posterior benchmark for the short posterior-proxy vignette
#
# Author      : Max Moldovan
# Version date: 26 Sep 2026
# Usage       : Rscript data-raw/vignette_results/posterior_proxy.R
#               (run from the package root)
#
# Data        : vignettes/articles/extended/posterior_proxy_sim.rds, written by
#               data-raw/articles/posterior_proxy_sim.R
# Requirements: none beyond base R; the simulation records its own run times
#
# Outputs     : vignettes/results/posterior_proxy.rds, read by vignettes/posterior_proxy.Rmd
#
# Structure   : 1. Read the stored simulation
#               2. Summarise the comparison
#               3. Write the summary
# =============================================================================


rm(list = ls())

ptm1 <- proc.time()	# total execution time


# Configuration: paths
########################################################

file_name_sim <- paste0(getwd(), "/vignettes/articles/extended/posterior_proxy_sim.rds")
dir_out <- paste0(getwd(), "/vignettes/results/")
########################################################


# =============================================================================
# 1. Read the stored simulation
# =============================================================================
print("==== 1. Read the stored simulation ====")
ptm_local <- proc.time()

sim <- readRDS(file_name_sim)

# the version stamp is the proxymix version the simulation ran under
stopifnot(all(c("sim_tab", "coverage", "secs_ensemble", "ref_ess_min", "n",
                "n_rep", "p_vec", "prior_sd", "beta_pop", "versions",
                "run_date", "elapsed_secs", "cpu_secs") %in% names(sim)))
stopifnot("proxymix" %in% names(sim$versions))
print(sim$versions)

# the simulation times each method inside every replicate (mean seconds)
stopifnot(!anyNA(sim$sim_tab$secs))
per_ds <- attr(sim$sim_tab, "per_dataset")
stopifnot(!is.null(per_ds))


time_section <- proc.time() - ptm_local
print("Section: 1. Read the stored simulation")
print(time_section)


# =============================================================================
# 2. Summarise the comparison
# =============================================================================
print("==== 2. Summarise the comparison ====")
ptm_local <- proc.time()

method_vec <- c("NUTS, 4000 draws", "proxymix", "Laplace", "INLA", "DEzs")

# one row per p and method; the interval-end error is the mean of the two ends
tab_sim <- sim$sim_tab
tab_sim <- tab_sim[order(tab_sim$p, match(tab_sim$method, method_vec)), ]
tab_sim$endpoints <- (tab_sim$err_lower + tab_sim$err_upper) / 2
tab_sim <- tab_sim[, c("p", "method", "err_mean", "endpoints", "err_p_pos",
                       "err_log_z", "secs")]
rownames(tab_sim) <- NULL
print(tab_sim)

per_ds$endpoints <- (per_ds$err_lower + per_ds$err_upper) / 2

# two methods on the same datasets: mean difference in dataset-level error
# divided by its paired standard error
paired_z <- function(p, m1, m2, what) {
	a <- per_ds[per_ds$p == p & per_ds$method == m1, ]
	b <- per_ds[per_ds$p == p & per_ds$method == m2, ]
	d_vec <- a[[what]][order(a$rep)] - b[[what]][order(b$rep)]
	mean(d_vec) / (sd(d_vec) / sqrt(length(d_vec)))
}

# the closest method on one error, the next closest, and the paired z
lead_tab <- data.frame()
for (what in c("err_mean", "endpoints", "err_p_pos", "err_log_z")) {
	for (p in sim$p_vec) {

		s1 <- tab_sim$p == p & !is.na(tab_sim[[what]])
		m_vec <- tab_sim$method[s1][order(tab_sim[[what]][s1])[1:2]]
		lead_tab <- rbind(lead_tab,
		                  data.frame(what = what, p = p, first = m_vec[1L],
		                             second = m_vec[2L],
		                             z = paired_z(p, m_vec[2L], m_vec[1L], what)))

	} # end of for (p in sim$p_vec)
} # end of for (what in ...)
print(lead_tab)

# proxymix against the 4000-draw NUTS run on the posterior means
z_px_nuts_mean <- vapply(sim$p_vec, paired_z, numeric(1L), m1 = "proxymix",
                         m2 = "NUTS, 4000 draws", what = "err_mean")
names(z_px_nuts_mean) <- sim$p_vec
print(z_px_nuts_mean)

# average error of the reference mean itself, in posterior standard deviations
ref_mean_noise <- sqrt(2 / pi) / sqrt(sim$ref_ess_min)


time_section <- proc.time() - ptm_local
print("Section: 2. Summarise the comparison")
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

res_out <- list(sim_tab = tab_sim,
                lead_tab = lead_tab,
                z_px_nuts_mean = z_px_nuts_mean,
                coverage = sim$coverage,
                secs_ensemble = sim$secs_ensemble,
                ref_ess_min = sim$ref_ess_min,
                ref_mean_noise = ref_mean_noise,
                n = sim$n,
                n_rep = sim$n_rep,
                p_vec = sim$p_vec,
                prior_sd = sim$prior_sd,
                beta_pop = sim$beta_pop,
                proxymix_version = sim$versions[["proxymix"]],
                versions = desc_version(sim$versions),
                run_date = sim$run_date,
                elapsed_secs = sim$elapsed_secs,
                cpu_secs = sim$cpu_secs)

dir.create(dir_out, showWarnings = FALSE)
file_name_out <- paste0(dir_out, "posterior_proxy.rds")
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
