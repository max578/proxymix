# =============================================================================
# Summary of the stored imputation benchmark for the short missing-data vignette
#
# Author      : Max Moldovan
# Version date: 26 Sep 2026
# Usage       : Rscript data-raw/vignette_results/missing_data.R
#               (run from the package root)
#
# Data        : vignettes/articles/extended/missing_data_sim.rds, written by
#               data-raw/articles/missing_data_sim.R
# Requirements: proxymix, mice, Amelia (section 2 only)
#
# Outputs     : vignettes/results/missing_data.rds, read by vignettes/missing_data.Rmd
#
# Structure   : 1. Read the stored simulation
#               2. Time one imputation per package
#               3. Write the summary
# =============================================================================


rm(list = ls())

ptm1 <- proc.time()	# total execution time


# Configuration: paths
########################################################

file_name_sim <- paste0(getwd(), "/vignettes/articles/extended/missing_data_sim.rds")
dir_out <- paste0(getwd(), "/vignettes/results/")
########################################################


# =============================================================================
# 1. Read the stored simulation
# =============================================================================
print("==== 1. Read the stored simulation ====")
ptm_local <- proc.time()

sim <- readRDS(file_name_sim)

# the version stamp is the proxymix version the simulation ran under
stopifnot(all(c("sim_tab", "n", "m", "n_rep", "population", "versions",
                "run_date", "elapsed_secs") %in% names(sim)))
stopifnot("proxymix" %in% names(sim$versions))
print(sim$versions)


time_section <- proc.time() - ptm_local
print("Section: 1. Read the stored simulation")
print(time_section)


# =============================================================================
# 2. Time one imputation per package
# =============================================================================
print("==== 2. Time one imputation per package ====")
ptm_local <- proc.time()

require(proxymix)
require(mice)
require(Amelia)

# one dataset from the two-group design, as in the vignette's competitor chunk
set.seed(1L)
n <- 500L
grp <- runif(n) < 0.5
rho <- ifelse(grp, -0.3, 0.6)
z1 <- rnorm(n)
full <- data.frame(x1 = ifelse(grp, 1.5, -1.5) + z1,
                   x2 = ifelse(grp, 2, -1.5) + rho * z1 + sqrt(1 - rho^2) * rnorm(n))
obs <- full
obs$x2[runif(n) < plogis(0.4 + 0.6 * full$x1)] <- NA

n_time <- 5L
time_mat <- matrix(NA_real_, nrow = n_time, ncol = 3L,
                   dimnames = list(NULL, c("proxymix", "mice", "Amelia")))
for (i1 in seq_len(n_time)) {
  time_mat[i1, "proxymix"] <- system.time(gmm_impute(obs, m = 20L, seed = 1L))[["elapsed"]]
  time_mat[i1, "mice"] <- system.time(mice(obs, m = 20L, seed = 1L, printFlag = FALSE))[["elapsed"]]
  time_mat[i1, "Amelia"] <- system.time(amelia(obs, m = 20L, p2s = 0L))[["elapsed"]]
} # end of for (i1 in seq_len(n_time))

# median seconds for 20 completions of one 500-row dataset
time_vec <- apply(time_mat, 2L, median)
print(time_vec)
time_versions <- c(proxymix = as.character(packageVersion("proxymix")),
                   mice = as.character(packageVersion("mice")),
                   Amelia = as.character(packageVersion("Amelia")))


time_section <- proc.time() - ptm_local
print("Section: 2. Time one imputation per package")
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
                m = sim$m,
                n_rep = sim$n_rep,
                population = sim$population,
                proxymix_version = sim$versions[["proxymix"]],
                versions = desc_version(sim$versions),
                run_date = sim$run_date,
                elapsed_secs = sim$elapsed_secs,
                time_secs = time_vec,
                time_versions = desc_version(time_versions),
                time_platform = R.version$platform)

dir.create(dir_out, showWarnings = FALSE)
file_name_out <- paste0(dir_out, "missing_data.rds")
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
