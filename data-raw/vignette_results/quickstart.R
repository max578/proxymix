# =============================================================================
# Summary of the stored banana benchmark for the short quickstart vignette
#
# Author      : Max Moldovan
# Version date: 26 Sep 2026
# Usage       : Rscript data-raw/vignette_results/quickstart.R
#               (run from the package root)
#
# Data        : vignettes/articles/extended/quickstart_sim.rds, written by
#               data-raw/articles/quickstart_sim.R
# Requirements: none beyond base R; the simulation records its own run times
#
# Outputs     : vignettes/results/quickstart.rds, read by vignettes/quickstart.Rmd
#
# Structure   : 1. Read the stored simulation
#               2. Write the summary
# =============================================================================


rm(list = ls())

ptm1 <- proc.time()	# total execution time


# Configuration: paths
########################################################

file_name_sim <- paste0(getwd(), "/vignettes/articles/extended/quickstart_sim.rds")
dir_out <- paste0(getwd(), "/vignettes/results/")
########################################################


# =============================================================================
# 1. Read the stored simulation
# =============================================================================
print("==== 1. Read the stored simulation ====")
ptm_local <- proc.time()

sim <- readRDS(file_name_sim)

# the version stamp is the proxymix version the simulation ran under
stopifnot(all(c("sim_tab", "n_rep", "h", "tail_ref", "tuning", "versions",
                "run_date", "elapsed_secs") %in% names(sim)))
stopifnot("proxymix" %in% names(sim$versions))
print(sim$versions)

# the simulation times each method inside every replicate (median seconds)
stopifnot(!anyNA(sim$sim_tab$secs))
print(sim$sim_tab[, c("method", "kl_mean", "tail_rmse", "secs")])


time_section <- proc.time() - ptm_local
print("Section: 1. Read the stored simulation")
print(time_section)


# =============================================================================
# 2. Write the summary
# =============================================================================
print("==== 2. Write the summary ====")
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
                n_rep = sim$n_rep,
                h = sim$h,
                tail_ref = sim$tail_ref,
                proxymix_version = sim$versions[["proxymix"]],
                versions = desc_version(sim$versions),
                run_date = sim$run_date,
                elapsed_secs = sim$elapsed_secs)

dir.create(dir_out, showWarnings = FALSE)
file_name_out <- paste0(dir_out, "quickstart.rds")
saveRDS(res_out, file_name_out)

print(paste0("written: ", file_name_out, " (", file.size(file_name_out), " bytes)"))
print(paste0("version stamp: ", res_out$proxymix_version))


time_section <- proc.time() - ptm_local
print("Section: 2. Write the summary")
print(time_section)


# =============================================================================
# Completion
# =============================================================================
exec_time1 <- proc.time() - ptm1
print(exec_time1)
