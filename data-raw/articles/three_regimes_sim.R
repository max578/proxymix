# =============================================================================
# Stored results for the simulation benchmark in the extended three-regimes vignette
#
# Author      : Max Moldovan
# Version date: 25 Sep 2026
# Usage       : Rscript data-raw/articles/three_regimes_sim.R
#               Rscript data-raw/articles/three_regimes_sim.R smoke
#               (run from the package root)
#
# Data        : simulated by the vignette chunk sim-code, which this script runs unchanged
# Requirements: proxymix, mclust, flexmix (assumed installed)
#
# Outputs     : vignettes/articles/extended/three_regimes_sim.rds, read by the vignette when it renders
#               three_regimes_sim_smoke.rds in tempdir(), from a smoke run of 3 datasets per sample size
#               data-raw/articles/three_regimes_sessionInfo.txt
#
# Structure   : 1. Extract the simulation chunk
#               2. Run the simulation
#               3. Write the artefacts
# =============================================================================


rm(list = ls())

ptm1 <- proc.time()	# total execution time


# Configuration: switches and paths
########################################################

# add "smoke" on the command line for a 3-dataset pipeline check
smoke_run <- "smoke" %in% commandArgs(trailingOnly = TRUE)

dir_art <- paste0(getwd(), "/vignettes/articles/extended/")
dir_gen <- paste0(getwd(), "/data-raw/articles/")
########################################################


# =============================================================================
# 1. Extract the simulation chunk
# =============================================================================
print("==== 1. Extract the simulation chunk ====")
ptm_local <- proc.time()

lines_vec <- readLines(paste0(dir_art, "three_regimes.Rmd"), warn = FALSE)

# the chunk runs from its opening fence to the next closing fence
i_start <- grep("^```\\{r sim-code", lines_vec)
stopifnot(length(i_start) == 1L)
i_end <- i_start + which(lines_vec[(i_start + 1L):length(lines_vec)] == "```")[1L]
code_vec <- lines_vec[(i_start + 1L):(i_end - 1L)]

# the smoke run changes one line, the number of datasets per design
code_run_vec <- code_vec
if (smoke_run) {

	code_run_vec <- sub("^n_rep <- 300L", "n_rep <- 3L", code_run_vec)
	stopifnot(sum(code_run_vec != code_vec) == 1L)

} # ends smoke_run

print(paste0("chunk lines: ", length(code_vec), ", smoke run: ", smoke_run))


time_section <- proc.time() - ptm_local
print("Section: 1. Extract the simulation chunk")
print(time_section)


# =============================================================================
# 2. Run the simulation
# =============================================================================
print("==== 2. Run the simulation ====")
ptm_local <- proc.time()

env_sim <- new.env(parent = globalenv())
warn_vec <- c()

# warnings are recorded and reported at the end of the chapter, not dropped
withCallingHandlers(
	eval(parse(text = code_run_vec), envir = env_sim),
	warning = function(w) {
		warn_vec <<- c(warn_vec, conditionMessage(w))
		invokeRestart("muffleWarning")
	}
)

print(paste0("warnings raised: ", length(warn_vec)))
if (length(warn_vec) > 0L) print(table(warn_vec))


time_section <- proc.time() - ptm_local
print("Section: 2. Run the simulation")
print(time_section)


# =============================================================================
# 3. Write the artefacts
# =============================================================================
print("==== 3. Write the artefacts ====")
ptm_local <- proc.time()

versions_vec <- c(R = as.character(getRversion()),
                  proxymix = as.character(packageVersion("proxymix")),
                  mclust = as.character(packageVersion("mclust")),
                  flexmix = as.character(packageVersion("flexmix")),
                  mixtools = as.character(packageVersion("mixtools")))

sim_out <- list(sim_tab = env_sim$sim_tab,
                code = code_vec,
                n_sizes = env_sim$n_sizes,
                n_rep = env_sim$n_rep,
                mixtools_secs = env_sim$mixtools_secs,
                res = env_sim$res,
                warnings = warn_vec,
                versions = versions_vec,
                run_date = format(Sys.Date(), "%d %B %Y"),
                elapsed_secs = round((proc.time() - ptm1)[[3L]]))

file_name_out <- paste0(dir_art, "three_regimes_sim.rds")
if (smoke_run) file_name_out <- file.path(tempdir(), "three_regimes_sim_smoke.rds")

saveRDS(sim_out, file_name_out)
print(paste0("written: ", file_name_out))
print(sim_out$sim_tab)


time_section <- proc.time() - ptm_local
print("Section: 3. Write the artefacts")
print(time_section)


# =============================================================================
# Completion
# =============================================================================
writeLines(capture.output(sessionInfo()), paste0(dir_gen, "three_regimes_sessionInfo.txt"))

exec_time1 <- proc.time() - ptm1
print(exec_time1)
