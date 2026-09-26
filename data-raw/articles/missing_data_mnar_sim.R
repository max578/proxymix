# =============================================================================
# Stored results for the simulation benchmark in the extended missing-not-at-random vignette
#
# Author      : Max Moldovan
# Version date: 25 Sep 2026
# Usage       : Rscript data-raw/articles/missing_data_mnar_sim.R
#               Rscript data-raw/articles/missing_data_mnar_sim.R smoke
#               (run from the package root)
#
# Data        : simulated by the vignette chunk sim-code, which this script runs unchanged
# Requirements: proxymix, mice, Amelia, AER, survival, sampleSelection (assumed installed)
#
# Outputs     : vignettes/articles/extended/missing_data_mnar_sim.rds, read by the vignette when it renders
#               missing_data_mnar_sim_smoke.rds in tempdir(), from a smoke run of 3 datasets per design
#               data-raw/articles/missing_data_mnar_sessionInfo.txt
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

lines_vec <- readLines(paste0(dir_art, "missing_data_mnar.Rmd"), warn = FALSE)

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
                  mice = as.character(packageVersion("mice")),
                  Amelia = as.character(packageVersion("Amelia")),
                  AER = as.character(packageVersion("AER")),
                  survival = as.character(packageVersion("survival")),
                  sampleSelection = as.character(packageVersion("sampleSelection")))

sim_out <- list(sim_tab = env_sim$sim_tab,
                code = code_vec,
                n = env_sim$n,
                m = env_sim$m,
                n_rep = env_sim$n_rep,
                beta_grid = env_sim$beta_grid,
                delta_grid = env_sim$delta_grid,
                truth = env_sim$truth,
                n_undefined = env_sim$n_undefined,
                warnings = warn_vec,
                versions = versions_vec,
                run_date = format(Sys.Date(), "%d %B %Y"),
                elapsed_secs = round((proc.time() - ptm1)[[3L]]))

file_name_out <- paste0(dir_art, "missing_data_mnar_sim.rds")
if (smoke_run) file_name_out <- file.path(tempdir(), "missing_data_mnar_sim_smoke.rds")

saveRDS(sim_out, file_name_out)
print(paste0("written: ", file_name_out))
print(sim_out$sim_tab)


time_section <- proc.time() - ptm_local
print("Section: 3. Write the artefacts")
print(time_section)


# =============================================================================
# Completion
# =============================================================================
writeLines(capture.output(sessionInfo()), paste0(dir_gen, "missing_data_mnar_sessionInfo.txt"))

exec_time1 <- proc.time() - ptm1
print(exec_time1)
