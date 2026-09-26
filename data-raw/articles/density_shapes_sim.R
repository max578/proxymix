# =============================================================================
# Stored results for the simulation benchmark in the extended density-shapes vignette
#
# Author      : Max Moldovan
# Version date: 25 Sep 2026
# Usage       : Rscript data-raw/articles/density_shapes_sim.R
#               Rscript data-raw/articles/density_shapes_sim.R smoke
#               (run from the package root)
#
# Data        : simulated by the vignette chunk sim-code, which this script runs unchanged
# Requirements: proxymix, cmdstanr with CmdStan installed, BayesianTools, mclust
#
# Outputs     : vignettes/articles/extended/density_shapes_sim.rds, read by the vignette when it renders
#               density_shapes_sim_smoke.rds in tempdir(), from a smoke run of 2 replicates per shape
#               data-raw/articles/density_shapes_stan/*.stan and their executables, written by the chunk
#               data-raw/articles/density_shapes_sessionInfo.txt
#
# Structure   : 1. Extract the simulation chunk
#               2. Run the simulation
#               3. Write the artefacts
# =============================================================================


rm(list = ls())

ptm1 <- proc.time()	# total execution time


# Configuration: switches and paths
########################################################

# add "smoke" on the command line for a 2-replicate pipeline check
smoke_run <- "smoke" %in% commandArgs(trailingOnly = TRUE)

dir_art <- paste0(getwd(), "/vignettes/articles/extended/")
dir_gen <- paste0(getwd(), "/data-raw/articles/")

# one chain at a time: the cores are shared
options(cmdstanr_no_ver_check = TRUE, mc.cores = 1L)
########################################################


# =============================================================================
# 1. Extract the simulation chunk
# =============================================================================
print("==== 1. Extract the simulation chunk ====")
ptm_local <- proc.time()

lines_vec <- readLines(paste0(dir_art, "density_shapes.Rmd"), warn = FALSE)

# the chunk runs from its opening fence to the next closing fence
i_start <- grep("^```\\{r sim-code", lines_vec)
stopifnot(length(i_start) == 1L)
i_end <- i_start + which(lines_vec[(i_start + 1L):length(lines_vec)] == "```")[1L]
code_vec <- lines_vec[(i_start + 1L):(i_end - 1L)]

# the smoke run changes one line, the number of replicates per shape
code_run_vec <- code_vec
if (smoke_run) {

	code_run_vec <- sub("^n_rep <- 60L", "n_rep <- 2L", code_run_vec)
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

# the chunk writes its Stan programs to density_shapes_stan/ under the
# working directory, which is data-raw/articles/ for the length of the run
setwd(dir_gen)

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
                  cmdstanr = as.character(packageVersion("cmdstanr")),
                  CmdStan = as.character(cmdstanr::cmdstan_version()),
                  BayesianTools = as.character(packageVersion("BayesianTools")),
                  mclust = as.character(packageVersion("mclust")))

# the event scored on each shape, and its mass under the target by quadrature
event_tab <- data.frame(
	shape = names(env_sim$shapes),
	event = vapply(env_sim$shapes, function(s) s$event_name, character(1L)),
	truth = vapply(env_sim$grids, function(g) sum(g$f[g$in_event]) * g$cell, numeric(1L)),
	stringsAsFactors = FALSE
)

sim_out <- list(sim_tab = env_sim$sim_tab,
                res = env_sim$res,
                event_tab = event_tab,
                code = code_vec,
                n_rep = env_sim$n_rep,
                budget = env_sim$budget,
                warnings = warn_vec,
                versions = versions_vec,
                run_date = format(Sys.Date(), "%d %B %Y"),
                elapsed_secs = round((proc.time() - ptm1)[[3L]]),
                # processor seconds of this process and its Stan child processes
                cpu_secs = round(sum((proc.time() - ptm1)[c(1L, 2L, 4L, 5L)])))

file_name_out <- paste0(dir_art, "density_shapes_sim.rds")
if (smoke_run) file_name_out <- file.path(tempdir(), "density_shapes_sim_smoke.rds")

saveRDS(sim_out, file_name_out)
print(paste0("written: ", file_name_out))
print(sim_out$event_tab)
print(sim_out$sim_tab, digits = 4L)


time_section <- proc.time() - ptm_local
print("Section: 3. Write the artefacts")
print(time_section)


# =============================================================================
# Completion
# =============================================================================
writeLines(capture.output(sessionInfo()), paste0(dir_gen, "density_shapes_sessionInfo.txt"))

exec_time1 <- proc.time() - ptm1
print(exec_time1)
