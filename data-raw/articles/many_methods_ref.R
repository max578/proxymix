# =============================================================================
# Stored tmle reference estimate for the extended many-methods vignette
#
# Author      : Max Moldovan
# Version date: 25 Sep 2026
# Usage       : Rscript data-raw/articles/many_methods_ref.R
#               (run from the package root)
#
# Data        : MatchIt::lalonde, analysed by the vignette chunk ref-code, which this script runs unchanged
# Requirements: tmle, SuperLearner, dbarts, glmnet, gam, MatchIt (assumed installed)
#
# Outputs     : vignettes/articles/extended/many_methods_ref.rds, read by the vignette when it renders
#               data-raw/articles/many_methods_ref_sessionInfo.txt
#
# Structure   : 1. Extract the reference chunk
#               2. Run the reference computation
#               3. Write the artefacts
# =============================================================================


rm(list = ls())

ptm1 <- proc.time()	# total execution time


# Configuration: switches and paths
########################################################

dir_art <- paste0(getwd(), "/vignettes/articles/extended/")
dir_gen <- paste0(getwd(), "/data-raw/articles/")
########################################################


# =============================================================================
# 1. Extract the reference chunk
# =============================================================================
print("==== 1. Extract the reference chunk ====")
ptm_local <- proc.time()

lines_vec <- readLines(paste0(dir_art, "many_methods.Rmd"), warn = FALSE)

# the chunk runs from its opening fence to the next closing fence
i_start <- grep("^```\\{r ref-code", lines_vec)
stopifnot(length(i_start) == 1L)
i_end <- i_start + which(lines_vec[(i_start + 1L):length(lines_vec)] == "```")[1L]
code_vec <- lines_vec[(i_start + 1L):(i_end - 1L)]

print(paste0("chunk lines: ", length(code_vec)))


time_section <- proc.time() - ptm_local
print("Section: 1. Extract the reference chunk")
print(time_section)


# =============================================================================
# 2. Run the reference computation
# =============================================================================
print("==== 2. Run the reference computation ====")
ptm_local <- proc.time()

env_sim <- new.env(parent = globalenv())
warn_vec <- c()

# warnings are recorded and reported at the end of the chapter, not dropped
withCallingHandlers(
	eval(parse(text = code_vec), envir = env_sim),
	warning = function(w) {
		warn_vec <<- c(warn_vec, conditionMessage(w))
		invokeRestart("muffleWarning")
	}
)

print(paste0("warnings raised: ", length(warn_vec)))
if (length(warn_vec) > 0L) print(table(warn_vec))


time_section <- proc.time() - ptm_local
print("Section: 2. Run the reference computation")
print(time_section)


# =============================================================================
# 3. Write the artefacts
# =============================================================================
print("==== 3. Write the artefacts ====")
ptm_local <- proc.time()

pkgs_vec <- c("tmle", "SuperLearner", "dbarts", "glmnet", "gam", "MatchIt")
versions_vec <- c(R = as.character(getRversion()),
                  vapply(pkgs_vec, function(p) as.character(packageVersion(p)),
                         character(1L)))

sim_out <- list(att_tmle = env_sim$att_tmle,
                code = code_vec,
                warnings = warn_vec,
                versions = versions_vec,
                run_date = format(Sys.Date(), "%d %B %Y"),
                elapsed_secs = round((proc.time() - ptm1)[[3L]]))

file_name_out <- paste0(dir_art, "many_methods_ref.rds")

saveRDS(sim_out, file_name_out)
print(paste0("written: ", file_name_out))
print(sim_out$att_tmle)


time_section <- proc.time() - ptm_local
print("Section: 3. Write the artefacts")
print(time_section)


# =============================================================================
# Completion
# =============================================================================
writeLines(capture.output(sessionInfo()), paste0(dir_gen, "many_methods_ref_sessionInfo.txt"))

exec_time1 <- proc.time() - ptm1
print(exec_time1)
