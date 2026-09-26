# =============================================================================
# Stored NUTS reference for the Old Faithful illustration in the extended three-regimes vignette
#
# Author      : Max Moldovan
# Version date: 25 Sep 2026
# Usage       : Rscript data-raw/articles/three_regimes_ref.R
#               (run from the package root)
#
# Data        : datasets::faithful, split by the vignette chunk ref-code, which this script runs unchanged
# Requirements: cmdstanr with CmdStan installed, ks, mclust (assumed installed)
#
# Outputs     : vignettes/articles/extended/three_regimes_ref.rds, read by the vignette when it renders
#
# Structure   : 1. Extract the reference chunk
#               2. Run the sampler
#               3. Write the artefacts
# =============================================================================


rm(list = ls())

ptm1 <- proc.time()	# total execution time


# Configuration: paths
########################################################

dir_art <- paste0(getwd(), "/vignettes/articles/extended/")
dir_gen <- paste0(getwd(), "/data-raw/articles/")
########################################################


# =============================================================================
# 1. Extract the reference chunk
# =============================================================================
print("==== 1. Extract the reference chunk ====")
ptm_local <- proc.time()

lines_vec <- readLines(paste0(dir_art, "three_regimes.Rmd"), warn = FALSE)

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
# 2. Run the sampler
# =============================================================================
print("==== 2. Run the sampler ====")
ptm_local <- proc.time()

# the chunk names the Stan file relative to the data-raw/articles/ folder
dir_root <- getwd()
setwd(dir_gen)

env_ref <- new.env(parent = globalenv())
warn_vec <- c()

# warnings are recorded and reported at the end of the chapter, not dropped
withCallingHandlers(
	eval(parse(text = code_vec), envir = env_ref),
	warning = function(w) {
		warn_vec <<- c(warn_vec, conditionMessage(w))
		invokeRestart("muffleWarning")
	}
)

setwd(dir_root)

print(paste0("warnings raised: ", length(warn_vec)))
if (length(warn_vec) > 0L) print(table(warn_vec))
print(env_ref$nuts_summary)
print(env_ref$nuts$diagnostic_summary())


time_section <- proc.time() - ptm_local
print("Section: 2. Run the sampler")
print(time_section)


# =============================================================================
# 3. Write the artefacts
# =============================================================================
print("==== 3. Write the artefacts ====")
ptm_local <- proc.time()

versions_vec <- c(R = as.character(getRversion()),
                  cmdstanr = as.character(packageVersion("cmdstanr")),
                  CmdStan = as.character(cmdstanr::cmdstan_version()),
                  ks = as.character(packageVersion("ks")),
                  mclust = as.character(packageVersion("mclust")))

ref_out <- list(code = code_vec,
                i_train = env_ref$i_train,
                H = env_ref$kde_fit$H,
                n_draws = nrow(env_ref$draws),
                nuts_summary = env_ref$nuts_summary,
                nuts_divergent = env_ref$nuts_divergent,
                two_stage_par = env_ref$two_stage_par,
                warnings = warn_vec,
                versions = versions_vec,
                run_date = format(Sys.Date(), "%d %B %Y"),
                elapsed_secs = round((proc.time() - ptm1)[[3L]]))

file_name_out <- paste0(dir_art, "three_regimes_ref.rds")
saveRDS(ref_out, file_name_out)
print(paste0("written: ", file_name_out))
str(ref_out$two_stage_par)


time_section <- proc.time() - ptm_local
print("Section: 3. Write the artefacts")
print(time_section)


# =============================================================================
# Completion
# =============================================================================
exec_time1 <- proc.time() - ptm1
print(exec_time1)
