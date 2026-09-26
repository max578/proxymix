# =============================================================================
# Stored NUTS reference for the Affairs illustration in the extended posterior-proxy vignette
#
# Author      : Max Moldovan
# Version date: 25 Sep 2026
# Usage       : Rscript data-raw/articles/posterior_proxy_ref.R
#               (run from the package root)
#
# Data        : AER::Affairs, read by the vignette chunk affairs-data
# Requirements: AER, cmdstanr (with CmdStan), bridgesampling (assumed installed)
#
# Outputs     : vignettes/articles/extended/posterior_proxy_ref.rds, read by the vignette when it renders
#
# Structure   : 1. Extract the chunks
#               2. Run the reference
#               3. Write the artefact
# =============================================================================


rm(list = ls())

ptm1 <- proc.time()	# total execution time


# Configuration: paths
########################################################

dir_art <- paste0(getwd(), "/vignettes/articles/extended/")
dir_gen <- paste0(getwd(), "/data-raw/articles/")
########################################################


# =============================================================================
# 1. Extract the chunks
# =============================================================================
print("==== 1. Extract the chunks ====")
ptm_local <- proc.time()

lines_vec <- readLines(paste0(dir_art, "posterior_proxy.Rmd"), warn = FALSE)

# each chunk runs from its opening fence to the next closing fence; the
# reference chunk uses the data and log posterior of the data chunk
chunk_lines <- function(label) {
	i_start <- grep(paste0("^```\\{r ", label, "[,}]"), lines_vec)
	stopifnot(length(i_start) == 1L)
	i_end <- i_start + which(lines_vec[(i_start + 1L):length(lines_vec)] == "```")[1L]
	lines_vec[(i_start + 1L):(i_end - 1L)]
}
code_vec <- c(chunk_lines("affairs-data"), chunk_lines("ref-code"))

print(paste0("chunk lines: ", length(code_vec)))


time_section <- proc.time() - ptm_local
print("Section: 1. Extract the chunks")
print(time_section)


# =============================================================================
# 2. Run the reference
# =============================================================================
print("==== 2. Run the reference ====")
ptm_local <- proc.time()

env_ref <- new.env(parent = globalenv())
warn_vec <- c()

withCallingHandlers(
	eval(parse(text = code_vec), envir = env_ref),
	warning = function(w) {
		warn_vec <<- c(warn_vec, conditionMessage(w))
		invokeRestart("muffleWarning")
	}
)

print(paste0("warnings raised: ", length(warn_vec)))
if (length(warn_vec) > 0L) print(table(warn_vec))

# the Stan program the chunk wrote must be the one kept beside this script
stan_chunk_vec <- readLines(env_ref$stan_file)
stan_kept_vec <- readLines(paste0(dir_gen, "posterior_proxy_stan/probit.stan"))
stan_kept_vec <- stan_kept_vec[!grepl("^//", stan_kept_vec)]
stopifnot(identical(stan_chunk_vec, stan_kept_vec))


time_section <- proc.time() - ptm_local
print("Section: 2. Run the reference")
print(time_section)


# =============================================================================
# 3. Write the artefact
# =============================================================================
print("==== 3. Write the artefact ====")
ptm_local <- proc.time()

versions_vec <- c(R = as.character(getRversion()),
                  AER = as.character(packageVersion("AER")),
                  cmdstanr = as.character(packageVersion("cmdstanr")),
                  CmdStan = as.character(cmdstanr::cmdstan_version()),
                  bridgesampling = as.character(packageVersion("bridgesampling")))

ref_out <- c(env_ref$ref,
             list(code = code_vec,
                  warnings = warn_vec,
                  versions = versions_vec,
                  run_date = format(Sys.Date(), "%d %B %Y"),
                  elapsed_secs = round((proc.time() - ptm1)[[3L]])))

file_name_out <- paste0(dir_art, "posterior_proxy_ref.rds")
saveRDS(ref_out, file_name_out)
print(paste0("written: ", file_name_out))
print(ref_out$long$summary)
print(ref_out$short$summary)
print(c(log_z_long = ref_out$long$log_z, log_z_short = ref_out$short$log_z))
print(c(ess_min_long = ref_out$long$ess_min, rhat_max_long = ref_out$long$rhat_max))
print(c(secs_long = ref_out$long$secs_sample, secs_short = ref_out$short$secs_sample))


time_section <- proc.time() - ptm_local
print("Section: 3. Write the artefact")
print(time_section)


# =============================================================================
# Completion
# =============================================================================
exec_time1 <- proc.time() - ptm1
print(exec_time1)
