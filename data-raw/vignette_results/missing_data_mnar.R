# =============================================================================
# Summary of the stored missing-not-at-random benchmark for the short vignette
#
# Author      : Max Moldovan
# Version date: 26 Sep 2026
# Usage       : Rscript data-raw/vignette_results/missing_data_mnar.R
#               (run from the package root)
#
# Data        : vignettes/articles/extended/missing_data_mnar_sim.rds, written by
#               data-raw/articles/missing_data_mnar_sim.R
# Requirements: proxymix, mice, Amelia, AER, survival, sampleSelection (section 2 only)
#
# Outputs     : vignettes/results/missing_data_mnar.rds, read by vignettes/missing_data_mnar.Rmd
#
# Structure   : 1. Read the stored simulation
#               2. Time one dataset per method
#               3. Write the summary
# =============================================================================


rm(list = ls())

ptm1 <- proc.time()	# total execution time


# Configuration: paths
########################################################

file_name_sim <- paste0(getwd(), "/vignettes/articles/extended/missing_data_mnar_sim.rds")
dir_out <- paste0(getwd(), "/vignettes/results/")
########################################################


# =============================================================================
# 1. Read the stored simulation
# =============================================================================
print("==== 1. Read the stored simulation ====")
ptm_local <- proc.time()

sim <- readRDS(file_name_sim)

# the version stamp is the proxymix version the simulation ran under
stopifnot(all(c("sim_tab", "n", "m", "n_rep", "beta_grid", "delta_grid",
                "truth", "n_undefined", "warnings", "versions", "run_date",
                "elapsed_secs") %in% names(sim)))
stopifnot("proxymix" %in% names(sim$versions))
print(sim$versions)
print(paste0("warnings in the stored run: ", length(sim$warnings)))


time_section <- proc.time() - ptm_local
print("Section: 1. Read the stored simulation")
print(time_section)


# =============================================================================
# 2. Time one dataset per method
# =============================================================================
print("==== 2. Time one dataset per method ====")
ptm_local <- proc.time()

require(proxymix)
require(mice)
require(Amelia)
require(AER)
require(survival)
require(sampleSelection)

# one dataset from each design, as in the vignette's competitor chunk
set.seed(1L)
n <- 300L
comp <- sample(1:2, n, replace = TRUE)
mu <- rbind(c(0, 0), c(1.5, 0.5))
chol_r <- chol(matrix(c(1, 0.6, 0.6, 1), 2L))
z <- matrix(rnorm(2 * n), n, 2L) %*% chol_r + mu[comp, ]
full <- data.frame(x1 = z[, 1L], y = z[, 2L])
obs <- full
obs$y[runif(n) < plogis(-0.5 + 0.7 * full$y)] <- NA
obs_sel <- obs
obs_sel$seen <- !is.na(obs$y)

set.seed(1L)
a_cens <- -sqrt(2) * qnorm(0.3)
x <- rnorm(n)
y_star <- a_cens + x + rnorm(n)
cens <- data.frame(x = x, y = pmax(y_star, 0))
holes <- cens
holes$y[cens$y == 0] <- NA

n_time <- 5L
method_vec <- c("proxymix sweep", "mice sweep", "Amelia", "heckit",
                "proxymix censored", "tobit")
time_mat <- matrix(NA_real_, nrow = n_time, ncol = length(method_vec),
                   dimnames = list(NULL, method_vec))
for (i1 in seq_len(n_time)) {
  time_mat[i1, "proxymix sweep"] <- system.time(
    proxy_mnar_sensitivity(obs, "y", beta_grid = c(0, 0.35, 0.7, 1.05),
                           m = 10L, seed = 1L))[["elapsed"]]
  time_mat[i1, "mice sweep"] <- system.time(
    lapply(c(0, 0.2, 0.4, 0.6), function(delta) {
      post <- make.post(obs)
      post["y"] <- paste0("imp[[j]][, i] <- imp[[j]][, i] + ", delta)
      imp <- mice(obs, m = 10L, method = "norm", post = post, seed = 1L,
                  printFlag = FALSE)
      summary(pool(with(imp, lm(y ~ 1))), conf.int = TRUE)
    }))[["elapsed"]]
  time_mat[i1, "Amelia"] <- system.time(
    summary(pool(lapply(amelia(obs, m = 10L, p2s = 0L)$imputations,
                        function(d) lm(y ~ 1, data = d))),
            conf.int = TRUE))[["elapsed"]]
  time_mat[i1, "heckit"] <- system.time(
    heckit(seen ~ x1, y ~ x1, data = obs_sel))[["elapsed"]]
  time_mat[i1, "proxymix censored"] <- system.time(
    summary(pool(lapply(complete(as_mids(gmm_impute(
      holes, m = 10L, mechanism = censored("y", upper = 0), seed = 1L)), "all"),
      function(d) lm(y ~ x, data = d))), conf.int = TRUE))[["elapsed"]]
  time_mat[i1, "tobit"] <- system.time(
    tobit(y ~ x, data = cens))[["elapsed"]]
} # end of for (i1 in seq_len(n_time))

# median seconds for one 300-row dataset, ten completions per imputer
time_vec <- apply(time_mat, 2L, median)
print(time_mat)
print(time_vec)
time_versions <- c(proxymix = as.character(packageVersion("proxymix")),
                   mice = as.character(packageVersion("mice")),
                   Amelia = as.character(packageVersion("Amelia")),
                   AER = as.character(packageVersion("AER")),
                   survival = as.character(packageVersion("survival")),
                   sampleSelection = as.character(packageVersion("sampleSelection")))


# the stored versions as each DESCRIPTION writes them (AER 1.2-17, not 1.2.17),
# taken from the installed package only when it is the version the simulation used
pkg_vec <- setdiff(names(sim$versions), "R")
versions_desc <- vapply(pkg_vec, function(s1) {
  if (identical(as.character(packageVersion(s1)), sim$versions[[s1]])) {
    utils::packageDescription(s1)$Version
  } else {
    sim$versions[[s1]]
  }
}, character(1L))
print(versions_desc)


time_section <- proc.time() - ptm_local
print("Section: 2. Time one dataset per method")
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
                beta_grid = sim$beta_grid,
                delta_grid = sim$delta_grid,
                truth = sim$truth,
                n_undefined = sim$n_undefined,
                n_warnings = length(sim$warnings),
                proxymix_version = sim$versions[["proxymix"]],
                versions = desc_version(sim$versions),
                versions_desc = versions_desc,
                run_date = sim$run_date,
                elapsed_secs = sim$elapsed_secs,
                time_secs = time_vec,
                time_versions = desc_version(time_versions),
                time_platform = R.version$platform)

dir.create(dir_out, showWarnings = FALSE)
file_name_out <- paste0(dir_out, "missing_data_mnar.rds")
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
