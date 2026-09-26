# =============================================================================
# Summary of the mixture-package comparison for the short three-regimes vignette
#
# Author      : Max Moldovan
# Version date: 26 Sep 2026
# Usage       : Rscript data-raw/vignette_results/three_regimes.R
#               (run from the package root)
#
# Data        : vignettes/articles/extended/three_regimes_sim.rds, written by
#               data-raw/articles/three_regimes_sim.R
#               datasets::faithful (section 2)
# Requirements: proxymix, mclust, mixtools, flexmix (section 2 only)
#
# Outputs     : vignettes/results/three_regimes.rds, read by vignettes/three_regimes.Rmd
#
# Structure   : 1. Read the stored simulation
#               2. Fit the Old Faithful split with four packages
#               3. Write the summary
# =============================================================================


rm(list = ls())

ptm1 <- proc.time()	# total execution time


# Configuration: paths
########################################################

file_name_sim <- paste0(getwd(), "/vignettes/articles/extended/three_regimes_sim.rds")
dir_out <- paste0(getwd(), "/vignettes/results/")
########################################################


# =============================================================================
# 1. Read the stored simulation
# =============================================================================
print("==== 1. Read the stored simulation ====")
ptm_local <- proc.time()

sim <- readRDS(file_name_sim)

# the version stamp is the proxymix version the simulation ran under
stopifnot(all(c("sim_tab", "res", "n_sizes", "n_rep", "mixtools_secs",
                "versions", "run_date", "elapsed_secs") %in% names(sim)))
stopifnot("proxymix" %in% names(sim$versions))
print(sim$versions)

# a misfit is a fit whose divergence from the true mixture exceeds 0.1
tab_res <- sim$res
tab_res$misfit <- tab_res$kl > 0.1
tab_misfit <- aggregate(misfit ~ n + method, data = tab_res, FUN = sum)
print(tab_misfit)

# paired differences, proxymix regime (ii) minus mclust, over the same datasets
paired_diff <- function(n_local, what) {
	x_local <- tab_res[tab_res$n == n_local & tab_res$method == "proxymix, regime (ii)", ]
	y_local <- tab_res[tab_res$n == n_local & tab_res$method == "mclust", ]
	d_vec <- x_local[[what]][order(x_local$rep)] - y_local[[what]][order(y_local$rep)]
	c(diff = mean(d_vec), se = sd(d_vec) / sqrt(length(d_vec)))
} # end of paired_diff
tab_paired <- rbind(kl_small = paired_diff(sim$n_sizes[1L], "kl"),
                    kl_large = paired_diff(sim$n_sizes[2L], "kl"))
print(tab_paired)

# misfits at the smaller size for one method and not the other, and the exact McNemar test
misfit_of <- function(method_local) {
	x_local <- tab_res[tab_res$n == sim$n_sizes[1L] & tab_res$method == method_local, ]
	x_local$misfit[order(x_local$rep)]
} # end of misfit_of
mis_ii_vec <- misfit_of("proxymix, regime (ii)")
mis_mc_vec <- misfit_of("mclust")
mcnemar <- c(proxymix_only = sum(mis_ii_vec & !mis_mc_vec),
             mclust_only = sum(mis_mc_vec & !mis_ii_vec),
             p_value = binom.test(sum(mis_ii_vec & !mis_mc_vec),
                                  sum(xor(mis_ii_vec, mis_mc_vec)))$p.value)
print(mcnemar)


time_section <- proc.time() - ptm_local
print("Section: 1. Read the stored simulation")
print(time_section)


# =============================================================================
# 2. Fit the Old Faithful split with four packages
# =============================================================================
print("==== 2. Fit the Old Faithful split with four packages ====")
ptm_local <- proc.time()

# the code below is the vignette's compare-code chunk, unchanged
library(proxymix)
library(mclust)
library(mixtools)
library(flexmix)

# split the Old Faithful eruptions into a training half and a held-out half
faithful_mat <- as.matrix(datasets::faithful)
set.seed(20260925)
i_train <- sort(sample.int(nrow(faithful_mat), nrow(faithful_mat) / 2L))
train <- faithful_mat[i_train, ]
test <- faithful_mat[-i_train, ]
train_df <- data.frame(eruptions = train[, 1L], waiting = train[, 2L])

# convert each package's fit to a proxymix mixture, so dgmm() scores all four
as_gmm_mclust <- function(fit) {
  gmm(
    weights = fit$parameters$pro,
    means = lapply(seq_len(fit$G), function(k) fit$parameters$mean[, k]),
    covariances = lapply(seq_len(fit$G), function(k) {
      fit$parameters$variance$sigma[, , k]
    })
  )
}
as_gmm_mixtools <- function(fit) {
  gmm(weights = fit$lambda, means = fit$mu, covariances = fit$sigma)
}
as_gmm_flexmix <- function(fit) {
  comps <- lapply(fit@components, function(cc) cc[[1L]]@parameters)
  gmm(
    weights = prior(fit),
    means = lapply(comps, function(p) unname(p$center)),
    covariances = lapply(comps, function(p) unname(p$cov))
  )
}

# two-component fits to the training half
set.seed(1L)
fits <- list(
  proxymix = fit_proxymix(gmm_target_from_samples(train), N = 2L,
                          regime = "sample"),
  mclust = as_gmm_mclust(Mclust(train, G = 2L, verbose = FALSE)),
  mixtools = as_gmm_mixtools(mvnormalmixEM(train, k = 2L, verb = FALSE)),
  flexmix = as_gmm_flexmix(stepFlexmix(
    cbind(eruptions, waiting) ~ 1, data = train_df, k = 2L, nrep = 5L,
    model = FLXMCmvnorm(diagonal = FALSE), verbose = FALSE
  ))
)

# mean log-likelihood of the held-out eruptions under each fit
# (the vignette chunk prints this vector without assigning it)
loglik_vec <- vapply(fits, function(g) mean(dgmm(test, g, log = TRUE)), numeric(1L))
print(loglik_vec)
print(paste0("range: ", format(max(loglik_vec) - min(loglik_vec))))
faithful_versions <- c(proxymix = as.character(packageVersion("proxymix")),
                       mclust = as.character(packageVersion("mclust")),
                       mixtools = as.character(packageVersion("mixtools")),
                       flexmix = as.character(packageVersion("flexmix")))


time_section <- proc.time() - ptm_local
print("Section: 2. Fit the Old Faithful split with four packages")
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
                misfit_tab = tab_misfit,
                paired_kl = tab_paired,
                mcnemar = mcnemar,
                n_sizes = sim$n_sizes,
                n_rep = sim$n_rep,
                mixtools_secs = sim$mixtools_secs,
                proxymix_version = sim$versions[["proxymix"]],
                versions = desc_version(sim$versions),
                run_date = sim$run_date,
                elapsed_secs = sim$elapsed_secs,
                faithful_loglik = loglik_vec,
                faithful_n_train = nrow(train),
                faithful_n_test = nrow(test),
                faithful_versions = desc_version(faithful_versions),
                faithful_platform = R.version$platform)

dir.create(dir_out, showWarnings = FALSE)
file_name_out <- paste0(dir_out, "three_regimes.rds")
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
