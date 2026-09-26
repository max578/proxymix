# =============================================================================
# Summary of the stored kernel-density benchmark for the short from_kde vignette
#
# Author      : Max Moldovan
# Version date: 26 Sep 2026
# Usage       : Rscript data-raw/vignette_results/from_kde.R
#               (run from the package root)
#
# Data        : vignettes/articles/extended/from_kde_sim.rds, written by
#               data-raw/articles/from_kde_sim.R; the faithful and penguins
#               data shipped with R (penguins from R 4.5.0)
# Requirements: proxymix, ks, np, mclust, mixtools, KernSmooth (section 3 only)
#
# Outputs     : vignettes/results/from_kde.rds, read by vignettes/from_kde.Rmd
#
# Structure   : 1. Read the stored simulation
#               2. Paired differences in integrated squared error
#               3. Time conditional queries per package
#               4. Write the summary
# =============================================================================


rm(list = ls())

ptm1 <- proc.time()	# total execution time


# Configuration: paths
########################################################

file_name_sim <- paste0(getwd(), "/vignettes/articles/extended/from_kde_sim.rds")
dir_out <- paste0(getwd(), "/vignettes/results/")
########################################################


# =============================================================================
# 1. Read the stored simulation
# =============================================================================
print("==== 1. Read the stored simulation ====")
ptm_local <- proc.time()

sim <- readRDS(file_name_sim)

# the version stamp is the proxymix version the simulation ran under
stopifnot(all(c("sim_tab", "res", "n", "n_rep", "grid_side", "versions",
                "run_date", "elapsed_secs") %in% names(sim)))
stopifnot("proxymix" %in% names(sim$versions))
print(sim$versions)
print(sim$sim_tab)


time_section <- proc.time() - ptm_local
print("Section: 1. Read the stored simulation")
print(time_section)


# =============================================================================
# 2. Paired differences in integrated squared error
# =============================================================================
print("==== 2. Paired differences in integrated squared error ====")
ptm_local <- proc.time()

# each method minus ks on the same datasets, in thousandths, with its standard error
tab_res <- sim$res
method_vec <- setdiff(unique(tab_res$method), "ks, kde")
tab_ks <- tab_res[tab_res$method == "ks, kde", ]
tab_ks <- tab_ks[order(tab_ks$rep), ]

s1 <- length(method_vec)
tab_paired <- data.frame(method = method_vec, mean = NA_real_, se = NA_real_)
for (i1 in seq_len(s1)) {
  tab_local <- tab_res[tab_res$method == method_vec[i1], ]
  tab_local <- tab_local[order(tab_local$rep), ]
  stopifnot(identical(tab_local$rep, tab_ks$rep))
  d_vec <- 1000 * (tab_local$ise - tab_ks$ise)
  tab_paired$mean[i1] <- mean(d_vec)
  tab_paired$se[i1] <- sd(d_vec) / sqrt(length(d_vec))
} # end of for (i1 in seq_len(s1))
print(tab_paired)

# number of datasets in which mclust chose each number of components
mclust_G <- table(tab_res$G[tab_res$method == "mclust, G by BIC"])
print(mclust_G)


time_section <- proc.time() - ptm_local
print("Section: 2. Paired differences in integrated squared error")
print(time_section)


# =============================================================================
# 3. Time conditional queries per package
# =============================================================================
print("==== 3. Time conditional queries per package ====")
ptm_local <- proc.time()

require(proxymix)
require(mclust)
require(mixtools)
require(np)
options(np.messages = FALSE)

# the same conditional queries as the extended article's real-data illustration
as_gmm <- function(w, mu, sigma) {
  gmm(weights = w, means = mu, covariances = sigma)
}

# conditional mean and 90% quantile from a density tabulated on y_grid
grid_summary <- function(y_grid, f) {
  w <- pmax(f, 0) / sum(pmax(f, 0))
  c(mean = sum(w * y_grid),
    q90 = approx(cumsum(w), y_grid, xout = 0.9, ties = "ordered")$y)
}

gmm_summary <- function(g, x0) {
  s <- gmm_conditionalise(g, given = c(x0, NA))
  c(mean = gmm_mean(s), q90 = qgmm(0.9, s))
}

# median milliseconds per query over n_time timings, per method; each timing
# repeats the sweep of 20 queries n_inner times, so that the fastest methods
# take long enough to register on the system clock
query_times <- function(xy, n_time, n_inner) {
  h <- sqrt(diag(ks::Hpi.diag(xy)))
  x_eval <- quantile(xy[, 1L], probs = seq(0.05, 0.95, length.out = 20L),
                     names = FALSE)
  y_grid <- seq(min(xy[, 2L]) - 3 * h[2L], max(xy[, 2L]) + 3 * h[2L],
                length.out = 401L)

  kde_cond <- function(x0) {
    f <- ks::kde(xy, H = diag(h^2), eval.points = cbind(x0, y_grid))
    grid_summary(y_grid, f$estimate)
  }

  set.seed(20260925)
  fit_k3 <- from_kde(xy, N = 3L, bandwidth = h, is_size = 20000L, seed = 1L)
  mc <- Mclust(xy, verbose = FALSE)
  g_mc <- as_gmm(
    mc$parameters$pro,
    lapply(seq_len(mc$G), function(j) mc$parameters$mean[, j]),
    lapply(seq_len(mc$G), function(j) mc$parameters$variance$sigma[, , j])
  )
  invisible(capture.output(mt <- mvnormalmixEM(xy)))
  g_mt <- as_gmm(mt$lambda, mt$mu, mt$sigma)

  cd <- npcdens(npcdensbw(xdat = xy[, 1L], ydat = xy[, 2L]))
  np_cond <- function(x0) {
    f <- predict(cd, newdata = data.frame(xdat = x0, ydat = y_grid))
    grid_summary(y_grid, f)
  }

  bk <- KernSmooth::bkde2D(xy, bandwidth = h, gridsize = c(401L, 401L),
                           range.x = list(range(x_eval), range(y_grid)))
  bk_cond <- function(x0) {
    i <- findInterval(x0, bk$x1, all.inside = TRUE)
    w <- (x0 - bk$x1[i]) / (bk$x1[i + 1L] - bk$x1[i])
    grid_summary(bk$x2, (1 - w) * bk$fhat[i, ] + w * bk$fhat[i + 1L, ])
  }

  query_list <- list(
    "proxymix, K = 3" = function(x0) gmm_summary(fit_k3, x0),
    "mclust, G by BIC" = function(x0) gmm_summary(g_mc, x0),
    "mixtools, k = 2" = function(x0) gmm_summary(g_mt, x0),
    "np, npcdens" = np_cond,
    "KernSmooth, bkde2D" = bk_cond,
    "ks, kde" = kde_cond
  )
  time_mat <- matrix(NA_real_, nrow = n_time, ncol = length(query_list),
                     dimnames = list(NULL, names(query_list)))
  for (i1 in seq_len(n_time)) {
    for (s2 in names(query_list)) {
      time_mat[i1, s2] <- system.time(
        for (i2 in seq_len(n_inner)) vapply(x_eval, query_list[[s2]], numeric(2L))
      )[["elapsed"]]
    } # end of for (s2 in names(query_list))
  } # end of for (i1 in seq_len(n_time))
  1000 * apply(time_mat, 2L, median) / (n_inner * length(x_eval))
}

faithful_xy <- as.matrix(faithful[, c("eruptions", "waiting")])
penguins_xy <- as.matrix(na.omit(penguins[, c("bill_len", "flipper_len")]))

n_time <- 5L
n_inner <- 10L
query_mat <- rbind(faithful = query_times(faithful_xy, n_time, n_inner),
                   penguins = query_times(penguins_xy, n_time, n_inner))
print(query_mat)

# milliseconds per query, averaged over the two datasets
query_ms <- colMeans(query_mat)
print(round(query_ms, 3))
query_versions <- c(proxymix = as.character(packageVersion("proxymix")),
                    ks = as.character(packageVersion("ks")),
                    np = as.character(packageVersion("np")),
                    mclust = as.character(packageVersion("mclust")),
                    mixtools = as.character(packageVersion("mixtools")),
                    KernSmooth = as.character(packageVersion("KernSmooth")))


time_section <- proc.time() - ptm_local
print("Section: 3. Time conditional queries per package")
print(time_section)


# =============================================================================
# 4. Write the summary
# =============================================================================
print("==== 4. Write the summary ====")
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
                paired = tab_paired,
                mclust_G = mclust_G,
                mixtools_secs = c(median = median(tab_res$secs[tab_res$method == "mixtools, k = 3"]),
                                  max = max(tab_res$secs[tab_res$method == "mixtools, k = 3"])),
                n = sim$n,
                n_rep = sim$n_rep,
                grid_side = sim$grid_side,
                proxymix_version = sim$versions[["proxymix"]],
                versions = desc_version(sim$versions),
                run_date = sim$run_date,
                elapsed_secs = sim$elapsed_secs,
                query_ms = query_ms,
                query_mat = query_mat,
                query_n = c(faithful = nrow(faithful_xy), penguins = nrow(penguins_xy)),
                query_versions = desc_version(query_versions),
                query_platform = R.version$platform)

dir.create(dir_out, showWarnings = FALSE)
file_name_out <- paste0(dir_out, "from_kde.rds")
saveRDS(res_out, file_name_out)

print(paste0("written: ", file_name_out, " (", file.size(file_name_out), " bytes)"))
print(paste0("version stamp: ", res_out$proxymix_version))


time_section <- proc.time() - ptm_local
print("Section: 4. Write the summary")
print(time_section)


# =============================================================================
# Completion
# =============================================================================
exec_time1 <- proc.time() - ptm1
print(exec_time1)
