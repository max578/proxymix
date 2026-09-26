# =============================================================================
# Summary of the stored entropy benchmark for the short entropy vignette
#
# Author      : Max Moldovan
# Version date: 26 Sep 2026
# Usage       : Rscript data-raw/vignette_results/entropy.R
#               (run from the package root)
#
# Data        : vignettes/articles/extended/entropy_sim.rds, written by
#               data-raw/articles/entropy_sim.R
# Requirements: proxymix, FNN (sections 2 and 3), mclust, pcalg (section 3)
#
# Outputs     : vignettes/results/entropy.rds, read by vignettes/entropy.Rmd
#
# Structure   : 1. Read the stored simulation
#               2. Check FNN against the stored estimates and pair the errors
#               3. Time one dataset per package
#               3b. Graph recovery of the energy field over seeds
#               4. Write the summary
# =============================================================================


rm(list = ls())

ptm1 <- proc.time()	# total execution time


# Configuration: paths
########################################################

file_name_sim <- paste0(getwd(), "/vignettes/articles/extended/entropy_sim.rds")
dir_out <- paste0(getwd(), "/vignettes/results/")
########################################################


# =============================================================================
# 1. Read the stored simulation
# =============================================================================
print("==== 1. Read the stored simulation ====")
ptm_local <- proc.time()

sim <- readRDS(file_name_sim)

# the version stamp is the proxymix version the simulation ran under
stopifnot(all(c("entropy_tab", "k_tab", "graph_tab", "truth", "est", "k_sel", "code",
                "n", "n_rep", "k_nn", "versions", "run_date",
                "elapsed_secs") %in% names(sim)))
stopifnot("proxymix" %in% names(sim$versions))
print(sim$versions)

# the designs, the blocks and the sampler, taken from the stored code
env_sim <- new.env()
for (e in parse(text = sim$code)) {
  if (is.call(e) && identical(e[[1L]], as.name("<-")) && is.name(e[[2L]]) &&
      as.character(e[[2L]]) %in% c("mixture", "rmix", "designs", "blocks")) {
    eval(e, env_sim)
  }
} # end of for (e in parse(text = sim$code))


time_section <- proc.time() - ptm_local
print("Section: 1. Read the stored simulation")
print(time_section)


# =============================================================================
# 2. Check FNN against the stored estimates and pair the errors
# =============================================================================
print("==== 2. Check FNN against the stored estimates and pair the errors ====")
ptm_local <- proc.time()

require(FNN)

# the simulation's mutual-information estimator is written out in R;
# FNN::mutinfo() on the same datasets must give the same values
s1 <- length(env_sim$designs)
gap_mat <- matrix(NA_real_, nrow = sim$n_rep * s1, ncol = 2L,
                  dimnames = list(NULL, c("entropy", "mutinfo")))
i_row <- 0L
for (i1 in seq_len(s1)) {

  m_local <- env_sim$designs[[i1]]
  a_vec <- env_sim$blocks[[i1]][[1L]]
  b_vec <- env_sim$blocks[[i1]][[2L]]

  for (r in seq_len(sim$n_rep)) {

    set.seed(r)
    x <- env_sim$rmix(sim$n, m_local)
    tab_local <- sim$est[sim$est$design == names(env_sim$designs)[i1] &
                           sim$est$rep == r, ]
    i_row <- i_row + 1L
    gap_mat[i_row, "entropy"] <- abs(entropy(x, k = sim$k_nn)[sim$k_nn] -
                                       tab_local$estimate[tab_local$method == "FNN"])
    gap_mat[i_row, "mutinfo"] <- abs(mutinfo(x[, a_vec, drop = FALSE],
                                             x[, b_vec, drop = FALSE],
                                             k = sim$k_nn) -
                                       tab_local$estimate[tab_local$method == "KSG"])

  } # end of for (r in seq_len(sim$n_rep))

} # end of for (i1 in seq_len(s1))

fnn_gap_vec <- apply(gap_mat, 2L, max)
print(fnn_gap_vec)
stopifnot(all(fnn_gap_vec < 1e-8))

# proxymix minus FNN absolute error, paired over the datasets
pair_tab <- expand.grid(quantity = c("Shannon entropy", "Shannon mutual information"),
                        design = names(env_sim$designs), stringsAsFactors = FALSE)
pair_tab$rival <- ifelse(pair_tab$quantity == "Shannon entropy", "FNN", "KSG")
pair_tab$diff <- NA_real_
pair_tab$diff_se <- NA_real_
for (i1 in seq_len(nrow(pair_tab))) {

  s_vec <- sim$est$design == pair_tab$design[i1] &
    sim$est$quantity == pair_tab$quantity[i1]
  a_tab <- sim$est[s_vec & sim$est$method == "proxymix", ]
  b_tab <- sim$est[s_vec & sim$est$method == pair_tab$rival[i1], ]
  d_vec <- abs(a_tab$error[order(a_tab$rep)]) - abs(b_tab$error[order(b_tab$rep)])
  pair_tab$diff[i1] <- mean(d_vec)
  pair_tab$diff_se[i1] <- sd(d_vec) / sqrt(length(d_vec))

} # end of for (i1 in seq_len(nrow(pair_tab)))
print(pair_tab)


time_section <- proc.time() - ptm_local
print("Section: 2. Check FNN against the stored estimates and pair the errors")
print(time_section)


# =============================================================================
# 3. Time one dataset per package
# =============================================================================
print("==== 3. Time one dataset per package ====")
ptm_local <- proc.time()

require(proxymix)
require(mclust)
require(pcalg)

# one dataset from the 2-d design, as in the vignette's competitor chunk
w <- c(0.4, 0.35, 0.25)
mu <- list(c(-2, 0), c(2, 1), c(0, 3))
sig <- list(matrix(c(1, 0.5, 0.5, 1), 2L),
            matrix(c(1.5, -0.6, -0.6, 0.8), 2L),
            diag(c(0.5, 1.2)))
set.seed(1L)
z <- sample.int(3L, 500L, replace = TRUE, prob = w)
x <- matrix(rnorm(500L * 2L), 500L, 2L)
for (k in 1:3) {
  x[z == k, ] <- x[z == k, , drop = FALSE] %*% chol(sig[[k]]) +
    matrix(mu[[k]], sum(z == k), 2L, byrow = TRUE)
} # end of for (k in 1:3)

# one dataset from the graph design, as in the vignette's competitor chunk
omega <- diag(6L)
omega[abs(row(omega) - col(omega)) == 1L] <- -0.4
set.seed(1L)
z_g <- sample.int(2L, 500L, replace = TRUE, prob = c(0.5, 0.5))
x_g <- matrix(rnorm(500L * 6L), 500L, 6L) %*% chol(solve(omega))
x_g[z_g == 2L, 1L] <- x_g[z_g == 2L, 1L] + 4

# the same draws through the simulation's own sampler
m_g <- env_sim$mixture(w = c(0.5, 0.5), mu = list(rep(0, 6L), c(4, rep(0, 5L))),
                       sig = list(solve(omega), solve(omega)))
set.seed(1L)
stopifnot(isTRUE(all.equal(x_g, env_sim$rmix(500L, m_g))))
set.seed(1L)
stopifnot(isTRUE(all.equal(x, env_sim$rmix(500L, env_sim$designs[[1L]]))))

h_mc <- function(g) gmm_entropy(g, order = "shannon", n_mc = 5000L, seed = 1L)$mc

run_list <- list(
  entropy_proxymix = function() {
    fit <- fit_em_samples(gmm_target_from_samples(x), N = 3L, seed = 1L)
    c(h_mc(fit),
      h_mc(gmm_marginalise(fit, 1L)) + h_mc(gmm_marginalise(fit, 2L)) - h_mc(fit))
  },
  entropy_FNN = function() {
    c(entropy(x, k = 10L)[10L],
      mutinfo(x[, 1L, drop = FALSE], x[, 2L, drop = FALSE], k = 10L))
  },
  count_proxymix = function() {
    tgt <- gmm_target_from_samples(x)
    lapply(1:5, function(k) bic_aic(fit_em_samples(tgt, N = k, seed = 1L)))
  },
  count_mclust = function() mclustICL(x, G = 1:5, verbose = FALSE),
  graph_proxymix = function() {
    sel <- select_N(gmm_target_from_samples(x_g), candidates = 1:4, seed = 1L)
    gmm_independence_graph(sel$best_fit)
  },
  graph_pcalg = function() {
    pc(suffStat = list(C = cor(x_g), n = nrow(x_g)),
       indepTest = gaussCItest, alpha = 0.01, p = ncol(x_g))
  }
)

n_time <- 5L
time_mat <- matrix(NA_real_, nrow = n_time, ncol = length(run_list),
                   dimnames = list(NULL, names(run_list)))
for (i1 in seq_len(n_time)) {
  for (s2 in names(run_list)) {
    time_mat[i1, s2] <- system.time(run_list[[s2]]())[["elapsed"]]
  } # end of for (s2 in names(run_list))
} # end of for (i1 in seq_len(n_time))

# median seconds for one 500-row dataset
time_vec <- apply(time_mat, 2L, median)
print(time_vec)
time_versions <- c(proxymix = as.character(packageVersion("proxymix")),
                   FNN = as.character(packageVersion("FNN")),
                   mclust = as.character(packageVersion("mclust")),
                   pcalg = as.character(packageVersion("pcalg")))


time_section <- proc.time() - ptm_local
print("Section: 3. Time one dataset per package")
print(time_section)


# =============================================================================
# 3b. Graph recovery of the energy field over seeds
# =============================================================================
print("==== 3b. Graph recovery of the energy field over seeds ====")
ptm_local <- proc.time()

# the vignette's regime (iii) example, refitted over ten seeds at its own
# trial-point count and at twice that count
energy <- function(x_mat) {
  x_mat <- matrix(x_mat, ncol = 3L)
  rowSums((x_mat^2 - 1)^2) -
    0.7 * (x_mat[, 1L] * x_mat[, 2L] + x_mat[, 2L] * x_mat[, 3L])
}
field <- gmm_target(n_dim = 3L, log_density = function(x_mat) -energy(x_mat))

field_tab <- expand.grid(seed = 1:10, is_size = c(6000L, 12000L))
field_tab$chain <- NA
field_tab$pcor13 <- NA_real_
field_tab$converged <- NA
field_tab$degenerate <- NA
field_tab$kld_approx <- NA_real_
field_tab$flagged <- NA
for (i1 in seq_len(nrow(field_tab))) {

  fit_local <- fit_kld_em(field, N = 8L, proposal = proposal_uniform(3L, -3, 3),
                          is_size = field_tab$is_size[i1], anneal = TRUE,
                          seed = field_tab$seed[i1], support_warn = FALSE)
  adj_local <- suppressMessages(gmm_independence_graph(fit_local))
  q_local <- gmm_fit_quality(fit_local)
  field_tab$chain[i1] <- sum(adj_local) == 4L && adj_local[1L, 3L] == 0L
  field_tab$pcor13[i1] <- attr(adj_local, "pcor")[1L, 3L]
  field_tab$converged[i1] <- isTRUE(q_local$converged)
  field_tab$degenerate[i1] <- isTRUE(q_local$degenerate)
  field_tab$kld_approx[i1] <- q_local$kld_approx
  ## the package's flag: not converged, degenerate, or approximate KL above 0.3
  field_tab$flagged[i1] <- isTRUE(q_local$degenerate) || isFALSE(q_local$converged) ||
    q_local$kld_approx > 0.3

} # end of for (i1 in seq_len(nrow(field_tab)))
print(field_tab)
print(aggregate(cbind(chain, flagged) ~ is_size, data = field_tab, FUN = sum))


time_section <- proc.time() - ptm_local
print("Section: 3b. Graph recovery of the energy field over seeds")
print(time_section)


# =============================================================================
# 4. Write the summary
# =============================================================================
print("==== 4. Write the summary ====")
ptm_local <- proc.time()

# share of datasets on which each selector chose two components
k_two_tab <- aggregate(chosen ~ design + selector, data = sim$k_sel,
                       FUN = function(v) mean(v == 2L))
names(k_two_tab)[3L] <- "share_two"
print(k_two_tab)

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

res_out <- list(entropy_tab = sim$entropy_tab,
                truth = sim$truth,
                pair_tab = pair_tab,
                k_tab = sim$k_tab,
                k_two_tab = k_two_tab,
                graph_tab = sim$graph_tab,
                field_tab = field_tab,
                fnn_gap = fnn_gap_vec,
                n = sim$n,
                n_rep = sim$n_rep,
                k_nn = sim$k_nn,
                proxymix_version = sim$versions[["proxymix"]],
                versions = desc_version(sim$versions),
                run_date = sim$run_date,
                elapsed_secs = sim$elapsed_secs,
                time_secs = time_vec,
                time_versions = desc_version(time_versions),
                time_platform = R.version$platform)

dir.create(dir_out, showWarnings = FALSE)
file_name_out <- paste0(dir_out, "entropy.rds")
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
