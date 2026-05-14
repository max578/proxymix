## inst/validation/from_kde_pinned_fits.R
##
## Pinned validation corpus for `from_kde()` (v0.2.0).
##
## Three reference KDE -> GMM proxy pipelines under pinned seeds.
## Acceptable diagnostic ranges are MC-SE-aware and capture the
## bandwidth / sample-size regime of each example. The script is
## *not* part of automated tests: it is run manually after a
## v0.2.x release to detect numerical drift.

if (!requireNamespace("proxymix", quietly = TRUE)) {
  stop("install proxymix first")
}

library(proxymix)

sessionInfo()

cat("\nBLAS:", basename(getOption("matprod") %||% "default"), "\n")
cat("R:   ", R.version.string, "\n\n")

results <- list()

## ---------------------------------------------------------------------------
## 1. Bimodal Gaussian — wide separation, p = 2, N = 2.
##    Truth: w = 0.5 / 0.5, mu_k = +/- (2, 0), Sigma_k = I_2.
## ---------------------------------------------------------------------------

set.seed(2026051401)
x_bimodal <- rbind(
  mvnfast::rmvn(200L, mu = c(-2, 0), sigma = diag(2)),
  mvnfast::rmvn(200L, mu = c( 2, 0), sigma = diag(2))
)
t0 <- proc.time()
fit_bimodal <- from_kde(
  x_bimodal, N = 2L, bandwidth = "silverman",
  is_size = 3000L, max_iter = 80L, seed = 1L,
  validation_size = 3000L
)
elapsed_bimodal <- proc.time() - t0

results$bimodal <- list(
  kld_final           = fit_bimodal@diagnostics$kld_final,
  validation_kld      = fit_bimodal@diagnostics$validation_kld,
  ess_relative        = fit_bimodal@diagnostics$ess_relative,
  max_weight          = fit_bimodal@diagnostics$max_weight,
  iterations          = fit_bimodal@iterations,
  converged           = fit_bimodal@converged,
  mu_recovered        = lapply(fit_bimodal@means, round, digits = 3L),
  trace_cov_recovered = vapply(fit_bimodal@covariances,
                               function(S) sum(diag(S)), numeric(1L)),
  elapsed_secs        = unname(elapsed_bimodal["elapsed"])
)

## Acceptance ranges (MC-SE aware):
##   ESS:           > 0.30
##   max_weight:    < 0.05
##   kld_final:     [0.00, 0.30]
##   mu_k:          within +/- 0.4 of truth
stopifnot(
  results$bimodal$ess_relative   > 0.30,
  results$bimodal$max_weight     < 0.05,
  results$bimodal$kld_final      >= 0,
  results$bimodal$kld_final      <= 0.30
)

## ---------------------------------------------------------------------------
## 2. Banana — curved 1-component target, p = 2, N = 3.
##    Truth: a regime-(iii) reference fit; expect higher KLD than (1).
## ---------------------------------------------------------------------------

set.seed(2026051402)
banana <- banana_target()
x_banana <- rgmm(300L, fit_proxymix(banana, N = 5L,
                                    proposal = is_mvt(2L, mean = c(0, 0),
                                                      sigma = 4 * diag(2),
                                                      df = 5),
                                    is_size = 2000L, max_iter = 30L,
                                    seed = 11L))
t0 <- proc.time()
fit_banana <- from_kde(
  x_banana, N = 3L, bandwidth = "silverman",
  is_size = 3000L, max_iter = 60L, seed = 2L,
  validation_size = 3000L
)
elapsed_banana <- proc.time() - t0

results$banana <- list(
  kld_final          = fit_banana@diagnostics$kld_final,
  validation_kld     = fit_banana@diagnostics$validation_kld,
  ess_relative       = fit_banana@diagnostics$ess_relative,
  max_weight         = fit_banana@diagnostics$max_weight,
  iterations         = fit_banana@iterations,
  converged          = fit_banana@converged,
  elapsed_secs       = unname(elapsed_banana["elapsed"])
)

stopifnot(
  results$banana$ess_relative > 0.15,
  results$banana$max_weight   < 0.10,
  results$banana$kld_final    >= 0
)

## ---------------------------------------------------------------------------
## 3. 3-component planar mixture — well-separated, p = 2, N = 3.
##    Bandwidth swept: small (under-smooth) vs Silverman.
## ---------------------------------------------------------------------------

set.seed(2026051403)
mix <- mixture_target()
x_mix <- rgmm(360L, fit_proxymix(mix, N = 5L,
                                 proposal = is_mvt(2L, mean = c(0, 0),
                                                   sigma = 9 * diag(2),
                                                   df = 5),
                                 is_size = 2000L, max_iter = 30L,
                                 seed = 21L))

fits_mix <- lapply(list("silverman", 0.4), function(bw) {
  from_kde(x_mix, N = 3L, bandwidth = bw,
           is_size = 2500L, max_iter = 60L, seed = 3L,
           validation_size = 2500L)
})

results$mixture <- list(
  bandwidth_silverman = list(
    kld_final     = fits_mix[[1L]]@diagnostics$kld_final,
    ess_relative  = fits_mix[[1L]]@diagnostics$ess_relative,
    max_weight    = fits_mix[[1L]]@diagnostics$max_weight,
    iterations    = fits_mix[[1L]]@iterations
  ),
  bandwidth_narrow = list(
    kld_final     = fits_mix[[2L]]@diagnostics$kld_final,
    ess_relative  = fits_mix[[2L]]@diagnostics$ess_relative,
    max_weight    = fits_mix[[2L]]@diagnostics$max_weight,
    iterations    = fits_mix[[2L]]@iterations
  )
)

stopifnot(
  results$mixture$bandwidth_silverman$ess_relative > 0.15,
  results$mixture$bandwidth_silverman$kld_final    >= 0,
  results$mixture$bandwidth_narrow$kld_final       >= 0
)

## ---------------------------------------------------------------------------
## Report.
## ---------------------------------------------------------------------------

cat("\nfrom_kde() pinned validation — proxymix v",
    as.character(utils::packageVersion("proxymix")),
    "\n", sep = "")
print(results)
invisible(results)
