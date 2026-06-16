## inst/validation/regime_iii_pinned_fits.R
##
## proxymix v0.1.1 - regime (iii) numerical validation corpus (starter).
##
## Run with:
##   Rscript inst/validation/regime_iii_pinned_fits.R
## or interactively:
##   source(system.file("validation", "regime_iii_pinned_fits.R",
##                      package = "proxymix"))
##
## This script fits the three built-in 2-D targets by regime (iii) with
## pinned seeds, captures the headline numerical diagnostics, and asserts
## that each one lies within an expected range. The expected ranges were
## chosen with deliberate slack (MC SE-aware) and are intended to fail
## loudly on a future regression rather than to pin a single number.
##
## The script exits with status 0 on success, 1 on failure.

suppressPackageStartupMessages({
  library(proxymix)
})

writeLines("=== proxymix validation: regime (iii) pinned fits ===")
writeLines(sprintf("R version : %s", R.version.string))
writeLines(sprintf("proxymix  : %s", as.character(utils::packageVersion("proxymix"))))
writeLines(sprintf("Platform  : %s", R.version$platform))
writeLines(sprintf("BLAS      : %s", utils::sessionInfo()$BLAS %||% "default"))

`%||%` <- function(a, b) if (is.null(a)) b else a

stamp <- function(name, fit, ok = TRUE) {
  d <- fit@diagnostics
  cat(sprintf(
    "%-14s | kld = %+7.4f (se %5.4f) | ESS = %5.1f / %4d (%.0f%%) | maxW = %.4f | valKLD = %s | %s\n",
    name,
    d$kld_final, d$mc_se_kld,
    d$ess, d$is_size, 100 * d$ess_relative,
    d$max_weight,
    if (is.na(d$validation_kld)) "    NA" else sprintf("%+7.4f", d$validation_kld),
    if (ok) "OK" else "FAIL"
  ))
}

assert_in_range <- function(label, value, lo, hi) {
  if (!is.finite(value) || value < lo || value > hi) {
    stop(sprintf("VALIDATION FAILURE: %s = %s outside [%s, %s]",
                 label, format(value), format(lo), format(hi)),
         call. = FALSE)
  }
  invisible(TRUE)
}

# -------------------------------------------------------------------------
# Banana target
# -------------------------------------------------------------------------
banana <- banana_target()
fit_b <- fit_proxymix(
  banana, N = 4L, regime = "kld",
  proposal = is_mvt(n_dim = 2L, sigma = 4 * diag(2), df = 5),
  is_size = 3000L, max_iter = 50L, seed = 1L,
  validation_size = 3000L, validation_seed = 2L
)
stamp("banana(N=4)", fit_b)
assert_in_range("banana.kld_final",  fit_b@diagnostics$kld_final, -0.05, 0.20)
assert_in_range("banana.ess_relative", fit_b@diagnostics$ess_relative, 0.20, 1.00)
assert_in_range("banana.validation_kld", fit_b@diagnostics$validation_kld,
                -0.05, 0.20)

# -------------------------------------------------------------------------
# Donut target
# -------------------------------------------------------------------------
donut <- donut_target()
fit_d <- fit_proxymix(
  donut, N = 6L, regime = "kld",
  proposal = is_mvt(n_dim = 2L, sigma = 9 * diag(2), df = 5),
  is_size = 3500L, max_iter = 60L, seed = 1L,
  validation_size = 3500L, validation_seed = 2L
)
stamp("donut(N=6)", fit_d)
assert_in_range("donut.kld_final",  fit_d@diagnostics$kld_final, -0.10, 1.50)
assert_in_range("donut.ess_relative", fit_d@diagnostics$ess_relative, 0.05, 1.00)
assert_in_range("donut.validation_kld", fit_d@diagnostics$validation_kld,
                -0.10, 1.50)

# -------------------------------------------------------------------------
# Three-mixture target
# -------------------------------------------------------------------------
mt <- mixture_target()
fit_m <- fit_proxymix(
  mt, N = 3L, regime = "kld",
  proposal = is_mvt(n_dim = 2L, sigma = 6 * diag(2), df = 5),
  is_size = 3000L, max_iter = 50L, seed = 1L,
  validation_size = 3000L, validation_seed = 2L
)
stamp("mixture(N=3)", fit_m)
assert_in_range("mixture.kld_final",  fit_m@diagnostics$kld_final, -0.05, 0.10)
assert_in_range("mixture.ess_relative", fit_m@diagnostics$ess_relative, 0.15, 1.00)
assert_in_range("mixture.validation_kld", fit_m@diagnostics$validation_kld,
                -0.05, 0.10)

# -------------------------------------------------------------------------
# Epanechnikov target (compact support; support-matched proposal, v0.4.1)
# -------------------------------------------------------------------------
# No explicit proposal: the declared compact support drives automatic
# is_uniform selection. Locks the no-NaN-weight guarantee (criterion S10).
epan <- epanechnikov_target(n_dim = 1L)
fit_e <- fit_proxymix(
  epan, N = 3L, regime = "kld",
  is_size = 4000L, max_iter = 60L, seed = 1L,
  validation_size = 4000L, validation_seed = 2L
)
stamp("epanechnikov(N=3)", fit_e)
assert_in_range("epanechnikov.kld_final",  fit_e@diagnostics$kld_final, -0.05, 1.00)
assert_in_range("epanechnikov.ess_relative", fit_e@diagnostics$ess_relative,
                0.30, 1.00)
assert_in_range("epanechnikov.support_fraction",
                fit_e@diagnostics$support_fraction, 0.999, 1.001)

# -------------------------------------------------------------------------
# All-targets summary
# -------------------------------------------------------------------------
writeLines("\nSummary table (all-targets):")
tbl <- rbind(
  data.frame(target = "banana",  N = 4L,
             kld = fit_b@diagnostics$kld_final,
             validation_kld = fit_b@diagnostics$validation_kld,
             ess_relative = fit_b@diagnostics$ess_relative,
             max_weight = fit_b@diagnostics$max_weight),
  data.frame(target = "donut",   N = 6L,
             kld = fit_d@diagnostics$kld_final,
             validation_kld = fit_d@diagnostics$validation_kld,
             ess_relative = fit_d@diagnostics$ess_relative,
             max_weight = fit_d@diagnostics$max_weight),
  data.frame(target = "mixture", N = 3L,
             kld = fit_m@diagnostics$kld_final,
             validation_kld = fit_m@diagnostics$validation_kld,
             ess_relative = fit_m@diagnostics$ess_relative,
             max_weight = fit_m@diagnostics$max_weight),
  data.frame(target = "epanechnikov", N = 3L,
             kld = fit_e@diagnostics$kld_final,
             validation_kld = fit_e@diagnostics$validation_kld,
             ess_relative = fit_e@diagnostics$ess_relative,
             max_weight = fit_e@diagnostics$max_weight)
)
print(tbl, row.names = FALSE, digits = 4L)

writeLines("\n=== Validation corpus PASSED ===")
invisible(NULL)
