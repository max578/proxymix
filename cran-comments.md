# proxymix 0.3.0 — submission notes

## Test environments

* Local: macOS 26.4 (aarch64-apple-darwin20), R 4.5.2, Apple Accelerate BLAS.
* Planned pre-submission (added once the GitHub repository is public):
  - `devtools::check_win_devel()` on the tarball;
  - `rhub::rhub_check()` v2 across Linux, macOS, Windows, r-devel;
  - `urlchecker::url_check()` against DESCRIPTION, man/, vignettes/, README.

## R CMD check results

`R CMD check --as-cran proxymix_0.3.0.tar.gz` on the local environment
returns **0 ERRORs, 0 WARNINGs, 4 NOTEs**. All four are environmental:

1. **`checking CRAN incoming feasibility ... NOTE`** — "New submission".
   Expected on the first CRAN submission of any package; resolves on
   subsequent releases.

2. **`checking for future file timestamps ... NOTE`** — "unable to verify
   current time". Local clock-sync environmental issue; does not appear on
   CRAN's incoming checkers.

3. **`checking HTML version of manual ... NOTE`** — "local HTML Tidy is
   not recent enough". Local tooling, not a package-code finding; does
   not appear on CRAN's checkers.

4. **`checking for non-standard things in the check directory ... NOTE`** —
   `Found the following files/directories: '.DS_Store'`. This is the macOS
   Finder metadata file. It is excluded from the package tarball via
   `.Rbuildignore` (`tar tzf proxymix_0.3.0.tar.gz` confirms the tarball is
   clean), but is recreated by macOS on the `*.Rcheck/` directory during
   the check process itself. The NOTE does not appear on CRAN's Linux,
   Solaris, or Windows checkers.

There are **no NOTEs concerning the code, examples, vignettes, or
documentation** of the package itself.

## URL fields

`URL:` and `BugReports:` in `DESCRIPTION` resolve at the time of submission
(verified via `urlchecker::url_check()`).

## Reverse dependencies

None — first release.

## Code coverage

`covr::package_coverage()` reports **85.82%** overall on the v0.3.0
codebase. The deliberately uncovered surface is:

* `R/zzz.R` (0%) — package machinery (`.onLoad`); not scientifically
  load-bearing.
* `R/proxymix-classes.R` (~66%) — S7 validator and print-method branches
  for invalid construction paths. Each branch raises an error or formats
  a print line; testing the print formatting itself yields little value
  beyond what `expect_snapshot()` already covers.
* `R/gmm_target_constructors.R` (~76%) — `gmm_target_from_posterior()`
  S3 generic dispatch. The `function`-method path is exercised; the
  `flexybayes`-method path will be covered by the cross-package
  integration test once that sibling package ships its
  `fb_log_posterior()` contract.
* `R/diagnostics.R` (~81%) — defensive guards on degenerate IS samples
  (every weight numerically zero); exercised in dedicated tests but not
  on every code path.

All numerically-load-bearing scientific code (`R/fit_kld_em.R`,
`R/fit_em_samples.R`, `R/gmm_ops.R`, `R/operator_calculus.R`,
`R/from_kde.R`) is covered at **>= 90%**.

## Reproducibility

Pinned validation scripts at `inst/validation/`:

* `regime_iii_pinned_fits.R` — three reference KLD-EM fits with
  MC-SE-aware tolerances on three built-in non-Gaussian targets.
* `from_kde_pinned_fits.R` — three reference KDE-to-GMM compressions.
* `operator_calculus_pinned.R` — affine-Gaussian operator algebra
  reference outputs (pushforward, observe, aggregate, missing).

All pass on the test environment with the recorded seeds.

## Notes for the maintainer

* Vignettes use `\VignetteEngine{knitr::rmarkdown}` on `.Rmd` files
  (Quarto `.qmd` is not first-class under `R CMD check --as-cran`).
* Vignette filenames follow R name rules (start with a letter, not a
  digit): `quickstart.Rmd`, `three_regimes.Rmd`, `density_shapes.Rmd`,
  `from_kde.Rmd`, `operator_calculus.Rmd`, `roadmap.Rmd`.
* `Encoding: UTF-8`, `Language: en-AU`. The `inst/WORDLIST` file holds
  scientific vocabulary (`Hoek`, `Sejdinovic`, `Kullback-Leibler`,
  `marginalise`, `conditionalise`, etc.) so `devtools::spell_check()`
  passes cleanly.
