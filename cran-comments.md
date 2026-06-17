# proxymix 0.6.0 — submission notes

First CRAN submission. The package implements the three-regime
Gaussian-mixture proxy hierarchy of Hoek and Elliott (2024), including
the importance-sampled KLD-EM regime for targets that can be evaluated but
not sampled, a closed-form affine-Gaussian operator calculus, support-aware
importance proposals for bounded and one-sided targets, and (new in 0.6.0) a
closed-form decision layer (heterogeneous treatment effects, next-best-action,
off-line policy value and an identification audit) read off a single
Gaussian-mixture fit, with the `do`-operator and the counterfactual as
first-class operators.

## Test environments

* Local: macOS 26.4 (aarch64-apple-darwin20), R 4.5.2, Apple Accelerate BLAS.
* Planned before submission:
  - `devtools::check_win_devel()` on the tarball;
  - `rhub::rhub_check()` v2 across Linux, macOS, Windows, r-release and r-devel;
  - `urlchecker::url_check()` (see the URL note below).

## R CMD check results

`R CMD check --as-cran proxymix_0.6.0.tar.gz` on the local environment
returns **0 ERRORs, 0 WARNINGs, 2 NOTEs**, both environmental:

1. **`checking for future file timestamps ... NOTE`** — "unable to verify
   current time". Local clock-sync / offline environmental issue; does not
   appear on CRAN's checkers.

2. **`checking HTML version of manual ... NOTE`** — local HTML Tidy is not
   recent enough. Local tooling, not a package-code finding; does not appear
   on CRAN's checkers.

On the actual first CRAN submission, the **`checking CRAN incoming
feasibility ... NOTE` ("New submission")** is also expected; it resolves on
subsequent releases. The macOS `.DS_Store` NOTE seen on earlier local runs
is a Finder artefact recreated on the `*.Rcheck/` directory during the check;
the tarball itself is clean (`.DS_Store` is in `.Rbuildignore`).

There are **no NOTEs concerning the code, examples, vignettes, or
documentation** of the package itself.

## URL / DOI note

`urlchecker::url_check()` reports HTTP **403 Forbidden** on the anchor-paper
DOI `https://doi.org/10.1080/07362994.2024.2372605` (cited in `README.md`
and four vignettes). This is a false positive: the DOI resolves correctly
(HTTP 302 to `https://www.tandfonline.com/doi/full/10.1080/07362994.2024.2372605`),
and the 403 is Taylor & Francis bot-protection that blocks every automated
agent (including a browser user-agent), not a broken link. The DOI is open
access (CC-BY) and resolves normally in an interactive browser. `R CMD check
--as-cran` does not flag it: `DESCRIPTION` cites the paper with the
`<doi:10.1080/07362994.2024.2372605>` token form, which CRAN does not fetch.

`URL:` and `BugReports:` (`https://github.com/max578/proxymix`) resolve.

## Reverse dependencies

None — first release.

## Code coverage

`covr::package_coverage()` reports **91.7%** overall; the new decision-layer
sources (`R/intervene.R`, `R/uplift.R`, `R/decide.R`, `R/identification.R`,
`R/uplift_internal.R`) are each covered at **>= 89%**. The largest remaining
uncovered surface is:

* `R/zzz.R` (0%) — package machinery (`.onLoad`, S7 method registration and
  the conditional `ggplot2` `autoplot()` registration); it runs at load time
  and is not creditable by `covr`.
* `R/init.R` (~82%) and `R/from_fb_posterior.R` (~84%) — multi-start and
  seam-dispatch branches, exercised indirectly through the fitters.

Every numerically load-bearing path (`R/fit_kld_em.R`, `R/fit_em_samples.R`,
`R/gmm_ops.R`, `R/operator_calculus.R`, `R/from_kde.R`, `R/proposals.R`) and
the S7 class / target-constructor / diagnostic surfaces are covered at
**>= 90%** (most at >= 96%).

## Reproducibility

Pinned validation scripts at `inst/validation/`:

* `regime_iii_pinned_fits.R` — reference KLD-EM fits on the built-in
  non-Gaussian targets, now including a compact-support Epanechnikov fit
  under the auto-selected support-matched proposal (MC-SE-aware tolerances;
  asserts no NaN weights).
* `from_kde_pinned_fits.R` — reference KDE-to-GMM compressions.
* `operator_calculus_pinned.R` — affine-Gaussian operator algebra reference
  outputs (pushforward, observe, aggregate, missing).

All pass on the test environment with the recorded seeds.

## Notes for the maintainer

* Vignettes use `\VignetteEngine{knitr::rmarkdown}` on `.Rmd` files
  (Quarto `.qmd` is not first-class under `R CMD check --as-cran`).
* Vignette filenames follow R name rules (start with a letter, not a digit):
  `quickstart.Rmd`, `three_regimes.Rmd`, `density_shapes.Rmd`, `from_kde.Rmd`,
  `operator_calculus.Rmd`, `uplift.Rmd`, `roadmap.Rmd`.
* `inst/benchmarks/` ships a reproducible benchmark script and its results;
  it is not run at check time and depends only on Suggests/external packages
  (`grf`, `DoubleML`) that are not package dependencies.
* `Encoding: UTF-8`, `Language: en-AU`. `inst/WORDLIST` holds the scientific
  vocabulary so `devtools::spell_check()` passes cleanly.
