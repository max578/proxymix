# `inst/validation/` — `proxymix` numerical validation corpus

This directory holds **runnable validation scripts** for the package's
methods, separate from the test suite. The audit recommendation
(2026-05-14) was for a small validation corpus with pinned reference fits
and slack-aware tolerances; this is the seed of that corpus.

Scripts live at the top level. Each script:

* loads `proxymix` from the installed library;
* prints the R, BLAS, and `proxymix` versions;
* fits a reference scenario with a pinned seed;
* asserts each headline number lies in a documented range;
* exits with a non-zero status on assertion failure.

## Running

```r
source(system.file("validation", "regime_iii_pinned_fits.R",
                   package = "proxymix"))
```

or, from a shell:

```sh
Rscript inst/validation/regime_iii_pinned_fits.R
```

## Current corpus

| Script | What it validates |
|---|---|
| `regime_iii_pinned_fits.R` | KLD-EM on banana / donut / 3-mixture with pinned seeds. Asserts that final KLD, validation KLD, relative ESS, and max weight all fall in pre-specified ranges. |

## Why not just unit tests?

Unit tests (`tests/testthat/`) protect API contracts and ensure new code
does not regress against old code. The validation corpus protects
**numerical claims**: that regime (iii) fits the donut to KLD ≲ 1.5 on
3500 IS draws with a `df = 5` Student-t proposal, that the banana fits
to KLD ≲ 0.2, and so on. Unit tests can move when the API moves; the
validation corpus moves only when the *science* moves.

The corpus is intentionally a small seed in v0.1.1. Subsequent releases
will grow it to cover more shapes, more dimensions, and a comparison
against `mclust` parity on regime (ii).
