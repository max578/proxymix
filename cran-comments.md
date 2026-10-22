## Update: proxymix 0.17.0

This version adds a fit to a distribution known only through its
characteristic function, with bounds on the Kolmogorov distance, changes
the default bandwidth of `from_kde()` to a plug-in rule, and corrects four
defects (see NEWS.md).

## Test environments

* local: macOS 26.7 (Apple silicon), R 4.6.1
* win-builder: Windows, R-devel (2026-10-21 r90733 ucrt)
* R-hub: Ubuntu 24.04, R-devel (2026-10-20 r90728)
* R-hub: Ubuntu 26.04 (ubuntu-next), R 4.6.2 RC (2026-10-21 r90732)

## R CMD check results

0 errors | 0 warnings | 0 notes

## Comments

* Examples run in about 5 s in total. The slowest single example takes
  about 1.3 s.
* The thirteen vignettes build in about 32 s. Longer versions that compare
  the package with other methods are published only on the package website
  and are excluded from the build.
* The tests, the characteristic-function examples and their vignette use
  one BLAS thread (through RhpcBLASctl, now suggested), and the tests run on
  at most two workers. On the two R-hub Linux checks the tests used CPU time
  about 2.0 times elapsed time.
* No package on CRAN depends on, imports or suggests proxymix.
