## Update: proxymix 0.17.0

This version adds a fit to a distribution known only through its
characteristic function, with bounds on the Kolmogorov distance, changes
the default bandwidth of `from_kde()` to a plug-in rule, and corrects four
defects (see NEWS.md).

## Test environments

* local: macOS 26.7 (Apple silicon), R 4.6.1

## R CMD check results

0 errors | 0 warnings | 0 notes

## Comments

* Examples run in about 4 s in total. The slowest single example takes
  about 1 s.
* The thirteen vignettes build in about 32 s. Longer versions that compare
  the package with other methods are published only on the package website
  and are excluded from the build.
* No package on CRAN depends on, imports or suggests proxymix.
