## Submission: proxymix 0.16.0 (first submission)

An earlier attempt to submit version 0.15.1 was not completed, so this is
the first submission of the package to CRAN.

## Test environments

* local: macOS 26.6 (Apple silicon), R 4.6.1
* win-builder: Windows, R 4.6.1 (release)
* win-builder: Windows, R-devel (2026-09-25 r90590)

## R CMD check results

0 errors | 0 warnings | 1 note

* This is a new submission.
* The note also lists Hoek, Kullback, Leibler and KLD as possibly
  misspelled. Hoek, Kullback and Leibler are surnames, and KLD is the
  usual abbreviation for Kullback-Leibler divergence.

## Comments

* The package fits Gaussian mixtures that approximate a target density in
  Kullback-Leibler divergence, including targets that can be evaluated but
  not sampled, following van der Hoek and Elliott (2024)
  <doi:10.1080/07362994.2024.2372605>.
* Examples run in about 20 s in total. The slowest single example takes
  about 4 s.
* The twelve vignettes build in about 60 s. Longer versions that compare
  the package with other methods are published only on the package website
  and are excluded from the build.
