## Resubmission: proxymix 0.16.0 (first submission)

An earlier attempt to submit version 0.15.1 was not completed, so this is
the first submission of the package to CRAN.

This resubmission fixes the two NOTEs raised by the incoming pretest on
Debian (r-devel):

* Rd files without \usage (autoplot.gmm_fit, glance.gmm_fit, tidy.gmm):
  these document methods registered on generics from suggested packages.
  Their argument descriptions now sit in a section of their own, so the
  pages no longer have an \arguments section without \usage.
* The example of proxy_mnar_sensitivity() took 6.3 s. It now uses a
  smaller data set and runs in under 1 s locally.

It also fixes a defect found after the pretest: the small ridge added to
fitted covariances is now scaled per coordinate, so coordinates on a small
scale are no longer over-smoothed.

The second pretest reported an overall check time of 11 min on Windows
(r-devel). Two slow tests now skip on CRAN and two vignettes compute less,
and the check time on win-builder (R-devel, 2026-09-29) is now 7.1 min.

## Test environments

* local: macOS 26.6 (Apple silicon), R 4.6.1
* win-builder: Windows, R 4.6.1 (release)
* win-builder: Windows, R-devel (2026-09-25 r90590)
* R-hub: Ubuntu 24.04, R-devel (2026-09-25 r90590)
* R-hub: Ubuntu 24.04 (ubuntu-next), R 4.6.1 patched (2026-09-17 r90590)

## R CMD check results

0 errors | 0 warnings | 1 note

* This is a new submission.
* The note may also list Hoek, der, Kullback and Leibler as possibly
  misspelled. They come from the surnames van der Hoek, Kullback and
  Leibler.

## Comments

* The package fits Gaussian mixtures that approximate a target density in
  Kullback-Leibler divergence, including targets that can be evaluated but
  not sampled, following van der Hoek and Elliott (2024)
  <doi:10.1080/07362994.2024.2372605>.
* Examples run in about 4 s in total. The slowest single example takes
  about 1.3 s.
* The twelve vignettes build in about 40 s. Longer versions that compare
  the package with other methods are published only on the package website
  and are excluded from the build.
