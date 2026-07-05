## Submission: proxymix 0.15.1 (first submission)

## Test environments

* local: macOS (Apple silicon), R release
* win-builder (devel): submitted alongside

## R CMD check results

0 errors | 0 warnings | 0 notes on the local platform.

Expected on incoming checks: the "New submission" NOTE.

## Comments

* First CRAN release of a package that has been developed publicly at
  https://github.com/max578/proxymix through eighteen tagged releases,
  with an independent-oracle validation battery shipped under
  inst/validation.
* All vignettes build with knitr::rmarkdown; the heavier simulation
  studies are pre-computed and the vignette code paths are fast.
