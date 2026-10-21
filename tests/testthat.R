library(testthat)
library(proxymix)

## CRAN allows at most two cores: at most two test workers, each with
## one BLAS thread and one data.table thread (tests/testthat/setup.R).
options(Ncpus = min(2L, as.integer(getOption("Ncpus", 2L))))

test_check("proxymix")
