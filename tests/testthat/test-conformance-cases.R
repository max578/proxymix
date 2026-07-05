## Run the shipped two-sided conformance case set end-to-end with the
## built-in serial driver, so the contract sweep gates `R CMD check`-time
## test runs rather than living only in an external harness. Skipped on
## CRAN for runtime.

test_that("the shipped conformance case set has no findings", {
  skip_on_cran()
  path <- system.file("validation", "contract-conformance.R",
                      package = "proxymix")
  skip_if(!nzchar(path), "conformance script not installed")
  withr::local_envvar(PROXYMIX_CONFORMANCE_OUT = "",
                      PROXYMIX_CONFORMANCE_ENGINE = "")
  env <- new.env(parent = globalenv())
  withr::with_preserve_seed(
    suppressWarnings(suppressMessages(
      capture.output(sys.source(path, envir = env))
    ))
  )
  out <- env$out
  expect_gt(out$n_cells, 200L)
  expect_equal(out$n_findings, 0L)
  expect_equal(out$n_meta_fail, 0L)
})
