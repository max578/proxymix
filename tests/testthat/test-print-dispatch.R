## The S3 print methods must be found by print() called from outside the
## package, which looks them up through base::print.

test_that("the namespace has no copy of print", {
  expect_false(exists("print", envir = asNamespace("proxymix"),
                      inherits = FALSE))
})

test_that("S3 print methods are found from the global environment", {
  for (cl in c("gmm_eos_test", "gmm_ensemble", "proxymix_evidence",
               "proxymix_gate", "proxymix_selection", "proxymix_cf_sources")) {
    expect_false(
      is.null(utils::getS3method("print", cl, optional = TRUE,
                                 envir = globalenv())),
      label = cl
    )
  }
  prior <- gmm(weights = 1, means = list(0), covariances = list(matrix(10)))
  y <- withr::with_seed(1L, c(stats::rnorm(119), 6))
  r <- gmm_eos_test(prior, list(A = matrix(1), Q = matrix(0.04)),
                    list(C = matrix(1), R = matrix(1)), y, m = 1L,
                    method = "andrews")
  expect_output(evalq(print(r), list(r = r), globalenv()),
                "End-of-sample instability test")
})
