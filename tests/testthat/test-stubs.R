## Signature-stability tests for Tier-2 stubs.
##
## These keep the public API frozen across releases: when a stub graduates
## to a real implementation, the *function signature* must not change, only
## the body. Update these expectations only when the graduation deliberately
## changes the public API.

stub_args <- function(fn) names(formals(fn))

test_that("from_aggregate_likelihood has stable arguments", {
  expect_equal(stub_args(from_aggregate_likelihood),
               c("y", "latent_aggregator", "N", "..."))
})

test_that("fit_kld_em_collider has stable arguments", {
  expect_equal(stub_args(fit_kld_em_collider),
               c("target", "dag", "N", "..."))
})

test_that("to_apsim_scenarios has stable arguments", {
  expect_equal(stub_args(to_apsim_scenarios),
               c("fit", "n", "schema", "..."))
})

test_that("from_simulator has stable arguments", {
  expect_equal(stub_args(from_simulator),
               c("simulator", "design", "..."))
})

test_that("stubs signal a proxymix_not_yet_implemented condition", {
  for (fn in list(
    function() from_aggregate_likelihood(matrix(0, 1, 1),
                                         latent_aggregator = identity),
    function() fit_kld_em_collider(banana_target(),
                                   dag = matrix(0, 2, 2)),
    function() to_apsim_scenarios(
      fit_proxymix(
        gmm_target_from_samples(matrix(stats::rnorm(200), ncol = 2)),
        N = 2L, regime = "sample", max_iter = 10L
      ),
      n = 5L, schema = list()
    ),
    function() from_simulator(identity, design = matrix(0, 1, 1))
  )) {
    expect_error(fn(), class = "proxymix_not_yet_implemented")
  }
})
