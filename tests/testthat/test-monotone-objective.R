## Stronger monotonicity test for regime (iii).
##
## EM minimises the fixed-IS-weighted objective
##   Q(theta) = - sum_n W_n log g_theta(x_n)
## or equivalently *maximises*
##   sum_n W_n log g_theta(x_n).
## With the IS sample and weights held fixed across iterations, the weighted
## objective must be monotone-up under exact EM (modulo ridge regularisation
## and empty-component re-seeding). This is a tighter property than "the
## reported KLD trace decreases", because the KLD trace omits the log f(x)
## term that EM does not control.

test_that("KLD-EM weighted Q-objective is monotone-up across iterations", {
  withr::with_seed(2026, {
    fit <- fit_kld_em(banana_target(), N = 3L,
                      is_size = 3000L, max_iter = 40L, seed = 7L)
    obj <- fit@diagnostics$weighted_obj_trace
    expect_true(length(obj) >= 2L)
    ## Allow tiny non-monotonicities (<= 1e-6) due to ridge regularisation.
    deltas <- diff(obj)
    expect_true(all(deltas >= -1e-6))
  })
})

test_that("weighted_obj_trace and kld_trace are consistent (kld = const - obj)", {
  withr::with_seed(2026, {
    fit <- fit_kld_em(banana_target(), N = 2L,
                      is_size = 2000L, max_iter = 15L, seed = 11L)
    obj <- fit@diagnostics$weighted_obj_trace
    kld <- fit@diagnostics$kld_trace
    ## kld_n = sum_n W_n (log f - log g) = (sum_n W_n log f) - obj_n.
    ## The first term is constant across iterations, so kld + obj should be
    ## constant up to floating point.
    sums <- obj + kld
    expect_lt(max(abs(sums - sums[1L])), 1e-10)
  })
})
