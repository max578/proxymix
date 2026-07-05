# Glance at a fitted Gaussian-mixture proxy

A broom-style `glance()` method: a one-row summary of a
[gmm_fit](https://max578.github.io/proxymix/reference/gmm_fit.md) with
the regime, the component count and dimension, convergence, iteration
count, and the regime's headline fit statistics. Available as
`generics::glance(fit)` when the `generics` package is installed.

## Arguments

- x:

  A [gmm_fit](https://max578.github.io/proxymix/reference/gmm_fit.md).

- ...:

  Ignored, for generic compatibility.

## Value

A one-row data frame.

## See also

Other classes:
[`autoplot.gmm_fit`](https://max578.github.io/proxymix/reference/autoplot.gmm_fit.md),
[`gmm()`](https://max578.github.io/proxymix/reference/gmm.md),
[`gmm_counterfactual_law()`](https://max578.github.io/proxymix/reference/gmm_counterfactual_law.md),
[`gmm_dim()`](https://max578.github.io/proxymix/reference/gmm_dim.md),
[`gmm_fit()`](https://max578.github.io/proxymix/reference/gmm_fit.md),
[`gmm_n_components()`](https://max578.github.io/proxymix/reference/gmm_n_components.md),
[`gmm_target()`](https://max578.github.io/proxymix/reference/gmm_target.md),
[`gmm_weights()`](https://max578.github.io/proxymix/reference/gmm_weights.md),
[`is_proposal()`](https://max578.github.io/proxymix/reference/is_proposal.md),
[`tidy.gmm`](https://max578.github.io/proxymix/reference/tidy.gmm.md)

## Examples

``` r
fit <- fit_proxymix(banana_target(), N = 2L, regime = "kld",
                    is_size = 1000L, max_iter = 10L, seed = 1L)
generics::glance(fit)
#>   regime n_components dim converged iterations      ess kld_final loglik_final
#> 1    kld            2   2      TRUE         10 706.9948 0.1123726           NA
#>   bic
#> 1  NA
```
