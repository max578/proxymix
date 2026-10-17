# Tidy a Gaussian mixture into a component table

A broom-style `tidy()` method: one row per component, with the weight,
the mean coordinates (`mean_1`, ...), and the marginal variances
(`var_1`, ...). Available as `generics::tidy(g)` (or `broom::tidy(g)`)
when the `generics` package is installed.

## Value

A data frame with `K` rows.

## Arguments

Called as `generics::tidy(x, ...)`, where `x` is a
[gmm](https://max578.github.io/proxymix/dev/reference/gmm.md) (or
[gmm_fit](https://max578.github.io/proxymix/dev/reference/gmm_fit.md))
and `...` is ignored.

## See also

Other classes:
[`autoplot.gmm_fit`](https://max578.github.io/proxymix/dev/reference/autoplot.gmm_fit.md),
[`glance.gmm_fit`](https://max578.github.io/proxymix/dev/reference/glance.gmm_fit.md),
[`gmm()`](https://max578.github.io/proxymix/dev/reference/gmm.md),
[`gmm_counterfactual_law()`](https://max578.github.io/proxymix/dev/reference/gmm_counterfactual_law.md),
[`gmm_dim()`](https://max578.github.io/proxymix/dev/reference/gmm_dim.md),
[`gmm_fit()`](https://max578.github.io/proxymix/dev/reference/gmm_fit.md),
[`gmm_n_components()`](https://max578.github.io/proxymix/dev/reference/gmm_n_components.md),
[`gmm_target()`](https://max578.github.io/proxymix/dev/reference/gmm_target.md),
[`gmm_weights()`](https://max578.github.io/proxymix/dev/reference/gmm_weights.md),
[`is_proposal()`](https://max578.github.io/proxymix/dev/reference/is_proposal.md)

## Examples

``` r
g <- gmm(weights = c(0.3, 0.7), means = list(c(-1, 0), c(2, 1)),
         covariances = list(diag(2), 0.5 * diag(2)))
generics::tidy(g)
#>   component weight mean_1 mean_2 var_1 var_2
#> 1         1    0.3     -1      0   1.0   1.0
#> 2         2    0.7      2      1   0.5   0.5
```
