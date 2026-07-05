# Number of components in a Gaussian mixture

Number of components in a Gaussian mixture

## Usage

``` r
gmm_n_components(x)
```

## Arguments

- x:

  A [gmm](https://max578.github.io/proxymix/reference/gmm.md) (or
  [gmm_fit](https://max578.github.io/proxymix/reference/gmm_fit.md))
  object.

## Value

Integer scalar.

## See also

Other classes:
[`autoplot.gmm_fit`](https://max578.github.io/proxymix/reference/autoplot.gmm_fit.md),
[`glance.gmm_fit`](https://max578.github.io/proxymix/reference/glance.gmm_fit.md),
[`gmm()`](https://max578.github.io/proxymix/reference/gmm.md),
[`gmm_counterfactual_law()`](https://max578.github.io/proxymix/reference/gmm_counterfactual_law.md),
[`gmm_dim()`](https://max578.github.io/proxymix/reference/gmm_dim.md),
[`gmm_fit()`](https://max578.github.io/proxymix/reference/gmm_fit.md),
[`gmm_target()`](https://max578.github.io/proxymix/reference/gmm_target.md),
[`gmm_weights()`](https://max578.github.io/proxymix/reference/gmm_weights.md),
[`is_proposal()`](https://max578.github.io/proxymix/reference/is_proposal.md),
[`tidy.gmm`](https://max578.github.io/proxymix/reference/tidy.gmm.md)

## Examples

``` r
g <- gmm(weights = c(0.5, 0.5), means = list(c(0, 0), c(1, 1)),
         covariances = list(diag(2), diag(2)))
gmm_n_components(g)
#> [1] 2
```
