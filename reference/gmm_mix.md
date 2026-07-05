# Mix Gaussian mixtures into one mixture

Flattens a list of Gaussian mixtures into a single mixture whose density
is the weighted average \\\sum_b \alpha_b\\ g_b(x)\\ – model averaging,
prior pooling, or the mixture-of-mixtures construction.

## Usage

``` r
gmm_mix(gmms, weights = NULL)
```

## Arguments

- gmms:

  A non-empty list of
  [gmm](https://max578.github.io/proxymix/reference/gmm.md) (or
  [gmm_fit](https://max578.github.io/proxymix/reference/gmm_fit.md))
  objects sharing one ambient dimension.

- weights:

  Optional non-negative mixing weights, one per element of `gmms`
  (normalised internally). Default: equal weights.

## Value

A [gmm](https://max578.github.io/proxymix/reference/gmm.md) with
`sum(K_b)` components.

## See also

Other operators:
[`gmm_affine()`](https://max578.github.io/proxymix/reference/gmm_affine.md),
[`gmm_aggregate()`](https://max578.github.io/proxymix/reference/gmm_aggregate.md),
[`gmm_convolve()`](https://max578.github.io/proxymix/reference/gmm_convolve.md),
[`gmm_counterfactual()`](https://max578.github.io/proxymix/reference/gmm_counterfactual.md),
[`gmm_filter()`](https://max578.github.io/proxymix/reference/gmm_filter.md),
[`gmm_intervene()`](https://max578.github.io/proxymix/reference/gmm_intervene.md),
[`gmm_missing()`](https://max578.github.io/proxymix/reference/gmm_missing.md),
[`gmm_observe()`](https://max578.github.io/proxymix/reference/gmm_observe.md),
[`gmm_product()`](https://max578.github.io/proxymix/reference/gmm_product.md),
[`gmm_reduce()`](https://max578.github.io/proxymix/reference/gmm_reduce.md)

## Examples

``` r
g1 <- gmm(weights = 1, means = list(-1), covariances = list(matrix(1)))
g2 <- gmm(weights = 1, means = list(2), covariances = list(matrix(0.5)))
gmm_mix(list(g1, g2), weights = c(0.7, 0.3))
#> <mix[2]>: K = 2 components in p = 1 dimensions
#>   [1] w = 0.7000, |mu| = 1.0000, tr(Sigma) = 1.0000
#>   [2] w = 0.3000, |mu| = 2.0000, tr(Sigma) = 0.5000
```
