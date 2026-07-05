# Aggregation pushforward of a Gaussian mixture

A named alias for
[`gmm_affine()`](https://max578.github.io/proxymix/reference/gmm_affine.md)
when `A` is a (row-wise) aggregation matrix — e.g. a block-sum,
block-average, or unequal-weight aggregation used in downscaling
pipelines. The mathematics is identical to
[`gmm_affine()`](https://max578.github.io/proxymix/reference/gmm_affine.md);
the alias gives the public API a clearer hook for aggregation-specific
diagnostics in later releases.

## Usage

``` r
gmm_aggregate(g, A, noise_cov = NULL, ridge_eps = 1e-06)
```

## Arguments

- g:

  A [gmm](https://max578.github.io/proxymix/reference/gmm.md) (or
  [gmm_fit](https://max578.github.io/proxymix/reference/gmm_fit.md)) in
  `R^p`.

- A:

  An `m` by `p` numeric matrix.

- noise_cov:

  Optional `m` by `m` SPD numeric matrix. Default `NULL` (deterministic
  aggregation).

- ridge_eps:

  Tiny ridge added to the output covariances for numerical hygiene.

## Value

A [gmm](https://max578.github.io/proxymix/reference/gmm.md) in `R^m`.

## See also

Other operators:
[`gmm_affine()`](https://max578.github.io/proxymix/reference/gmm_affine.md),
[`gmm_convolve()`](https://max578.github.io/proxymix/reference/gmm_convolve.md),
[`gmm_counterfactual()`](https://max578.github.io/proxymix/reference/gmm_counterfactual.md),
[`gmm_filter()`](https://max578.github.io/proxymix/reference/gmm_filter.md),
[`gmm_intervene()`](https://max578.github.io/proxymix/reference/gmm_intervene.md),
[`gmm_missing()`](https://max578.github.io/proxymix/reference/gmm_missing.md),
[`gmm_mix()`](https://max578.github.io/proxymix/reference/gmm_mix.md),
[`gmm_observe()`](https://max578.github.io/proxymix/reference/gmm_observe.md),
[`gmm_product()`](https://max578.github.io/proxymix/reference/gmm_product.md),
[`gmm_reduce()`](https://max578.github.io/proxymix/reference/gmm_reduce.md)

## Examples

``` r
g <- gmm(weights = c(0.5, 0.5),
         means = list(c(-1, 0, 1), c(1, 0, -1)),
         covariances = list(diag(3), diag(3)))
# Sum coordinates 1 and 2 into a single aggregate; pass coord 3 through.
A <- matrix(c(1, 1, 0,
              0, 0, 1), nrow = 2L, byrow = TRUE)
gmm_aggregate(g, A)
#> <aggregate(gmm)>: K = 2 components in p = 2 dimensions
#>   [1] w = 0.5000, |mu| = 1.4142, tr(Sigma) = 3.0000
#>   [2] w = 0.5000, |mu| = 1.4142, tr(Sigma) = 3.0000
```
