# Condition a Gaussian mixture on the exact values of some coordinates

A structured wrapper around
[`gmm_conditionalise()`](https://max578.github.io/proxymix/reference/gmm_conditionalise.md)
for the common case where the observed coordinates are specified by
integer index rather than `NA`-padded vector. Equivalent to
[`gmm_observe()`](https://max578.github.io/proxymix/reference/gmm_observe.md)
with a selection matrix `A` and zero noise covariance, but routed
through the Schur-complement path for efficiency.

## Usage

``` r
gmm_missing(g, observed, values)
```

## Arguments

- g:

  A [gmm](https://max578.github.io/proxymix/reference/gmm.md) (or
  [gmm_fit](https://max578.github.io/proxymix/reference/gmm_fit.md)) in
  `R^p`.

- observed:

  Integer vector of indices in `seq_len(p)`. The coordinates to
  condition on (fully observed).

- values:

  Numeric vector of length `length(observed)`. The observed values, in
  the same order as `observed`.

## Value

A [gmm](https://max578.github.io/proxymix/reference/gmm.md) in
`R^(p - length(observed))`.

## See also

Other operators:
[`gmm_affine()`](https://max578.github.io/proxymix/reference/gmm_affine.md),
[`gmm_aggregate()`](https://max578.github.io/proxymix/reference/gmm_aggregate.md),
[`gmm_convolve()`](https://max578.github.io/proxymix/reference/gmm_convolve.md),
[`gmm_counterfactual()`](https://max578.github.io/proxymix/reference/gmm_counterfactual.md),
[`gmm_filter()`](https://max578.github.io/proxymix/reference/gmm_filter.md),
[`gmm_intervene()`](https://max578.github.io/proxymix/reference/gmm_intervene.md),
[`gmm_mix()`](https://max578.github.io/proxymix/reference/gmm_mix.md),
[`gmm_observe()`](https://max578.github.io/proxymix/reference/gmm_observe.md),
[`gmm_product()`](https://max578.github.io/proxymix/reference/gmm_product.md),
[`gmm_reduce()`](https://max578.github.io/proxymix/reference/gmm_reduce.md)

## Examples

``` r
g <- gmm(weights = c(0.4, 0.6),
         means = list(c(-1, 0), c(1, 0)),
         covariances = list(diag(2), diag(2)))
# Condition coord 2 on the value 0.5; keep coord 1.
gmm_missing(g, observed = 2L, values = 0.5)
#> <missing(gmm)>: K = 2 components in p = 1 dimensions
#>   [1] w = 0.4000, |mu| = 1.0000, tr(Sigma) = 1.0000
#>   [2] w = 0.6000, |mu| = 1.0000, tr(Sigma) = 1.0000
```
