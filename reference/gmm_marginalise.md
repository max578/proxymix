# Marginal of a Gaussian mixture

Computes the marginal distribution of a Gaussian mixture over a subset
of coordinates. The marginal of a Gaussian mixture is itself a Gaussian
mixture with the same weights.

## Usage

``` r
gmm_marginalise(g, keep)
```

## Arguments

- g:

  A [gmm](https://max578.github.io/proxymix/reference/gmm.md) (or
  [gmm_fit](https://max578.github.io/proxymix/reference/gmm_fit.md))
  object.

- keep:

  Integer vector of coordinate indices to retain (in `1..p`).

## Value

A [gmm](https://max578.github.io/proxymix/reference/gmm.md) object in
dimension `length(keep)`.

## See also

Other ops:
[`dgmm()`](https://max578.github.io/proxymix/reference/dgmm.md),
[`gmm_canonicalise()`](https://max578.github.io/proxymix/reference/gmm_canonicalise.md),
[`gmm_conditionalise()`](https://max578.github.io/proxymix/reference/gmm_conditionalise.md),
[`gmm_divergence()`](https://max578.github.io/proxymix/reference/gmm_divergence.md),
[`gmm_kld()`](https://max578.github.io/proxymix/reference/gmm_kld.md),
[`gmm_mean()`](https://max578.github.io/proxymix/reference/gmm_mean.md),
[`gmm_modes()`](https://max578.github.io/proxymix/reference/gmm_modes.md),
[`pgmm()`](https://max578.github.io/proxymix/reference/pgmm.md),
[`rgmm()`](https://max578.github.io/proxymix/reference/rgmm.md)

## Examples

``` r
g <- gmm(weights = c(0.5, 0.5),
         means = list(c(-1, 0, 2), c(1, 0, -2)),
         covariances = list(diag(3), diag(3)))
gmm_marginalise(g, keep = c(1L, 3L))
#> <marginalise(gmm)>: K = 2 components in p = 2 dimensions
#>   [1] w = 0.5000, |mu| = 2.2361, tr(Sigma) = 2.0000
#>   [2] w = 0.5000, |mu| = 2.2361, tr(Sigma) = 2.0000
```
