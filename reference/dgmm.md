# Density of a Gaussian mixture

Evaluates the density (or log-density) of a Gaussian mixture at one or
more points.

## Usage

``` r
dgmm(x, g, log = FALSE)
```

## Arguments

- x:

  A numeric matrix with one observation per row, or a length-`p` numeric
  vector (treated as a single observation).

- g:

  A [gmm](https://max578.github.io/proxymix/reference/gmm.md) (or
  [gmm_fit](https://max578.github.io/proxymix/reference/gmm_fit.md))
  object.

- log:

  Logical. If `TRUE`, return log-densities.

## Value

A numeric vector of length `nrow(x)`.

## See also

Other ops:
[`gmm_canonicalise()`](https://max578.github.io/proxymix/reference/gmm_canonicalise.md),
[`gmm_conditionalise()`](https://max578.github.io/proxymix/reference/gmm_conditionalise.md),
[`gmm_divergence()`](https://max578.github.io/proxymix/reference/gmm_divergence.md),
[`gmm_kld()`](https://max578.github.io/proxymix/reference/gmm_kld.md),
[`gmm_marginalise()`](https://max578.github.io/proxymix/reference/gmm_marginalise.md),
[`gmm_mean()`](https://max578.github.io/proxymix/reference/gmm_mean.md),
[`gmm_modes()`](https://max578.github.io/proxymix/reference/gmm_modes.md),
[`pgmm()`](https://max578.github.io/proxymix/reference/pgmm.md),
[`rgmm()`](https://max578.github.io/proxymix/reference/rgmm.md)

## Examples

``` r
g <- gmm(weights = c(0.5, 0.5),
         means = list(c(-1, 0), c(1, 0)),
         covariances = list(diag(2), diag(2)))
dgmm(c(0, 0), g)
#> [1] 0.09653235
dgmm(c(0, 0), g, log = TRUE)
#> [1] -2.337877
```
