# Conditional of a Gaussian mixture

Computes the conditional distribution of a Gaussian mixture given fixed
values of a subset of coordinates, by the Schur-complement formula
applied component-wise and re-weighted by the marginal evidence
\\p(\textit{x}\_b)\\ of each component.

## Usage

``` r
gmm_conditionalise(g, given)
```

## Arguments

- g:

  A [gmm](https://max578.github.io/proxymix/reference/gmm.md) (or
  [gmm_fit](https://max578.github.io/proxymix/reference/gmm_fit.md))
  object.

- given:

  A length-`p` numeric vector. Coordinates to *condition on* take their
  numeric value; coordinates left *free* are `NA`.

## Value

A [gmm](https://max578.github.io/proxymix/reference/gmm.md) object in
dimension equal to the number of free coordinates.

## See also

Other ops:
[`dgmm()`](https://max578.github.io/proxymix/reference/dgmm.md),
[`gmm_canonicalise()`](https://max578.github.io/proxymix/reference/gmm_canonicalise.md),
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
gmm_conditionalise(g, given = c(NA, 0.5))
#> <conditionalise(gmm)>: K = 2 components in p = 1 dimensions
#>   [1] w = 0.5000, |mu| = 1.0000, tr(Sigma) = 1.0000
#>   [2] w = 0.5000, |mu| = 1.0000, tr(Sigma) = 1.0000
```
