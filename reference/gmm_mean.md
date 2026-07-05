# Mean and covariance of a Gaussian mixture

The exact first two moments of a mixture: \$\$\mu = \sum_k w_k \mu_k,
\qquad \Sigma = \sum_k w_k \left(\Sigma_k + \mu_k \mu_k^\top\right) -
\mu \mu^\top.\$\$

## Usage

``` r
gmm_mean(g)

gmm_cov(g)
```

## Arguments

- g:

  A [gmm](https://max578.github.io/proxymix/reference/gmm.md) (or
  [gmm_fit](https://max578.github.io/proxymix/reference/gmm_fit.md)).

## Value

`gmm_mean()` returns a length-`p` numeric vector; `gmm_cov()` returns a
`p` by `p` numeric matrix.

## See also

Other ops:
[`dgmm()`](https://max578.github.io/proxymix/reference/dgmm.md),
[`gmm_canonicalise()`](https://max578.github.io/proxymix/reference/gmm_canonicalise.md),
[`gmm_conditionalise()`](https://max578.github.io/proxymix/reference/gmm_conditionalise.md),
[`gmm_divergence()`](https://max578.github.io/proxymix/reference/gmm_divergence.md),
[`gmm_kld()`](https://max578.github.io/proxymix/reference/gmm_kld.md),
[`gmm_marginalise()`](https://max578.github.io/proxymix/reference/gmm_marginalise.md),
[`gmm_modes()`](https://max578.github.io/proxymix/reference/gmm_modes.md),
[`pgmm()`](https://max578.github.io/proxymix/reference/pgmm.md),
[`rgmm()`](https://max578.github.io/proxymix/reference/rgmm.md)

## Examples

``` r
g <- gmm(weights = c(0.3, 0.7), means = list(c(-1, 0), c(1, 2)),
         covariances = list(diag(2), 0.5 * diag(2)))
gmm_mean(g)
#> [1] 0.4 1.4
gmm_cov(g)
#>      [,1] [,2]
#> [1,] 1.49 0.84
#> [2,] 0.84 1.49
```
