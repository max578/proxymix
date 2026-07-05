# Convolution of two independent Gaussian mixtures

The exact distribution of \\X + Y\\ for independent \\X \sim g_1\\ and
\\Y \sim g_2\\: a Gaussian mixture with \\K_1 K_2\\ components, \$\$g_1
\* g_2 = \sum\_{ij} w_i v_j\\ \mathcal{N}\\\left(\mu_i + m_j,\\
\Sigma_i + S_j\right).\$\$

## Usage

``` r
gmm_convolve(g1, g2)
```

## Arguments

- g1, g2:

  Two [gmm](https://max578.github.io/proxymix/reference/gmm.md) (or
  [gmm_fit](https://max578.github.io/proxymix/reference/gmm_fit.md))
  objects of the same ambient dimension.

## Value

A [gmm](https://max578.github.io/proxymix/reference/gmm.md) with
`K1 * K2` components.

## Details

For the affine special case \\X + c\\ with a constant `c`, or \\A X +
\epsilon\\ with Gaussian \\\epsilon\\, use
[`gmm_affine()`](https://max578.github.io/proxymix/reference/gmm_affine.md);
the convolution operator is the general mixture-plus-mixture case that
[`gmm_affine()`](https://max578.github.io/proxymix/reference/gmm_affine.md)
cannot express.

## See also

Other operators:
[`gmm_affine()`](https://max578.github.io/proxymix/reference/gmm_affine.md),
[`gmm_aggregate()`](https://max578.github.io/proxymix/reference/gmm_aggregate.md),
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
g1 <- gmm(weights = c(0.5, 0.5), means = list(-1, 1),
          covariances = list(matrix(0.5), matrix(0.5)))
g2 <- gmm(weights = 1, means = list(2), covariances = list(matrix(1)))
gmm_convolve(g1, g2)
#> <convolve(gmm, gmm)>: K = 2 components in p = 1 dimensions
#>   [1] w = 0.5000, |mu| = 1.0000, tr(Sigma) = 1.5000
#>   [2] w = 0.5000, |mu| = 3.0000, tr(Sigma) = 1.5000
```
