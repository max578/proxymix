# Pointwise product of two Gaussian mixtures

The normalised pointwise product of two mixture densities – the
conjugate Bayes update when one mixture plays the prior and the other
the likelihood. The product of two Gaussian mixtures is again a Gaussian
mixture with \\K_1 K_2\\ components: \$\$p(x)\\ q(x) \propto \sum\_{ij}
w_i v_j z\_{ij}\\ \mathcal{N}\\\left(x;\\ \mu\_{ij},\\
\Sigma\_{ij}\right),\$\$ with \\\Sigma\_{ij} = (\Sigma_i^{-1} +
S_j^{-1})^{-1}\\, \\\mu\_{ij} = \Sigma\_{ij} (\Sigma_i^{-1} \mu_i +
S_j^{-1} m_j)\\, and \\z\_{ij} = \mathcal{N}(\mu_i;\\ m_j,\\ \Sigma_i +
S_j)\\ evaluated in the log domain. The overall normaliser \\\int p\\ q
= \sum\_{ij} w_i v_j z\_{ij}\\ is returned in the result's metadata as
`log_integral` (it is the marginal evidence of the update).

## Usage

``` r
gmm_product(g1, g2, reduce = NULL)
```

## Arguments

- g1, g2:

  Two [gmm](https://max578.github.io/proxymix/reference/gmm.md) (or
  [gmm_fit](https://max578.github.io/proxymix/reference/gmm_fit.md))
  objects of the same ambient dimension.

- reduce:

  Optional positive integer: reduce the product to at most this many
  components via
  [`gmm_reduce()`](https://max578.github.io/proxymix/reference/gmm_reduce.md).
  Default `NULL` (no reduction).

## Value

A [gmm](https://max578.github.io/proxymix/reference/gmm.md) with
`K1 * K2` components (or at most `reduce`), with `metadata$log_integral`
recording \\\log \int p\\ q\\.

## Details

Because the component count multiplies, chained products grow
geometrically; pass `reduce` to cap the count via
[`gmm_reduce()`](https://max578.github.io/proxymix/reference/gmm_reduce.md)
after the product (the standard assumed-density trick), or manage the
budget yourself.

## See also

Other operators:
[`gmm_affine()`](https://max578.github.io/proxymix/reference/gmm_affine.md),
[`gmm_aggregate()`](https://max578.github.io/proxymix/reference/gmm_aggregate.md),
[`gmm_convolve()`](https://max578.github.io/proxymix/reference/gmm_convolve.md),
[`gmm_counterfactual()`](https://max578.github.io/proxymix/reference/gmm_counterfactual.md),
[`gmm_filter()`](https://max578.github.io/proxymix/reference/gmm_filter.md),
[`gmm_intervene()`](https://max578.github.io/proxymix/reference/gmm_intervene.md),
[`gmm_missing()`](https://max578.github.io/proxymix/reference/gmm_missing.md),
[`gmm_mix()`](https://max578.github.io/proxymix/reference/gmm_mix.md),
[`gmm_observe()`](https://max578.github.io/proxymix/reference/gmm_observe.md),
[`gmm_reduce()`](https://max578.github.io/proxymix/reference/gmm_reduce.md)

## Examples

``` r
prior <- gmm(weights = c(0.5, 0.5), means = list(-2, 2),
             covariances = list(matrix(1), matrix(1)))
lik <- gmm(weights = 1, means = list(0.5), covariances = list(matrix(0.5)))
post <- gmm_product(prior, lik)
post@weights
#> [1] 0.2086085 0.7913915
```
