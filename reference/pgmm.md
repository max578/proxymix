# Distribution and quantile functions of a one-dimensional mixture

The exact cumulative distribution function \\F(x) = \sum_k w_k
\Phi\\\left((x - \mu_k) / \sigma_k\right)\\ of a one-dimensional
Gaussian mixture, and its inverse by monotone root-finding. Together
with [`dgmm()`](https://max578.github.io/proxymix/reference/dgmm.md) and
[`rgmm()`](https://max578.github.io/proxymix/reference/rgmm.md) these
complete the usual d/p/q/r quartet for the one-dimensional case; for
tail probabilities of a multivariate mixture, marginalise first
([`gmm_marginalise()`](https://max578.github.io/proxymix/reference/gmm_marginalise.md))
or push through the relevant linear functional
([`gmm_affine()`](https://max578.github.io/proxymix/reference/gmm_affine.md)).

## Usage

``` r
pgmm(q, g, lower.tail = TRUE)

qgmm(p, g)
```

## Arguments

- q:

  Numeric vector of quantiles.

- g:

  A one-dimensional
  [gmm](https://max578.github.io/proxymix/reference/gmm.md) (or
  [gmm_fit](https://max578.github.io/proxymix/reference/gmm_fit.md)).

- lower.tail:

  Logical; if `TRUE` (default), probabilities are \\P(X \le x)\\.

- p:

  Numeric vector of probabilities in `(0, 1)`.

## Value

A numeric vector the length of the first argument.

## See also

Other ops:
[`dgmm()`](https://max578.github.io/proxymix/reference/dgmm.md),
[`gmm_canonicalise()`](https://max578.github.io/proxymix/reference/gmm_canonicalise.md),
[`gmm_conditionalise()`](https://max578.github.io/proxymix/reference/gmm_conditionalise.md),
[`gmm_divergence()`](https://max578.github.io/proxymix/reference/gmm_divergence.md),
[`gmm_kld()`](https://max578.github.io/proxymix/reference/gmm_kld.md),
[`gmm_marginalise()`](https://max578.github.io/proxymix/reference/gmm_marginalise.md),
[`gmm_mean()`](https://max578.github.io/proxymix/reference/gmm_mean.md),
[`gmm_modes()`](https://max578.github.io/proxymix/reference/gmm_modes.md),
[`rgmm()`](https://max578.github.io/proxymix/reference/rgmm.md)

## Examples

``` r
g <- gmm(weights = c(0.4, 0.6), means = list(-2, 1),
         covariances = list(matrix(0.5), matrix(1)))
pgmm(c(-2, 0, 2), g)
#> [1] 0.2008099 0.4942576 0.9048068
qgmm(c(0.1, 0.5, 0.9), g)
#> [1] -2.47778037  0.03783803  1.96742159
```
