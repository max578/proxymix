# Characteristic function of a Gaussian mixture

Evaluates \\\varphi_g(t) = \sum_k w_k \exp(i t^\top \mu_k - t^\top
\Sigma_k t / 2)\\ at each row of `t`.

## Usage

``` r
gmm_cf(g, t)
```

## Arguments

- g:

  A [gmm](https://max578.github.io/proxymix/dev/reference/gmm.md) or
  [gmm_fit](https://max578.github.io/proxymix/dev/reference/gmm_fit.md).

- t:

  Numeric matrix of frequencies with one frequency per row and as many
  columns as the mixture has dimensions. A vector is read as one
  frequency per element for a one-dimensional mixture and as a single
  frequency otherwise.

## Value

A complex vector with one value for each row of `t`.

## See also

[`from_cf()`](https://max578.github.io/proxymix/dev/reference/from_cf.md),
[`gmm_cf_distance()`](https://max578.github.io/proxymix/dev/reference/gmm_cf_distance.md).

Other ops:
[`dgmm()`](https://max578.github.io/proxymix/dev/reference/dgmm.md),
[`gmm_canonicalise()`](https://max578.github.io/proxymix/dev/reference/gmm_canonicalise.md),
[`gmm_conditionalise()`](https://max578.github.io/proxymix/dev/reference/gmm_conditionalise.md),
[`gmm_divergence()`](https://max578.github.io/proxymix/dev/reference/gmm_divergence.md),
[`gmm_kld()`](https://max578.github.io/proxymix/dev/reference/gmm_kld.md),
[`gmm_marginalise()`](https://max578.github.io/proxymix/dev/reference/gmm_marginalise.md),
[`gmm_mean()`](https://max578.github.io/proxymix/dev/reference/gmm_mean.md),
[`gmm_modes()`](https://max578.github.io/proxymix/dev/reference/gmm_modes.md),
[`pgmm()`](https://max578.github.io/proxymix/dev/reference/pgmm.md),
[`rgmm()`](https://max578.github.io/proxymix/dev/reference/rgmm.md)

## Examples

``` r
g <- gmm(
  weights = c(0.3, 0.7),
  means = list(c(0, 0), c(2, 1)),
  covariances = list(diag(2), 0.5 * diag(2))
)
gmm_cf(g, rbind(c(0, 0), c(0.5, -1)))
#> [1] 1.0000000+0i 0.6727094+0i
```
