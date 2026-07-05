# Canonicalise the component ordering of a Gaussian mixture

Returns a new [gmm](https://max578.github.io/proxymix/reference/gmm.md)
(or [gmm_fit](https://max578.github.io/proxymix/reference/gmm_fit.md))
with the components permuted into a canonical order: weight descending,
then `||mu||` descending as a tiebreaker. The mixture distribution is
unchanged — only the bookkeeping order is — but the canonical ordering
removes the EM label-switching nuisance from snapshot tests, cross-run
comparisons, and printed summaries.

## Usage

``` r
gmm_canonicalise(g)
```

## Arguments

- g:

  A [gmm](https://max578.github.io/proxymix/reference/gmm.md) (or
  [gmm_fit](https://max578.github.io/proxymix/reference/gmm_fit.md))
  object.

## Value

A [gmm](https://max578.github.io/proxymix/reference/gmm.md) (or
[gmm_fit](https://max578.github.io/proxymix/reference/gmm_fit.md)) of
the same subclass as `g`, with the components permuted into canonical
order.

## Details

Applied automatically by the regime-specific fitters
([`fit_moment_match()`](https://max578.github.io/proxymix/reference/fit_moment_match.md),
[`fit_em_samples()`](https://max578.github.io/proxymix/reference/fit_em_samples.md),
[`fit_kld_em()`](https://max578.github.io/proxymix/reference/fit_kld_em.md))
and by the top-level dispatcher
[`fit_proxymix()`](https://max578.github.io/proxymix/reference/fit_proxymix.md)
when `canonicalise = TRUE` (the default).

## See also

Other ops:
[`dgmm()`](https://max578.github.io/proxymix/reference/dgmm.md),
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
g <- gmm(weights = c(0.1, 0.6, 0.3),
         means = list(c(0, 0), c(3, 0), c(-1, 1)),
         covariances = list(diag(2), diag(2), diag(2)))
gmm_canonicalise(g)
#> <gmm>: K = 3 components in p = 2 dimensions
#>   [1] w = 0.6000, |mu| = 3.0000, tr(Sigma) = 2.0000
#>   [2] w = 0.3000, |mu| = 1.4142, tr(Sigma) = 2.0000
#>   [3] w = 0.1000, |mu| = 0.0000, tr(Sigma) = 2.0000
```
