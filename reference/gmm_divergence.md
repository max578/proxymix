# Divergence between two Gaussian mixtures

Computes a divergence between two Gaussian mixtures of the same ambient
dimension. The Cauchy-Schwarz divergence \$\$D\_{\mathrm{CS}}(p, q) =
\tfrac{1}{2}\log V(p, p) + \tfrac{1}{2}\log V(q, q) - \log V(p, q),\$\$
with \\V(p, q) = \int p(x) q(x)\\ dx\\, is closed-form, symmetric,
non-negative, and zero exactly when \\p \propto q\\. The `"kl"` option
delegates to
[`gmm_kld()`](https://max578.github.io/proxymix/reference/gmm_kld.md), a
Monte-Carlo estimate of the asymmetric Kullback-Leibler divergence
\\\mathrm{KL}(p \Vert q)\\.

## Usage

``` r
gmm_divergence(p, q, type = c("cs", "kl"), n_mc = 5000L)
```

## Arguments

- p, q:

  Two [gmm](https://max578.github.io/proxymix/reference/gmm.md) (or
  [gmm_fit](https://max578.github.io/proxymix/reference/gmm_fit.md))
  objects of the same ambient dimension.

- type:

  `"cs"` (closed-form symmetric Cauchy-Schwarz divergence, the default)
  or `"kl"` (delegates to
  [`gmm_kld()`](https://max578.github.io/proxymix/reference/gmm_kld.md)).

- n_mc:

  Number of Monte Carlo samples used when `type = "kl"`.

## Value

For `type = "cs"`, a non-negative numeric scalar. For `type = "kl"`, the
list returned by
[`gmm_kld()`](https://max578.github.io/proxymix/reference/gmm_kld.md).

## See also

[`gmm_entropy()`](https://max578.github.io/proxymix/reference/gmm_entropy.md),
[`gmm_kld()`](https://max578.github.io/proxymix/reference/gmm_kld.md)

Other ops:
[`dgmm()`](https://max578.github.io/proxymix/reference/dgmm.md),
[`gmm_canonicalise()`](https://max578.github.io/proxymix/reference/gmm_canonicalise.md),
[`gmm_conditionalise()`](https://max578.github.io/proxymix/reference/gmm_conditionalise.md),
[`gmm_kld()`](https://max578.github.io/proxymix/reference/gmm_kld.md),
[`gmm_marginalise()`](https://max578.github.io/proxymix/reference/gmm_marginalise.md),
[`gmm_mean()`](https://max578.github.io/proxymix/reference/gmm_mean.md),
[`gmm_modes()`](https://max578.github.io/proxymix/reference/gmm_modes.md),
[`pgmm()`](https://max578.github.io/proxymix/reference/pgmm.md),
[`rgmm()`](https://max578.github.io/proxymix/reference/rgmm.md)

## Examples

``` r
p <- gmm(weights = c(0.5, 0.5),
         means = list(c(-1, 0), c(1, 0)),
         covariances = list(diag(2), diag(2)))
q <- gmm(weights = 1, means = list(c(0, 0)),
         covariances = list(diag(2) * 2))
gmm_divergence(p, q)
#> [1] 0.03561544
gmm_divergence(p, p)
#> [1] 0
```
