# Kullback-Leibler divergence between two Gaussian mixtures

Estimates `KL(p || q)` between two Gaussian mixtures by Monte Carlo and,
optionally, evaluates the Hershey–Olsen variational approximation as a
deterministic sanity check.

## Usage

``` r
gmm_kld(p, q, n_mc = 5000L, variational = TRUE)
```

## Arguments

- p, q:

  Two [gmm](https://max578.github.io/proxymix/reference/gmm.md) (or
  [gmm_fit](https://max578.github.io/proxymix/reference/gmm_fit.md))
  objects of the same ambient dimension.

- n_mc:

  Number of Monte Carlo samples drawn from `p`.

- variational:

  If `TRUE`, also return the Hershey–Olsen variational approximation.

## Value

A list with components

- `mc` - the Monte Carlo estimate of `KL(p || q)`,

- `mc_se` - its Monte Carlo standard error,

- `variational` - the variational approximation (`NA` if
  `variational = FALSE`),

- `n_mc` - the number of Monte Carlo samples used.

## Details

The Monte Carlo estimator draws `n_mc` samples from `p` and returns the
empirical mean of `log p(x) - log q(x)`, together with a Monte Carlo
standard error.

The variational approximation is \$\$\widehat{D}\_{\mathrm{var}}(p \Vert
q) = \sum_a \pi_a \log\\\left(\frac{\sum\_{a'} \pi\_{a'} \\
e^{-\mathrm{KL}(p_a \Vert p\_{a'})}}{\sum_b \omega_b \\
e^{-\mathrm{KL}(p_a \Vert q_b)}}\right),\$\$ which is exact when
`p == q` and tends to be a usable lower bound when the components of `p`
and `q` are well-separated. The closed-form Gaussian–Gaussian KL
`KL(p_a || q_b)` is used internally.

## See also

Other ops:
[`dgmm()`](https://max578.github.io/proxymix/reference/dgmm.md),
[`gmm_canonicalise()`](https://max578.github.io/proxymix/reference/gmm_canonicalise.md),
[`gmm_conditionalise()`](https://max578.github.io/proxymix/reference/gmm_conditionalise.md),
[`gmm_divergence()`](https://max578.github.io/proxymix/reference/gmm_divergence.md),
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
q <- gmm(weights = 1,
         means = list(c(0, 0)),
         covariances = list(diag(2) * 2))
gmm_kld(p, q, n_mc = 500L)
#> $mc
#> [1] 0.1067818
#> 
#> $mc_se
#> [1] 0.01585305
#> 
#> $variational
#> [1] -0.123072
#> 
#> $n_mc
#> [1] 500
#> 
```
