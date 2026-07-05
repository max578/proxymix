# Moment-seed initialisation

Computes the global mean and covariance of the supplied samples and
spreads `N` components along the leading principal direction. Useful as
a deterministic starting point that survives multi-modal targets better
than a single-Gaussian fit.

## Usage

``` r
init_moment_seed(samples, N = 2L, spread = 1.5)
```

## Arguments

- samples:

  An `n` by `p` numeric matrix.

- N:

  Number of components.

- spread:

  Multiplier on the principal-direction standard deviation used to place
  the component means symmetrically about the global mean.

## Value

A [gmm](https://max578.github.io/proxymix/reference/gmm.md) of `N`
components in dimension `ncol(samples)`.

## See also

Other init:
[`init_kmeans()`](https://max578.github.io/proxymix/reference/init_kmeans.md),
[`init_random()`](https://max578.github.io/proxymix/reference/init_random.md),
[`init_warm_start()`](https://max578.github.io/proxymix/reference/init_warm_start.md),
[`multi_start_best_of()`](https://max578.github.io/proxymix/reference/multi_start_best_of.md)

## Examples

``` r
x <- matrix(stats::rnorm(200), ncol = 2)
init_moment_seed(x, N = 3L)
#> <init_moment_seed>: K = 3 components in p = 2 dimensions
#>   [1] w = 0.3333, |mu| = 1.7233, tr(Sigma) = 2.0251
#>   [2] w = 0.3333, |mu| = 0.1768, tr(Sigma) = 2.0251
#>   [3] w = 0.3333, |mu| = 1.5184, tr(Sigma) = 2.0251
```
