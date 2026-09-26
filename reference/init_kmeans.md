# k-means initialisation

Runs [`stats::kmeans()`](https://rdrr.io/r/stats/kmeans.html) on the
supplied samples and uses the resulting cluster centres and
within-cluster covariances to seed an EM-style fitter.

## Usage

``` r
init_kmeans(samples, N = 2L, ridge_eps = 1e-06, nstart = 10L)
```

## Arguments

- samples:

  An `n` by `p` numeric matrix of samples from (or close to) the target.

- N:

  Number of components.

- ridge_eps:

  Ridge added to each cluster covariance for numerical stability when a
  cluster has fewer than two points.

- nstart:

  [`stats::kmeans`](https://rdrr.io/r/stats/kmeans.html) `nstart`
  argument.

## Value

A [gmm](https://max578.github.io/proxymix/reference/gmm.md) of `N`
components in dimension `ncol(samples)`.

## See also

Other init:
[`init_moment_seed()`](https://max578.github.io/proxymix/reference/init_moment_seed.md),
[`init_random()`](https://max578.github.io/proxymix/reference/init_random.md),
[`init_warm_start()`](https://max578.github.io/proxymix/reference/init_warm_start.md),
[`multi_start_best_of()`](https://max578.github.io/proxymix/reference/multi_start_best_of.md)

## Examples

``` r
x <- matrix(stats::rnorm(200), ncol = 2)
init_kmeans(x, N = 3L)
#> <init_kmeans>: K = 3 components in p = 2 dimensions
#>   [1] w = 0.3200, |mu| = 1.0904, tr(Sigma) = 1.0259
#>   [2] w = 0.4500, |mu| = 0.7372, tr(Sigma) = 0.7710
#>   [3] w = 0.2300, |mu| = 1.5437, tr(Sigma) = 1.0441
```
