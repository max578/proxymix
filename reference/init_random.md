# Random initialisation

Generates an `N`-component random initialisation by perturbing isotropic
means around a centre. Useful as one starting point in a multi-start
best-of strategy.

## Usage

``` r
init_random(
  N = 1L,
  p = 2L,
  centre = rep(0, p),
  scale = 1,
  sigma_diag = 1,
  seed = NULL
)
```

## Arguments

- N:

  Number of components.

- p:

  Ambient dimension.

- centre:

  Length-`p` numeric vector — the location around which means are drawn.

- scale:

  Standard deviation of the mean perturbation.

- sigma_diag:

  Diagonal value used for the initial component covariances.

- seed:

  Optional integer seed.

## Value

A [gmm](https://max578.github.io/proxymix/reference/gmm.md) of `N`
components in dimension `p`.

## See also

Other init:
[`init_kmeans()`](https://max578.github.io/proxymix/reference/init_kmeans.md),
[`init_moment_seed()`](https://max578.github.io/proxymix/reference/init_moment_seed.md),
[`init_warm_start()`](https://max578.github.io/proxymix/reference/init_warm_start.md),
[`multi_start_best_of()`](https://max578.github.io/proxymix/reference/multi_start_best_of.md)

## Examples

``` r
init_random(N = 3L, p = 2L, seed = 1L)
#> <init_random>: K = 3 components in p = 2 dimensions
#>   [1] w = 0.3333, |mu| = 0.6528, tr(Sigma) = 2.0000
#>   [2] w = 0.3333, |mu| = 1.8009, tr(Sigma) = 2.0000
#>   [3] w = 0.3333, |mu| = 0.8842, tr(Sigma) = 2.0000
```
