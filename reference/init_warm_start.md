# Warm-start initialisation from an existing fit

Returns the input as-is. Provided as a name so that the multi-start
driver can include "warm starts" by symbolic name.

## Usage

``` r
init_warm_start(g)
```

## Arguments

- g:

  A [gmm](https://max578.github.io/proxymix/reference/gmm.md) (or
  [gmm_fit](https://max578.github.io/proxymix/reference/gmm_fit.md)).

## Value

The input `g`, validated.

## See also

Other init:
[`init_kmeans()`](https://max578.github.io/proxymix/reference/init_kmeans.md),
[`init_moment_seed()`](https://max578.github.io/proxymix/reference/init_moment_seed.md),
[`init_random()`](https://max578.github.io/proxymix/reference/init_random.md),
[`multi_start_best_of()`](https://max578.github.io/proxymix/reference/multi_start_best_of.md)

## Examples

``` r
g <- gmm(weights = 1, means = list(c(0, 0)),
         covariances = list(diag(2)))
init_warm_start(g)
#> <gmm>: K = 1 components in p = 2 dimensions
#>   [1] w = 1.0000, |mu| = 0.0000, tr(Sigma) = 2.0000
```
