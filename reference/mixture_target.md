# Three-component Gaussian-mixture target

A toy three-component planar mixture target where everything is known
exactly: the `log_density` matches the GMM density formula and the
attached `samples` are drawn from that same mixture. Useful for sanity-
checking the three fitting regimes against ground truth.

## Usage

``` r
mixture_target(with_samples = FALSE, n = 2000L, seed = 1L)
```

## Arguments

- with_samples:

  If `TRUE`, attach `n` exact mixture draws.

- n:

  Number of samples to attach when `with_samples = TRUE`.

- seed:

  Optional integer seed used when drawing the samples.

## Value

A
[gmm_target](https://max578.github.io/proxymix/reference/gmm_target.md)
in dimension 2.

## See also

Other targets:
[`banana_target()`](https://max578.github.io/proxymix/reference/banana_target.md),
[`donut_target()`](https://max578.github.io/proxymix/reference/donut_target.md),
[`epanechnikov_target()`](https://max578.github.io/proxymix/reference/epanechnikov_target.md),
[`gmm_target_from_samples()`](https://max578.github.io/proxymix/reference/gmm_target_from_samples.md),
[`maxent_target()`](https://max578.github.io/proxymix/reference/maxent_target.md)

## Examples

``` r
m <- mixture_target(with_samples = TRUE, n = 100L)
m
#> <gmm_target>: "three_mixture" in p = 2 dimensions
#>   log_density : supplied
#>   samples     : 100 x 2 matrix
#>   normalised  : TRUE
#>   log Z(f)    : 0
```
