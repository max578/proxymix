# Build a target from samples alone

Wraps a numeric matrix of i.i.d. samples as a
[gmm_target](https://max578.github.io/proxymix/reference/gmm_target.md).
The resulting target carries no `log_density`, so it can only feed
regimes `"moment"` (via empirical moments) and `"sample"` (classical
EM).

## Usage

``` r
gmm_target_from_samples(samples, name = "target_from_samples")
```

## Arguments

- samples:

  An `n` by `p` numeric matrix.

- name:

  Optional human-readable name. Defaults to `"target_from_samples"`.

## Value

A
[gmm_target](https://max578.github.io/proxymix/reference/gmm_target.md)
object.

## See also

Other targets:
[`banana_target()`](https://max578.github.io/proxymix/reference/banana_target.md),
[`donut_target()`](https://max578.github.io/proxymix/reference/donut_target.md),
[`epanechnikov_target()`](https://max578.github.io/proxymix/reference/epanechnikov_target.md),
[`maxent_target()`](https://max578.github.io/proxymix/reference/maxent_target.md),
[`mixture_target()`](https://max578.github.io/proxymix/reference/mixture_target.md)

## Examples

``` r
x <- matrix(stats::rnorm(200), ncol = 2)
tgt <- gmm_target_from_samples(x)
tgt
#> <gmm_target>: "target_from_samples" in p = 2 dimensions
#>   log_density : <absent>
#>   samples     : 100 x 2 matrix
#>   normalised  : unknown
```
