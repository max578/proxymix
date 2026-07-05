# Banana-shaped 2-D target

A 2-D "banana" density obtained by warping an isotropic Gaussian through
the map \\(z_1, z_2) \mapsto (z_1, z_2 + \tfrac{1}{2}(z_1^2 - 1))\\. The
map has unit Jacobian, so the resulting density is exactly normalised.

## Usage

``` r
banana_target(with_samples = FALSE, n = 2000L, seed = 1L)
```

## Arguments

- with_samples:

  If `TRUE`, attach `n` exact samples drawn by the change-of-variables
  trick. Default `FALSE` — the target then exposes only its
  `log_density`, which is the regime-(iii) use case.

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
[`donut_target()`](https://max578.github.io/proxymix/reference/donut_target.md),
[`epanechnikov_target()`](https://max578.github.io/proxymix/reference/epanechnikov_target.md),
[`gmm_target_from_samples()`](https://max578.github.io/proxymix/reference/gmm_target_from_samples.md),
[`maxent_target()`](https://max578.github.io/proxymix/reference/maxent_target.md),
[`mixture_target()`](https://max578.github.io/proxymix/reference/mixture_target.md)

## Examples

``` r
b <- banana_target()
b
#> <gmm_target>: "banana" in p = 2 dimensions
#>   log_density : supplied
#>   samples     : <absent>
#>   normalised  : TRUE
#>   log Z(f)    : 0
b@log_density(matrix(c(0, 0, 1, 0), ncol = 2, byrow = TRUE))
#> [1] -1.962877 -2.337877
```
