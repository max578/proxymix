# Donut-shaped 2-D target

A rotationally symmetric annulus on \\\mathbb{R}^2\\, with density
\$\$f(x) \propto \exp\\\left(-\tfrac{(\Vert x \Vert - r_0)^2}{2
\sigma^2}\right).\$\$ Numerical integration in polar coordinates fixes
the normaliser; the returned target exposes a normalised `log_density`.

## Usage

``` r
donut_target(r0 = 2.5, sigma = 0.5, with_samples = FALSE, n = 2000L, seed = 1L)
```

## Arguments

- r0:

  Centre radius of the annulus.

- sigma:

  Annulus width.

- with_samples:

  If `TRUE`, attach `n` exact samples via polar change-of-variables and
  a one-dimensional rejection step.

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
[`epanechnikov_target()`](https://max578.github.io/proxymix/reference/epanechnikov_target.md),
[`gmm_target_from_samples()`](https://max578.github.io/proxymix/reference/gmm_target_from_samples.md),
[`maxent_target()`](https://max578.github.io/proxymix/reference/maxent_target.md),
[`mixture_target()`](https://max578.github.io/proxymix/reference/mixture_target.md)

## Examples

``` r
d <- donut_target()
d
#> <gmm_target>: "donut" in p = 2 dimensions
#>   log_density : supplied
#>   samples     : <absent>
#>   normalised  : TRUE
#>   log Z(f)    : 0
```
