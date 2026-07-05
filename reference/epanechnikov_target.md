# Compact-support Epanechnikov target

The canonical compact-support target: a product of one-dimensional
Epanechnikov densities \\K(u) = \tfrac{3}{4}(1 -
u^2)\\\mathbf{1}\_{\|u\| \le 1}\\, rescaled to
`[center - half_width, center + half_width]` in each coordinate. No
mixture of full-support Gaussians can have compact support, so this
target is the clean case where regime (iii) is the only viable fitting
route. It declares its `support`, which makes
[`fit_kld_em()`](https://max578.github.io/proxymix/reference/fit_kld_em.md)
(and
[`fit_proxymix()`](https://max578.github.io/proxymix/reference/fit_proxymix.md)
with `regime = "kld"`) select a support-matched uniform proposal
automatically instead of the default multivariate-t, which would place
importance mass where the log-density is `-Inf`.

## Usage

``` r
epanechnikov_target(
  n_dim = 1L,
  center = 0,
  half_width = 1,
  with_samples = FALSE,
  n = 2000L,
  seed = 1L
)
```

## Arguments

- n_dim:

  Ambient dimension `p`.

- center:

  Length-1 (recycled) or length-`p` numeric centre per coordinate.

- half_width:

  Length-1 (recycled) or length-`p` positive numeric half-width per
  coordinate.

- with_samples:

  If `TRUE`, attach `n` exact samples drawn by the Devroye (1986)
  three-uniform method. Default `FALSE`.

- n:

  Number of samples to attach when `with_samples = TRUE`.

- seed:

  Optional integer seed used when drawing the samples.

## Value

A
[gmm_target](https://max578.github.io/proxymix/reference/gmm_target.md)
in dimension `n_dim` with a declared compact `support`.

## See also

Other targets:
[`banana_target()`](https://max578.github.io/proxymix/reference/banana_target.md),
[`donut_target()`](https://max578.github.io/proxymix/reference/donut_target.md),
[`gmm_target_from_samples()`](https://max578.github.io/proxymix/reference/gmm_target_from_samples.md),
[`maxent_target()`](https://max578.github.io/proxymix/reference/maxent_target.md),
[`mixture_target()`](https://max578.github.io/proxymix/reference/mixture_target.md)

## Examples

``` r
e <- epanechnikov_target()
e
#> <gmm_target>: "epanechnikov" in p = 1 dimensions
#>   log_density : supplied
#>   samples     : <absent>
#>   normalised  : TRUE
#>   log Z(f)    : 0
#>   support     : [-1, 1]
e@log_density(matrix(c(0, 0.5, 1.5), ncol = 1L)) # finite, finite, -Inf
#> [1] -0.2876821 -0.5753641       -Inf
```
