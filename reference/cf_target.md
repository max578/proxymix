# A target given by its characteristic function

Builds a
[gmm_target](https://max578.github.io/proxymix/reference/gmm_target.md)
that carries only a characteristic function \\\varphi(t) = E \exp(i
t^\top X)\\.
[`fit_proxymix()`](https://max578.github.io/proxymix/reference/fit_proxymix.md)
fits such a target with regime `"cf"`, which calls
[`from_cf()`](https://max578.github.io/proxymix/reference/from_cf.md).

## Usage

``` r
cf_target(cf, n_dim, name = "cf_target", metadata = list())
```

## Arguments

- cf:

  A function that takes an `n` by `n_dim` matrix of frequencies, one per
  row, and returns the characteristic function at each row as a complex
  vector of length `n`, as in
  [`from_cf()`](https://max578.github.io/proxymix/reference/from_cf.md).

- n_dim:

  Dimension of the target.

- name:

  Human-readable name.

- metadata:

  Optional list of additional descriptors.

## Value

A
[gmm_target](https://max578.github.io/proxymix/reference/gmm_target.md)
with `cf` set and no `log_density` or `samples`.

## Details

The function is checked before the target is built: it must return one
finite value for each row of its argument, equal 1 at zero, have modulus
at most 1, satisfy \\\varphi(-t) = \overline{\varphi(t)}\\, and fall to
one half along each coordinate axis between frequencies `1e-4` and
`1e4`. The checks use random frequencies, and the state of the random
number generator is left unchanged.

## See also

[`fit_proxymix()`](https://max578.github.io/proxymix/reference/fit_proxymix.md),
[`from_cf()`](https://max578.github.io/proxymix/reference/from_cf.md),
[`from_cf_elliptical()`](https://max578.github.io/proxymix/reference/from_cf_elliptical.md),
[`from_cf_sources()`](https://max578.github.io/proxymix/reference/from_cf_sources.md),
[`gmm_cf_distance()`](https://max578.github.io/proxymix/reference/gmm_cf_distance.md).

Other targets:
[`banana_target()`](https://max578.github.io/proxymix/reference/banana_target.md),
[`donut_target()`](https://max578.github.io/proxymix/reference/donut_target.md),
[`epanechnikov_target()`](https://max578.github.io/proxymix/reference/epanechnikov_target.md),
[`gmm_target_from_samples()`](https://max578.github.io/proxymix/reference/gmm_target_from_samples.md),
[`maxent_target()`](https://max578.github.io/proxymix/reference/maxent_target.md),
[`mixture_target()`](https://max578.github.io/proxymix/reference/mixture_target.md)

## Examples

``` r
## The sum of a Gamma(3) variable and independent N(0, 0.25) noise.
tgt <- cf_target(
  function(t) (1 - 1i * t[, 1])^(-3) * exp(-0.125 * t[, 1]^2),
  n_dim = 1L
)
tgt
#> <gmm_target>: "cf_target" in p = 1 dimensions
#>   log_density : <absent>
#>   samples     : <absent>
#>   cf          : supplied
#>   normalised  : unknown
fit <- fit_proxymix(tgt, N = 3L, n_freq = 1024L, seed = 1L)
#> Auto-selected regime: "cf".
gmm_fit_quality(fit)$cf_upper
#> [1] 0.00859195
```
