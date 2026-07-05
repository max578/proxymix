# Maximum-entropy target under moment and support constraints

Builds the maximum-entropy
[gmm_target](https://max578.github.io/proxymix/reference/gmm_target.md)
consistent with the supplied constraints, the least-committal density
given what is known. The maximum- entropy density under linear
constraints is an exponential family \\p(x) \propto \exp(\eta^\top
T(x))\\ on the support, and three cases admit an exact closed form:

## Usage

``` r
maxent_target(moments = NULL, support = NULL, name = "maxent_target")
```

## Arguments

- moments:

  Either `NULL` (only a support constraint, giving the uniform) or a
  list with `mean` (length-`p` numeric) and `cov` (a `p`-by-`p`
  symmetric positive-definite matrix). Supplying `mean` without `cov`
  (or vice versa) is an error.

- support:

  `NULL` for the full \\\mathbb{R}^p\\ (only valid with second-moment
  `moments`, giving the Gaussian), or a list `list(lower = , upper = )`
  of per-coordinate box bounds (each length 1, recycled, or length `p`).
  The uniform case requires finite bounds.

- name:

  Human-readable name.

## Value

A
[gmm_target](https://max578.github.io/proxymix/reference/gmm_target.md).

## Details

- **First- and second-moment constraints on full support** – the
  Gaussian \\\mathcal{N}(\mathrm{mean}, \mathrm{cov})\\. The moment
  constraints are realised exactly, so a regime-(i) moment match
  recovers the target, and the target carries its `moments` in
  `metadata` for that purpose.

- **First- and second-moment constraints on a box support** – the same
  canonical Gaussian form restricted to the box, a *truncated Gaussian*.
  The normaliser is closed-form when `cov` is diagonal (a product of
  univariate Gaussian box probabilities) and the target is then exactly
  normalised; otherwise the target is declared unnormalised and
  regime (iii) fits it up to the unknown constant. Truncation shifts the
  realised moments inward, so the truncated density's mean and
  covariance are not `mean` and `cov`; it is the canonical-form
  maximum-entropy density on the box, not the moment-matched one.

- **A support constraint alone (no moments) on a finite box** – the
  uniform density, the maximum-entropy density on a compact support. Its
  differential entropy is exactly \\\log \mathrm{vol}(\mathrm{box})\\,
  the largest attainable on that support.

The bounded-support cases declare their `support`, so
[`fit_kld_em()`](https://max578.github.io/proxymix/reference/fit_kld_em.md)
(and
[`fit_proxymix()`](https://max578.github.io/proxymix/reference/fit_proxymix.md)
with `regime = "kld"`) selects a support-matched uniform importance
proposal automatically. Together with the Gaussian as the least-
committal full-support density and the
[`epanechnikov_target()`](https://max578.github.io/proxymix/reference/epanechnikov_target.md)
as a compact- support obstruction, these complete a family of principled
test targets.

## References

Jaynes, E. T. (1957) Information theory and statistical mechanics.
*Physical Review* 106(4), 620–630.
[doi:10.1103/PhysRev.106.620](https://doi.org/10.1103/PhysRev.106.620)

## See also

Other targets:
[`banana_target()`](https://max578.github.io/proxymix/reference/banana_target.md),
[`donut_target()`](https://max578.github.io/proxymix/reference/donut_target.md),
[`epanechnikov_target()`](https://max578.github.io/proxymix/reference/epanechnikov_target.md),
[`gmm_target_from_samples()`](https://max578.github.io/proxymix/reference/gmm_target_from_samples.md),
[`mixture_target()`](https://max578.github.io/proxymix/reference/mixture_target.md)

## Examples

``` r
## Full-support second-moment maximum entropy is the Gaussian.
g <- maxent_target(moments = list(mean = c(0, 0), cov = diag(2)))
g@log_density(matrix(c(0, 0), nrow = 1L))
#> [1] -1.837877

## Support alone on a box is the uniform.
u <- maxent_target(support = list(lower = 0, upper = 1))
exp(u@log_density(matrix(c(0.5), nrow = 1L)))
#> [1] 1

## Second moments on a box is a truncated Gaussian, fit via regime (iii).
tg <- maxent_target(moments = list(mean = 0, cov = matrix(1)),
                    support = list(lower = -2, upper = 2))
tg
#> <gmm_target>: "maxent_target" in p = 1 dimensions
#>   log_density : supplied
#>   samples     : <absent>
#>   normalised  : TRUE
#>   log Z(f)    : 0
#>   support     : [-2, 2]
```
