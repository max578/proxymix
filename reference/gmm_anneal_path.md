# Phase-transition component discovery by deterministic annealing

Tracks the number of distinct mixture centroids as a function of
temperature under mass-constrained deterministic annealing (Rose,
Gurewitz and Fox 1990), a physics-derived alternative to
information-criterion model selection. The system starts at a high
temperature where all `k_max` centroids collapse to the data centroid (a
single effective component) and is cooled along a geometric schedule; at
each critical temperature a centroid bifurcates, so the number of
distinct centroids grows in steps. The temperatures at which it grows
are the phase transitions, and the count occupying the widest
temperature range is the discovered component number.

## Usage

``` r
gmm_anneal_path(
  x,
  k_max = 8L,
  sigma = NULL,
  t_high = NULL,
  t_low = NULL,
  n_steps = 80L,
  n_inner = 30L,
  w = NULL,
  perturb = 0.02,
  merge_tol = 0.1,
  ridge_eps = 1e-06,
  seed = 1L
)
```

## Arguments

- x:

  A numeric `n` by `p` matrix of samples, or a
  [gmm_target](https://max578.github.io/proxymix/reference/gmm_target.md)
  carrying a `samples` matrix. For regime (iii) targets, pass an
  importance-resampled draw.

- k_max:

  Maximum number of centroids tracked (the discovered count is at most
  `k_max`).

- sigma:

  Reference scale: the shared covariance is `sigma^2 * I`. When `NULL`
  (the default) `sigma` is `1`, so the first critical temperature is the
  largest eigenvalue of the data covariance.

- t_high, t_low:

  Top and bottom of the temperature schedule. When `NULL` they default
  to `3 * t_critical_analytic` and `0.05 * t_critical_analytic`,
  bracketing the bifurcation cascade.

- n_steps:

  Number of temperatures on the geometric schedule.

- n_inner:

  Fixed-point iterations run at each temperature.

- w:

  Optional length-`n` vector of non-negative observation weights (e.g.
  importance weights). Defaults to uniform.

- perturb:

  Symmetry-breaking perturbation, as a fraction of the data scale,
  applied to the centroids at each temperature.

- merge_tol:

  Two centroids count as distinct when their distance exceeds
  `merge_tol` times the data scale.

- ridge_eps:

  Ridge added to the reference covariance for stability.

- seed:

  Optional integer seed for the perturbations (the result is
  deterministic given a seed).

## Value

A list with elements `path` (a data frame of `temperature`,
`n_effective` and `free_energy`), `critical_temperatures` (the
temperatures at which the count increased), `first_critical_temperature`
(the first such, or `NA` if none was detected), `t_critical_analytic`
(\\\lambda\_{\max}(\Sigma^{-1} C)\\), `k_selected` (the widest-plateau
component count), `lambda_max` and `sigma`.

## Details

The first bifurcation has a closed-form critical temperature \\T_c =
\lambda\_{\max}(\Sigma^{-1} C)\\, where \\C\\ is the (weighted) data
covariance and \\\Sigma = \sigma^2 I\\ the shared reference covariance.
This value is returned as `t_critical_analytic` and serves as an
independent analytic check on the empirically detected first transition.
Subsequent transitions have no comparably simple closed form, and the
count is a diagnostic rather than a guarantee.

Annealing fixes the component covariance to the reference \\\Sigma\\ so
the temperature is the only scale; this is the clean isotropic regime in
which the critical temperature is exact. For robust *fitting* under free
covariances, use `anneal = TRUE` on
[`fit_em_samples()`](https://max578.github.io/proxymix/reference/fit_em_samples.md)
or
[`fit_kld_em()`](https://max578.github.io/proxymix/reference/fit_kld_em.md)
instead.

## References

Rose, K., Gurewitz, E. and Fox, G. C. (1990) Statistical mechanics and
phase transitions in clustering. *Physical Review Letters* 65(8),
945–948.
[doi:10.1103/PhysRevLett.65.945](https://doi.org/10.1103/PhysRevLett.65.945)

## See also

Other diagnostics:
[`bic_aic()`](https://max578.github.io/proxymix/reference/bic_aic.md),
[`ess_summary()`](https://max578.github.io/proxymix/reference/ess_summary.md),
[`ess_trace()`](https://max578.github.io/proxymix/reference/ess_trace.md),
[`gmm_conditional_entropy()`](https://max578.github.io/proxymix/reference/gmm_conditional_entropy.md),
[`gmm_entropy()`](https://max578.github.io/proxymix/reference/gmm_entropy.md),
[`gmm_evidence()`](https://max578.github.io/proxymix/reference/gmm_evidence.md),
[`gmm_fit_ensemble()`](https://max578.github.io/proxymix/reference/gmm_fit_ensemble.md),
[`gmm_fit_quality()`](https://max578.github.io/proxymix/reference/gmm_fit_quality.md),
[`gmm_independence_graph()`](https://max578.github.io/proxymix/reference/gmm_independence_graph.md),
[`gmm_mutual_information()`](https://max578.github.io/proxymix/reference/gmm_mutual_information.md),
[`hellinger_mc()`](https://max578.github.io/proxymix/reference/hellinger_mc.md),
[`kld_trace()`](https://max578.github.io/proxymix/reference/kld_trace.md),
[`proxy_functional_ci()`](https://max578.github.io/proxymix/reference/proxy_functional_ci.md)

## Examples

``` r
set.seed(1)
x <- rbind(
  matrix(stats::rnorm(120, mean = -4), ncol = 2),
  matrix(stats::rnorm(120, mean =  4), ncol = 2)
)
path <- gmm_anneal_path(x, k_max = 4L, n_steps = 40L)
path$k_selected
#> [1] 2
path$first_critical_temperature
#> [1] 26.82229
```
