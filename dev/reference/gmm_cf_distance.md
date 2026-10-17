# Kolmogorov distance between a characteristic-function target and a mixture

For a direction \\a\\, compares the distribution of \\a^\top X\\ under
the target with its distribution under the mixture and returns a lower
and an upper bound for the largest absolute difference between the two
distribution functions (the Kolmogorov distance). The differences are
obtained by numerical inversion of the two characteristic functions
(Gil-Pelaez, 1951), so neither distribution is sampled.

## Usage

``` r
gmm_cf_distance(
  g,
  cf,
  direction = NULL,
  n_grid = 201L,
  probs = c(0.005, 0.995)
)
```

## Arguments

- g:

  A [gmm](https://max578.github.io/proxymix/dev/reference/gmm.md) or
  [gmm_fit](https://max578.github.io/proxymix/dev/reference/gmm_fit.md).

- cf:

  The target's characteristic function, as in
  [`from_cf()`](https://max578.github.io/proxymix/dev/reference/from_cf.md).

- direction:

  Numeric vector with one entry for each dimension of `g`, not all zero.
  Default: the first coordinate axis.

- n_grid:

  Number of grid points. Default `201L`.

- probs:

  The lowest and highest probability of the grid. Default
  `c(0.005, 0.995)`.

## Value

A list with `distance` (the largest difference on the grid), `at` (the
grid point where it is attained), `upper` (the upper bound described
above), `upper_inside` (the bound between the two quantiles in `probs`),
`truncation` (a bound on the omitted part of the inversion integral from
the highest frequency integrated to 20 times that frequency, added to
`upper` and `upper_inside`; the part beyond is ignored), `x` (the grid),
`difference` (target minus mixture distribution function at each grid
point), and `cutoff_gap` (the modulus of the difference between the two
characteristic functions at the highest frequency integrated). A
`cutoff_gap` or `truncation` that is not small means the integral was
cut off early, which happens when the target has a point mass or a very
slowly decaying characteristic function. The bounds do not hold then.

## Details

If `upper` is at most \\\varepsilon\\, every quantile of the target's
projection at probability \\q\\ lies between the mixture's quantiles at
probabilities \\q - \varepsilon\\ and \\q + \varepsilon\\.

`distance` is the largest difference at `n_grid` points, spaced equally
in mixture probability between the two values in `probs`. It is a lower
bound.

`upper` is an upper bound, valid when the characteristic function of the
projected target is integrable and `truncation` and `cutoff_gap` are
small. At a point between two grid points, the difference exceeds the
larger end value by at most the probability step, and exceeds the nearer
end value by at most half the spacing times a bound on the largest
density difference. Outside the grid it adds the tail probability
`probs[1]` or `1 - probs[2]` to the difference at the nearest end. Move
`probs` towards 0 and 1 to shrink the outside term.

## References

Gil-Pelaez, J. (1951). Note on the inversion theorem. *Biometrika*,
38(3–4), 481–482.
[doi:10.1093/biomet/38.3-4.481](https://doi.org/10.1093/biomet/38.3-4.481)

## See also

[`from_cf()`](https://max578.github.io/proxymix/dev/reference/from_cf.md),
[`gmm_cf()`](https://max578.github.io/proxymix/dev/reference/gmm_cf.md).

Other diagnostics:
[`bic_aic()`](https://max578.github.io/proxymix/dev/reference/bic_aic.md),
[`ess_summary()`](https://max578.github.io/proxymix/dev/reference/ess_summary.md),
[`ess_trace()`](https://max578.github.io/proxymix/dev/reference/ess_trace.md),
[`gmm_anneal_path()`](https://max578.github.io/proxymix/dev/reference/gmm_anneal_path.md),
[`gmm_conditional_entropy()`](https://max578.github.io/proxymix/dev/reference/gmm_conditional_entropy.md),
[`gmm_entropy()`](https://max578.github.io/proxymix/dev/reference/gmm_entropy.md),
[`gmm_evidence()`](https://max578.github.io/proxymix/dev/reference/gmm_evidence.md),
[`gmm_fit_ensemble()`](https://max578.github.io/proxymix/dev/reference/gmm_fit_ensemble.md),
[`gmm_fit_quality()`](https://max578.github.io/proxymix/dev/reference/gmm_fit_quality.md),
[`gmm_independence_graph()`](https://max578.github.io/proxymix/dev/reference/gmm_independence_graph.md),
[`gmm_mutual_information()`](https://max578.github.io/proxymix/dev/reference/gmm_mutual_information.md),
[`hellinger_mc()`](https://max578.github.io/proxymix/dev/reference/hellinger_mc.md),
[`kld_trace()`](https://max578.github.io/proxymix/dev/reference/kld_trace.md),
[`proxy_functional_ci()`](https://max578.github.io/proxymix/dev/reference/proxy_functional_ci.md)

## Examples

``` r
## A normal target against a normal with a shifted mean.
cf <- function(t) exp(-0.5 * t[, 1]^2)
g <- gmm(weights = 1, means = list(0.1), covariances = list(matrix(1)))
gmm_cf_distance(g, cf)$distance
#> [1] 0.03987815
2 * pnorm(0.05) - 1
#> [1] 0.03987761
```
