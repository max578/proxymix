# Plot a fitted Gaussian-mixture proxy

An `autoplot()` method for
[gmm_fit](https://max578.github.io/proxymix/reference/gmm_fit.md)
objects, rendering the fitted mixture with `ggplot2`. The displayed
coordinates are reduced to the requested one or two dimensions through
the closed-form marginal
[`gmm_marginalise()`](https://max578.github.io/proxymix/reference/gmm_marginalise.md),
so the method works for a proxy of any ambient dimension `p`.

A one-dimensional request draws the marginal mixture density, optionally
with the per-component densities underneath and a rug of the target's
samples. A two-dimensional request draws the marginal density as a
viridis raster with white contour lines, optionally overlaying each
component's mean and a probability-contour ellipse.

## Arguments

- object:

  A [gmm_fit](https://max578.github.io/proxymix/reference/gmm_fit.md),
  typically from
  [`fit_proxymix()`](https://max578.github.io/proxymix/reference/fit_proxymix.md).

- dims:

  Integer vector of length one or two giving the coordinate(s) to
  display, in `1:p`. Defaults to the first two coordinates (or the only
  coordinate when `p == 1`).

- n_grid:

  Integer scalar — the number of grid points per axis at which the
  density is evaluated. A two-dimensional plot evaluates `n_grid^2`
  points.

- n_sd:

  Numeric scalar — how many component standard deviations beyond the
  extreme component means the plotting window extends.

- level:

  Numeric scalar in `(0, 1)` — the probability level of the
  per-component ellipse drawn on a two-dimensional plot.

- show_components:

  Logical scalar — whether to overlay the per-component densities (one
  dimension) or mean-and-ellipse glyphs (two dimensions).

- show_data:

  Logical scalar — whether to overlay the target's samples, when the
  fitted target carries any.

- ...:

  Currently ignored, present for generic compatibility.

## Value

A `ggplot` object.

## Details

The method is registered against the
[`ggplot2::autoplot()`](https://ggplot2.tidyverse.org/reference/autoplot.html)
generic only when `ggplot2` is installed; call it as
`ggplot2::autoplot(fit)` or load `ggplot2` first. It returns the
`ggplot` object, so the usual `+` layering applies for further
customisation.

## See also

Other classes:
[`glance.gmm_fit`](https://max578.github.io/proxymix/reference/glance.gmm_fit.md),
[`gmm()`](https://max578.github.io/proxymix/reference/gmm.md),
[`gmm_counterfactual_law()`](https://max578.github.io/proxymix/reference/gmm_counterfactual_law.md),
[`gmm_dim()`](https://max578.github.io/proxymix/reference/gmm_dim.md),
[`gmm_fit()`](https://max578.github.io/proxymix/reference/gmm_fit.md),
[`gmm_n_components()`](https://max578.github.io/proxymix/reference/gmm_n_components.md),
[`gmm_target()`](https://max578.github.io/proxymix/reference/gmm_target.md),
[`gmm_weights()`](https://max578.github.io/proxymix/reference/gmm_weights.md),
[`is_proposal()`](https://max578.github.io/proxymix/reference/is_proposal.md),
[`tidy.gmm`](https://max578.github.io/proxymix/reference/tidy.gmm.md)

## Examples

``` r
samples <- matrix(stats::rnorm(200), ncol = 2)
tgt <- gmm_target_from_samples(samples)
fit <- fit_proxymix(tgt, N = 2L, regime = "sample", max_iter = 25L)
ggplot2::autoplot(fit)

ggplot2::autoplot(fit, dims = 1L)
```
