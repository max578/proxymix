# Bounded Gaussian-sum filtering over an observation series

Runs a Gaussian-sum filter over a series of `n` observations by
alternating the predict operator
[`gmm_affine()`](https://max578.github.io/proxymix/reference/gmm_affine.md),
the update operator
[`gmm_observe()`](https://max578.github.io/proxymix/reference/gmm_observe.md)
and an optional reduction
[`gmm_reduce()`](https://max578.github.io/proxymix/reference/gmm_reduce.md).
At one component this is the classical Kalman filter; at several
components, with Gaussian-mixture process or measurement noise, it is
the Gaussian-sum filter of Alspach and Sorenson (1972). The filter is a
pure composition of the existing closed-form operators; nothing here
leaves the affine-Gaussian world.

## Usage

``` r
gmm_filter(
  prior,
  dynamics,
  measurement,
  y,
  k_max = NULL,
  reduce = c("merge", "anneal"),
  ridge_eps = 1e-08
)
```

## Arguments

- prior:

  A [gmm](https://max578.github.io/proxymix/reference/gmm.md) (or
  [gmm_fit](https://max578.github.io/proxymix/reference/gmm_fit.md)) in
  \\\mathbb{R}^p\\: the initial belief.

- dynamics:

  A list `list(A =, b =, Q =)` describing the predict step, or a
  function `function(t)` returning such a list for time-varying
  dynamics. `A` is the `p`-by-`p` state-transition matrix; `b` is an
  optional length-`p` offset (default `0`); `Q` is the process noise – a
  `p`-by-`p` covariance matrix, a
  [gmm](https://max578.github.io/proxymix/reference/gmm.md) on
  \\\mathbb{R}^p\\ (a Gaussian-sum process noise), or `NULL` (a
  deterministic predict).

- measurement:

  A list `list(C =, R =, d =)` describing the update step, or a function
  `function(t)` returning such a list. `C` is the `m`-by-`p` observation
  matrix; `R` is the measurement noise – an `m`-by-`m` covariance matrix
  or a [gmm](https://max578.github.io/proxymix/reference/gmm.md) on
  \\\mathbb{R}^m\\ (a Gaussian-sum measurement noise); `d` is an
  optional length-`m` offset (default `0`).

- y:

  The observations: an `n`-by-`m` numeric matrix, a length-`n` list of
  length-`m` numeric vectors, or (when `m = 1`) a numeric vector of
  length `n`.

- k_max:

  Optional component cap. `NULL` (the default) disables reduction and
  runs the exact Kalman / Gaussian-sum filter; a positive integer caps
  the component count after each step via
  [`gmm_reduce()`](https://max578.github.io/proxymix/reference/gmm_reduce.md).

- reduce:

  The reduction method passed to
  [`gmm_reduce()`](https://max578.github.io/proxymix/reference/gmm_reduce.md)
  when `k_max` is set: `"merge"` (the default, a deterministic
  moment-preserving merge) or `"anneal"`.

- ridge_eps:

  Tiny ridge passed to
  [`gmm_affine()`](https://max578.github.io/proxymix/reference/gmm_affine.md)
  and
  [`gmm_observe()`](https://max578.github.io/proxymix/reference/gmm_observe.md)
  for numerical hygiene. Set to `0` for an exact (ridge-free) recursion.

## Value

A list with elements

- `filtered`:

  a length-`n` list of the filtered
  [gmm](https://max578.github.io/proxymix/reference/gmm.md) beliefs, one
  per step (after predict, update and any reduction).

- `mean`:

  an `n`-by-`p` matrix of the filtered mixture means.

- `cov`:

  a length-`n` list of the `p`-by-`p` filtered mixture covariances.

- `summary`:

  a data frame with one row per step: the step index, the per-coordinate
  filtered mean (`mean_1`, ...) and standard deviation (`sd_1`, ...),
  the component count `n_components`, and the step log marginal evidence
  `log_evidence` (whose sum over steps is the log likelihood of the
  series).

## Details

At step \\t\\ the belief is propagated through the linear-Gaussian
dynamics \\x_t = A x\_{t-1} + b + w\\, then conditioned on the
observation \\y_t = C x_t + d + v\\, then (if `k_max` is set) reduced
back to at most `k_max` components. With Gaussian noise the component
count is constant and the recursion is the Kalman filter. Non-Gaussian
noise is supplied as a Gaussian sum – a
[gmm](https://max578.github.io/proxymix/reference/gmm.md) in place of
the covariance matrix `Q` or `R`. Each such mixture multiplies the
component count by its own component count every step, so a long horizon
is only runnable with reduction enabled. Reduction is moment-preserving
(see
[`gmm_reduce()`](https://max578.github.io/proxymix/reference/gmm_reduce.md)),
so the filtered mean and covariance are unaffected to the
moment-matching order while the count stays bounded.

## References

Alspach, D. L. and Sorenson, H. W. (1972) Nonlinear Bayesian estimation
using Gaussian sum approximations. *IEEE Transactions on Automatic
Control* 17(4), 439–448.
[doi:10.1109/TAC.1972.1100034](https://doi.org/10.1109/TAC.1972.1100034)

## See also

Other operators:
[`gmm_affine()`](https://max578.github.io/proxymix/reference/gmm_affine.md),
[`gmm_aggregate()`](https://max578.github.io/proxymix/reference/gmm_aggregate.md),
[`gmm_convolve()`](https://max578.github.io/proxymix/reference/gmm_convolve.md),
[`gmm_counterfactual()`](https://max578.github.io/proxymix/reference/gmm_counterfactual.md),
[`gmm_intervene()`](https://max578.github.io/proxymix/reference/gmm_intervene.md),
[`gmm_missing()`](https://max578.github.io/proxymix/reference/gmm_missing.md),
[`gmm_mix()`](https://max578.github.io/proxymix/reference/gmm_mix.md),
[`gmm_observe()`](https://max578.github.io/proxymix/reference/gmm_observe.md),
[`gmm_product()`](https://max578.github.io/proxymix/reference/gmm_product.md),
[`gmm_reduce()`](https://max578.github.io/proxymix/reference/gmm_reduce.md)

## Examples

``` r
## A one-dimensional local-level filter (random walk observed in noise).
prior <- gmm(weights = 1, means = list(0), covariances = list(matrix(1)))
truth <- cumsum(stats::rnorm(20, sd = 0.3))
y <- truth + stats::rnorm(20, sd = 0.5)
out <- gmm_filter(
  prior,
  dynamics    = list(A = matrix(1), Q = matrix(0.09)),
  measurement = list(C = matrix(1), R = matrix(0.25)),
  y = y
)
head(out$summary)
#>   step      mean_1      sd_1 n_components log_evidence
#> 1    1  0.80627732 0.4509526            1   -1.4318718
#> 2    2  0.10428756 0.3673889            1   -2.1696288
#> 3    3  0.07995021 0.3441134            1   -0.5494708
#> 4    4 -0.09431233 0.3371356            1   -0.6891914
#> 5    5 -0.01022096 0.3350101            1   -0.5624068
#> 6    6  0.66062459 0.3343599            1   -3.0103132

## Heavy-tailed process noise as a two-component Gaussian sum, capped at
## four components per step.
q_noise <- gmm(weights = c(0.9, 0.1),
               means = list(0, 0),
               covariances = list(matrix(0.05), matrix(1)))
out2 <- gmm_filter(
  prior,
  dynamics    = list(A = matrix(1), Q = q_noise),
  measurement = list(C = matrix(1), R = matrix(0.25)),
  y = y, k_max = 4L
)
max(out2$summary$n_components)
#> [1] 4
```
