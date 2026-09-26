# Map the optima of an objective with a Gaussian-mixture proxy

Fits a Gaussian-mixture proxy to the Gibbs measure \\\exp(-f(x) / T)\\
of a user-supplied objective `f` over a bounded box, by cooling a short
temperature ladder through regime-(iii) importance-sampled KLD-EM
([`fit_kld_em()`](https://max578.github.io/proxymix/reference/fit_kld_em.md)).
As the temperature falls the mixture mass concentrates on the low
regions of `f`, so the fitted mixture is a closed-form *map* over the
optima rather than a single point estimate. Pair it with
[`gmm_modes()`](https://max578.github.io/proxymix/reference/gmm_modes.md)
to read off the distinct optima.

## Usage

``` r
from_objective(
  objective,
  lower,
  upper,
  N = NULL,
  minimise = TRUE,
  temperature = NULL,
  n_steps = 6L,
  exploration = 0.5,
  inflate = 1.8,
  is_size = 10000L,
  max_iter = 70L,
  ridge_eps = 1e-04,
  seed = NULL,
  scale = c("objective", "loglik")
)
```

## Arguments

- objective:

  Function taking a length-`p` numeric vector and returning a finite
  numeric scalar – the objective to minimise (or maximise; see
  `minimise`).

- lower, upper:

  Numeric vectors of equal length `p` giving the box over which the
  optima are sought. Every `upper` must exceed its `lower`. `objective`
  is only evaluated inside the box; points outside it are treated like
  points where `objective` is not finite and receive a strongly
  unattractive value.

- N:

  Number of mixture components. Default `max(10L, 5L * p)`. Raise it for
  objectives with many optima or strong symmetry.

- minimise:

  Logical. If `TRUE` (the default) lower `objective` values are better
  (the proxy concentrates on the minima); if `FALSE` the proxy
  concentrates on the maxima.

- temperature:

  Optional control of the cooling ladder. `NULL` (the default) derives a
  ladder from the range `r` of `objective` over a uniform probe of the
  box: from `r / 2` down to `r / 80` when `scale = "objective"`, and
  from `max(r / 2, 1)` down to 1 when `scale = "loglik"`. A positive
  scalar sets the final (lowest) temperature; a length-2 numeric
  `c(high, low)` sets both ends explicitly. A supplied `temperature` is
  used as given, whatever `scale` is.

- n_steps:

  Number of temperatures in the cooling ladder. Default `6L`. When
  `scale = "loglik"`, `temperature` is `NULL` and `n_steps` is not
  supplied, the ladder takes as many steps as it needs to cool by no
  more than the default factor of \\40^{1/5} \approx 2.1\\ per step, and
  never fewer than 6.

- exploration:

  Probability mass the importance proposal places on uniform exploration
  of the box at each step, in `[0, 1]`. Default `0.5`. Larger values
  explore more (and starve no basin); smaller values exploit the basins
  found so far.

- inflate:

  Factor by which the current mixture covariances are inflated when used
  as the exploitation part of the proposal. Default `1.8`.

- is_size:

  Importance-sample size per cooling step. Default `1e4L`.

- max_iter:

  Maximum EM iterations per cooling step. Default `70L`.

- ridge_eps:

  Ridge added to each component covariance at every M-step. Forwarded to
  [`fit_kld_em()`](https://max578.github.io/proxymix/reference/fit_kld_em.md).

- seed:

  Optional integer seed for reproducibility.

- scale:

  What the objective's values mean, which sets the default temperature
  ladder. `"objective"` (the default) is any objective to optimise: the
  ladder scales with its range, so the result does not depend on the
  objective's units. `"loglik"` is a negative log-likelihood (or a
  log-likelihood with `minimise = FALSE`): the ladder ends at
  temperature 1, where the Gibbs measure is the likelihood itself.

## Value

A [gmm_fit](https://max578.github.io/proxymix/reference/gmm_fit.md) (the
fitted proxy) carrying a `from_objective` metadata record with the box,
the temperature ladder, `scale`, and the effective sample size at each
cooling step (`ess`). Pass it to
[`gmm_modes()`](https://max578.github.io/proxymix/reference/gmm_modes.md)
to extract the distinct optima.

## Details

The Gibbs measure can be evaluated point-wise but not directly sampled,
which is precisely the setting of regime (iii): minimising the
Kullback-Leibler divergence from a Gaussian mixture to a peaked target
is the rank-weighted Gaussian update at the heart of
[`fit_kld_em()`](https://max578.github.io/proxymix/reference/fit_kld_em.md).
This function is that fit, driven against a sequence of cooling Gibbs
targets and warm-started from the previous fit at each step. Because a
multimodal `f` produces a multimodal target, the components spread
across the basins and recover the optima together.

Recovery is most reliable with component headroom – a number of
components `N` comfortably larger than the number of optima you expect,
so the defensive proposal can keep a component on each basin. Symmetric
landscapes (where several optima are exchangeable) need the most
headroom.

Dimensional scope. As with the rest of regime (iii), the importance-
sampling effective sample size falls sharply with dimension; the guard
is `p <= 5` (recommended), `p <= 10` (allowed with a warning), `p > 10`
(rejected).

Effective sample size. Each cooling step records its importance-sampling
effective sample size in the `from_objective` metadata. A classed
`proxymix_low_ess` warning is raised at most once per call, and only
when the final step, the one that produced the returned fit, falls below
50; its message gives the final effective sample size, how many steps
were low, and the lowest value with its step.

Likelihoods. With `scale = "loglik"` the objective is read as a negative
log-likelihood (or, with `minimise = FALSE`, a log-likelihood) and the
ladder ends at temperature 1. The fitted map is then a proxy for the
likelihood itself, restricted to the box, which is the posterior under a
uniform prior on the box.

## See also

[`gmm_modes()`](https://max578.github.io/proxymix/reference/gmm_modes.md)
to resolve the fitted map into distinct optima.

Other fitting:
[`fit_em_samples()`](https://max578.github.io/proxymix/reference/fit_em_samples.md),
[`fit_kld_em()`](https://max578.github.io/proxymix/reference/fit_kld_em.md),
[`fit_moment_match()`](https://max578.github.io/proxymix/reference/fit_moment_match.md),
[`from_kde()`](https://max578.github.io/proxymix/reference/from_kde.md),
[`select_N()`](https://max578.github.io/proxymix/reference/select_N.md)

## Examples

``` r
## A bimodal 1-D objective with minima at +/- 2.
f <- function(v) (v[1]^2 - 4)^2
fit <- from_objective(f, lower = -5, upper = 5, N = 6L,
                      is_size = 2000L, n_steps = 5L, seed = 1L)
gmm_modes(fit)$modes
#>           [,1]
#> [1,] -2.059751
#> [2,]  2.009906

## The likelihood surface of a logistic-regression slope.
x <- c(-1.2, -0.6, -0.1, 0.3, 0.8, 1.4)
y <- c(0, 0, 1, 0, 1, 1)
nll <- function(b) sum(log1p(exp(b[1] * x))) - sum(y * b[1] * x)
lik <- from_objective(nll, lower = -5, upper = 10, N = 4L,
                      is_size = 2000L, seed = 1L, scale = "loglik")
lik@metadata$from_objective$temperatures
#> [1] 9.172354 5.888228 3.779971 2.426567 1.557744 1.000000
```
