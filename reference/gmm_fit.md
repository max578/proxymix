# A fitted Gaussian-mixture proxy

A `gmm_fit` is the result of
[`fit_proxymix()`](https://max578.github.io/proxymix/reference/fit_proxymix.md)
(or one of the regime-specific fitters). It inherits the mixture
parameters of [gmm](https://max578.github.io/proxymix/reference/gmm.md)
and adds a record of the target it was fitted to, the regime used, and
the iteration diagnostics.

## Usage

``` r
gmm_fit(
  weights = numeric(0),
  means = list(),
  covariances = list(),
  name = "gmm",
  metadata = list(),
  target = NULL,
  regime = NA_character_,
  diagnostics = list(),
  converged = NA,
  iterations = NA_integer_,
  call = NULL
)
```

## Arguments

- weights, means, covariances, name, metadata:

  See [gmm](https://max578.github.io/proxymix/reference/gmm.md).

- target:

  The
  [gmm_target](https://max578.github.io/proxymix/reference/gmm_target.md)
  the mixture was fitted to.

- regime:

  One of `"moment"`, `"sample"`, `"kld"`.

- diagnostics:

  A list of regime-specific diagnostics (see
  [`kld_trace()`](https://max578.github.io/proxymix/reference/kld_trace.md),
  [`ess_trace()`](https://max578.github.io/proxymix/reference/ess_trace.md)).

- converged:

  Logical scalar.

- iterations:

  Integer scalar.

- call:

  The matched call.

## Value

An S7 object inheriting from `gmm_fit` (and `gmm`).

## See also

Other classes:
[`autoplot.gmm_fit`](https://max578.github.io/proxymix/reference/autoplot.gmm_fit.md),
[`glance.gmm_fit`](https://max578.github.io/proxymix/reference/glance.gmm_fit.md),
[`gmm()`](https://max578.github.io/proxymix/reference/gmm.md),
[`gmm_counterfactual_law()`](https://max578.github.io/proxymix/reference/gmm_counterfactual_law.md),
[`gmm_dim()`](https://max578.github.io/proxymix/reference/gmm_dim.md),
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
inherits(fit, "proxymix::gmm_fit")
#> [1] TRUE
```
