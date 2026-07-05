# Fit an uplift / next-best-action model from a data frame

Assembles a joint Gaussian-mixture proxy over the outcome, the
treatment, and the covariates, and returns an
[uplift_model](https://max578.github.io/proxymix/reference/uplift_model.md)
that the decision verbs read in closed form. One fit yields prediction,
heterogeneous treatment effects, optimal actions, off-line policy value,
and an identification audit – see
[`proxy_cate()`](https://max578.github.io/proxymix/reference/proxy_cate.md),
[`proxy_decide()`](https://max578.github.io/proxymix/reference/proxy_decide.md),
[`proxy_policy_value()`](https://max578.github.io/proxymix/reference/proxy_policy_value.md)
and
[`proxy_identification_report()`](https://max578.github.io/proxymix/reference/proxy_identification_report.md).

## Usage

``` r
fit_uplift(
  data,
  outcome,
  treatment,
  covariates,
  N = "auto",
  regime = "auto",
  assume = c("ignorability", "latent_confounder"),
  outcome_type = c("continuous", "binary", "count"),
  n_grid = 1:4,
  seed = NULL,
  ...
)
```

## Arguments

- data:

  A data frame holding the outcome, treatment and covariate columns.

- outcome:

  A single column name – the outcome `Y`.

- treatment:

  A single column name – the binary treatment `T`.

- covariates:

  A character vector of one or more column names – the pre-treatment
  covariates `X`.

- N:

  Number of mixture components, or `"auto"` (the default) to select by
  BIC over `n_grid`.

- regime:

  One of `"auto"`, `"moment"`, `"sample"`, `"kld"`, forwarded to
  [`fit_proxymix()`](https://max578.github.io/proxymix/reference/fit_proxymix.md).
  The default `"auto"` uses classical EM on the supplied rows.

- assume:

  One of `"ignorability"` (the default) or `"latent_confounder"` – the
  identification regime the effects are read under. See
  [`proxy_cate()`](https://max578.github.io/proxymix/reference/proxy_cate.md)
  and
  [`proxy_confounding_gap()`](https://max578.github.io/proxymix/reference/proxy_confounding_gap.md).

- outcome_type:

  One of `"continuous"` (the default), `"binary"` or `"count"`. Effects
  are reported on the response scale via a discretised predictive for
  the non-continuous types; see
  [`proxy_cate()`](https://max578.github.io/proxymix/reference/proxy_cate.md).

- n_grid:

  Integer vector of candidate component counts used when `N = "auto"`.
  Default `1:4`.

- seed:

  Optional integer. When supplied, the fitting (including any random EM
  starts) runs under a fixed seed and the global RNG state is restored
  on exit, so the fit is reproducible without disturbing the caller's
  stream.

- ...:

  Additional arguments forwarded to
  [`fit_proxymix()`](https://max578.github.io/proxymix/reference/fit_proxymix.md)
  (e.g. `max_iter`, `n_starts`, `ridge_eps`).

## Value

An
[uplift_model](https://max578.github.io/proxymix/reference/uplift_model.md).

## Details

The component count `N` may be fixed or chosen automatically. With
`N = "auto"` the function sweeps `n_grid` and keeps the K that minimises
the joint BIC; the full BIC trace is stored in the model's metadata. The
treatment is binary at this version (`{t0, t1}`); a continuous dose is a
future extension.

## See also

[`proxy_cate()`](https://max578.github.io/proxymix/reference/proxy_cate.md),
[`proxy_decide()`](https://max578.github.io/proxymix/reference/proxy_decide.md),
[`proxy_identification_report()`](https://max578.github.io/proxymix/reference/proxy_identification_report.md)

Other decision:
[`gmm_cf_mean()`](https://max578.github.io/proxymix/reference/gmm_cf_mean.md),
[`gmm_cf_tail_prob()`](https://max578.github.io/proxymix/reference/gmm_cf_tail_prob.md),
[`gmm_cf_variance()`](https://max578.github.io/proxymix/reference/gmm_cf_variance.md),
[`gmm_counterfactual()`](https://max578.github.io/proxymix/reference/gmm_counterfactual.md),
[`gmm_intervene()`](https://max578.github.io/proxymix/reference/gmm_intervene.md),
[`proxy_cate()`](https://max578.github.io/proxymix/reference/proxy_cate.md),
[`proxy_confounding_gap()`](https://max578.github.io/proxymix/reference/proxy_confounding_gap.md),
[`proxy_decide()`](https://max578.github.io/proxymix/reference/proxy_decide.md),
[`proxy_identification_report()`](https://max578.github.io/proxymix/reference/proxy_identification_report.md),
[`proxy_overlap()`](https://max578.github.io/proxymix/reference/proxy_overlap.md),
[`proxy_policy_value()`](https://max578.github.io/proxymix/reference/proxy_policy_value.md),
[`proxy_predict()`](https://max578.github.io/proxymix/reference/proxy_predict.md),
[`proxy_regime_segments()`](https://max578.github.io/proxymix/reference/proxy_regime_segments.md),
[`proxy_retrospective_uplift()`](https://max578.github.io/proxymix/reference/proxy_retrospective_uplift.md),
[`proxy_uplift()`](https://max578.github.io/proxymix/reference/proxy_uplift.md),
[`uplift_identification()`](https://max578.github.io/proxymix/reference/uplift_identification.md),
[`uplift_model()`](https://max578.github.io/proxymix/reference/uplift_model.md)

## Examples

``` r
set.seed(1)
n <- 400L
x <- stats::rnorm(n)
t <- stats::rbinom(n, 1L, 0.5)
y <- 1 + 0.5 * x + (1 + x) * t + stats::rnorm(n, sd = 0.5)
dat <- data.frame(y = y, t = t, x = x)
m <- fit_uplift(dat, outcome = "y", treatment = "t", covariates = "x",
                N = 2L, regime = "sample", max_iter = 50L, seed = 1L)
m
#> <uplift_model>: K = 2 regimes, assume = "ignorability"
#>   outcome   : y (continuous)
#>   treatment : t  (arms 0 vs 1)
#>   covariates: x
#>   trained on: 400 units
```
