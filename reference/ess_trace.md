# Effective sample size of the importance-sampling weights

Returns the effective sample size (`1 / sum(W^2)`) of the
self-normalised importance weights used by a regime-(iii) fit. `NA` for
regimes that do not use importance sampling.

## Usage

``` r
ess_trace(fit)
```

## Arguments

- fit:

  A [gmm_fit](https://max578.github.io/proxymix/reference/gmm_fit.md).

## Value

Numeric scalar (or `NA_real_`).

## See also

Other diagnostics:
[`bic_aic()`](https://max578.github.io/proxymix/reference/bic_aic.md),
[`ess_summary()`](https://max578.github.io/proxymix/reference/ess_summary.md),
[`gmm_anneal_path()`](https://max578.github.io/proxymix/reference/gmm_anneal_path.md),
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
fit <- fit_proxymix(banana_target(), N = 2L, regime = "kld",
                    is_size = 1000L, max_iter = 15L, seed = 1L)
ess_trace(fit)
#> [1] 706.9948
```
