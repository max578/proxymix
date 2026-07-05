# Per-iteration KLD trace of a fit

Returns the per-iteration estimate of `KL(f || g_theta)` produced during
a regime-(iii) fit, or `NA` for regimes that do not estimate the KLD
internally.

## Usage

``` r
kld_trace(fit)
```

## Arguments

- fit:

  A [gmm_fit](https://max578.github.io/proxymix/reference/gmm_fit.md).

## Value

Numeric vector (or `NA_real_`).

## See also

Other diagnostics:
[`bic_aic()`](https://max578.github.io/proxymix/reference/bic_aic.md),
[`ess_summary()`](https://max578.github.io/proxymix/reference/ess_summary.md),
[`ess_trace()`](https://max578.github.io/proxymix/reference/ess_trace.md),
[`gmm_anneal_path()`](https://max578.github.io/proxymix/reference/gmm_anneal_path.md),
[`gmm_conditional_entropy()`](https://max578.github.io/proxymix/reference/gmm_conditional_entropy.md),
[`gmm_entropy()`](https://max578.github.io/proxymix/reference/gmm_entropy.md),
[`gmm_evidence()`](https://max578.github.io/proxymix/reference/gmm_evidence.md),
[`gmm_fit_ensemble()`](https://max578.github.io/proxymix/reference/gmm_fit_ensemble.md),
[`gmm_fit_quality()`](https://max578.github.io/proxymix/reference/gmm_fit_quality.md),
[`gmm_independence_graph()`](https://max578.github.io/proxymix/reference/gmm_independence_graph.md),
[`gmm_mutual_information()`](https://max578.github.io/proxymix/reference/gmm_mutual_information.md),
[`hellinger_mc()`](https://max578.github.io/proxymix/reference/hellinger_mc.md),
[`proxy_functional_ci()`](https://max578.github.io/proxymix/reference/proxy_functional_ci.md)

## Examples

``` r
fit <- fit_proxymix(banana_target(), N = 2L, regime = "kld",
                    is_size = 1000L, max_iter = 15L, seed = 1L)
kld_trace(fit)
#>  [1] 0.1532797 0.1184025 0.1149076 0.1135196 0.1129067 0.1126276 0.1124956
#>  [8] 0.1124295 0.1123938 0.1123726
```
