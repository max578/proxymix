# Summary of importance-sampling diagnostics

Convenience accessor returning the headline IS-quality numbers for a
regime-(iii) fit: effective sample size and its ratio to `is_size`, the
largest self-normalised weight, the support fraction (proportion of
draws that received a finite weight), and the Monte-Carlo standard error
of the final KLD estimate. Returns `NA` fields for regimes that do not
use importance sampling.

## Usage

``` r
ess_summary(fit)
```

## Arguments

- fit:

  A [gmm_fit](https://max578.github.io/proxymix/reference/gmm_fit.md).

## Value

A list of numeric scalars (or `NA`s where not applicable).

## Details

Validation-side numbers (`validation_*`) are populated only when the fit
was called with `validation_size > 0`.

## See also

Other diagnostics:
[`bic_aic()`](https://max578.github.io/proxymix/reference/bic_aic.md),
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
[`kld_trace()`](https://max578.github.io/proxymix/reference/kld_trace.md),
[`proxy_functional_ci()`](https://max578.github.io/proxymix/reference/proxy_functional_ci.md)

## Examples

``` r
fit <- fit_proxymix(banana_target(), N = 3L, regime = "kld",
                    is_size = 1500L, max_iter = 20L, seed = 1L,
                    validation_size = 1500L)
ess_summary(fit)
#> $is_size
#> [1] 1500
#> 
#> $ess
#> [1] 1048.118
#> 
#> $ess_relative
#> [1] 0.6987454
#> 
#> $max_weight
#> [1] 0.005472202
#> 
#> $support_fraction
#> [1] 1
#> 
#> $mc_se_kld
#> [1] 0.009224739
#> 
#> $validation_size
#> [1] 1500
#> 
#> $validation_ess
#> [1] 1111.561
#> 
#> $validation_ess_relative
#> [1] 0.7410408
#> 
#> $validation_max_weight
#> [1] 0.003698124
#> 
#> $validation_support_fraction
#> [1] 1
#> 
#> $validation_kld
#> [1] 0.02868117
#> 
```
