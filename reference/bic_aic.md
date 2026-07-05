# Information criteria: BIC, AIC, and ICL

Returns the Bayesian and Akaike information criteria of a regime-(ii)
fit, together with the integrated completed likelihood (ICL). All three
are computed against the *empirical* log-likelihood of the samples used
to fit the model and are reported on the same scale (smaller is better).
They are `NA` for regimes that do not have an empirical likelihood
(`"moment"`, `"kld"`).

## Usage

``` r
bic_aic(fit)
```

## Arguments

- fit:

  A [gmm_fit](https://max578.github.io/proxymix/reference/gmm_fit.md).

## Value

A list with `bic`, `aic`, `icl`, `classification_entropy`, and
`n_params`.

## Details

The ICL of Biernacki, Celeux and Govaert (2000) adds to the BIC twice
the entropy of the fitted classification, \\\mathrm{ICL} =
\mathrm{BIC} + 2 E_N\\, where \\E_N = -\sum\_{i,k} \gamma\_{ik} \log
\gamma\_{ik} \ge 0\\ is the entropy of the responsibilities
\\\gamma\_{ik}\\. It therefore penalises mixtures whose components
overlap (uncertain assignments), and favours well-separated clustering
solutions over the merely best-fitting ones. Because \\E_N \ge 0\\, the
ICL is never smaller than the BIC, and the two coincide for a single
component (\\K = 1\\), where every responsibility is one. The
classification entropy itself is returned as `classification_entropy`.

## References

Biernacki, C., Celeux, G. and Govaert, G. (2000) Assessing a mixture
model for clustering with the integrated completed likelihood. *IEEE
Transactions on Pattern Analysis and Machine Intelligence* 22(7),
719–725. [doi:10.1109/34.865189](https://doi.org/10.1109/34.865189)

## See also

Other diagnostics:
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
[`kld_trace()`](https://max578.github.io/proxymix/reference/kld_trace.md),
[`proxy_functional_ci()`](https://max578.github.io/proxymix/reference/proxy_functional_ci.md)

## Examples

``` r
x <- matrix(stats::rnorm(200), ncol = 2)
tgt <- gmm_target_from_samples(x)
fit <- fit_proxymix(tgt, N = 2L, regime = "sample", max_iter = 25L)
bic_aic(fit)
#> $bic
#> [1] 637.4143
#> 
#> $aic
#> [1] 608.7574
#> 
#> $icl
#> [1] 718.8163
#> 
#> $classification_entropy
#> [1] 40.70102
#> 
#> $n_params
#> [1] 11
#> 
```
