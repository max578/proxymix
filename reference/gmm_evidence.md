# Estimate the target's normalising constant from a fitted proxy

Importance-sampling estimate of \\Z = \int f(x)\\ dx\\ for the fit's
target \\f\\, using the fitted mixture \\\hat g\\ as the proposal:
\$\$\widehat{Z} = \frac{1}{n} \sum\_{i=1}^n \frac{f(x_i)}{\hat g(x_i)},
\qquad x_i \sim \hat g,\$\$ computed in the log domain. For a Bayesian
posterior handed over as `likelihood x prior`, \\\log Z\\ is the log
marginal likelihood, so a fitted proxy doubles as a model-comparison
device.

## Usage

``` r
gmm_evidence(fit, n = 4000L, seed = NULL)
```

## Arguments

- fit:

  A [gmm_fit](https://max578.github.io/proxymix/reference/gmm_fit.md)
  whose target carries a `log_density`.

- n:

  Number of evidence draws from the fitted proxy.

- seed:

  Optional integer seed for the evidence draw.

## Value

A list of class `proxymix_evidence` with elements `log_z`, `se_log_z`,
`n`, `ess`, `max_weight_share`, `top_decile_share`, and `flagged` (the
heavy-tail indicator).

## Details

The estimator is exact in expectation for any proposal that dominates
\\f\\, and its Monte Carlo error is driven by how well \\\hat g\\
matches \\f\\ – which is precisely what the fit optimised. The variance
is finite only when \\\hat g\\ has tails at least as heavy as \\f\\; a
right-tail diagnostic is returned (the effective sample size and the
share of the estimate carried by the largest ten percent of weights,
which sits near 0.10 for a well-matched proxy), and a classed warning
(`proxymix_heavy_tail`) is raised when it indicates an untrustworthy
tail. Results also report the delta-method standard error of \\\log
\widehat{Z}\\.

When the target declares itself normalised (`normalised = TRUE`), the
true value is \\\log Z = 0\\ and the function still estimates it – a
useful end-to-end diagnostic of the fit.

## See also

Other diagnostics:
[`bic_aic()`](https://max578.github.io/proxymix/reference/bic_aic.md),
[`ess_summary()`](https://max578.github.io/proxymix/reference/ess_summary.md),
[`ess_trace()`](https://max578.github.io/proxymix/reference/ess_trace.md),
[`gmm_anneal_path()`](https://max578.github.io/proxymix/reference/gmm_anneal_path.md),
[`gmm_conditional_entropy()`](https://max578.github.io/proxymix/reference/gmm_conditional_entropy.md),
[`gmm_entropy()`](https://max578.github.io/proxymix/reference/gmm_entropy.md),
[`gmm_fit_ensemble()`](https://max578.github.io/proxymix/reference/gmm_fit_ensemble.md),
[`gmm_fit_quality()`](https://max578.github.io/proxymix/reference/gmm_fit_quality.md),
[`gmm_independence_graph()`](https://max578.github.io/proxymix/reference/gmm_independence_graph.md),
[`gmm_mutual_information()`](https://max578.github.io/proxymix/reference/gmm_mutual_information.md),
[`hellinger_mc()`](https://max578.github.io/proxymix/reference/hellinger_mc.md),
[`kld_trace()`](https://max578.github.io/proxymix/reference/kld_trace.md),
[`proxy_functional_ci()`](https://max578.github.io/proxymix/reference/proxy_functional_ci.md)

## Examples

``` r
## An unnormalised target with a known constant: log f = log N(., 0, I) + 3.
tgt <- gmm_target(
  n_dim = 2L,
  log_density = function(x) {
    if (is.null(dim(x))) x <- matrix(x, ncol = 2L)
    -0.5 * rowSums(x^2) - log(2 * pi) + 3
  },
  normalised = FALSE, name = "shifted_gaussian"
)
fit <- fit_kld_em(tgt, N = 1L, is_size = 2000L, max_iter = 40L, seed = 1L)
ev <- gmm_evidence(fit, n = 2000L, seed = 2L)
ev$log_z   # close to 3
#> [1] 3.00043
```
