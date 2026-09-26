# Conditional-independence (Gaussian graphical model) structure of a mixture

Returns the undirected second-order conditional-independence graph of a
fitted Gaussian mixture: the partial-correlation (Gaussian graphical
model) structure of the mixture's overall covariance. An edge \\i - j\\
is present when the partial correlation \\r\_{ij}\\ of coordinates \\i\\
and \\j\\ given all the others is judged nonzero, and absent otherwise –
the latter is the Markov statement \\x_i \perp x_j \mid
x\_{\mathrm{rest}}\\ at second order. The overall covariance
\$\$\mathrm{Cov}(X) = \sum_k w_k (\Sigma_k + \mu_k \mu_k^\top) -
\bar\mu\\\bar\mu^\top, \qquad \bar\mu = \sum_k w_k \mu_k,\$\$ is
closed-form in the mixture parameters, so no sampling is required.

## Usage

``` r
gmm_independence_graph(g, threshold = NULL, alpha = 0.05, n = NULL)
```

## Arguments

- g:

  A [gmm](https://max578.github.io/proxymix/reference/gmm.md) (or
  [gmm_fit](https://max578.github.io/proxymix/reference/gmm_fit.md))
  mixture.

- threshold:

  Non-negative partial-correlation magnitude above which an edge is
  drawn. `NULL` (the default) uses the Fisher-z test when a sample size
  is available and `0.05` otherwise. Cannot be combined with `n`.

- alpha:

  Per-edge level of the Fisher-z test, in `(0, 1)`. Used only when a
  sample size is available and `threshold` is not given.

- n:

  Number of data points the mixture was fitted to. `NULL` (the default)
  takes it from a regime `"sample"` or `"moment"` fit and otherwise
  leaves it unknown.

## Value

A symmetric integer adjacency matrix (`1` = edge, `0` = none) with the
coordinate names of `g`, carrying the partial-correlation matrix as the
`"pcor"` attribute.

## Details

When the mixture was fitted to \\n\\ data points, each partial
correlation is tested with Fisher's z: an edge is drawn when
\$\$\|\operatorname{atanh}(r\_{ij})\| \sqrt{n - p - 1} \> \Phi^{-1}(1 -
\alpha / 2).\$\$ A
[gmm_fit](https://max578.github.io/proxymix/reference/gmm_fit.md) from
regime `"sample"` or `"moment"` supplies \\n\\ as the number of rows of
its target's samples; `n` can also be given directly. The test is per
edge, and a graph on \\p\\ coordinates carries \\p(p-1)/2\\ of them:
`alpha = 0.05 / choose(p, 2)` keeps the chance of any false edge below
0.05. Without a sample size (a regime `"kld"` fit, or a plain
[gmm](https://max578.github.io/proxymix/reference/gmm.md)), an edge is
drawn when \\\|r\_{ij}\|\\ exceeds `threshold`, which defaults to `0.05`
there. Giving `threshold` always selects the magnitude rule.

This is a graphical-model (dependency-structure) diagnostic, not a
causal discovery method: it recovers the undirected Markov skeleton, not
edge directions. Its distinctive use is **regime (iii)**: composed with
[`fit_kld_em()`](https://max578.github.io/proxymix/reference/fit_kld_em.md),
it recovers the dependency structure of a target you can only *evaluate*
(an unnormalised energy / Gibbs density), where no sample exists to hand
a sampling-based estimator. Being second order, it sees dependencies
that enter the covariance; a purely higher-order coupling (zero
correlation, nonzero dependence) is not detected – raise the mixture's
component count, or read the coordinate-block dependence with
[`gmm_mutual_information()`](https://max578.github.io/proxymix/reference/gmm_mutual_information.md)
instead.

## See also

[`gmm_mutual_information()`](https://max578.github.io/proxymix/reference/gmm_mutual_information.md),
[`fit_kld_em()`](https://max578.github.io/proxymix/reference/fit_kld_em.md),
[`gmm_conditionalise()`](https://max578.github.io/proxymix/reference/gmm_conditionalise.md)

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
[`gmm_mutual_information()`](https://max578.github.io/proxymix/reference/gmm_mutual_information.md),
[`hellinger_mc()`](https://max578.github.io/proxymix/reference/hellinger_mc.md),
[`kld_trace()`](https://max578.github.io/proxymix/reference/kld_trace.md),
[`proxy_functional_ci()`](https://max578.github.io/proxymix/reference/proxy_functional_ci.md)

## Examples

``` r
## Regime (iii): the Markov structure of an evaluable-but-unsampleable density.
## A continuous chain field x1 - x2 - x3 (couplings only between neighbours).
energy <- function(X) {
  X <- matrix(X, ncol = 3)
  rowSums((X^2 - 1)^2) - 0.7 * (X[, 1] * X[, 2] + X[, 2] * X[, 3])
}
target <- gmm_target(n_dim = 3L, log_density = function(X) -energy(X))
g <- fit_kld_em(target, N = 8L, proposal = is_uniform(3L, -3, 3),
                is_size = 8000L, anneal = TRUE, seed = 1L, support_warn = FALSE)
gmm_independence_graph(g)            # recovers x1 - x2 - x3 (no x1 - x3 edge)
#>    x1 x2 x3
#> x1  0  1  0
#> x2  1  0  1
#> x3  0  1  0
#> attr(,"pcor")
#>            x1        x2         x3
#> x1 1.00000000 0.4793587 0.01689514
#> x2 0.47935874 1.0000000 0.49064195
#> x3 0.01689514 0.4906420 1.00000000
```
