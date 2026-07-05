# A target density on R^p

An S7 representation of a target density that `proxymix` is asked to
approximate. A target may carry an evaluable `log_density`, a matrix of
i.i.d. `samples`, or both. Each of the three fitting regimes consumes a
different subset:

## Usage

``` r
gmm_target(
  n_dim = integer(0),
  log_density = NULL,
  samples = NULL,
  support = NULL,
  normalised = NA,
  log_normalizer = NA_real_,
  name = "gmm_target",
  metadata = list()
)
```

## Arguments

- n_dim:

  Integer scalar — the ambient dimension `p` (the property is called
  `n_dim` rather than `dim` because S7 reserves `dim` as an attribute
  name).

- log_density:

  Optional function: `function(x)` taking a numeric matrix `n` by `p`
  and returning a length-`n` numeric vector of `log f(x)`.

- samples:

  Optional `n` by `p` numeric matrix of i.i.d. samples from the target.

- support:

  Optional declaration of the target's support. `NULL` (the default)
  means the full \\\mathbb{R}^p\\. Otherwise a list
  `list(lower = , upper = )` of per-coordinate bounds, each of length 1
  (recycled) or `n_dim`, with `-Inf` / `Inf` permitted for unbounded
  coordinates. A bounded or one-sided support drives automatic
  support-matched importance proposal selection in regime (iii); see
  [`fit_kld_em()`](https://max578.github.io/proxymix/reference/fit_kld_em.md).

- normalised:

  Logical scalar declaring whether `log_density` integrates to one.
  `TRUE`, `FALSE`, or `NA` (unknown). Defaults to `NA`. Downstream
  diagnostics treat `NA` and `FALSE` identically and label any KLD
  estimate as shifted.

- log_normalizer:

  Numeric scalar `log Z(f)` of the supplied `log_density`, if known.
  Default `NA_real_`. When `normalised = FALSE` and `log_normalizer` is
  finite, downstream diagnostics can correct shifted KLD estimates by
  `+ log_normalizer`.

- name:

  Human-readable name.

- metadata:

  Optional list of additional descriptors.

## Value

An S7 object of class `gmm_target`.

## Details

- regime `"moment"` needs samples (or both moments via metadata);

- regime `"sample"` needs samples;

- regime `"kld"` needs the log-density.

Use `gmm_target()` or
[`gmm_target_from_samples()`](https://max578.github.io/proxymix/reference/gmm_target_from_samples.md)
to construct.

Importance-sampled KLD-EM (regime `"kld"`) only requires `log_density`
to be specified up to an unknown additive constant — the self-normalised
weights are invariant to scaling. The package's *diagnostics*
downstream, however, do depend on normalisation: an importance-sampled
KLD estimate against an unnormalised log-density measures
\\\widehat{KL}(f \Vert g) - \log Z(f)\\ rather than \\\widehat{KL}(f
\Vert g)\\, and a squared-Hellinger Monte Carlo estimate is only
meaningful when both densities integrate to one. Declare the target's
normalisation explicitly via `normalised` (and, where possible, supply
`log_normalizer`) so that the package can label shifted KLDs as shifted
and refuse misleading Hellinger reports.

## See also

Other classes:
[`autoplot.gmm_fit`](https://max578.github.io/proxymix/reference/autoplot.gmm_fit.md),
[`glance.gmm_fit`](https://max578.github.io/proxymix/reference/glance.gmm_fit.md),
[`gmm()`](https://max578.github.io/proxymix/reference/gmm.md),
[`gmm_counterfactual_law()`](https://max578.github.io/proxymix/reference/gmm_counterfactual_law.md),
[`gmm_dim()`](https://max578.github.io/proxymix/reference/gmm_dim.md),
[`gmm_fit()`](https://max578.github.io/proxymix/reference/gmm_fit.md),
[`gmm_n_components()`](https://max578.github.io/proxymix/reference/gmm_n_components.md),
[`gmm_weights()`](https://max578.github.io/proxymix/reference/gmm_weights.md),
[`is_proposal()`](https://max578.github.io/proxymix/reference/is_proposal.md),
[`tidy.gmm`](https://max578.github.io/proxymix/reference/tidy.gmm.md)

## Examples

``` r
tgt <- banana_target()
tgt
#> <gmm_target>: "banana" in p = 2 dimensions
#>   log_density : supplied
#>   samples     : <absent>
#>   normalised  : TRUE
#>   log Z(f)    : 0
```
