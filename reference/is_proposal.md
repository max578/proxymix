# An importance-sampling proposal

An `is_proposal` packages a sampler with its corresponding log-density.
Pass one to
[`fit_kld_em()`](https://max578.github.io/proxymix/reference/fit_kld_em.md)
(or
[`fit_proxymix()`](https://max578.github.io/proxymix/reference/fit_proxymix.md)
with `regime = "kld"`) to plug an alternative proposal into the
regime-(iii) loop.

## Usage

``` r
is_proposal(
  n_dim = integer(0),
  sample = NULL,
  log_density = NULL,
  name = "is_proposal",
  metadata = list()
)
```

## Arguments

- n_dim:

  Integer scalar — the ambient dimension `p` (the property is called
  `n_dim` rather than `dim` because S7 reserves `dim` as an attribute
  name).

- sample:

  Function: `function(n)` returning an `n` by `p` numeric matrix of
  independent draws.

- log_density:

  Function: `function(x)` taking a numeric matrix `n` by `p` and
  returning a length-`n` numeric vector of `log q(x)`.

- name:

  Human-readable name.

- metadata:

  Optional list of additional descriptors.

## Value

An S7 object of class `is_proposal`.

## See also

Other classes:
[`autoplot.gmm_fit`](https://max578.github.io/proxymix/reference/autoplot.gmm_fit.md),
[`glance.gmm_fit`](https://max578.github.io/proxymix/reference/glance.gmm_fit.md),
[`gmm()`](https://max578.github.io/proxymix/reference/gmm.md),
[`gmm_counterfactual_law()`](https://max578.github.io/proxymix/reference/gmm_counterfactual_law.md),
[`gmm_dim()`](https://max578.github.io/proxymix/reference/gmm_dim.md),
[`gmm_fit()`](https://max578.github.io/proxymix/reference/gmm_fit.md),
[`gmm_n_components()`](https://max578.github.io/proxymix/reference/gmm_n_components.md),
[`gmm_target()`](https://max578.github.io/proxymix/reference/gmm_target.md),
[`gmm_weights()`](https://max578.github.io/proxymix/reference/gmm_weights.md),
[`tidy.gmm`](https://max578.github.io/proxymix/reference/tidy.gmm.md)

## Examples

``` r
q <- is_mvn(n_dim = 2L, mean = c(0, 0), cov = diag(2))
q
#> <is_proposal>: "is_mvn" in p = 2 dimensions
```
