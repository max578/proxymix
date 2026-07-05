# Compile an unnormalised Bayesian posterior into a `gmm_target`

Generic S3 constructor that turns a Bayesian posterior — represented
either by a fitted model object (e.g. from `brms` or `Stan`) or by a
bare callable — into a `gmm_target` suitable for regime (iii) of
[`fit_proxymix()`](https://max578.github.io/proxymix/reference/fit_proxymix.md)
/
[`fit_kld_em()`](https://max578.github.io/proxymix/reference/fit_kld_em.md).

## Usage

``` r
gmm_target_from_posterior(model, ...)

# Default S3 method
gmm_target_from_posterior(model, ...)

# S3 method for class '`function`'
gmm_target_from_posterior(
  model,
  ...,
  parameter_names = NULL,
  log_normalizer = NA_real_,
  name = NULL
)
```

## Arguments

- model:

  One of:

  - a function — a bare callable satisfying the contract above;

  - a fitted model object whose class registers a
    `gmm_target_from_posterior.<class>` method in its own package.

- ...:

  Forwarded to method-specific implementations.

- parameter_names:

  Character vector of parameter names. Required for the `function`
  method (or attached as `attr(model, "parameter_names")`). The length
  determines `n_dim`.

- log_normalizer:

  Numeric scalar `log Z` of the posterior, if known. `NA_real_` (the
  default) otherwise; downstream diagnostics will label any KLD estimate
  as shifted.

- name:

  Optional human-readable target name. Defaults to `"posterior"`.

## Value

A
[gmm_target](https://max578.github.io/proxymix/reference/gmm_target.md)
with `normalised = FALSE` and the user-supplied `log_normalizer` (or
`NA_real_`).

## Details

The contract for the underlying callable is:

- **Vectorised**: accepts a numeric matrix with rows indexing
  independent parameter draws and columns indexing parameters; returns a
  length-`nrow(theta)` numeric vector of `log p(theta | data) + const`.

- **Unnormalised is fine**: the marginal likelihood `log Z` is not
  required. Where the source package can supply it, pass
  `log_normalizer`.

- **Side-effect free**: no plotting, no mutable state. Pure function.

- **Domain-safe**: returns `-Inf` outside support rather than raising an
  error.

The default method errors with a hint pointing the user at either (a) a
Bayesian package that registers a method, or (b) the `function` method
below.

## Examples

``` r
# A trivial unnormalised log-posterior: a 2D banana centred near (1, 0).
log_post <- function(theta) {
  x <- theta[, 1L]
  y <- theta[, 2L]
  -0.5 * (x^2 + (y - 0.1 * x^2 + 1)^2)
}
tgt <- gmm_target_from_posterior(
  log_post,
  parameter_names = c("x", "y")
)
tgt
#> <gmm_target>: "posterior" in p = 2 dimensions
#>   log_density : supplied
#>   samples     : <absent>
#>   normalised  : FALSE
```
