# Missingness mechanisms for multiple imputation

Specifications passed to the `mechanism` argument of
[`gmm_impute()`](https://max578.github.io/proxymix/reference/gmm_impute.md).
They describe how an entry came to be missing, which sets the
conditional the missing value is drawn from.

## Usage

``` r
mar()

mnar(coord, beta, link = c("logit", "probit"))

censored(coord, lower = -Inf, upper = Inf)
```

## Arguments

- coord:

  Name or index of the single coordinate the mechanism acts on.

- beta:

  Sensitivity slope of the log-odds (or probit score) of being missing
  in the missing value itself. Positive `beta` makes larger values more
  likely to be missing.

- link:

  Selection link, `"logit"` (the default) or `"probit"`.

- lower, upper:

  Bounds of the interval a censored missing entry is known to lie in. At
  least one must be finite.

## Value

A `proxymix_gate` object for
[`gmm_impute()`](https://max578.github.io/proxymix/reference/gmm_impute.md).

## Details

`mar()` is missing at random: the probability that an entry is missing
may depend on the observed entries but not on the missing value, so the
imputation conditional is the plain mixture conditional. This is the
default.

`censored()` is a known interval. A missing entry of `coord` is known
only to lie in `[lower, upper]` – a detection limit (`upper = LOD`), a
ceiling (`lower = cap`), or interval censoring. The imputation
conditional is the mixture conditional truncated to that interval, which
proxymix evaluates in closed form, so the imputations respect the bound
instead of substituting a constant such as half the detection limit.

`mnar()` is missing not at random through a selection model: an entry of
`coord` is missing with probability \\g(\alpha + \beta\\ y)\\ in its own
unobserved value \\y\\, where \\g\\ is the logistic or normal link. The
slope `beta` is the sensitivity parameter and is supplied, not estimated
– missing-not-at-random departures are not identified from the observed
data, so the appropriate use is to posit `beta`, propagate it, and
report how conclusions move with it (see
[`proxy_mnar_sensitivity()`](https://max578.github.io/proxymix/reference/proxy_mnar_sensitivity.md)).
The intercept is calibrated to the observed missingness rate. `beta = 0`
is missing at random.

## See also

[`gmm_impute()`](https://max578.github.io/proxymix/reference/gmm_impute.md),
[`proxy_mnar_sensitivity()`](https://max578.github.io/proxymix/reference/proxy_mnar_sensitivity.md).

Other imputation:
[`as_mids()`](https://max578.github.io/proxymix/reference/as_mids.md),
[`gmm_complete()`](https://max578.github.io/proxymix/reference/gmm_complete.md),
[`gmm_imputation()`](https://max578.github.io/proxymix/reference/gmm_imputation.md),
[`gmm_impute()`](https://max578.github.io/proxymix/reference/gmm_impute.md),
[`proxy_fmi()`](https://max578.github.io/proxymix/reference/proxy_fmi.md),
[`proxy_mnar_sensitivity()`](https://max578.github.io/proxymix/reference/proxy_mnar_sensitivity.md),
[`proxy_pool()`](https://max578.github.io/proxymix/reference/proxy_pool.md)

## Examples

``` r
mar()
#> $type
#> [1] "mar"
#> 
#> attr(,"class")
#> [1] "proxymix_mar"  "proxymix_gate"
censored("y", upper = 0.5)             # a lower detection limit at 0.5
#> $type
#> [1] "censored"
#> 
#> $coord
#> [1] "y"
#> 
#> $lower
#> [1] -Inf
#> 
#> $upper
#> [1] 0.5
#> 
#> attr(,"class")
#> [1] "proxymix_censored" "proxymix_gate"    
mnar("y", beta = 0.8)                  # larger y more likely missing
#> $type
#> [1] "mnar"
#> 
#> $coord
#> [1] "y"
#> 
#> $beta
#> [1] 0.8
#> 
#> $link
#> [1] "logit"
#> 
#> attr(,"class")
#> [1] "proxymix_mnar" "proxymix_gate"
```
