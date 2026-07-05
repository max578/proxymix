# Multiple imputation by Gaussian-mixture conditioning

Fits a Gaussian mixture to a numeric dataset that contains missing
values and draws `m` completed datasets from the mixture conditional
\\p(x\_{\mathrm{missing}} \mid x\_{\mathrm{observed}})\\. Because the
mixture can be multimodal and heteroscedastic, the imputations follow
the shape of the joint distribution rather than a single Gaussian, which
keeps downstream inference valid on data that a single-Gaussian or
linear-Gaussian imputer mis-specifies.

## Usage

``` r
gmm_impute(
  data,
  N = NULL,
  m = 20L,
  mechanism = mar(),
  seed = NULL,
  max_iter = 100L,
  tol = 1e-06,
  ridge_eps = 1e-06
)
```

## Arguments

- data:

  A numeric matrix or data frame with `NA` for missing entries.

- N:

  Number of mixture components. `NULL` (the default) selects it by the
  Bayesian information criterion over `1:6`.

- m:

  Number of completed datasets to draw. Default `20L`.

- mechanism:

  A missingness mechanism:
  [`mar()`](https://max578.github.io/proxymix/reference/mechanism.md),
  [`censored()`](https://max578.github.io/proxymix/reference/mechanism.md),
  or
  [`mnar()`](https://max578.github.io/proxymix/reference/mechanism.md).
  The string `"mar"` is also accepted. Default
  [`mar()`](https://max578.github.io/proxymix/reference/mechanism.md).

- seed:

  Optional integer seed. When supplied the result is reproducible and
  the ambient random-number state is restored on exit.

- max_iter:

  Maximum EM iterations per fit. Default `100L`.

- tol:

  Relative log-likelihood tolerance for EM convergence. Default `1e-6`.

- ridge_eps:

  Ridge added to each component covariance at every M-step. Default
  `1e-6`.

## Value

A
[gmm_imputation](https://max578.github.io/proxymix/reference/gmm_imputation.md)
object.

## Details

Imputation is conditioning. For a row with observed coordinates the
missing coordinates follow the closed-form mixture conditional (the same
Schur-complement algebra as
[`gmm_conditionalise()`](https://max578.github.io/proxymix/reference/gmm_conditionalise.md)).
The mixture is fitted to the incomplete data by expectation-maximisation
whose E-step uses each row's observed margin and whose M-step restores
the conditional covariance of the filled entries, so component variances
are not under-estimated.

Proper multiple imputation requires the fitting parameters themselves to
carry uncertainty, otherwise the pooled intervals are too narrow. Each
of the `m` imputations is therefore drawn under a mixture fitted to an
independent bootstrap resample of the rows, so
[`proxy_pool()`](https://max578.github.io/proxymix/reference/proxy_pool.md)
reflects both imputation and parameter uncertainty.

The `mechanism` says how an entry came to be missing, which sets the
conditional the missing value is drawn from:
[`mar()`](https://max578.github.io/proxymix/reference/mechanism.md) (the
default) for missing at random,
[`censored()`](https://max578.github.io/proxymix/reference/mechanism.md)
for a known interval such as a detection limit, or
[`mnar()`](https://max578.github.io/proxymix/reference/mechanism.md) for
a value-dependent selection model. The interval and value-dependent
gates act on a single coordinate, and a row missing that coordinate must
have its other coordinates observed. Numeric data only; categorical
variables are out of scope.

## See also

[`gmm_complete()`](https://max578.github.io/proxymix/reference/gmm_complete.md)
to extract completions,
[`proxy_pool()`](https://max578.github.io/proxymix/reference/proxy_pool.md)
to pool an estimand across them,
[`gmm_conditionalise()`](https://max578.github.io/proxymix/reference/gmm_conditionalise.md)
for the conditioning algebra.

Other imputation:
[`as_mids()`](https://max578.github.io/proxymix/reference/as_mids.md),
[`gmm_complete()`](https://max578.github.io/proxymix/reference/gmm_complete.md),
[`gmm_imputation()`](https://max578.github.io/proxymix/reference/gmm_imputation.md),
[`mechanism`](https://max578.github.io/proxymix/reference/mechanism.md),
[`proxy_fmi()`](https://max578.github.io/proxymix/reference/proxy_fmi.md),
[`proxy_mnar_sensitivity()`](https://max578.github.io/proxymix/reference/proxy_mnar_sensitivity.md),
[`proxy_pool()`](https://max578.github.io/proxymix/reference/proxy_pool.md)

## Examples

``` r
set.seed(1)
x1 <- rnorm(200)
x2 <- x1 + rnorm(200)
x2[runif(200) < plogis(x1)] <- NA          # missing at random on x1
imp <- gmm_impute(cbind(x1, x2), N = 1L, m = 10L, seed = 1L)
proxy_pool(imp, "x2")$estimate             # pooled mean of x2
#> [1] -0.06410561
```
