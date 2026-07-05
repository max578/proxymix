# A Gaussian-mixture multiple-imputation result

The object returned by
[`gmm_impute()`](https://max578.github.io/proxymix/reference/gmm_impute.md).
It carries the `m` completed data matrices, the bootstrap-fitted
mixtures behind them (used by the analytic pooling in
[`proxy_pool()`](https://max578.github.io/proxymix/reference/proxy_pool.md)),
the mixture fitted to the full data, and a record of the missingness.
Pass it to
[`gmm_complete()`](https://max578.github.io/proxymix/reference/gmm_complete.md)
to extract the completed datasets and to
[`proxy_pool()`](https://max578.github.io/proxymix/reference/proxy_pool.md)
/
[`proxy_fmi()`](https://max578.github.io/proxymix/reference/proxy_fmi.md)
for inference.

## Usage

``` r
gmm_imputation(
  data = NULL,
  completions = list(),
  fits = list(),
  point_fit = NULL,
  n_components = integer(0),
  m = integer(0),
  mechanism = "mar",
  observed = NULL,
  var_names = character(0),
  is_data_frame = FALSE,
  diagnostics = list(),
  call = NULL
)
```

## Arguments

- data:

  The numeric data matrix supplied to
  [`gmm_impute()`](https://max578.github.io/proxymix/reference/gmm_impute.md),
  with `NA` for missing entries.

- completions:

  List of `m` completed data matrices.

- fits:

  List of `m` bootstrap-fitted
  [gmm](https://max578.github.io/proxymix/reference/gmm.md) objects
  behind the completions.

- point_fit:

  The [gmm](https://max578.github.io/proxymix/reference/gmm.md) fitted
  to the full data.

- n_components:

  Integer number of mixture components.

- m:

  Integer number of completions.

- mechanism:

  Missingness mechanism (currently `"mar"`).

- observed:

  Logical matrix marking the observed entries.

- var_names:

  Character vector of column names.

- is_data_frame:

  Logical; whether the input was a data frame.

- diagnostics:

  List of fit diagnostics (per-column missing rates, convergence,
  iterations).

- call:

  The matched call.

## Value

An S7 object of class `gmm_imputation`.

## See also

Other imputation:
[`as_mids()`](https://max578.github.io/proxymix/reference/as_mids.md),
[`gmm_complete()`](https://max578.github.io/proxymix/reference/gmm_complete.md),
[`gmm_impute()`](https://max578.github.io/proxymix/reference/gmm_impute.md),
[`mechanism`](https://max578.github.io/proxymix/reference/mechanism.md),
[`proxy_fmi()`](https://max578.github.io/proxymix/reference/proxy_fmi.md),
[`proxy_mnar_sensitivity()`](https://max578.github.io/proxymix/reference/proxy_mnar_sensitivity.md),
[`proxy_pool()`](https://max578.github.io/proxymix/reference/proxy_pool.md)
