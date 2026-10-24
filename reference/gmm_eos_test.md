# End-of-sample instability test on a Gaussian state-space filter

Tests whether the last `m` observations of a series are consistent with
a linear-Gaussian state-space model fitted on the rest, in the regime
where `m` is small (even `m = 1`) and smaller than the parameter count –
where ordinary structural-break tests (Chow, sup-Wald) are undefined
because the post-break parameters cannot be estimated. The statistic is
the sum of the last `m` squared standardised one-step innovations from
the filter; a break inflates it.

## Usage

``` r
gmm_eos_test(
  prior,
  dynamics,
  measurement,
  y,
  m = 1L,
  method = c("chisq", "andrews"),
  alpha = 0.05
)
```

## Arguments

- prior:

  A single-component
  [gmm](https://max578.github.io/proxymix/reference/gmm.md) giving the
  state prior (the test is defined for a fitted linear-Gaussian model;
  multi-component priors are not yet supported).

- dynamics:

  A list with `A` (the state-transition matrix), `Q` (the process-noise
  covariance) and an optional offset `b`, or a function `function(t)`
  returning such a list, exactly as for
  [`gmm_filter()`](https://max578.github.io/proxymix/reference/gmm_filter.md).
  Gaussian-sum (mixture) process noise is rejected: the calibrations
  below are defined for Gaussian innovations only.

- measurement:

  A list with `C` (the observation matrix), `R` (the observation-noise
  covariance) and an optional offset `d`, or a function `function(t)`
  returning such a list, exactly as for
  [`gmm_filter()`](https://max578.github.io/proxymix/reference/gmm_filter.md).
  Gaussian-sum measurement noise is rejected for the same reason.

- y:

  A numeric vector or an `n x d` matrix of observations.

- m:

  Integer; the number of end-of-sample observations to test,
  `1 <= m < nrow(y)`. The tiny-`m` regime (`m = 1, 2, 3`) is the point
  of the test.

- method:

  Either `"chisq"` (parametric) or `"andrews"` (subsampling calibration
  after Andrews' P-test, asymptotic in `n`).

- alpha:

  The nominal level used to set `reject`.

## Value

An object of class `gmm_eos_test`: a list with the `statistic`, the
`p_value`, the logical `reject`, the `method`, `m`, `alpha`, and the
in-sample block statistics used for the subsampling calibration.

## Details

Two calibrations are offered. `method = "chisq"` refers the statistic to
a chi-square distribution on `m * ncol(y)` degrees of freedom, which is
exact when the standardised innovations are Gaussian.
`method = "andrews"` is a subsampling calibration after the P-test of
Andrews (2003): the statistic's rank among the \\n - 2m + 1\\
overlapping `m`-blocks of the earlier innovations gives the p-value
\\(1 + k) / (n - 2m + 2)\\, where \\k\\ is the number of blocks at least
as large as the statistic. It differs from Andrews' test in two ways:
the p-value counts the tested block itself, and the block statistics are
computed at the supplied model instead of at estimates that leave out
each block. It does not assume Gaussian innovations, so it suits
heavy-tailed observation noise, but its validity is asymptotic: `n`
large with `m` fixed, and innovations that are stationary and ergodic
before the tested block. Its size equals the nominal level only when the
end block is exchangeable with the earlier blocks, and in short series
it can over-reject. The model is supplied exactly as for
[`gmm_filter()`](https://max578.github.io/proxymix/reference/gmm_filter.md).

Two finite-sample cautions. The chi-square calibration is exact when the
model is given; with parameters estimated on a short series it
over-rejects (the end-of-sample vignette reports its size by series
length). The subsampling calibration can also over-reject in short
series, so neither method holds its nominal size when the model is
estimated on little data. The subsampling p-value has a floor of
`1 / (n - 2m + 2)`, so it can reject at level 0.05 only when
`n > 2m + 18`.

## References

Andrews, D. W. K. (2003). End-of-Sample Instability Tests.
*Econometrica*, 71(6), 1661–1694.

## See also

[`gmm_filter()`](https://max578.github.io/proxymix/reference/gmm_filter.md)

## Examples

``` r
prior <- gmm(weights = 1, means = list(0), covariances = list(matrix(10)))
dyn   <- list(A = matrix(1), Q = matrix(0.04))
meas  <- list(C = matrix(1), R = matrix(1))
set.seed(1)
y <- c(rnorm(119), 6)                # a stable series with a jump at the end
gmm_eos_test(prior, dyn, meas, y, m = 1L, method = "andrews")
#> End-of-sample instability test (andrews), m = 1
#>   statistic = 29.1038   p-value = 0.0083   REJECT (instability) at alpha = 0.05
```
