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

  Either `"chisq"` (parametric) or `"andrews"` (distribution-free
  subsampling P-test).

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
`method = "andrews"` is the distribution-free Andrews (2003) subsampling
P-test: the statistic's rank among the in-sample overlapping `m`-blocks
of the same innovations gives the p-value, so it stays calibrated when
the innovations are non-Gaussian (heavy-tailed observation noise, say).
The model is supplied exactly as for
[`gmm_filter()`](https://max578.github.io/proxymix/reference/gmm_filter.md).

Two finite-sample cautions. The chi-square calibration is exact when the
model is given; with parameters estimated on a short series it
over-rejects (size 0.068 at `n = 30` against a nominal 0.05 in the
validation study, settling to 0.044 by `n = 120`), so prefer
`method = "andrews"` when the model is estimated on little data. The
subsampling p-value has a floor of `1 / (n - 2m + 2)`, so it can reject
at level 0.05 only when `n > 2m + 18`.

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
#> $statistic
#> [1] 29.10382
#> 
#> $p_value
#> [1] 0.008333333
#> 
#> $reject
#> [1] TRUE
#> 
#> $method
#> [1] "andrews"
#> 
#> $m
#> [1] 1
#> 
#> $alpha
#> [1] 0.05
#> 
#> $df
#> [1] 1
#> 
#> $in_sample_blocks
#>   [1] 3.554750e-02 2.911332e-01 2.622507e-01 2.935641e+00 2.740067e-02
#>   [6] 7.895244e-01 2.170989e-01 3.560326e-01 1.118222e-01 2.736059e-01
#>  [11] 1.474811e+00 5.579593e-04 8.706621e-01 4.854571e+00 1.485184e+00
#>  [16] 3.742538e-03 5.791188e-04 7.207274e-01 3.412181e-01 7.431643e-02
#>  [21] 2.677158e-01 8.993352e-02 1.558138e-01 4.800755e+00 3.212741e-01
#>  [26] 2.176074e-02 4.452547e-02 1.857312e+00 4.746573e-02 4.000792e-01
#>  [31] 1.875224e+00 4.043256e-02 7.793913e-02 2.920177e-02 1.788863e+00
#>  [36] 5.051060e-02 2.733394e-02 2.813809e-02 1.407948e+00 4.448343e-01
#>  [41] 8.603388e-02 1.027993e-01 3.569358e-01 1.312872e-01 6.894557e-01
#>  [46] 4.858162e-01 1.594944e-01 4.797690e-01 5.285386e-02 5.052012e-01
#>  [51] 2.102961e-02 6.326433e-01 4.458853e-02 1.340566e+00 1.878720e+00
#>  [56] 2.617749e+00 6.391465e-01 1.606209e+00 1.785440e-01 8.511620e-02
#>  [61] 4.230079e+00 2.750905e-01 5.297272e-02 1.683986e-01 1.069329e+00
#>  [66] 1.162422e-05 3.265654e+00 2.189644e+00 5.903947e-04 3.412866e+00
#>  [71] 5.214728e-04 1.191428e+00 9.073974e-02 1.325579e+00 1.518131e+00
#>  [76] 1.514492e-01 1.198558e-01 1.407395e-02 2.671228e-02 2.180238e-01
#>  [81] 1.321657e-01 8.941230e-03 1.602585e+00 1.982918e+00 5.822228e-01
#>  [86] 1.511109e-01 9.587286e-01 1.896176e-01 6.426711e-02 1.310650e-02
#>  [91] 4.082348e-01 1.125279e+00 6.820047e-01 6.754606e-02 1.030703e+00
#>  [96] 9.832867e-03 3.034342e+00 6.243341e-01 1.529166e+00 1.108467e-01
#> [101] 1.645748e-01 7.144235e-02 4.141931e-01 1.938593e-01 1.404771e-01
#> [106] 3.552470e+00 3.515232e-01 4.364804e-01 4.234652e-03 1.507889e+00
#> [111] 1.192417e+00 5.428482e-01 1.233330e+00 9.516459e-01 1.582138e-01
#> [116] 2.436178e-01 1.144764e-01 5.765224e-02 2.531887e-01
#> 
#> attr(,"class")
#> [1] "gmm_eos_test"
```
