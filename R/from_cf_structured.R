## Two fits for characteristic functions with known structure; see the help
## pages of from_cf_elliptical() and from_cf_sources().

# ---------------------------------------------------------------------------
# internal helpers: shared
# ---------------------------------------------------------------------------

## Non-negative least squares by the active-set method of Lawson and
## Hanson (1974): minimise |A x - b| subject to x >= 0.
.cf_nnls <- function(A, b) {
  n <- ncol(A)
  x <- numeric(n)
  inside <- rep(FALSE, n)
  tol <- 10 * .Machine$double.eps * max(colSums(abs(A))) * nrow(A)
  ls_on <- function(set) {
    s <- numeric(n)
    co <- qr.coef(qr(A[, set, drop = FALSE], tol = 1e-12), b)
    co[is.na(co)] <- 0
    s[set] <- co
    s
  }
  w <- as.vector(crossprod(A, b))
  it <- 0L
  while (any(!inside) && max(w[!inside]) > tol && it < 3L * n) {
    it <- it + 1L
    inside[which(!inside)[which.max(w[!inside])]] <- TRUE
    s <- ls_on(inside)
    while (min(s[inside]) <= 0) {
      neg <- inside & s <= 0
      step <- min(x[neg] / (x[neg] - s[neg]))
      x <- x + step * (s - x)
      inside <- inside & x > tol
      if (!any(inside)) break
      s <- ls_on(inside)
    }
    x <- s
    w <- as.vector(crossprod(A, b - A %*% x))
  }
  x
}

## Hessian of the log characteristic function at `t0` by central
## differences with step `h[d]` along axis `d`. Complex and symmetric.
.cf_log_hessian <- function(phi, p, t0, h) {
  lp <- function(v) log(phi(matrix(v, 1L, p)))
  e <- function(d) {
    v <- numeric(p)
    v[d] <- h[d]
    v
  }
  f0 <- lp(t0)
  H <- matrix(0i, p, p)
  for (d in seq_len(p)) {
    H[d, d] <- (lp(t0 + e(d)) - 2 * f0 + lp(t0 - e(d))) / h[d]^2
    for (f in seq_len(d - 1L)) {
      H[d, f] <- (
        lp(t0 + e(d) + e(f)) - lp(t0 + e(d) - e(f)) -
          lp(t0 - e(d) + e(f)) + lp(t0 - e(d) - e(f))
      ) / (4 * h[d] * h[f])
      H[f, d] <- H[d, f]
    }
  }
  H
}

## Rows of unit length, uniform on the sphere.
.cf_unit_directions <- function(n, p) {
  D <- matrix(stats::rnorm(n * p), n, p)
  D / sqrt(rowSums(D^2))
}

# ---------------------------------------------------------------------------
# internal helpers: elliptical targets
# ---------------------------------------------------------------------------

## Location from the phase of the characteristic function: a first value
## at a low frequency, where the phase cannot wrap, refined at `th`.
.cf_phase_location <- function(phi, p, th) {
  vapply(seq_len(p), function(d) {
    v <- numeric(p)
    v[d] <- 0.01 * th[d]
    m0 <- Arg(phi(matrix(v, 1L, p))) / v[d]
    v[d] <- th[d]
    m0 + Arg(phi(matrix(v, 1L, p)) * exp(-1i * th[d] * m0)) / th[d]
  }, numeric(1L))
}

## Frequency at which the modulus falls to one half along each row of `D`.
## `guess` is a first value for each row.
.cf_half_frequency <- function(phi, D, guess) {
  vapply(seq_len(nrow(D)), function(i) {
    s <- guess[i] * 10^seq(-1.5, 1.5, length.out = 61L)
    below <- which(Mod(phi(outer(s, D[i, ]))) < 0.5)[1L]
    if (is.na(below) || below == 1L) return(NA_real_)
    s <- exp(seq(log(s[below - 1L]), log(s[below]), length.out = 41L))
    m <- Mod(phi(outer(s, D[i, ])))
    below <- which(m < 0.5)[1L]
    if (is.na(below) || below == 1L) return(NA_real_)
    lo <- below - 1L
    exp(
      log(s[lo]) + (log(s[below]) - log(s[lo])) * (m[lo] - 0.5) /
        (m[lo] - m[below])
    )
  }, numeric(1L))
}

## Symmetric matrix S with a' S a closest to `y` over the rows a of `D`.
.cf_quadratic_form <- function(D, y) {
  p <- ncol(D)
  ix <- which(lower.tri(diag(p), diag = TRUE), arr.ind = TRUE)
  X <- vapply(seq_len(nrow(ix)), function(k) {
    (if (ix[k, 1L] == ix[k, 2L]) 1 else 2) * D[, ix[k, 1L]] * D[, ix[k, 2L]]
  }, numeric(nrow(D)))
  S <- matrix(0, p, p)
  S[ix] <- qr.solve(X, y)
  S + t(S) - diag(diag(S), p)
}

## Tail probability constant of a symmetric stable law of index `alpha` and
## unit scale: P(X > x) is close to this constant times x^(-alpha).
.cf_stable_tail_constant <- function(alpha) {
  sin(pi * alpha / 2) * gamma(alpha) / pi
}

.check_cf_elliptical_arguments <- function(cf, n_dim, N, scale_range) {
  if (!is.function(cf)) {
    cli::cli_abort("{.arg cf} must be a function of a matrix of frequencies.")
  }
  count <- function(x) {
    is.numeric(x) && length(x) == 1L && !is.na(x) && x >= 1 && x == round(x)
  }
  if (!count(n_dim)) {
    cli::cli_abort("{.arg n_dim} must be a positive integer scalar.")
  }
  if (!count(N) || N < 2L) {
    cli::cli_abort("{.arg N} must be an integer of at least 2.")
  }
  if (!is.null(scale_range) &&
      (!is.numeric(scale_range) || length(scale_range) != 2L ||
         anyNA(scale_range) || scale_range[1L] <= 0 ||
         scale_range[1L] >= scale_range[2L])) {
    cli::cli_abort(
      "{.arg scale_range} must be two increasing positive numbers or {.code NULL}."
    )
  }
  invisible(NULL)
}

# ---------------------------------------------------------------------------
# from_cf_elliptical
# ---------------------------------------------------------------------------

#' Fit concentric normals to an elliptical characteristic function
#'
#' Fits a Gaussian mixture to a symmetric distribution whose characteristic
#' function has the form
#' \eqn{\varphi(t) = \exp(i t^\top \mu)\, \rho(t^\top \Sigma t)}, where
#' \eqn{\rho(u)} is a mixture of \eqn{\exp(-s u / 2)} over \eqn{s > 0} (a
#' Gaussian scale mixture). Elliptical stable laws, the multivariate Cauchy,
#' Student and Laplace distributions are of this form. Every component has
#' mean \eqn{\mu} and a covariance proportional to \eqn{\Sigma}, and the fit
#' reaches further into heavy tails than [from_cf()] does.
#'
#' The centre \eqn{\mu} is read from the phase of `cf`. The matrix
#' \eqn{\Sigma} is found by least squares from the frequency at which
#' \eqn{|\varphi|} falls to one half along each of a set of random
#' directions. The standard deviations of the components are `N` values on
#' a geometric grid. Their weights are found by non-negative least squares
#' on the relative error of \eqn{1 - \rho}, over frequencies from very small
#' to well past the half-height. The smallest frequencies fix the tails.
#' Grid values that receive no weight are dropped, so the result can have
#' fewer than `N` components.
#'
#' `diagnostics$shape_ratio` should be 1 for an elliptical target. It is the
#' smallest and largest value, over 200 random directions, of fitted scale
#' divided by target scale. A warning of class `proxymix_cf_not_elliptical`
#' is given if either value is outside 0.9 to 1.1, or if the imaginary part
#' of the centred characteristic function exceeds 0.001 of its modulus.
#'
#' @param cf,n_dim As in [from_cf()].
#' @param N Number of values on the grid of component standard deviations.
#'   Default `20L`.
#' @param scale_range Smallest and largest multiplier \eqn{s}. A component
#'   has covariance \eqn{s^2 \Sigma}, and \eqn{s = 1} is the normal whose
#'   half-height frequency matches the target's. The default `NULL` reads the
#'   tail index from `cf` and runs from 0.05 to 30 times the 0.999 quantile
#'   of a stable law with that index (to 20 when the index is 1.9 or more),
#'   both in units of that law's scale.
#' @param seed Optional integer seed for the random directions. When
#'   supplied, the random number generator state is restored on exit.
#' @param bracket Logical. Compute Kolmogorov bounds against `cf` along
#'   eight directions, as in [from_cf()]. Default `TRUE`.
#'
#' @returns A [gmm_fit] with `regime` set to `"cf"`. Its `diagnostics` hold
#'   `shape` (the matrix \eqn{\Sigma}), `shape_ratio`, `asymmetry` (the
#'   largest imaginary part of the centred characteristic function, relative
#'   to its modulus), `scales` (all `N` grid values, including those
#'   dropped),
#'   `tail_index`, `radial_error` (the largest absolute difference between
#'   \eqn{\rho} and its fit on the frequencies used), `heldout_rms` and
#'   `heldout_max`, and `bracket` and `bracket_directions`, as in
#'   [from_cf()].
#' @references
#' Lawson, C. L. and Hanson, R. J. (1974). *Solving Least Squares
#' Problems*. Prentice-Hall.
#' @family fitting
#' @seealso [cf_target()], [from_cf()], [from_cf_sources()],
#'   [gmm_cf_distance()].
#' @export
#' @examples
#' \dontshow{if (requireNamespace("RhpcBLASctl", quietly = TRUE)) {
#'   blas_threads <- RhpcBLASctl::blas_get_num_procs()
#'   RhpcBLASctl::blas_set_num_threads(1L)
#' }}
#' ## A bivariate Cauchy distribution.
#' S <- matrix(c(2, 0.5, 0.5, 1), 2L)
#' cf <- function(t) exp(1i * t[, 1] - sqrt(rowSums((t %*% S) * t)))
#' fit <- from_cf_elliptical(cf, n_dim = 2L, seed = 1L)
#' ## the print shows the Kolmogorov bounds along eight directions
#' fit
#' fit@diagnostics$shape_ratio
#' ## the 0.99 quantile of the first coordinate against the exact value
#' qgmm(0.99, gmm_marginalise(fit, 1L))
#' 1 + sqrt(2) * qcauchy(0.99)
#' \dontshow{if (exists("blas_threads")) {
#'   RhpcBLASctl::blas_set_num_threads(blas_threads)
#'   rm(blas_threads)
#' }}
from_cf_elliptical <- function(
  cf,
  n_dim,
  N = 20L,
  scale_range = NULL,
  seed = NULL,
  bracket = TRUE
) {
  .check_cf_elliptical_arguments(cf, n_dim, N, scale_range)
  .check_cf_bracket_argument(bracket)
  p <- as.integer(n_dim)
  K <- as.integer(N)
  phi <- .cf_checked(cf, p)

  run <- function() {
    # Centre -------------------------------------------------------------------
    th <- .cf_axis_scales(phi, p)
    .check_cf_properties(phi, p, th)
    mu <- .cf_phase_location(phi, p, th)
    phi0 <- function(Tm) exp(-1i * as.vector(Tm %*% mu)) * phi(Tm)

    # Shape matrix -------------------------------------------------------------
    n_dir <- max(200L, 2L * p * (p + 1L))
    half <- function(D) {
      s <- .cf_half_frequency(phi0, D, 1 / sqrt(rowSums((D / rep(th, each = nrow(D)))^2)))
      if (anyNA(s)) {
        cli::cli_abort(
          "The modulus of {.arg cf} does not fall to one half along every direction.",
          class = "proxymix_invalid_cf"
        )
      }
      s
    }
    D <- .cf_unit_directions(n_dir, p)
    Sigma <- .cf_quadratic_form(D, (.cf_half_height / half(D))^2)
    ev <- eigen(Sigma, symmetric = TRUE, only.values = TRUE)$values
    if (min(ev) <= 1e-10 * max(ev)) {
      cli::cli_abort(
        c(
          "The shape matrix read from {.arg cf} is not positive definite.",
          "i" = "The target is not close to elliptical. Use {.fn from_cf}."
        ),
        class = "proxymix_cf_not_elliptical"
      )
    }
    D2 <- .cf_unit_directions(200L, p)
    scale_fit <- sqrt(rowSums((D2 %*% Sigma) * D2))
    shape_ratio <- range(scale_fit * half(D2) / .cf_half_height)

    # Radial profile -----------------------------------------------------------
    ## rows of `R` have unit length in the metric of Sigma
    R <- D2[seq_len(8L), , drop = FALSE] / scale_fit[seq_len(8L)]
    r <- .cf_half_height * 10^seq(-9, 1.5, length.out = 900L)
    val <- vapply(seq_len(nrow(R)), function(i) phi0(outer(r, R[i, ])),
                  complex(length(r)))
    asymmetry <- max(abs(Im(val)) / pmax(Mod(val), 1e-12))
    y <- 1 - rowMeans(Re(val))
    ## below 1e-7 the subtraction has lost too many digits
    keep <- y > 1e-7
    r <- r[keep]
    y <- y[keep]
    tail_index <- unname(stats::coef(stats::lm.fit(
      cbind(1, log(r[seq_len(20L)])), log(y[seq_len(20L)])
    ))[2L])

    # Scales and weights -------------------------------------------------------
    alpha <- min(tail_index, 2)
    unit <- log(2)^(1 / alpha) / .cf_half_height
    if (is.null(scale_range)) {
      reach <- if (alpha >= 1.9) {
        20
      } else {
        30 * (.cf_stable_tail_constant(alpha) / 1e-3)^(1 / alpha)
      }
      scale_range <- unit * c(0.05, reach)
    }
    sds <- exp(seq(log(scale_range[1L]), log(scale_range[2L]), length.out = K))
    B <- -expm1(-0.5 * outer(r^2, sds^2))
    w <- .cf_nnls(rbind(B / y, rep(1e3, K)), c(rep(1, length(r)), 1e3))
    if (sum(w) <= 0) {
      cli::cli_abort(
        "No weights could be fitted; widen {.arg scale_range}.",
        class = "proxymix_invalid_cf"
      )
    }
    w <- w / sum(w)
    radial_error <- max(abs(y - as.vector(B %*% w)))
    used <- w > 1e-8
    w <- w[used] / sum(w[used])

    # Error on fresh frequencies -----------------------------------------------
    bw <- c(0.3, 0.6, 1.2, 2.5, 5)
    fresh <- .cf_frequencies(2048L, p, bw)$U %*% solve(t(chol(Sigma)))
    quad <- rowSums((fresh %*% Sigma) * fresh)
    gap <- Mod(phi0(fresh) - as.vector(exp(-0.5 * outer(quad, sds[used]^2)) %*% w))
    list(
      mu = mu, Sigma = Sigma, w = w, sds = sds, used = used,
      shape_ratio = shape_ratio, asymmetry = asymmetry,
      tail_index = tail_index, radial_error = radial_error,
      heldout_rms = sqrt(mean(gap^2)), heldout_max = max(gap)
    )
  }
  res <- if (is.null(seed)) run() else withr::with_seed(seed, run())

  if (res$shape_ratio[1L] < 0.9 || res$shape_ratio[2L] > 1.1 ||
      res$asymmetry > 1e-3) {
    cli::cli_warn(
      c(
        "The target does not look elliptical.",
        "i" = "The fitted scale is between {signif(res$shape_ratio[1L], 3)} and {signif(res$shape_ratio[2L], 3)} times the target's, depending on direction.",
        "i" = "After centring, the imaginary part of the characteristic function is {signif(res$asymmetry, 2)} of its modulus.",
        "i" = "{.fn from_cf} makes no assumption on the shape."
      ),
      class = "proxymix_cf_not_elliptical"
    )
  }

  sds <- res$sds[res$used]
  means <- rep(list(res$mu), length(sds))
  covs <- lapply(sds, function(s) s^2 * res$Sigma)
  kb <- .cf_bracket(
    gmm(weights = res$w, means = means, covariances = covs), phi, p,
    bracket, seed
  )
  fit <- gmm_fit(
    weights = res$w,
    means = means,
    covariances = covs,
    target = NULL,
    regime = "cf",
    diagnostics = list(
      shape = res$Sigma,
      shape_ratio = res$shape_ratio,
      asymmetry = res$asymmetry,
      scales = res$sds,
      tail_index = res$tail_index,
      radial_error = res$radial_error,
      heldout_rms = res$heldout_rms,
      heldout_max = res$heldout_max,
      bracket = kb$table,
      bracket_directions = kb$directions
    ),
    converged = TRUE,
    iterations = NA_integer_,
    call = match.call(),
    name = sprintf("from_cf_elliptical[N=%d]", length(sds)),
    metadata = list(quality = list(
      regime = "cf",
      converged = TRUE,
      degenerate = FALSE,
      ess = NA_real_,
      ess_relative = NA_real_,
      min_component_ess = NA_real_,
      max_weight = NA_real_,
      support_fraction = NA_real_,
      kld_final = NA_real_,
      kld_approx = NA_real_,
      validation_gap = NA_real_,
      heldout_kld = NA_real_,
      cf_heldout_rms = res$heldout_rms,
      cf_distance = kb$distance,
      cf_upper = kb$upper
    ))
  )
  gmm_canonicalise(fit)
}

# ---------------------------------------------------------------------------
# internal helpers: independent sources
# ---------------------------------------------------------------------------

## Mixing directions of x = A s + Gaussian noise with independent sources
## s. The Hessian of the log characteristic function at t is
## A D(t) A' minus the noise covariance, with D(t) diagonal, so differences
## of Hessians are diagonalised together by A (Yeredor, 2000). `pts` holds
## three frequencies as rows: two give the directions and the third checks
## them.
.cf_mixing_directions <- function(phi, p, pts, h) {
  H0 <- .cf_log_hessian(phi, p, numeric(p), h)
  M <- lapply(seq_len(3L), function(i) {
    .cf_log_hessian(phi, p, pts[i, ], h) - H0
  })
  if (max(Mod(M[[2L]])) < 1e-4 * max(Mod(H0)) || rcond(M[[2L]]) < 1e-8) {
    cli::cli_abort(
      c(
        "The mixing directions cannot be read from {.arg cf}.",
        "i" = "This happens when a source is Gaussian or there are fewer independent sources than dimensions."
      ),
      class = "proxymix_cf_not_sources"
    )
  }
  V <- eigen(M[[1L]] %*% solve(M[[2L]]))$vectors
  for (j in seq_len(p)) V[, j] <- V[, j] / V[which.max(Mod(V[, j])), j]
  A <- Re(V)
  B <- solve(A)
  check <- B %*% M[[3L]] %*% t(B)
  list(
    A = A,
    imaginary = max(abs(Im(V))),
    off_diagonal = max(Mod(check[row(check) != col(check)])) /
      max(Mod(diag(check)))
  )
}

## Gauss-Hermite nodes and weights for the standard normal density.
.cf_normal_nodes <- function(n) {
  J <- matrix(0, n, n)
  J[cbind(seq_len(n - 1L), seq_len(n - 1L) + 1L)] <- sqrt(seq_len(n - 1L))
  e <- eigen(J + t(J), symmetric = TRUE)
  list(x = e$values, w = e$vectors[1L, ]^2)
}

## Largest pairwise dependence of z = B x (Eriksson and Koivunen, 2003):
## the root of the integral of |joint cf - product of marginal cfs|^2
## against the N(0, I) density, over pairs of coordinates. Gaussian noise
## multiplies the product by exp(-u_j u_k noise_cov[j, k]), so the value is
## 0 for x = A s + noise with independent s and B the inverse of A.
.cf_pair_dependence <- function(phi, B, noise_cov, n_node = 32L) {
  p <- nrow(B)
  gh <- .cf_normal_nodes(n_node)
  marginal <- vapply(
    seq_len(p), function(j) phi(outer(gh$x, B[j, ])), complex(n_node)
  )
  g <- expand.grid(a = seq_len(n_node), b = seq_len(n_node))
  w <- gh$w[g$a] * gh$w[g$b]
  ua <- gh$x[g$a]
  ub <- gh$x[g$b]
  worst <- 0
  for (j in seq_len(p - 1L)) {
    for (k in (j + 1L):p) {
      joint <- phi(outer(ua, B[j, ]) + outer(ub, B[k, ]))
      product <- marginal[g$a, j] * marginal[g$b, k] *
        exp(-ua * ub * noise_cov[j, k])
      worst <- max(worst, sum(w * Mod(joint - product)^2))
    }
  }
  sqrt(worst)
}

.check_cf_sources_arguments <- function(cf, n_dim, N) {
  if (!is.function(cf)) {
    cli::cli_abort("{.arg cf} must be a function of a matrix of frequencies.")
  }
  count <- function(x) {
    is.numeric(x) && length(x) == 1L && !is.na(x) && x >= 1 && x == round(x)
  }
  if (!count(n_dim) || n_dim < 2L) {
    cli::cli_abort("{.arg n_dim} must be an integer of at least 2.")
  }
  if (!count(N)) {
    cli::cli_abort("{.arg N} must be a positive integer scalar.")
  }
  invisible(NULL)
}

# ---------------------------------------------------------------------------
# from_cf_sources
# ---------------------------------------------------------------------------

#' Fit a characteristic function that is a linear mix of independent sources
#'
#' Fits a distribution of the form \eqn{X = A S + G}, where the coordinates
#' of \eqn{S} are independent and not Gaussian, \eqn{A} is an invertible
#' matrix and \eqn{G} is Gaussian noise or absent. Only the characteristic
#' function of \eqn{X} is needed. The mixing matrix is read from `cf`, and
#' each source is then fitted in one dimension with [from_cf()]. A joint fit
#' in `n_dim` dimensions would be less accurate.
#'
#' The fitted family, independent sources that are each a one-dimensional
#' Gaussian mixture, is the independent factor model of Attias (1999), who
#' fits it to a sample by maximum likelihood.
#'
#' The result is a list, not a [gmm]. As a single mixture it would need
#' `N^n_dim` components. It has three functions: `cf()`, `draw()` (sampling
#' from the fit) and `project()` (the one-dimensional mixture along a
#' direction). `draw()` gives an error when the noise covariance estimated
#' from `cf` is too far from positive semi-definite to be sampled.
#'
#' The mixing matrix is found from second derivatives of the logarithm of
#' `cf` at two random frequencies (Yeredor, 2000). The coordinates it
#' separates are then tested for pairwise independence with the criterion of
#' Eriksson and Koivunen (2003), which compares the joint characteristic
#' function of each pair with the product of the two marginal ones at all
#' frequencies. If the criterion is above 0.01 the frequencies are drawn
#' again, up to four times, and the directions with the smallest value are
#' kept. A draw at which the second derivatives are too ill-conditioned to
#' give directions is drawn again in the same way; the function stops only
#' when all five draws fail.
#'
#' The target must have a finite, positive definite covariance. A warning of class
#' `proxymix_cf_not_sources` is raised when the criterion is still above
#' 0.01, or when the directions do not separate the sources at a third
#' frequency. Either is the sign that the target is not of the form above.
#'
#' @param cf,n_dim As in [from_cf()]. `n_dim` must be at least 2.
#' @param N Number of mixture components for each source.
#' @param n_freq,bandwidths,max_iter,bracket Passed to [from_cf()] for each
#'   source. With `bracket = TRUE` each source fit carries Kolmogorov
#'   bounds for its coordinate of \eqn{Z = \hat A^{-1} X}, not for \eqn{X}.
#'   For a direction of \eqn{X}, apply [gmm_cf_distance()] to the result of
#'   `project()`.
#' @param seed Optional integer seed. When supplied, the random number
#'   generator state is restored on exit.
#'
#' @returns A list of class `proxymix_cf_sources` with
#'   * `mixing`, the matrix \eqn{\hat A}. Its columns are the mixing
#'     directions, scaled so that each coordinate of
#'     \eqn{Z = \hat A^{-1} X} has unit variance. Column order and sign are
#'     not identified.
#'   * `sources`, a list of `n_dim` one-dimensional [gmm_fit] objects, the
#'     fits to the coordinates of \eqn{Z}.
#'   * `noise_cov`, the off-diagonal covariance of \eqn{Z}, with the
#'     diagonal set to zero. It is zero when there is no noise. The noise
#'     variance on the diagonal stays inside `sources`.
#'   * `cf(t)`, the characteristic function of the fit at the rows of `t`.
#'   * `draw(n)`, an `n` by `n_dim` matrix of draws from the fit.
#'   * `project(direction, max_components = 1e5)`, the [gmm] of
#'     \eqn{a^\top X} under the fit for a direction \eqn{a}. Sources whose
#'     coefficient is below `1e-3` of the largest are replaced by a single
#'     normal with the same mean and variance. It gives an error if the
#'     projection needs more than `max_components` components.
#'   * `diagnostics`, with `off_diagonal` (the largest off-diagonal entry
#'     divided by the largest diagonal entry of the Hessian difference at the
#'     third frequency, transformed by the inverse of `mixing`, 0 if the sources are separated
#'     exactly, with a warning above 0.05), `dependence` (the largest over
#'     pairs of coordinates of \eqn{Z} of the root of the integral, against
#'     the standard normal density in two dimensions, of the squared modulus
#'     of the joint characteristic function minus the product of the
#'     marginal ones, with the Gaussian noise term allowed for; 0 for
#'     independent sources, with a warning above 0.01), `imaginary` (the largest imaginary part of the
#'     directions before it is dropped) and `heldout_rms` (the root mean
#'     square of the difference between `cf` and the fit's characteristic
#'     function on fresh frequencies).
#' @references
#' Attias, H. (1999). Independent factor analysis. *Neural Computation*,
#' 11(4), 803--851. \doi{10.1162/089976699300016458}
#'
#' Eriksson, J., and Koivunen, V. (2003). Characteristic-function-based
#' independent component analysis. *Signal Processing*, 83(10), 2195--2208.
#' \doi{10.1016/S0165-1684(03)00162-2}
#'
#' Yeredor, A. (2000). Blind source separation via the second
#' characteristic function. *Signal Processing*, 80(5), 897--902.
#' \doi{10.1016/S0165-1684(00)00062-1}
#' @family fitting
#' @seealso [cf_target()], [from_cf()], [from_cf_elliptical()],
#'   [gmm_cf_distance()].
#' @export
#' @examples
#' \dontshow{if (requireNamespace("RhpcBLASctl", quietly = TRUE)) {
#'   blas_threads <- RhpcBLASctl::blas_get_num_procs()
#'   RhpcBLASctl::blas_set_num_threads(1L)
#' }}
#' ## Two independent Gamma(3) sources, mixed, with N(0, 0.25 I) noise.
#' A <- matrix(c(0.7, 0.3, 0, -0.7), 2L)
#' cf <- function(t) {
#'   u <- t %*% A
#'   (1 - 1i * u[, 1])^(-3) * (1 - 1i * u[, 2])^(-3) *
#'     exp(-0.125 * rowSums(t^2))
#' }
#' fit <- from_cf_sources(cf, n_dim = 2L, N = 3L, n_freq = 1024L, seed = 1L)
#' fit$diagnostics$off_diagonal
#' first <- fit$project(c(1, 0))
#' gmm_cf_distance(first, function(s) cf(cbind(s[, 1], 0)))$distance
#' \dontshow{if (exists("blas_threads")) {
#'   RhpcBLASctl::blas_set_num_threads(blas_threads)
#'   rm(blas_threads)
#' }}
from_cf_sources <- function(
  cf,
  n_dim,
  N,
  n_freq = 8192L,
  bandwidths = c(0.3, 0.6, 1.2, 2.5, 5),
  max_iter = 300L,
  seed = NULL,
  bracket = TRUE
) {
  .check_cf_sources_arguments(cf, n_dim, N)
  .check_cf_bracket_argument(bracket)
  .check_cf_arguments(cf, 1L, N, n_freq, bandwidths, 1L, max_iter)
  p <- as.integer(n_dim)
  K <- as.integer(N)
  phi <- .cf_checked(cf, p)

  run <- function() {
    # Scale and moments --------------------------------------------------------
    th <- .cf_axis_scales(phi, p)
    .check_cf_properties(phi, p, th)
    co <- .cf_standardise(phi, p, th, min(.cf_tail_index(phi, p, th)) < 1.9)
    if (is.null(co$cov)) {
      cli::cli_abort(
        "{.fn from_cf_sources} needs a target with a finite, positive definite covariance.",
        class = "proxymix_invalid_cf"
      )
    }

    # Mixing directions --------------------------------------------------------
    directions <- function() {
      pts <- (0.8 * matrix(stats::rnorm(3L * p), 3L, p) / sqrt(p)) %*%
        solve(co$A)
      dirs <- .cf_mixing_directions(phi, p, pts, 0.02 * th)
      B <- solve(dirs$A)
      sd_z <- sqrt(diag(B %*% co$cov %*% t(B)))
      B <- B / sd_z
      noise_cov <- B %*% co$cov %*% t(B)
      diag(noise_cov) <- 0
      c(dirs, list(
        B = B, A_hat = dirs$A * rep(sd_z, each = p), noise_cov = noise_cov,
        dependence = .cf_pair_dependence(phi, B, noise_cov)
      ))
    }
    ## Up to five draws of the frequencies, while the separated coordinates
    ## are still dependent or a draw cannot separate them.
    dirs <- NULL
    for (i in seq_len(5L)) {
      if (!is.null(dirs) && dirs$dependence <= 0.01) break
      again <- tryCatch(
        directions(), proxymix_cf_not_sources = function(e) e
      )
      if (inherits(again, "proxymix_cf_not_sources")) {
        failed <- again
      } else if (is.null(dirs) || again$dependence < dirs$dependence) {
        dirs <- again
      }
    }
    if (is.null(dirs)) stop(failed)
    B <- dirs$B
    A_hat <- dirs$A_hat
    noise_cov <- dirs$noise_cov

    # One-dimensional fits -----------------------------------------------------
    sources <- lapply(seq_len(p), function(j) {
      from_cf(
        function(t) phi(outer(t[, 1L], B[j, ])), 1L, K,
        n_freq = n_freq, bandwidths = bandwidths, max_iter = max_iter,
        bracket = bracket
      )
    })
    list(
      A = A_hat, sources = sources, noise_cov = noise_cov, cov = co$cov,
      imaginary = dirs$imaginary, off_diagonal = dirs$off_diagonal,
      dependence = dirs$dependence,
      fresh = .cf_frequencies(2048L, p, bandwidths)$U %*% solve(co$A)
    )
  }
  res <- if (is.null(seed)) run() else withr::with_seed(seed, run())

  if (res$off_diagonal > 0.05 || res$dependence > 0.01) {
    cli::cli_warn(
      c(
        "The target does not look like a linear mix of independent sources.",
        "i" = "At a third test frequency, off-diagonal entries are {signif(100 * res$off_diagonal, 2)}% of diagonal ones; for independent sources they would be 0.",
        "i" = "The pairwise dependence of the separated coordinates is {signif(res$dependence, 2)}; for independent sources it would be 0.",
        "i" = "{.fn from_cf} makes no assumption on the structure."
      ),
      class = "proxymix_cf_not_sources"
    )
  }

  A <- res$A
  sources <- res$sources
  noise_cov <- res$noise_cov
  w_l <- lapply(sources, function(f) f@weights)
  m_l <- lapply(sources, function(f) unlist(f@means))
  v_l <- lapply(sources, function(f) unlist(f@covariances))

  # The fit as functions --------------------------------------------------------
  cf_fit <- function(t) {
    if (is.null(dim(t))) t <- matrix(t, nrow = 1L)
    U <- t %*% A
    out <- exp(-0.5 * rowSums((U %*% noise_cov) * U))
    for (j in seq_len(p)) out <- out * gmm_cf(sources[[j]], U[, j])
    out
  }
  ## noise variance on the diagonal sits inside each source's components;
  ## moving the smallest of them into the noise makes its matrix sampleable
  spare <- vapply(v_l, min, numeric(1L))
  noise_eigen <- eigen(noise_cov + diag(spare, p), symmetric = TRUE)
  draw <- function(n) {
    if (min(noise_eigen$values) < -1e-8) {
      cli::cli_abort(
        c(
          "The fit cannot be sampled: the noise covariance estimated from {.arg cf} is not positive semi-definite.",
          "i" = "Its smallest eigenvalue is {signif(min(noise_eigen$values), 3)}."
        )
      )
    }
    Z <- matrix(0, n, p)
    for (j in seq_len(p)) {
      k <- sample.int(length(w_l[[j]]), n, replace = TRUE, prob = w_l[[j]])
      Z[, j] <- m_l[[j]][k] +
        sqrt(v_l[[j]][k] - spare[j]) * stats::rnorm(n)
    }
    Z <- Z + matrix(stats::rnorm(n * p), n, p) %*%
      (t(noise_eigen$vectors) * sqrt(pmax(noise_eigen$values, 0)))
    Z %*% t(A)
  }
  project <- function(direction, max_components = 1e5) {
    if (!is.numeric(direction) || length(direction) != p ||
        any(!is.finite(direction)) || all(direction == 0)) {
      cli::cli_abort(
        "{.arg direction} must be {p} finite numbers, not all zero."
      )
    }
    co <- as.vector(direction %*% A)
    big <- which(abs(co) > 1e-3 * max(abs(co)))
    small <- setdiff(seq_len(p), big)
    n_comp <- prod(lengths(w_l[big]))
    if (n_comp > max_components) {
      cli::cli_abort(
        "The projection needs {n_comp} components, more than {.arg max_components}."
      )
    }
    idx <- as.matrix(expand.grid(lapply(w_l[big], seq_along)))
    w <- rep(1, nrow(idx))
    m <- rep(0, nrow(idx))
    v <- rep(0, nrow(idx))
    for (i in seq_along(big)) {
      j <- big[i]
      w <- w * w_l[[j]][idx[, i]]
      m <- m + co[j] * m_l[[j]][idx[, i]]
      v <- v + co[j]^2 * v_l[[j]][idx[, i]]
    }
    for (j in small) {
      mean_j <- sum(w_l[[j]] * m_l[[j]])
      m <- m + co[j] * mean_j
      v <- v + co[j]^2 * sum(w_l[[j]] * (v_l[[j]] + (m_l[[j]] - mean_j)^2))
    }
    v <- v + sum(outer(co, co) * noise_cov)
    if (min(v) <= 0) {
      cli::cli_abort(
        "A component of the projection has zero or negative variance."
      )
    }
    gmm(weights = w / sum(w), means = as.list(m),
        covariances = lapply(v, matrix))
  }

  structure(
    list(
      mixing = A,
      sources = sources,
      noise_cov = noise_cov,
      cf = cf_fit,
      draw = draw,
      project = project,
      diagnostics = list(
        off_diagonal = res$off_diagonal,
        dependence = res$dependence,
        imaginary = res$imaginary,
        heldout_rms = sqrt(mean(Mod(phi(res$fresh) - cf_fit(res$fresh))^2))
      )
    ),
    class = "proxymix_cf_sources"
  )
}

#' @export
print.proxymix_cf_sources <- function(x, ...) {
  p <- nrow(x$mixing)
  cli::cli_text(
    "Fit to a linear mix of {p} independent sources, {length(x$sources[[1L]]@weights)} component{?s} each"
  )
  cli::cli_text(
    "Root mean square error of the characteristic function on fresh frequencies: {signif(x$diagnostics$heldout_rms, 3)}"
  )
  cli::cli_text("Elements: {.field {names(x)}}")
  invisible(x)
}
