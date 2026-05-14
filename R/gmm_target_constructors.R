## gmm_target constructors and the three built-in toy targets used in the
## vignettes and the test surface.

#' Build a target from samples alone
#'
#' Wraps a numeric matrix of i.i.d. samples as a [gmm_target]. The resulting
#' target carries no `log_density`, so it can only feed regimes `"moment"`
#' (via empirical moments) and `"sample"` (classical EM).
#'
#' @param samples An `n` by `p` numeric matrix.
#' @param name Optional human-readable name. Defaults to
#'   `"target_from_samples"`.
#'
#' @returns A [gmm_target] object.
#' @family targets
#' @export
#' @examples
#' x <- matrix(stats::rnorm(200), ncol = 2)
#' tgt <- gmm_target_from_samples(x)
#' tgt
gmm_target_from_samples <- function(samples, name = "target_from_samples") {
  if (!is.matrix(samples) || !is.numeric(samples)) {
    cli::cli_abort("`samples` must be a numeric matrix.")
  }
  gmm_target(
    n_dim = ncol(samples),
    samples = samples,
    name = name
  )
}

#' Banana-shaped 2-D target
#'
#' A 2-D "banana" density obtained by warping an isotropic Gaussian
#' through the map \eqn{(z_1, z_2) \mapsto (z_1, z_2 + \tfrac{1}{2}(z_1^2 - 1))}.
#' The map has unit Jacobian, so the resulting density is exactly normalised.
#'
#' @param with_samples If `TRUE`, attach `n` exact samples drawn by the
#'   change-of-variables trick. Default `FALSE` — the target then exposes
#'   only its `log_density`, which is the regime-(iii) use case.
#' @param n Number of samples to attach when `with_samples = TRUE`.
#' @param seed Optional integer seed used when drawing the samples.
#'
#' @returns A [gmm_target] in dimension 2.
#' @family targets
#' @export
#' @examples
#' b <- banana_target()
#' b
#' b@log_density(matrix(c(0, 0, 1, 0), ncol = 2, byrow = TRUE))
banana_target <- function(with_samples = FALSE, n = 2000L, seed = 1L) {
  log_density <- function(x) {
    if (is.null(dim(x))) x <- matrix(x, nrow = 1L)
    z1 <- x[, 1L]
    z2 <- x[, 2L] - 0.5 * (z1^2 - 1)
    -0.5 * (z1^2 + z2^2) - log(2 * pi)
  }

  samples <- NULL
  if (isTRUE(with_samples)) {
    samples <- withr::with_seed(seed, {
      z <- matrix(stats::rnorm(n * 2L), ncol = 2L)
      cbind(z[, 1L], z[, 2L] + 0.5 * (z[, 1L]^2 - 1))
    })
  }

  gmm_target(
    n_dim = 2L,
    log_density = log_density,
    samples = samples,
    normalised = TRUE,
    log_normalizer = 0,
    name = "banana",
    metadata = list(family = "banana")
  )
}

#' Donut-shaped 2-D target
#'
#' A rotationally symmetric annulus on \eqn{\mathbb{R}^2}, with density
#' \deqn{f(x) \propto \exp\!\left(-\tfrac{(\Vert x \Vert - r_0)^2}{2 \sigma^2}\right).}
#' Numerical integration in polar coordinates fixes the normaliser; the
#' returned target exposes a normalised `log_density`.
#'
#' @param r0 Centre radius of the annulus.
#' @param sigma Annulus width.
#' @param with_samples If `TRUE`, attach `n` exact samples via polar
#'   change-of-variables and a one-dimensional rejection step.
#' @param n Number of samples to attach when `with_samples = TRUE`.
#' @param seed Optional integer seed used when drawing the samples.
#'
#' @returns A [gmm_target] in dimension 2.
#' @family targets
#' @export
#' @examples
#' d <- donut_target()
#' d
donut_target <- function(r0 = 2.5, sigma = 0.5,
                         with_samples = FALSE, n = 2000L, seed = 1L) {
  ## Normaliser: Z = integral over R^2 of exp(-(|x| - r0)^2 / (2 sigma^2)) dx
  ##   = 2 pi * integral_0^Inf r * exp(-(r - r0)^2 / (2 sigma^2)) dr
  integrand <- function(r) r * exp(-(r - r0)^2 / (2 * sigma^2))
  z_int <- stats::integrate(integrand, lower = 0, upper = r0 + 10 * sigma,
                            rel.tol = 1e-10)$value
  log_z <- log(2 * pi) + log(z_int)

  log_density <- function(x) {
    if (is.null(dim(x))) x <- matrix(x, nrow = 1L)
    r <- sqrt(rowSums(x^2))
    -((r - r0)^2) / (2 * sigma^2) - log_z
  }

  samples <- NULL
  if (isTRUE(with_samples)) {
    samples <- withr::with_seed(seed, {
      ## Sample radius via rejection on the marginal r * exp(-(r-r0)^2/(2 sigma^2)).
      r_max <- r0 + 8 * sigma
      r_lo <- max(0, r0 - 8 * sigma)
      m <- r_max * exp(0) ## envelope
      r_out <- numeric(0L)
      while (length(r_out) < n) {
        candidates <- stats::runif(2 * n, min = r_lo, max = r_max)
        u <- stats::runif(length(candidates))
        keep <- u <
          candidates * exp(-(candidates - r0)^2 / (2 * sigma^2)) / m
        r_out <- c(r_out, candidates[keep])
      }
      r_out <- r_out[seq_len(n)]
      theta <- stats::runif(n, min = 0, max = 2 * pi)
      cbind(r_out * cos(theta), r_out * sin(theta))
    })
  }

  gmm_target(
    n_dim = 2L,
    log_density = log_density,
    samples = samples,
    normalised = TRUE,
    log_normalizer = 0,
    name = "donut",
    metadata = list(family = "donut", r0 = r0, sigma = sigma, log_z = log_z)
  )
}

#' Three-component Gaussian-mixture target
#'
#' A toy three-component planar mixture target where everything is known
#' exactly: the `log_density` matches the GMM density formula and the
#' attached `samples` are drawn from that same mixture. Useful for sanity-
#' checking the three fitting regimes against ground truth.
#'
#' @param with_samples If `TRUE`, attach `n` exact mixture draws.
#' @param n Number of samples to attach when `with_samples = TRUE`.
#' @param seed Optional integer seed used when drawing the samples.
#'
#' @returns A [gmm_target] in dimension 2.
#' @family targets
#' @export
#' @examples
#' m <- mixture_target(with_samples = TRUE, n = 100L)
#' m
mixture_target <- function(with_samples = FALSE, n = 2000L, seed = 1L) {
  w <- c(0.3, 0.4, 0.3)
  mus <- list(c(-2, -2), c(0, 0), c(2, 2))
  S1 <- matrix(c(0.6, 0.0, 0.0, 0.6), 2, 2)
  S2 <- matrix(c(0.5, 0.2, 0.2, 0.5), 2, 2)
  S3 <- matrix(c(0.4, -0.1, -0.1, 0.4), 2, 2)
  Ss <- list(S1, S2, S3)

  log_density <- function(x) {
    if (is.null(dim(x))) x <- matrix(x, nrow = 1L)
    parts <- gmm_log_unnorm(x, w, mus, Ss)
    logsumexp_rows(parts)
  }

  samples <- NULL
  if (isTRUE(with_samples)) {
    samples <- withr::with_seed(seed, {
      z <- sample.int(3L, size = n, replace = TRUE, prob = w)
      out <- matrix(0, nrow = n, ncol = 2L)
      for (k in 1:3) {
        idx <- which(z == k)
        if (length(idx) > 0) {
          out[idx, ] <- mvnfast::rmvn(length(idx), mu = mus[[k]], sigma = Ss[[k]])
        }
      }
      out
    })
  }

  gmm_target(
    n_dim = 2L,
    log_density = log_density,
    samples = samples,
    normalised = TRUE,
    log_normalizer = 0,
    name = "three_mixture",
    metadata = list(
      family = "mixture",
      weights = w,
      means = mus,
      covariances = Ss
    )
  )
}

