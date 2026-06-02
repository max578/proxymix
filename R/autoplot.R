## autoplot.R -- ggplot2 visualisation of a fitted Gaussian-mixture proxy.
##
## Provides a single S7 method registering `gmm_fit` against the
## `ggplot2::autoplot()` generic. ggplot2 lives in Suggests; the method is
## an S7 external-generic registration, wired through `.onLoad()`'s
## `S7::methods_register()`, so it activates only when ggplot2 is present and
## the package keeps its standalone-functional `R CMD check`. Any requested
## subset of coordinates is reduced to one or two dimensions through the
## package's own closed-form `gmm_marginalise()` before rendering.

#' Plot a fitted Gaussian-mixture proxy
#'
#' @description
#' An `autoplot()` method for [gmm_fit] objects, rendering the fitted mixture
#' with `ggplot2`. The displayed coordinates are reduced to the requested one
#' or two dimensions through the closed-form marginal [gmm_marginalise()], so
#' the method works for a proxy of any ambient dimension `p`.
#'
#' A one-dimensional request draws the marginal mixture density, optionally
#' with the per-component densities underneath and a rug of the target's
#' samples. A two-dimensional request draws the marginal density as a viridis
#' raster with white contour lines, optionally overlaying each component's
#' mean and a probability-contour ellipse.
#'
#' @details
#' The method is registered against the `ggplot2::autoplot()` generic only
#' when `ggplot2` is installed; call it as `ggplot2::autoplot(fit)` or load
#' `ggplot2` first. It returns the `ggplot` object, so the usual `+` layering
#' applies for further customisation.
#'
#' @param object A [gmm_fit], typically from [fit_proxymix()].
#' @param dims Integer vector of length one or two giving the coordinate(s) to
#'   display, in `1:p`. Defaults to the first two coordinates (or the only
#'   coordinate when `p == 1`).
#' @param n_grid Integer scalar — the number of grid points per axis at which
#'   the density is evaluated. A two-dimensional plot evaluates `n_grid^2`
#'   points.
#' @param n_sd Numeric scalar — how many component standard deviations beyond
#'   the extreme component means the plotting window extends.
#' @param level Numeric scalar in `(0, 1)` — the probability level of the
#'   per-component ellipse drawn on a two-dimensional plot.
#' @param show_components Logical scalar — whether to overlay the per-component
#'   densities (one dimension) or mean-and-ellipse glyphs (two dimensions).
#' @param show_data Logical scalar — whether to overlay the target's samples,
#'   when the fitted target carries any.
#' @param ... Currently ignored, present for generic compatibility.
#'
#' @returns A `ggplot` object.
#' @family classes
#' @name autoplot.gmm_fit
#' @examplesIf requireNamespace("ggplot2", quietly = TRUE)
#' samples <- matrix(stats::rnorm(200), ncol = 2)
#' tgt <- gmm_target_from_samples(samples)
#' fit <- fit_proxymix(tgt, N = 2L, regime = "sample", max_iter = 25L)
#' ggplot2::autoplot(fit)
#' ggplot2::autoplot(fit, dims = 1L)
NULL

# The method body. ggplot2's `autoplot()` is an S3 generic in a Suggested
# package, and S7 objects carry a package-prefixed S3 class, so this worker is
# registered under "proxymix::gmm_fit" against ggplot2's generic in `.onLoad()`
# (see zzz.R) -- and only when ggplot2 is installed, so ggplot2 stays in
# Suggests and the package keeps a clean standalone `R CMD check`.
.autoplot_gmm_fit <- function(
  object,
  dims = c(1L, 2L),
  n_grid = 120L,
  n_sd = 3.5,
  level = 0.95,
  show_components = TRUE,
  show_data = TRUE,
  ...
) {
  rlang::check_installed("ggplot2", reason = "to `autoplot()` a `gmm_fit`.")

  # Resolve and validate the coordinates to display ---------------------------

  p_dim <- gmm_dim(object)
  dims <- if (p_dim == 1L) 1L else as.integer(dims)
  if (length(dims) < 1L || length(dims) > 2L) {
    cli::cli_abort("`dims` must name one or two coordinates.")
  }
  if (any(dims < 1L) || any(dims > p_dim) || anyDuplicated(dims) > 0L) {
    cli::cli_abort(
      "`dims` must be unique indices in {.val 1}:{.val {p_dim}}."
    )
  }

  g_disp <- gmm_marginalise(object, keep = dims)
  data_marginal <- .target_marginal_samples(object, dims, show_data)

  if (length(dims) == 1L) {
    .autoplot_gmm_1d(object, g_disp, dims, data_marginal,
                     n_grid, n_sd, show_components)
  } else {
    .autoplot_gmm_2d(object, g_disp, dims, data_marginal,
                     n_grid, n_sd, level, show_components)
  }
}

# One-dimensional density plot ------------------------------------------------

#' Render the one-dimensional marginal density of a fit
#'
#' @param object The original [gmm_fit] (used for the title only).
#' @param g_disp The one-dimensional marginal [gmm] to render.
#' @param dims The single retained coordinate index.
#' @param data_marginal A one-column matrix of target samples, or `NULL`.
#' @param n_grid,show_components,n_sd See [autoplot.gmm_fit].
#'
#' @returns A `ggplot` object.
#' @noRd
#' @keywords internal
.autoplot_gmm_1d <- function(object, g_disp, dims, data_marginal,
                             n_grid, n_sd, show_components) {
  mu <- vapply(g_disp@means, function(m) m[[1L]], numeric(1L))
  sdev <- vapply(g_disp@covariances, function(s) sqrt(s[[1L]]), numeric(1L))
  grid <- seq(min(mu - n_sd * sdev), max(mu + n_sd * sdev),
              length.out = as.integer(n_grid))
  mix_df <- data.frame(x = grid, density = dgmm(matrix(grid, ncol = 1L), g_disp))

  plt <- ggplot2::ggplot(
    mix_df,
    ggplot2::aes(x = .data$x, y = .data$density)
  )

  if (isTRUE(show_components)) {
    comp_df <- do.call(rbind, lapply(seq_along(mu), function(k) {
      data.frame(
        x = grid,
        density = g_disp@weights[[k]] * stats::dnorm(grid, mu[k], sdev[k]),
        component = factor(k)
      )
    }))
    plt <- plt + ggplot2::geom_line(
      data = comp_df,
      mapping = ggplot2::aes(group = .data$component),
      colour = "grey55", linetype = "dashed", linewidth = 0.3
    )
  }

  if (!is.null(data_marginal)) {
    plt <- plt + ggplot2::geom_rug(
      data = data.frame(x = data_marginal[, 1L]),
      mapping = ggplot2::aes(x = .data$x),
      inherit.aes = FALSE, sides = "b", alpha = 0.3
    )
  }

  plt +
    ggplot2::geom_line(colour = "#21908C", linewidth = 0.8) +
    ggplot2::labs(
      x = .dim_label(dims[1L]),
      y = "density",
      title = .fit_title(object)
    )
}

# Two-dimensional density plot ------------------------------------------------

#' Render the two-dimensional marginal density of a fit
#'
#' @param object The original [gmm_fit] (used for the title only).
#' @param g_disp The two-dimensional marginal [gmm] to render.
#' @param dims The two retained coordinate indices.
#' @param data_marginal A two-column matrix of target samples, or `NULL`.
#' @param n_grid,n_sd,level,show_components See [autoplot.gmm_fit].
#'
#' @returns A `ggplot` object.
#' @noRd
#' @keywords internal
.autoplot_gmm_2d <- function(object, g_disp, dims, data_marginal,
                             n_grid, n_sd, level, show_components) {
  centres <- do.call(rbind, g_disp@means)
  sdev <- vapply(g_disp@covariances, function(s) sqrt(diag(s)), numeric(2L))
  axis_grid <- function(j) {
    seq(min(centres[, j] - n_sd * sdev[j, ]),
        max(centres[, j] + n_sd * sdev[j, ]),
        length.out = as.integer(n_grid))
  }
  grid <- expand.grid(x = axis_grid(1L), y = axis_grid(2L))
  grid$density <- dgmm(as.matrix(grid[, c("x", "y")]), g_disp)

  plt <- ggplot2::ggplot(grid, ggplot2::aes(x = .data$x, y = .data$y)) +
    ggplot2::geom_raster(ggplot2::aes(fill = .data$density)) +
    ggplot2::scale_fill_viridis_c(option = "viridis") +
    ggplot2::geom_contour(
      ggplot2::aes(z = .data$density),
      colour = "white", alpha = 0.4, linewidth = 0.2
    )

  if (isTRUE(show_components)) {
    ellipses <- do.call(rbind, lapply(seq_len(nrow(centres)), function(k) {
      ell <- .gmm_ellipse(g_disp@means[[k]], g_disp@covariances[[k]], level)
      ell$component <- factor(k)
      ell
    }))
    plt <- plt +
      ggplot2::geom_path(
        data = ellipses,
        mapping = ggplot2::aes(group = .data$component),
        colour = "white", linewidth = 0.4
      ) +
      ggplot2::geom_point(
        data = data.frame(x = centres[, 1L], y = centres[, 2L]),
        colour = "white", shape = 3L, size = 2
      )
  }

  if (!is.null(data_marginal)) {
    plt <- plt + ggplot2::geom_point(
      data = data.frame(x = data_marginal[, 1L], y = data_marginal[, 2L]),
      colour = "white", alpha = 0.25, size = 0.5
    )
  }

  plt + ggplot2::labs(
    x = .dim_label(dims[1L]),
    y = .dim_label(dims[2L]),
    fill = "density",
    title = .fit_title(object)
  )
}

# Internal helpers ------------------------------------------------------------

#' Probability-contour ellipse of a bivariate Gaussian
#'
#' Returns the points tracing the `level`-probability contour of a bivariate
#' Gaussian, by scaling the unit circle along the covariance's principal axes
#' by the chi-squared radius.
#'
#' @param mu Length-two numeric vector — the component mean.
#' @param sigma A two-by-two covariance matrix.
#' @param level Numeric scalar in `(0, 1)` — the probability level.
#' @param n_seg Integer scalar — the number of segments tracing the ellipse.
#'
#' @returns A `data.frame` with columns `x` and `y`.
#' @noRd
#' @keywords internal
.gmm_ellipse <- function(mu, sigma, level = 0.95, n_seg = 100L) {
  radius <- sqrt(stats::qchisq(level, df = 2L))
  theta <- seq(0, 2 * pi, length.out = as.integer(n_seg))
  eig <- eigen(sigma, symmetric = TRUE)
  axes <- eig$vectors %*% diag(sqrt(pmax(eig$values, 0)), 2L)
  pts <- axes %*% (radius * rbind(cos(theta), sin(theta)))
  data.frame(x = mu[[1L]] + pts[1L, ], y = mu[[2L]] + pts[2L, ])
}

#' Marginal target samples for the displayed coordinates
#'
#' Extracts the requested columns of the fitted target's sample matrix, or
#' returns `NULL` when overlay is switched off or no samples are available.
#'
#' @param object A [gmm_fit].
#' @param dims The coordinate indices being displayed.
#' @param show_data Logical scalar — whether the overlay is requested.
#'
#' @returns A numeric matrix with `length(dims)` columns, or `NULL`.
#' @noRd
#' @keywords internal
.target_marginal_samples <- function(object, dims, show_data) {
  if (!isTRUE(show_data)) {
    return(NULL)
  }
  tgt <- object@target
  if (is.null(tgt) || !S7::S7_inherits(tgt, gmm_target)) {
    return(NULL)
  }
  samples <- tgt@samples
  if (is.null(samples) || !is.matrix(samples)) {
    return(NULL)
  }
  samples[, dims, drop = FALSE]
}

#' Axis label for a displayed coordinate
#'
#' @param d Integer scalar — the coordinate index.
#'
#' @returns A length-one character string.
#' @noRd
#' @keywords internal
.dim_label <- function(d) {
  sprintf("dimension %d", d)
}

#' Plot title summarising a fit
#'
#' @param object A [gmm_fit].
#'
#' @returns A length-one character string.
#' @noRd
#' @keywords internal
.fit_title <- function(object) {
  sprintf(
    "proxymix proxy: regime \"%s\", K = %d",
    object@regime,
    gmm_n_components(object)
  )
}
