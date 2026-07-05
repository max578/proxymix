.onLoad <- function(libname, pkgname) {
  S7::methods_register()

  # Register the autoplot() method only when ggplot2 is installed, so the
  # ggplot2 dependency stays in Suggests and `R CMD check` is clean without
  # it. S7 objects carry a package-prefixed S3 class, so the method registers
  # under "proxymix::gmm_fit" against ggplot2's S3 generic.
  if (requireNamespace("ggplot2", quietly = TRUE)) {
    registerS3method(
      "autoplot", "proxymix::gmm_fit", .autoplot_gmm_fit,
      envir = asNamespace("ggplot2")
    )
    ## Plain mixtures (operator-calculus results) plot too: the composed
    ## objects are exactly the ones users need to look at.
    registerS3method(
      "autoplot", "proxymix::gmm", .autoplot_gmm_fit,
      envir = asNamespace("ggplot2")
    )
  }
}
