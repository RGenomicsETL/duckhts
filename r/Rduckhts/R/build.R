# The package configure scripts are the sole extension-build implementation.
# Keep this exported compatibility entry point long enough for callers to get
# an explicit migration error instead of silently running a second compiler,
# source inventory, and metadata implementation.

#' Retired manual DuckHTS extension builder
#'
#' The former manual builder duplicated the package configure path and could
#' mutate an installed package tree. Build Rduckhts from a source tarball with
#' \code{R CMD build} and install that tarball instead.
#'
#' @param build_dir Retained for source compatibility; ignored.
#' @param make Retained for source compatibility; ignored.
#' @param force Retained for source compatibility; ignored.
#' @param verbose Retained for source compatibility; ignored.
#' @return This function always raises an error.
#' @export
duckhts_build <- function(build_dir = NULL, make = NULL, force = FALSE, verbose = TRUE) {
  .Deprecated(
    new = "R CMD build followed by R CMD INSTALL on the source tarball",
    package = "Rduckhts",
    msg = paste(
      "duckhts_build() is retired; use the package tarball workflow",
      "so configure owns the compiler, source inventory, and metadata receipt."
    )
  )
  stop(
    "duckhts_build() is retired. Run 'R CMD build .' and install the resulting ",
    "Rduckhts_*.tar.gz instead.",
    call. = FALSE
  )
}
