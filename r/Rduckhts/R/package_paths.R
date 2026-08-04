#' @keywords internal
duckhts_extension_dir <- function() {
  ext_path <- system.file(
    "duckhts_extension",
    package = "Rduckhts",
    mustWork = FALSE
  )
  if (nzchar(ext_path) && dir.exists(ext_path)) {
    return(ext_path)
  }
  return(NULL)
}
