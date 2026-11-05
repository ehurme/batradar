#' Copy the targets pipeline into a project folder
#'
#' @param path Destination folder.
#' @param overwrite Overwrite an existing `_targets.R`.
#' @export
use_pipeline <- function(path = ".", overwrite = FALSE) {
  src <- system.file("pipeline", "_targets.R", package = "batradar")
  dest <- file.path(path, "_targets.R")
  if (file.exists(dest) && !overwrite) cli::cli_abort("{.file {dest}} exists; use overwrite = TRUE.")
  file.copy(src, dest, overwrite = overwrite)
  cli::cli_inform("Pipeline written to {.file {dest}}. Edit the CONFIG block, then run targets::tar_make().")
  invisible(dest)
}
