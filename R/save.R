#' Save a volcano plot
#'
#' Thin wrapper around [ggplot2::ggsave()] that defaults to a white background
#' and uses `ragg` for PNG output when it is installed.
#'
#' @param plot A ggplot object.
#' @param filename Output path.
#' @param width,height Plot size in inches.
#' @param dpi Raster resolution.
#' @param bg Background color.
#' @param ... Additional arguments passed to [ggplot2::ggsave()].
#'
#' @return The normalized output filename, invisibly.
#' @export
save_volcano <- function(plot,
                         filename,
                         width = 7,
                         height = 5.5,
                         dpi = 320,
                         bg = "white",
                         ...) {
  filename <- check_scalar_string(filename, "filename")
  dir.create(dirname(filename), recursive = TRUE, showWarnings = FALSE)
  ext <- tolower(tools::file_ext(filename))
  if (!nzchar(ext)) {
    stop("filename must include a file extension", call. = FALSE)
  }

  device <- NULL
  if (identical(ext, "png") && requireNamespace("ragg", quietly = TRUE)) {
    device <- ragg::agg_png
  }

  ggplot2::ggsave(
    filename = filename,
    plot = plot,
    width = width,
    height = height,
    dpi = dpi,
    bg = bg,
    device = device,
    ...
  )
  invisible(normalizePath(filename, mustWork = FALSE))
}
