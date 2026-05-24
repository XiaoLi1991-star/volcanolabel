#' Volcano plot palettes
#'
#' @param name Palette name. Currently `"aurora"`, `"classic"`, and
#'   `"midnight"` are available.
#' @param up,down,normal Optional color overrides.
#'
#' @return A named character vector with `Up`, `Normal`, and `Down` colors.
#' @export
volcano_palette <- function(name = "aurora", up = NULL, down = NULL, normal = NULL) {
  palettes <- list(
    aurora = c(Up = "#E45756", Normal = "#B9C0CA", Down = "#2E86AB"),
    classic = c(Up = "#D95F59", Normal = "#B7B7B7", Down = "#4DBBD5"),
    midnight = c(Up = "#FF7A70", Normal = "#8B95A5", Down = "#58B8D8")
  )
  name <- tolower(check_scalar_string(name, "name"))
  if (!name %in% names(palettes)) {
    stop("unknown palette: ", name, call. = FALSE)
  }
  pal <- palettes[[name]]
  pal["Up"] <- up %||% pal["Up"]
  pal["Down"] <- down %||% pal["Down"]
  pal["Normal"] <- normal %||% pal["Normal"]
  pal
}

#' A refined volcano plot theme
#'
#' @param base_size Base font size.
#' @param base_family Base font family.
#' @param legend_position Legend position passed to [ggplot2::theme()].
#'
#' @return A ggplot2 theme.
#' @export
theme_volcano <- function(base_size = 12,
                          base_family = "",
                          legend_position = "top") {
  ggplot2::theme_minimal(base_size = base_size, base_family = base_family) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major.x = ggplot2::element_line(color = "#EDF0F3", linewidth = 0.35),
      panel.grid.major.y = ggplot2::element_line(color = "#EDF0F3", linewidth = 0.35),
      axis.line = ggplot2::element_line(color = "#2F3A48", linewidth = 0.35),
      axis.ticks = ggplot2::element_line(color = "#2F3A48", linewidth = 0.35),
      axis.title = ggplot2::element_text(color = "#202A36", face = "bold"),
      axis.text = ggplot2::element_text(color = "#475467"),
      plot.title = ggplot2::element_text(color = "#17202B", face = "bold", size = base_size * 1.18),
      plot.subtitle = ggplot2::element_text(color = "#667085", margin = ggplot2::margin(t = 3, b = 7)),
      plot.caption = ggplot2::element_text(color = "#667085"),
      plot.margin = ggplot2::margin(10, 28, 10, 28),
      legend.position = legend_position,
      legend.title = ggplot2::element_text(color = "#202A36", face = "bold"),
      legend.text = ggplot2::element_text(color = "#475467"),
      legend.key.height = grid::unit(0.45, "lines"),
      legend.key.width = grid::unit(1.1, "lines")
    )
}
