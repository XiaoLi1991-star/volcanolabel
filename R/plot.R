#' Create an elegant volcano plot
#'
#' Builds a publication-ready volcano plot from already prepared plotting
#' columns. The function does not read files, transform p-values, classify
#' genes, or choose labels. It focuses on the drawing layer: soft points,
#' threshold guides, counted legends, and adaptive outside label columns.
#' Every call prints resolved automatic parameters with [message()] so
#' agent-assisted tuning can start from the exact values used in the previous
#' plot.
#'
#' @param data A data frame containing plotting-ready columns.
#' @param x,y Column names for the x and y positions. For a typical volcano
#'   plot these are log2 fold change and `-log10(p)`.
#' @param label Optional column containing labels to draw. `NA` and empty
#'   strings are ignored.
#' @param color Optional column used for point groups and colors.
#' @param x_cutoff Optional positive x cutoff for symmetric vertical guide
#'   lines. This is only a guide; it does not classify rows.
#' @param y_cutoff Optional y cutoff for a horizontal guide line. This is only
#'   a guide; it does not classify rows.
#' @param ... Unused arguments are rejected with a clear error.
#' @param title,subtitle Plot title and subtitle.
#' @param palette Named color vector or palette name passed to
#'   [volcano_palette()]. Custom color groups are supported with named vectors.
#' @param point_size,point_alpha Main point size and alpha.
#' @param label_size,label_color Label text controls.
#' @param label_segment_color,label_segment_size,label_segment_alpha Connector
#'   controls.
#' @param label_point_ring_halo_color Color for the outer contrast halo around
#'   each labelled source point ring.
#' @param label_point_ring_halo_alpha,label_point_ring_halo_size,label_point_ring_halo_stroke
#'   Alpha, size, and stroke width for the outer contrast halo.
#' @param label_point_ring_color Color for the ring drawn around each labelled
#'   source point. Use `"auto"` to derive higher-contrast ring colors from
#'   each point group color, or `"group"` to follow group colors exactly.
#' @param label_point_ring_normal_color Ring color used for `Normal` labelled
#'   points when `label_point_ring_color = "auto"`. The default `"auto"`
#'   derives a higher-contrast color from the current `Normal` point color.
#'   Set a color value to override only the `Normal` ring color.
#' @param label_point_ring_alpha,label_point_ring_size,label_point_ring_stroke
#'   Alpha, size, and stroke width for the inner labelled source point rings.
#' @param label_wrap_width Maximum label line width. Use `"auto"` to choose
#'   from label length and density.
#' @param label_max_lines Maximum number of lines per label before capping with
#'   an ellipsis. Use `"auto"` to allow more detail only when labels are sparse.
#' @param label_anchor_x_left,label_anchor_x_right Optional x positions for the
#'   left and right outside label anchor columns. Leave as `NULL` for automatic
#'   placement.
#' @param plot_xmax Optional x-axis half-width.
#' @param ymax Optional y-axis maximum.
#' @param base_size Base font size.
#' @param legend_position Legend position.
#' @param xlab,ylab Axis labels.
#'
#' @return A ggplot object.
#' @export
volcano_plot <- function(data,
                         x,
                         y,
                         label = NULL,
                         color = NULL,
                         x_cutoff = NULL,
                         y_cutoff = NULL,
                         ...,
                         title = NULL,
                         subtitle = NULL,
                         palette = "aurora",
                         point_size = 1.55,
                         point_alpha = 0.78,
                         label_size = 3.4,
                         label_color = "#1F2937",
                         label_segment_color = "#9AA4B2",
                         label_segment_size = 0.32,
                         label_segment_alpha = 0.88,
                         label_point_ring_halo_color = "#E5E7EB",
                         label_point_ring_halo_alpha = 0.93,
                         label_point_ring_halo_size = 3.45,
                         label_point_ring_halo_stroke = 1.15,
                         label_point_ring_color = "auto",
                         label_point_ring_normal_color = "auto",
                         label_point_ring_alpha = 0.95,
                         label_point_ring_size = 2.75,
                         label_point_ring_stroke = 0.85,
                         label_wrap_width = "auto",
                         label_max_lines = "auto",
                         label_anchor_x_left = NULL,
                         label_anchor_x_right = NULL,
                         plot_xmax = NULL,
                         ymax = NULL,
                         base_size = 12,
                         legend_position = "top",
                         xlab = NULL,
                         ylab = NULL) {
  unused <- list(...)
  if (length(unused) > 0) {
    unused_names <- names(unused)
    unnamed <- !nzchar(unused_names %||% rep("", length(unused)))
    unused_names[unnamed] <- paste0("..", which(unnamed))
    stop("unused argument(s): ", paste(unused_names, collapse = ", "), call. = FALSE)
  }
  prepared <- build_volcano_plot_data(
    data = data,
    x = x,
    y = y,
    label = label,
    color = color
  )
  label_data <- prepared[!is.na(prepared$.volcano_label) & nzchar(prepared$.volcano_label), , drop = FALSE]

  x_cutoff_value <- normalize_plot_cutoff(x_cutoff, "x_cutoff")
  y_cutoff_value <- normalize_plot_cutoff(y_cutoff, "y_cutoff")
  label_point_ring_halo_color <- check_scalar_string(label_point_ring_halo_color, "label_point_ring_halo_color")
  label_point_ring_halo_alpha <- check_plot_number(label_point_ring_halo_alpha, "label_point_ring_halo_alpha", min = 0, max = 1)
  label_point_ring_halo_size <- check_plot_number(label_point_ring_halo_size, "label_point_ring_halo_size", min = .Machine$double.eps)
  label_point_ring_halo_stroke <- check_plot_number(label_point_ring_halo_stroke, "label_point_ring_halo_stroke", min = 0)
  label_point_ring_color <- normalize_label_ring_color(label_point_ring_color)
  label_point_ring_normal_color <- check_scalar_string(label_point_ring_normal_color, "label_point_ring_normal_color")
  label_point_ring_alpha <- check_plot_number(label_point_ring_alpha, "label_point_ring_alpha", min = 0, max = 1)
  label_point_ring_size <- check_plot_number(label_point_ring_size, "label_point_ring_size", min = .Machine$double.eps)
  label_point_ring_stroke <- check_plot_number(label_point_ring_stroke, "label_point_ring_stroke", min = 0)
  label_anchor_x_left <- normalize_optional_anchor_x(label_anchor_x_left, "label_anchor_x_left")
  label_anchor_x_right <- normalize_optional_anchor_x(label_anchor_x_right, "label_anchor_x_right")

  max_fc <- max(abs(prepared$.volcano_x), abs(x_cutoff_value %||% NA_real_), na.rm = TRUE)
  if (!is.finite(max_fc) || max_fc <= 0) {
    max_fc <- 1
  }
  max_y <- max(prepared$.volcano_y[is.finite(prepared$.volcano_y)], y_cutoff_value %||% NA_real_, na.rm = TRUE)
  if (!is.finite(max_y) || max_y <= 0) {
    max_y <- 1
  }
  plot_xmax_auto <- is.null(plot_xmax) || is.na(plot_xmax)
  if (plot_xmax_auto) {
    if (nrow(label_data) > 0) {
      label_plan <- plan_volcano_label_layout(
        add_label_sides(label_data),
        label_col = ".volcano_label",
        wrap_width = label_wrap_width,
        max_lines = label_max_lines
      )
      side_counts <- table(factor(add_label_sides(label_data)$.label_side, levels = c("left", "right")))
      max_side_count <- max(as.integer(side_counts), 0)
      count_pressure <- sqrt(max_side_count) * 0.11
      width_pressure <- min(label_plan$wrap_width, 34) * if (label_plan$max_lines > 1) 0.038 else 0.020
      label_padding <- if (label_plan$max_lines > 1 || label_plan$wrap_width > 18) 1.15 else 0.70
      relative_padding <- max_fc * (1.12 + min(max_side_count, 20) * 0.004)
      plot_xmax <- max(max_fc * 1.18, relative_padding, max_fc + label_padding + count_pressure + width_pressure)
    } else {
      plot_xmax <- max_fc * 1.10
    }
  }
  ymax_auto <- is.null(ymax) || is.na(ymax)
  if (ymax_auto) {
    ymax <- max_y * 1.06
    if (nrow(label_data) > 0) {
      ymax <- max(
        ymax,
        estimate_label_ymax(
          label_data,
          label_col = ".volcano_label",
          wrap_width = label_wrap_width,
          max_lines = label_max_lines
        )
      )
    }
  }
  label_layout <- NULL
  if (nrow(label_data) > 0) {
    for (i in seq_len(if (plot_xmax_auto) 3 else 1)) {
      label_layout <- compute_outside_label_layout(
        label_data,
        plot_xmax = plot_xmax,
        ymax = ymax,
        label_col = ".volcano_label",
        anchor_x_left = label_anchor_x_left,
        anchor_x_right = label_anchor_x_right,
        wrap_width = label_wrap_width,
        max_lines = label_max_lines,
        point_data = prepared
      )
      if (!plot_xmax_auto) {
        break
      }
      compact_xmax <- estimate_compact_plot_xmax(
        prepared,
        label_layout,
        plot_xmax = plot_xmax,
        ymax = ymax,
        label_col = ".volcano_label"
      )
      if (!is.finite(compact_xmax) || compact_xmax >= plot_xmax * 0.995) {
        break
      }
      plot_xmax <- compact_xmax
    }
  }

  group_levels <- levels(prepared$regulation)
  palette <- resolve_plot_palette(palette, group_levels)

  counts <- table(factor(prepared$regulation, levels = group_levels))
  legend_labels <- paste0(group_levels, " (", as.integer(counts[group_levels]), ")")
  names(legend_labels) <- group_levels
  legend_guide <- if (is.null(color)) {
    "none"
  } else {
    ggplot2::guide_legend(
      override.aes = list(
        shape = 21,
        fill = unname(palette[group_levels]),
        alpha = 1,
        size = point_size + 0.7
      )
    )
  }

  plot <- ggplot2::ggplot(prepared, ggplot2::aes(x = .volcano_x, y = .volcano_y)) +
    ggplot2::geom_point(
      ggplot2::aes(color = regulation, fill = regulation),
      shape = 21,
      size = point_size,
      alpha = point_alpha,
      stroke = 0.12
    )
  if (!is.null(x_cutoff_value) && is.finite(x_cutoff_value) && x_cutoff_value > 0) {
    plot <- plot +
      ggplot2::geom_vline(
        xintercept = c(-x_cutoff_value, x_cutoff_value),
        linetype = "22",
        color = "#A5ADBA",
        linewidth = 0.34
      )
  }
  if (!is.null(y_cutoff_value) && is.finite(y_cutoff_value)) {
    plot <- plot +
      ggplot2::geom_hline(
        yintercept = y_cutoff_value,
        linetype = "22",
        color = "#A5ADBA",
        linewidth = 0.34
      )
  }
  plot <- plot +
    ggplot2::scale_color_manual(
      name = color %||% "Group",
      values = palette,
      breaks = group_levels,
      labels = legend_labels,
      guide = legend_guide,
      drop = FALSE
    ) +
    ggplot2::scale_fill_manual(values = palette, breaks = group_levels, labels = legend_labels, guide = "none", drop = FALSE) +
    ggplot2::coord_cartesian(xlim = c(-plot_xmax, plot_xmax), ylim = c(0, ymax), clip = "off") +
    ggplot2::labs(
      x = xlab %||% x,
      y = ylab %||% y,
      title = normalize_optional_text(title),
      subtitle = normalize_optional_text(subtitle)
    ) +
    theme_volcano(base_size = base_size, legend_position = legend_position)

  if (nrow(label_data) > 0) {
    plot <- add_outside_label_layers(
      plot,
      label_layout,
      palette = palette,
      point_size = point_size,
      point_alpha = point_alpha,
      label_size = label_size,
      label_color = label_color,
      label_segment_color = label_segment_color,
      label_segment_size = label_segment_size,
      label_segment_alpha = label_segment_alpha,
      label_point_ring_halo_color = label_point_ring_halo_color,
      label_point_ring_halo_alpha = label_point_ring_halo_alpha,
      label_point_ring_halo_size = label_point_ring_halo_size,
      label_point_ring_halo_stroke = label_point_ring_halo_stroke,
      label_point_ring_color = label_point_ring_color,
      label_point_ring_normal_color = label_point_ring_normal_color,
      label_point_ring_alpha = label_point_ring_alpha,
      label_point_ring_size = label_point_ring_size,
      label_point_ring_stroke = label_point_ring_stroke
    )
  }

  message_resolved_volcano_params(
    x = x,
    y = y,
    label = label,
    color = color,
    palette = palette,
    x_cutoff = x_cutoff_value,
    y_cutoff = y_cutoff_value,
    plot_xmax = plot_xmax,
    ymax = ymax,
    plot_xmax_auto = plot_xmax_auto,
    ymax_auto = ymax_auto,
    label_layout = label_layout,
    point_size = point_size,
    point_alpha = point_alpha,
    label_size = label_size,
    label_segment_size = label_segment_size,
    label_point_ring_color = label_point_ring_color,
    label_point_ring_normal_color = label_point_ring_normal_color,
    label_point_ring_size = label_point_ring_size,
    label_point_ring_halo_size = label_point_ring_halo_size
  )

  plot
}

build_volcano_plot_data <- function(data,
                                    x,
                                    y,
                                    label = NULL,
                                    color = NULL) {
  data <- as_plain_data_frame(data)
  x <- check_scalar_string(x, "x")
  y <- check_scalar_string(y, "y")
  label <- normalize_optional_column(label, "label")
  color <- normalize_optional_column(color, "color")

  missing <- setdiff(c(x, y, label, color), names(data))
  if (length(missing) > 0) {
    stop("data is missing required column(s): ", paste(missing, collapse = ", "), call. = FALSE)
  }

  out <- data
  out$.volcano_x <- suppressWarnings(as.numeric(out[[x]]))
  out$.volcano_y <- suppressWarnings(as.numeric(out[[y]]))
  keep <- is.finite(out$.volcano_x) & is.finite(out$.volcano_y)
  if (any(!keep)) {
    warning(sum(!keep), " rows were removed because x or y was not finite", call. = FALSE)
    out <- out[keep, , drop = FALSE]
  }
  if (nrow(out) == 0) {
    stop("data has no rows with finite x and y values", call. = FALSE)
  }

  if (is.null(label)) {
    out$.volcano_label <- NA_character_
  } else {
    out$.volcano_label <- trimws(as.character(out[[label]]))
    out$.volcano_label[is.na(out[[label]]) | !nzchar(out$.volcano_label)] <- NA_character_
  }

  if (is.null(color)) {
    out$regulation <- factor("Point", levels = "Point")
  } else {
    groups <- trimws(as.character(out[[color]]))
    groups[is.na(out[[color]]) | !nzchar(groups)] <- "Point"
    levels <- volcano_group_levels(groups)
    out$regulation <- factor(groups, levels = levels)
  }

  out
}

normalize_optional_column <- function(value, name) {
  if (is.null(value)) {
    return(NULL)
  }
  check_scalar_string(value, name)
}

normalize_plot_cutoff <- function(value, name) {
  if (is.null(value) || length(value) == 0) {
    return(NULL)
  }
  if (!is.numeric(value) || length(value) != 1 || is.na(value)) {
    stop(name, " must be NULL or a single numeric value", call. = FALSE)
  }
  as.numeric(value)
}

normalize_optional_anchor_x <- function(value, name) {
  if (is.null(value) || length(value) == 0) {
    return(NULL)
  }
  if (!is.numeric(value) || length(value) != 1 || !is.finite(value)) {
    stop(name, " must be NULL or a single finite numeric value", call. = FALSE)
  }
  as.numeric(value)
}

normalize_label_ring_color <- function(value) {
  value <- check_scalar_string(value, "label_point_ring_color")
  lowered <- tolower(value)
  if (lowered %in% c("auto", "group")) {
    return(lowered)
  }
  value
}

check_plot_number <- function(value, name, min = -Inf, max = Inf) {
  if (!is.numeric(value) || length(value) != 1 || !is.finite(value) || value < min || value > max) {
    stop(name, " must be a single finite number", call. = FALSE)
  }
  as.numeric(value)
}

volcano_group_levels <- function(groups) {
  groups <- as.character(groups)
  canonical <- c("Up", "Normal", "Down")
  present <- canonical[canonical %in% groups]
  extras <- unique(groups[!groups %in% canonical])
  c(present, extras)
}

resolve_plot_palette <- function(palette, group_levels) {
  if (is.character(palette) && length(palette) == 1 && is.null(names(palette))) {
    base <- volcano_palette(palette)
    missing <- setdiff(group_levels, names(base))
    if (length(missing) == 0) {
      return(base[group_levels])
    }
    fallback <- fallback_volcano_palette(length(group_levels))
    names(fallback) <- group_levels
    shared <- intersect(group_levels, names(base))
    fallback[shared] <- base[shared]
    return(fallback)
  }

  if (!is.character(palette) || length(palette) == 0) {
    stop("palette must be a palette name or a character vector of colors", call. = FALSE)
  }
  if (is.null(names(palette)) || any(!nzchar(names(palette)))) {
    if (length(palette) < length(group_levels)) {
      stop("unnamed palette must provide at least one color per group", call. = FALSE)
    }
    palette <- palette[seq_along(group_levels)]
    names(palette) <- group_levels
    return(palette)
  }

  missing <- setdiff(group_levels, names(palette))
  if (length(missing) > 0) {
    stop("palette is missing colors for: ", paste(missing, collapse = ", "), call. = FALSE)
  }
  palette[group_levels]
}

fallback_volcano_palette <- function(n) {
  base <- c(
    "#E45756", "#4C78A8", "#72B7B2", "#F58518", "#54A24B",
    "#B279A2", "#FF9DA6", "#9D755D", "#BAB0AC", "#2E86AB"
  )
  if (n <= length(base)) {
    return(base[seq_len(n)])
  }
  grDevices::hcl.colors(n, palette = "Dark 3")
}

estimate_compact_plot_xmax <- function(prepared,
                                       label_layout,
                                       plot_xmax,
                                       ymax,
                                       label_col) {
  if (nrow(label_layout) == 0) {
    return(max(abs(prepared$.volcano_x), na.rm = TRUE) * 1.08)
  }
  wrap_width <- unique(label_layout$.label_wrap_width)
  max_lines <- unique(label_layout$.label_max_lines)
  layout_plan <- plan_volcano_label_layout(
    label_layout,
    label_col = label_col,
    wrap_width = wrap_width[1],
    max_lines = max_lines[1]
  )
  boxes <- label_text_boxes(
    label_layout,
    plot_xmax = plot_xmax,
    ymax = ymax,
    label_col = label_col,
    layout_plan = layout_plan
  )
  required <- max(
    abs(prepared$.volcano_x),
    abs(label_layout$.label_anchor_x),
    abs(c(boxes$xmin, boxes$xmax)),
    na.rm = TRUE
  )
  max(required * 1.045, max(abs(prepared$.volcano_x), na.rm = TRUE) * 1.035, 1)
}

estimate_label_ymax <- function(label_data,
                                label_col,
                                wrap_width,
                                max_lines) {
  label_data <- add_label_sides(label_data)
  label_plan <- plan_volcano_label_layout(
    label_data,
    label_col = label_col,
    wrap_width = wrap_width,
    max_lines = max_lines
  )
  labels <- wrap_volcano_labels(
    label_data[[label_col]],
    width = label_plan$wrap_width,
    max_lines = label_plan$max_lines
  )
  label_lines <- lengths(strsplit(labels, "\n", fixed = TRUE))
  side <- label_data$.label_side

  required <- vapply(
    split(label_lines, side),
    function(lines) {
      if (length(lines) <= 1) {
        return(0)
      }
      min_step <- 0.42 + 0.22 * (min(max(lines), 3) - 1)
      ((length(lines) - 1) * min_step) / 0.84
    },
    numeric(1)
  )
  max(required, 0)
}

add_label_sides <- function(label_data) {
  if (".label_side" %in% names(label_data)) {
    return(label_data)
  }
  direction <- if ("regulation" %in% names(label_data)) as.character(label_data$regulation) else NA_character_
  label_data$.label_side <- ifelse(
    direction == "Down" | (is.na(direction) & label_data$.volcano_x < 0),
    "left",
    ifelse(direction == "Up" | label_data$.volcano_x >= 0, "right", "left")
  )
  label_data
}

add_outside_label_layers <- function(plot,
                                     label_layout,
                                     palette,
                                     point_size,
                                     point_alpha,
                                     label_size,
                                     label_color,
                                     label_segment_color,
                                     label_segment_size,
                                     label_segment_alpha,
                                     label_point_ring_halo_color,
                                     label_point_ring_halo_alpha,
                                     label_point_ring_halo_size,
                                     label_point_ring_halo_stroke,
                                     label_point_ring_color,
                                     label_point_ring_normal_color,
                                     label_point_ring_alpha,
                                     label_point_ring_size,
                                     label_point_ring_stroke) {
  first_leg <- data.frame(
    x = label_layout$.label_anchor_x,
    y = label_layout$.label_y,
    xend = label_layout$.label_elbow_x,
    yend = label_layout$.label_y
  )
  second_leg <- data.frame(
    x = label_layout$.label_elbow_x,
    y = label_layout$.label_y,
    xend = label_layout$.volcano_x,
    yend = label_layout$.volcano_y
  )

  out <- plot +
    ggplot2::geom_segment(
      data = first_leg,
      ggplot2::aes(x = x, y = y, xend = xend, yend = yend),
      inherit.aes = FALSE,
      color = label_segment_color,
      linewidth = label_segment_size,
      alpha = label_segment_alpha,
      lineend = "round"
    ) +
    ggplot2::geom_segment(
      data = second_leg,
      ggplot2::aes(x = x, y = y, xend = xend, yend = yend),
      inherit.aes = FALSE,
      color = label_segment_color,
      linewidth = label_segment_size,
      alpha = label_segment_alpha,
      lineend = "round"
    ) +
    ggplot2::geom_point(
      data = label_layout,
      ggplot2::aes(x = .label_anchor_x, y = .label_y, fill = regulation),
      inherit.aes = FALSE,
      shape = 21,
      color = label_segment_color,
      size = 1.95,
      stroke = 0.24,
      alpha = 0.96,
      show.legend = FALSE
    ) +
    ggplot2::geom_point(
      data = label_layout,
      ggplot2::aes(x = .volcano_x, y = .volcano_y),
      inherit.aes = FALSE,
      shape = 21,
      color = label_point_ring_halo_color,
      fill = NA,
      size = label_point_ring_halo_size,
      stroke = label_point_ring_halo_stroke,
      alpha = label_point_ring_halo_alpha,
      show.legend = FALSE
    )
  if (label_point_ring_color %in% c("auto", "group")) {
    ring_layout <- label_layout
    ring_layout$.label_ring_color <- resolve_label_ring_colors(
      ring_layout,
      palette = palette,
      label_point_ring_color = label_point_ring_color,
      label_point_ring_normal_color = label_point_ring_normal_color
    )
    out <- out +
      ggplot2::geom_point(
        data = ring_layout,
        ggplot2::aes(x = .volcano_x, y = .volcano_y, color = I(.label_ring_color)),
        inherit.aes = FALSE,
        shape = 21,
        fill = NA,
        size = label_point_ring_size,
        stroke = label_point_ring_stroke,
        alpha = label_point_ring_alpha,
        show.legend = FALSE
      )
  } else {
    out <- out +
      ggplot2::geom_point(
        data = label_layout,
        ggplot2::aes(x = .volcano_x, y = .volcano_y),
        inherit.aes = FALSE,
        shape = 21,
        color = label_point_ring_color,
        fill = NA,
        size = label_point_ring_size,
        stroke = label_point_ring_stroke,
        alpha = label_point_ring_alpha,
        show.legend = FALSE
      )
  }
  out +
    ggplot2::geom_point(
      data = label_layout,
      ggplot2::aes(x = .volcano_x, y = .volcano_y, color = regulation, fill = regulation),
      inherit.aes = FALSE,
      shape = 21,
      size = point_size,
      stroke = 0.12,
      alpha = 1,
      show.legend = FALSE
    ) +
    ggplot2::geom_text(
      data = label_layout,
      ggplot2::aes(x = .label_text_x, y = .label_y, label = .volcano_label, hjust = .label_hjust),
      inherit.aes = FALSE,
      color = label_color,
      fontface = "bold",
      size = label_size,
      lineheight = 0.9,
      vjust = 0.5,
      show.legend = FALSE
    )
}

resolve_label_ring_colors <- function(label_layout,
                                      palette,
                                      label_point_ring_color,
                                      label_point_ring_normal_color) {
  groups <- as.character(label_layout$regulation)
  colors <- unname(palette[groups])
  colors[is.na(colors)] <- if (identical(tolower(label_point_ring_normal_color), "auto")) "#9CA3AF" else label_point_ring_normal_color
  if (identical(label_point_ring_color, "auto")) {
    colors <- vapply(colors, adaptive_ring_color, character(1), USE.NAMES = FALSE)
    if (!identical(tolower(label_point_ring_normal_color), "auto")) {
      normal <- groups %in% c("Normal", "Point")
      colors[normal] <- label_point_ring_normal_color
    }
  }
  colors
}

adaptive_ring_color <- function(color, min_contrast = 2.1) {
  color <- check_scalar_string(color, "color")
  luminance <- relative_color_luminance(color)
  target <- if (luminance >= 0.45) "#111827" else "#F9FAFB"
  amounts <- seq(0.25, 0.80, by = 0.05)
  candidates <- vapply(amounts, mix_hex_color, character(1), color = color, target = target, USE.NAMES = FALSE)
  contrasts <- vapply(candidates, color_contrast_ratio, numeric(1), color_b = color, USE.NAMES = FALSE)
  passing <- which(contrasts >= min_contrast)
  if (length(passing) > 0) {
    return(candidates[passing[1]])
  }
  candidates[which.max(contrasts)]
}

mix_hex_color <- function(amount, color, target) {
  rgb <- grDevices::col2rgb(color) / 255
  target_rgb <- grDevices::col2rgb(target) / 255
  mixed <- rgb * (1 - amount) + target_rgb * amount
  grDevices::rgb(mixed[1, 1], mixed[2, 1], mixed[3, 1])
}

relative_color_luminance <- function(color) {
  rgb <- grDevices::col2rgb(color) / 255
  linear <- ifelse(rgb <= 0.03928, rgb / 12.92, ((rgb + 0.055) / 1.055)^2.4)
  as.numeric(c(0.2126, 0.7152, 0.0722) %*% linear[, 1])
}

color_contrast_ratio <- function(color_a, color_b) {
  lum_a <- relative_color_luminance(color_a)
  lum_b <- relative_color_luminance(color_b)
  (max(lum_a, lum_b) + 0.05) / (min(lum_a, lum_b) + 0.05)
}

message_resolved_volcano_params <- function(x,
                                            y,
                                            label,
                                            color,
                                            palette,
                                            x_cutoff,
                                            y_cutoff,
                                            plot_xmax,
                                            ymax,
                                            plot_xmax_auto,
                                            ymax_auto,
                                            label_layout,
                                            point_size,
                                            point_alpha,
                                            label_size,
                                            label_segment_size,
                                            label_point_ring_color,
                                            label_point_ring_normal_color,
                                            label_point_ring_size,
                                            label_point_ring_halo_size) {
  label_counts <- resolved_label_counts(label_layout)
  lines <- c(
    "volcanolabel resolved parameters",
    "volcano_plot:",
    paste0("  x: ", x),
    paste0("  y: ", y),
    paste0("  label: ", format_resolved_value(label)),
    paste0("  color: ", format_resolved_value(color)),
    paste0("  palette: ", format_named_values(palette)),
    paste0("  x_cutoff: ", format_resolved_number(x_cutoff)),
    paste0("  y_cutoff: ", format_resolved_number(y_cutoff)),
    paste0("  plot_xmax: ", format_resolved_number(plot_xmax)),
    paste0("  plot_xmax_auto: ", format_resolved_value(plot_xmax_auto)),
    paste0("  ymax: ", format_resolved_number(ymax)),
    paste0("  ymax_auto: ", format_resolved_value(ymax_auto)),
    paste0("  point_size: ", format_resolved_number(point_size)),
    paste0("  point_alpha: ", format_resolved_number(point_alpha)),
    "label_layout:",
    paste0("  label_count: ", sum(label_counts)),
    paste0("  label_count_left: ", label_counts[["left"]]),
    paste0("  label_count_right: ", label_counts[["right"]]),
    paste0("  label_anchor_x_left: ", format_resolved_number(resolved_label_anchor(label_layout, "left"))),
    paste0("  label_anchor_x_right: ", format_resolved_number(resolved_label_anchor(label_layout, "right"))),
    paste0("  label_text_side_left: ", resolved_label_text_side(label_layout, "left")),
    paste0("  label_text_side_right: ", resolved_label_text_side(label_layout, "right")),
    paste0("  label_wrap_width: ", format_resolved_number(resolved_label_setting(label_layout, ".label_wrap_width"))),
    paste0("  label_max_lines: ", format_resolved_number(resolved_label_setting(label_layout, ".label_max_lines"))),
    paste0("  label_size: ", format_resolved_number(label_size)),
    paste0("  label_segment_size: ", format_resolved_number(label_segment_size)),
    paste0("  label_point_ring_color: ", label_point_ring_color),
    paste0("  label_point_ring_normal_color: ", label_point_ring_normal_color),
    paste0("  label_point_ring_size: ", format_resolved_number(label_point_ring_size)),
    paste0("  label_point_ring_halo_size: ", format_resolved_number(label_point_ring_halo_size)),
    paste0("  label_point_ring_colors: ", resolved_label_ring_color_summary(
      palette,
      label_point_ring_color,
      label_point_ring_normal_color
    ))
  )
  message(paste(lines, collapse = "\n"))
}

resolved_label_counts <- function(label_layout) {
  if (is.null(label_layout) || nrow(label_layout) == 0) {
    return(c(left = 0L, right = 0L))
  }
  counts <- table(factor(label_layout$.label_side, levels = c("left", "right")))
  c(left = as.integer(counts[["left"]]), right = as.integer(counts[["right"]]))
}

resolved_label_anchor <- function(label_layout, side) {
  if (is.null(label_layout) || nrow(label_layout) == 0) {
    return(NULL)
  }
  values <- unique(label_layout$.label_anchor_x[label_layout$.label_side == side])
  values <- values[is.finite(values)]
  if (length(values) == 0) {
    return(NULL)
  }
  values[1]
}

resolved_label_text_side <- function(label_layout, side) {
  if (is.null(label_layout) || nrow(label_layout) == 0) {
    return("none")
  }
  hjust <- unique(label_layout$.label_hjust[label_layout$.label_side == side])
  hjust <- hjust[is.finite(hjust)]
  if (length(hjust) == 0) {
    return("none")
  }
  if (all(hjust == 1)) {
    return("left")
  }
  if (all(hjust == 0)) {
    return("right")
  }
  "mixed"
}

resolved_label_setting <- function(label_layout, name) {
  if (is.null(label_layout) || nrow(label_layout) == 0 || !name %in% names(label_layout)) {
    return(NULL)
  }
  values <- unique(label_layout[[name]])
  values <- values[!is.na(values)]
  if (length(values) == 0) {
    return(NULL)
  }
  values[1]
}

resolved_label_ring_color_summary <- function(palette,
                                              label_point_ring_color,
                                              label_point_ring_normal_color) {
  if (!label_point_ring_color %in% c("auto", "group")) {
    return(label_point_ring_color)
  }
  ring_layout <- data.frame(regulation = names(palette), stringsAsFactors = FALSE)
  colors <- resolve_label_ring_colors(
    ring_layout,
    palette = palette,
    label_point_ring_color = label_point_ring_color,
    label_point_ring_normal_color = label_point_ring_normal_color
  )
  names(colors) <- names(palette)
  format_named_values(colors)
}

format_named_values <- function(values) {
  if (is.null(values) || length(values) == 0) {
    return("NULL")
  }
  names <- names(values)
  if (is.null(names) || any(!nzchar(names))) {
    return(paste(format_resolved_value(values), collapse = ", "))
  }
  paste(paste0(names, "=", unname(values)), collapse = ", ")
}

format_resolved_value <- function(value) {
  if (is.null(value) || length(value) == 0) {
    return("NULL")
  }
  if (is.logical(value)) {
    return(ifelse(value, "TRUE", "FALSE"))
  }
  if (is.numeric(value)) {
    return(format_resolved_number(value))
  }
  paste(as.character(value), collapse = ", ")
}

format_resolved_number <- function(value) {
  if (is.null(value) || length(value) == 0) {
    return("NULL")
  }
  if (!is.numeric(value)) {
    return(format_resolved_value(value))
  }
  if (any(!is.finite(value))) {
    value <- value[is.finite(value)]
  }
  if (length(value) == 0) {
    return("NULL")
  }
  paste(format(round(value, 4), trim = TRUE, scientific = FALSE), collapse = ", ")
}
