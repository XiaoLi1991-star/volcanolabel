#' Wrap volcano labels
#'
#' Wraps labels at separators and shortens very long uninterrupted tokens.
#'
#' @param labels Character vector of labels.
#' @param width Maximum line width. Use `0` to disable wrapping.
#' @param max_lines Maximum number of output lines. Longer labels are capped
#'   with `ellipsis`.
#' @param ellipsis Text appended when a label is capped.
#'
#' @return A character vector with embedded newlines where needed.
#' @export
wrap_volcano_labels <- function(labels,
                                width = 36,
                                max_lines = Inf,
                                ellipsis = "...") {
  if (is.null(labels)) {
    return(character())
  }
  labels <- as.character(labels)
  if (is.null(max_lines) || length(max_lines) != 1 || is.na(max_lines) || max_lines <= 0) {
    max_lines <- Inf
  }
  if (is.null(width) || is.na(width) || width <= 0) {
    return(labels)
  }

  split_long_part <- function(text, width) {
    if (nchar(text, type = "width") <= width) {
      return(text)
    }
    suffix_width <- nchar(ellipsis, type = "width")
    if (width <= suffix_width) {
      return(substr(text, 1, width))
    }
    paste0(substr(text, 1, width - suffix_width), ellipsis)
  }

  wrap_one <- function(x) {
    if (is.na(x)) {
      return(NA_character_)
    }
    source <- gsub("([,;/]|-)", "\\1 ", x, perl = TRUE)
    wrapped <- strwrap(source, width = width)
    wrapped <- gsub("([,;/]|-)\\s+", "\\1", wrapped, perl = TRUE)
    wrapped <- unlist(lapply(wrapped, split_long_part, width = width), use.names = FALSE)
    if (length(wrapped) > max_lines) {
      wrapped <- wrapped[seq_len(max_lines)]
      suffix_width <- nchar(ellipsis, type = "width")
      keep_width <- max(width - suffix_width, 1)
      wrapped[length(wrapped)] <- paste0(substr(wrapped[length(wrapped)], 1, keep_width), ellipsis)
    }
    paste(wrapped, collapse = "\n")
  }

  vapply(labels, wrap_one, character(1))
}

#' Compute outside label positions
#'
#' Creates left and right label columns, ordered by significance, with connector
#' elbows back to the original volcano points.
#'
#' @param label_data Prepared data containing labels.
#' @param plot_xmax Positive x-axis half-width.
#' @param ymax Positive y-axis maximum.
#' @param label_col Column containing label text.
#' @param spacing,spacing_left,spacing_right Vertical spacing multipliers.
#' @param anchor_x_left,anchor_x_right Optional anchor x positions.
#' @param text_side_left,text_side_right Whether text sits to the `"left"` or
#'   `"right"` of each side's anchor point. Use `"auto"` to choose one
#'   direction per side from the available width and point layout.
#' @param wrap_width Label wrapping width.
#' @param max_lines Maximum number of lines per label.
#' @param point_data Optional prepared full plotting data. When supplied,
#'   automatic side decisions avoid placing labels across dense point regions.
#'
#' @return `label_data` with `.label_*` layout columns.
#' @export
compute_outside_label_layout <- function(label_data,
                                         plot_xmax,
                                         ymax,
                                         label_col = ".volcano_label",
                                         spacing = 0.72,
                                         spacing_left = NULL,
                                         spacing_right = NULL,
                                         anchor_x_left = NULL,
                                         anchor_x_right = NULL,
                                         text_side_left = "auto",
                                         text_side_right = "auto",
                                         wrap_width = "auto",
                                         max_lines = "auto",
                                         point_data = NULL) {
  label_data <- as_plain_data_frame(label_data)
  label_col <- check_scalar_string(label_col, "label_col")
  text_side_left <- validate_text_side(text_side_left, "text_side_left")
  text_side_right <- validate_text_side(text_side_right, "text_side_right")

  required <- c(".volcano_x", ".volcano_y", label_col)
  missing <- setdiff(required, names(label_data))
  if (length(missing) > 0) {
    stop("label_data is missing required columns: ", paste(missing, collapse = ", "), call. = FALSE)
  }
  if (!is.numeric(plot_xmax) || length(plot_xmax) != 1 || !is.finite(plot_xmax) || plot_xmax <= 0) {
    stop("plot_xmax must be a positive finite number", call. = FALSE)
  }
  if (!is.numeric(ymax) || length(ymax) != 1 || !is.finite(ymax) || ymax <= 0) {
    stop("ymax must be a positive finite number", call. = FALSE)
  }

  if (nrow(label_data) == 0) {
    label_data$.label_side <- character()
    label_data$.label_anchor_x <- numeric()
    label_data$.label_text_x <- numeric()
    label_data$.label_y <- numeric()
    label_data$.label_hjust <- numeric()
    label_data$.label_elbow_x <- numeric()
    return(label_data)
  }
  if (is.null(point_data)) {
    point_data <- label_data
  } else {
    point_data <- as_plain_data_frame(point_data)
  }
  if (!all(c(".volcano_x", ".volcano_y") %in% names(point_data))) {
    stop("point_data must contain .volcano_x and .volcano_y columns", call. = FALSE)
  }
  point_data <- point_data[is.finite(point_data$.volcano_x) & is.finite(point_data$.volcano_y), , drop = FALSE]

  direction <- if ("regulation" %in% names(label_data)) as.character(label_data$regulation) else NA_character_
  label_data$.label_side <- ifelse(
    direction == "Down" | (is.na(direction) & label_data$.volcano_x < 0),
    "left",
    ifelse(direction == "Up" | label_data$.volcano_x >= 0, "right", "left")
  )
  layout_plan <- plan_volcano_label_layout(
    label_data,
    label_col = label_col,
    wrap_width = wrap_width,
    max_lines = max_lines
  )
  wrap_width <- layout_plan$wrap_width
  max_lines <- layout_plan$max_lines
  label_data[[label_col]] <- disambiguate_wrapped_labels(wrap_volcano_labels(
    label_data[[label_col]],
    width = wrap_width,
    max_lines = max_lines
  ), width = wrap_width)

  anchor_fraction <- layout_plan$anchor_fraction
  side_gap <- plot_xmax * layout_plan$point_gap_fraction
  text_offset <- plot_xmax * layout_plan$text_gap_fraction
  axis_gap <- plot_xmax * layout_plan$axis_gap_fraction
  outer_limit <- plot_xmax * layout_plan$anchor_outer_fraction
  left_labeled <- label_data[label_data$.label_side == "left", , drop = FALSE]
  right_labeled <- label_data[label_data$.label_side == "right", , drop = FALSE]
  spacing_left <- spacing_left %||% spacing
  spacing_right <- spacing_right %||% spacing
  spacing_left <- max(as.numeric(spacing_left) * layout_plan$spacing_multiplier, 0.1)
  spacing_right <- max(as.numeric(spacing_right) * layout_plan$spacing_multiplier, 0.1)

  anchor_for_side <- function(part, side, text_side) {
    if (side == "left") {
      if (!is.null(anchor_x_left)) {
        return(anchor_x_left)
      }
      left_edge <- if (nrow(part) > 0) min(part$.volcano_x, na.rm = TRUE) - side_gap else -plot_xmax * anchor_fraction
      anchor_x <- max(-outer_limit, min(-plot_xmax * anchor_fraction, left_edge))
      if (text_side == "left") {
        anchor_x <- min(anchor_x, -axis_gap + text_offset)
      }
      anchor_x
    } else {
      if (!is.null(anchor_x_right)) {
        return(anchor_x_right)
      }
      right_edge <- if (nrow(part) > 0) max(part$.volcano_x, na.rm = TRUE) + side_gap else plot_xmax * anchor_fraction
      anchor_x <- min(outer_limit, max(plot_xmax * anchor_fraction, right_edge))
      if (text_side == "right") {
        anchor_x <- max(anchor_x, axis_gap - text_offset)
      }
      anchor_x
    }
  }

  place_side <- function(part, side, text_side) {
    if (nrow(part) == 0) {
      return(part)
    }
    part <- part[order(-part$.volcano_y, abs(part$.volcano_x), na.last = TRUE), , drop = FALSE]
    spacing_side <- if (side == "left") spacing_left else spacing_right
    anchor_x <- anchor_for_side(part, side, text_side)

    top <- ymax * 0.92
    floor <- ymax * 0.08
    step <- max(ymax * 0.065, 0.60) * spacing_side
    bottom <- top - (nrow(part) - 1) * step
    if (bottom < floor) {
      step <- (top - floor) / max(nrow(part) - 1, 1)
    }

    part$.label_anchor_x <- anchor_x
    part$.label_y <- top - (seq_len(nrow(part)) - 1) * step
    if (text_side == "left") {
      part$.label_text_x <- anchor_x - text_offset
      part$.label_hjust <- 1
    } else {
      part$.label_text_x <- anchor_x + text_offset
      part$.label_hjust <- 0
    }
    part$.label_elbow_x <- anchor_x + (part$.volcano_x - anchor_x) * 0.68
    part$.label_wrap_width <- wrap_width
    part$.label_max_lines <- max_lines
    part$.label_line_count <- lengths(strsplit(part[[label_col]], "\n", fixed = TRUE))
    part$.label_anchor_fraction <- anchor_fraction
    part$.label_text_side <- text_side
    part
  }

  choose_side_layout <- function(part, side, requested) {
    if (nrow(part) == 0) {
      return(part)
    }
    candidate_sides <- if (requested == "auto") {
      if (side == "left") c("left", "right") else c("right", "left")
    } else {
      requested
    }
    candidates <- lapply(candidate_sides, function(text_side) place_side(part, side, text_side))
    if (length(candidates) == 1) {
      return(candidates[[1]])
    }
    scores <- vapply(
      candidates,
      score_label_side_layout,
      numeric(1),
      side = side,
      plot_xmax = plot_xmax,
      ymax = ymax,
      label_col = label_col,
      point_data = point_data,
      layout_plan = layout_plan
    )
    candidates[[which.min(scores)]]
  }

  out <- rbind(
    choose_side_layout(left_labeled, "left", text_side_left),
    choose_side_layout(right_labeled, "right", text_side_right)
  )
  rownames(out) <- NULL
  out
}

plan_volcano_label_layout <- function(label_data,
                                      label_col,
                                      wrap_width = "auto",
                                      max_lines = "auto") {
  labels <- as.character(label_data[[label_col]])
  label_widths <- nchar(labels, type = "width")
  label_widths <- label_widths[is.finite(label_widths)]
  if (length(label_widths) == 0) {
    label_widths <- 1
  }

  side_counts <- table(factor(label_data$.label_side, levels = c("left", "right")))
  max_side_count <- max(as.integer(side_counts), 0)
  q70 <- as.numeric(stats::quantile(label_widths, 0.70, names = FALSE, type = 7))
  q90 <- as.numeric(stats::quantile(label_widths, 0.90, names = FALSE, type = 7))
  longest <- max(label_widths, na.rm = TRUE)

  wrap_width_auto <- is_auto_value(wrap_width)
  max_lines_auto <- is_auto_value(max_lines)

  if (wrap_width_auto) {
    density_cap <- if (max_side_count >= 12) {
      14
    } else if (max_side_count >= 8) {
      16
    } else if (max_side_count >= 5) {
      18
    } else {
      20
    }
    detail_width <- if (max_side_count <= 4 && q90 > 30) {
      min(20, max(18, ceiling(q90 / 3)))
    } else {
      max(12, ceiling(min(q70, density_cap)))
    }
    wrap_width <- as.integer(max(10, min(density_cap, detail_width)))
  } else {
    wrap_width <- parse_positive_integer(wrap_width, "wrap_width")
  }

  if (max_lines_auto) {
    max_lines <- if (max_side_count <= 4 && longest > wrap_width * 1.15) 2L else 1L
  } else {
    max_lines <- parse_positive_integer(max_lines, "max_lines")
  }

  long_label_pressure <- min(max(q90 / max(wrap_width, 1), 1), 2)
  density_bonus <- if (max_side_count >= 16) {
    0.07
  } else if (max_side_count >= 10) {
    0.04
  } else if (max_side_count >= 6) {
    0.02
  } else {
    0
  }
  width_bonus <- if (long_label_pressure > 1.55) {
    0.07
  } else if (long_label_pressure > 1.25) {
    0.04
  } else {
    0
  }
  anchor_fraction <- min(0.60, max(0.46, 0.46 + density_bonus + width_bonus))
  point_gap_fraction <- if (max_side_count >= 10) {
    0.045
  } else if (long_label_pressure > 1.35) {
    0.050
  } else {
    0.035
  }
  list(
    wrap_width = wrap_width,
    max_lines = max_lines,
    anchor_fraction = anchor_fraction,
    point_gap_fraction = point_gap_fraction,
    text_gap_fraction = 0.030,
    axis_gap_fraction = 0.18,
    anchor_outer_fraction = 0.96,
    char_width_fraction = 0.0125,
    spacing_multiplier = if (max_lines > 1) 1.16 else 1
  )
}

score_label_side_layout <- function(layout,
                                    side,
                                    plot_xmax,
                                    ymax,
                                    label_col,
                                    point_data,
                                    layout_plan) {
  boxes <- label_text_boxes(
    layout,
    plot_xmax = plot_xmax,
    ymax = ymax,
    label_col = label_col,
    layout_plan = layout_plan
  )
  outer_limit <- plot_xmax * 0.985
  overflow <- sum(pmax(-outer_limit - boxes$xmin, 0) + pmax(boxes$xmax - outer_limit, 0))

  axis_gap <- plot_xmax * layout_plan$axis_gap_fraction
  axis_overlap <- if (side == "left") {
    sum(pmax(boxes$xmax - (-axis_gap), 0))
  } else {
    sum(pmax(axis_gap - boxes$xmin, 0))
  }

  point_overlap <- count_label_point_overlaps(boxes, point_data)
  side_density <- count_label_point_density(boxes, point_data)
  inward <- if (side == "left") {
    unique(layout$.label_hjust) == 0
  } else {
    unique(layout$.label_hjust) == 1
  }
  inward_penalty <- if (isTRUE(inward)) 1.15 else 0
  edge_reach <- max(abs(c(boxes$xmin, boxes$xmax, layout$.label_anchor_x)), na.rm = TRUE) / plot_xmax

  overflow * 1000 +
    axis_overlap * 90 +
    point_overlap * 35 +
    side_density * 0.07 +
    edge_reach +
    inward_penalty
}

label_text_boxes <- function(layout,
                             plot_xmax,
                             ymax,
                             label_col,
                             layout_plan) {
  labels <- strsplit(as.character(layout[[label_col]]), "\n", fixed = TRUE)
  line_widths <- vapply(
    labels,
    function(lines) max(nchar(lines, type = "width"), na.rm = TRUE),
    numeric(1)
  )
  text_width <- pmax(line_widths, 1) * plot_xmax * layout_plan$char_width_fraction
  line_count <- pmax(lengths(labels), 1)
  text_height <- ymax * (0.026 + 0.018 * pmin(line_count - 1, 3))

  xmin <- ifelse(layout$.label_hjust == 1, layout$.label_text_x - text_width, layout$.label_text_x)
  xmax <- ifelse(layout$.label_hjust == 1, layout$.label_text_x, layout$.label_text_x + text_width)
  data.frame(
    xmin = xmin,
    xmax = xmax,
    ymin = layout$.label_y - text_height,
    ymax = layout$.label_y + text_height
  )
}

count_label_point_overlaps <- function(boxes, point_data) {
  if (nrow(point_data) == 0 || nrow(boxes) == 0) {
    return(0)
  }
  sum(vapply(
    seq_len(nrow(boxes)),
    function(i) {
      sum(
        point_data$.volcano_x >= boxes$xmin[i] &
          point_data$.volcano_x <= boxes$xmax[i] &
          point_data$.volcano_y >= boxes$ymin[i] &
          point_data$.volcano_y <= boxes$ymax[i],
        na.rm = TRUE
      )
    },
    numeric(1)
  ))
}

count_label_point_density <- function(boxes, point_data) {
  if (nrow(point_data) == 0 || nrow(boxes) == 0) {
    return(0)
  }
  sum(vapply(
    seq_len(nrow(boxes)),
    function(i) {
      sum(
        point_data$.volcano_x >= boxes$xmin[i] &
          point_data$.volcano_x <= boxes$xmax[i],
        na.rm = TRUE
      )
    },
    numeric(1)
  ))
}

is_auto_value <- function(value) {
  is.character(value) && length(value) == 1 && identical(tolower(trimws(value)), "auto")
}

parse_positive_integer <- function(value, name) {
  value <- suppressWarnings(as.numeric(value[1]))
  if (!is.finite(value) || value <= 0) {
    stop(name, " must be a positive number or 'auto'", call. = FALSE)
  }
  as.integer(ceiling(value))
}

disambiguate_wrapped_labels <- function(labels, width = NULL) {
  labels <- as.character(labels)
  duplicated_groups <- split(seq_along(labels), labels)
  duplicated_groups <- duplicated_groups[lengths(duplicated_groups) > 1]
  if (length(duplicated_groups) == 0) {
    return(labels)
  }

  for (idx in duplicated_groups) {
    for (rank in seq_along(idx)) {
      labels[idx[rank]] <- append_label_suffix(labels[idx[rank]], paste0(" #", rank), width = width)
    }
  }
  labels
}

append_label_suffix <- function(label, suffix, width = NULL) {
  lines <- strsplit(label, "\n", fixed = TRUE)[[1]]
  last <- length(lines)
  width <- if (is.null(width) || length(width) == 0) NA_real_ else suppressWarnings(as.numeric(width[1]))
  if (is.finite(width) && width > nchar(suffix, type = "width") + 4) {
    ellipsis <- "..."
    line <- sub("\\.\\.\\.$", "", lines[last])
    keep_width <- max(width - nchar(suffix, type = "width") - nchar(ellipsis, type = "width"), 1)
    if (nchar(lines[last], type = "width") + nchar(suffix, type = "width") > width) {
      lines[last] <- paste0(substr(line, 1, keep_width), ellipsis, suffix)
      return(paste(lines, collapse = "\n"))
    }
  }
  if (grepl("\\.\\.\\.$", lines[last])) {
    lines[last] <- sub("\\.\\.\\.$", paste0("...", suffix), lines[last])
  } else {
    lines[last] <- paste0(lines[last], suffix)
  }
  paste(lines, collapse = "\n")
}
