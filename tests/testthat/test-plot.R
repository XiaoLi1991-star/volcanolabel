test_that("volcano_plot returns a ggplot with polished outside labels", {
  plot <- volcano_plot(
    mini_volcano_data(),
    x = "log2FC",
    y = "neg_log10_p",
    label = "plot_label",
    color = "direction",
    x_cutoff = log2(1.2),
    y_cutoff = -log10(0.05),
    title = "Demo volcano"
  )

  expect_s3_class(plot, "ggplot")
  expect_gte(length(plot$layers), 7)
  expect_equal(plot$labels$title, "Demo volcano")
  expect_equal(plot$coordinates$clip, "off")
  expect_false(any(vapply(plot$layers, function(layer) inherits(layer$geom, "GeomLabel"), logical(1))))
  layer_params <- unlist(lapply(plot$layers, function(layer) layer$aes_params), use.names = FALSE)
  expect_false(any(tolower(as.character(layer_params)) %in% c("white", "#ffffff", "#ffffffff")))
})

test_that("marked point rings can use independent fixed styling", {
  plot <- volcano_plot(
    mini_volcano_data(),
    x = "log2FC",
    y = "neg_log10_p",
    label = "plot_label",
    color = "direction",
    label_point_ring_halo_color = "#E5E7EB",
    label_point_ring_halo_alpha = 0.93,
    label_point_ring_halo_size = 3.7,
    label_point_ring_halo_stroke = 1.2,
    label_point_ring_color = "#111827",
    label_point_ring_alpha = 0.91,
    label_point_ring_size = 3.1,
    label_point_ring_stroke = 0.95
  )

  point_layers <- Filter(function(layer) inherits(layer$geom, "GeomPoint"), plot$layers)
  halo_layer <- point_layers[[length(point_layers) - 1]]
  ring_layer <- point_layers[[length(point_layers)]]

  expect_identical(halo_layer$aes_params$colour, "#E5E7EB")
  expect_equal(halo_layer$aes_params$alpha, 0.93)
  expect_equal(halo_layer$aes_params$size, 3.7)
  expect_equal(halo_layer$aes_params$stroke, 1.2)
  expect_identical(ring_layer$aes_params$colour, "#111827")
  expect_equal(ring_layer$aes_params$alpha, 0.91)
  expect_equal(ring_layer$aes_params$size, 3.1)
  expect_equal(ring_layer$aes_params$stroke, 0.95)
})

test_that("marked point rings use adaptive group-aware auto colors by default", {
  plot <- volcano_plot(
    mini_volcano_data(),
    x = "log2FC",
    y = "neg_log10_p",
    label = "plot_label",
    color = "direction"
  )

  point_layers <- Filter(function(layer) inherits(layer$geom, "GeomPoint"), plot$layers)
  ring_data <- point_layers[[length(point_layers)]]$data

  expect_true(".label_ring_color" %in% names(ring_data))
  expect_equal(ring_data$.label_ring_color[ring_data$.volcano_label == "Gad1"], volcanolabel:::adaptive_ring_color(volcano_palette()[["Normal"]]))
  expect_equal(ring_data$.label_ring_color[ring_data$.volcano_label == "Vip"], volcanolabel:::adaptive_ring_color(volcano_palette()[["Up"]]))
  expect_equal(ring_data$.label_ring_color[ring_data$.volcano_label == "Sst"], volcanolabel:::adaptive_ring_color(volcano_palette()[["Down"]]))
  expect_false(any(ring_data$.label_ring_color %in% unname(volcano_palette())))
  expect_false(any(ring_data$.label_ring_color == "#111827"))
})

test_that("auto rings adapt when point colors change", {
  palette <- c(Up = "#B42318", Normal = "#DDEAF7", Down = "#0077B6")
  plot <- volcano_plot(
    mini_volcano_data(),
    x = "log2FC",
    y = "neg_log10_p",
    label = "plot_label",
    color = "direction",
    palette = palette
  )

  point_layers <- Filter(function(layer) inherits(layer$geom, "GeomPoint"), plot$layers)
  ring_data <- point_layers[[length(point_layers)]]$data
  normal_ring <- ring_data$.label_ring_color[ring_data$.volcano_label == "Gad1"]

  expect_equal(normal_ring, volcanolabel:::adaptive_ring_color(palette[["Normal"]]))
  expect_false(identical(normal_ring, "#6B7280"))
  expect_false(identical(normal_ring, palette[["Normal"]]))
  expect_lt(volcanolabel:::relative_color_luminance(normal_ring), volcanolabel:::relative_color_luminance(palette[["Normal"]]))
  expect_equal(ring_data$.label_ring_color[ring_data$.volcano_label == "Vip"], volcanolabel:::adaptive_ring_color(palette[["Up"]]))
  expect_equal(ring_data$.label_ring_color[ring_data$.volcano_label == "Sst"], volcanolabel:::adaptive_ring_color(palette[["Down"]]))
  expect_false(any(ring_data$.label_ring_color %in% unname(palette)))
})

test_that("marked point rings can still follow group colors when requested", {
  plot <- volcano_plot(
    mini_volcano_data(),
    x = "log2FC",
    y = "neg_log10_p",
    label = "plot_label",
    color = "direction",
    label_point_ring_color = "group"
  )

  point_layers <- Filter(function(layer) inherits(layer$geom, "GeomPoint"), plot$layers)
  ring_data <- point_layers[[length(point_layers)]]$data

  expect_equal(ring_data$.label_ring_color[ring_data$.volcano_label == "Gad1"], volcano_palette()[["Normal"]])
})

test_that("rendered labels do not create a white box when placed directly over points", {
  skip_if_not_installed("png")

  count_white_pixels <- function(plot) {
    path <- tempfile(fileext = ".png")
    ggplot2::ggsave(
      path,
      plot = plot,
      width = 3,
      height = 2,
      dpi = 150,
      bg = "#D73027"
    )
    image <- png::readPNG(path)
    height <- dim(image)[1]
    width <- dim(image)[2]
    crop <- image[
      round(height * 0.36):round(height * 0.64),
      round(width * 0.20):round(width * 0.80),
      1:3,
      drop = FALSE
    ]
    sum(crop[, , 1] > 0.95 & crop[, , 2] > 0.95 & crop[, , 3] > 0.95)
  }

  red_panel <- ggplot2::theme_void() +
    ggplot2::theme(
      panel.background = ggplot2::element_rect(fill = "#D73027", color = NA),
      plot.background = ggplot2::element_rect(fill = "#D73027", color = NA),
      plot.margin = ggplot2::margin(0, 0, 0, 0),
      legend.position = "none"
    )
  base_plot <- ggplot2::ggplot() +
    ggplot2::coord_cartesian(xlim = c(-1, 1), ylim = c(-1, 1), expand = FALSE, clip = "off") +
    red_panel

  control_with_white_box <- base_plot +
    ggplot2::geom_label(
      ggplot2::aes(x = 0, y = 0, label = "WHITEBOXCHECK"),
      inherit.aes = FALSE,
      fill = "white",
      color = "black",
      label.padding = grid::unit(0.18, "lines"),
      label.r = grid::unit(0.05, "lines"),
      size = 10
    )
  expect_gt(count_white_pixels(control_with_white_box), 500)

  label_layout <- data.frame(
    regulation = "Up",
    .volcano_x = 0,
    .volcano_y = 0,
    .label_anchor_x = 0,
    .label_text_x = 0,
    .label_y = 0,
    .label_hjust = 0.5,
    .label_elbow_x = 0,
    .volcano_label = "WHITEBOXCHECK",
    stringsAsFactors = FALSE
  )
  actual_plot <- volcanolabel:::add_outside_label_layers(
    base_plot,
    label_layout,
    palette = volcano_palette(),
    label_size = 10,
    label_color = "#111111",
    label_segment_color = "#9AA4B2",
    label_segment_size = 0.32,
    label_segment_alpha = 0.88,
    label_point_ring_halo_color = "#E5E7EB",
    label_point_ring_halo_alpha = 0.93,
    label_point_ring_halo_size = 3.45,
    label_point_ring_halo_stroke = 1.15,
    label_point_ring_color = "#111827",
    label_point_ring_alpha = 0.95,
    label_point_ring_size = 2.75,
    label_point_ring_stroke = 0.85
  )
  expect_lte(count_white_pixels(actual_plot), 5)
})

test_that("volcano_plot keeps the demo layout compact when labels are short", {
  data <- utils::read.table(
    system.file("extdata", "gene_expression.txt", package = "volcanolabel"),
    header = TRUE,
    sep = "\t",
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  data <- make_plot_ready(
    data,
    labels = c("Vip", "Sst", "Gad1", "Gad2", "Slc17a7", "Pvalb", "Rps6ka2")
  )
  plot <- volcano_plot(
    data,
    x = "log2FC",
    y = "neg_log10_p",
    label = "plot_label",
    color = "direction",
    x_cutoff = log2(1.2),
    y_cutoff = -log10(0.05)
  )

  expect_lt(max(abs(plot$coordinates$limits$x)), 7.2)
})

test_that("volcano_plot tightens edge-hugging labels by allowing inside side text", {
  edge <- data.frame(
    gene = paste0("RightEdge", seq_len(8)),
    log2FC = c(8.7, 8.8, 8.9, 9.0, -8.7, -8.8, -8.9, -9.0),
    pvalue = 10^-seq(5.0, 5.7, length.out = 8),
    stringsAsFactors = FALSE
  )

  edge <- make_plot_ready(edge, labels = edge$gene)
  plot <- volcano_plot(edge, x = "log2FC", y = "neg_log10_p", label = "plot_label", color = "direction")

  expect_lt(max(abs(plot$coordinates$limits$x)), 10.5)
})

test_that("save_volcano writes a non-empty image file", {
  plot <- volcano_plot(mini_volcano_data(), x = "log2FC", y = "neg_log10_p", label = "plot_label", color = "direction")
  path <- tempfile(fileext = ".png")

  written <- save_volcano(plot, path, width = 4, height = 3, dpi = 96)

  expect_equal(normalizePath(written), normalizePath(path))
  expect_true(file.exists(path))
  expect_gt(file.info(path)$size, 1000)
})

test_that("bundled example data can be plotted after user prepares plotting columns", {
  data <- utils::read.table(
    system.file("extdata", "gene_expression.txt", package = "volcanolabel"),
    header = TRUE,
    sep = "\t",
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  data <- make_plot_ready(data, labels = c("Vip", "Sst"))
  plot <- volcano_plot(data, x = "log2FC", y = "neg_log10_p", label = "plot_label", color = "direction")

  expect_true(all(c("gene", "log2FC", "pvalue", "padjs") %in% names(data)))
  expect_gt(nrow(data), 1000)
  expect_s3_class(plot, "ggplot")
})
