test_that("wrap_volcano_labels wraps separators and long words cleanly", {
  wrapped <- wrap_volcano_labels("AlphaBetaGamma-Delta/Epsilon", width = 10)

  expect_match(wrapped, "\n", fixed = TRUE)
  expect_false(grepl(" \n", wrapped, fixed = TRUE))
})

test_that("wrap_volcano_labels can cap very long labels with ellipsis", {
  wrapped <- wrap_volcano_labels(
    "VeryLongGeneName-With/Many;Separators_19_AlphaBetaGammaDelta",
    width = 18,
    max_lines = 2
  )

  expect_lte(length(strsplit(wrapped, "\n", fixed = TRUE)[[1]]), 2)
  expect_match(wrapped, "...", fixed = TRUE)
})

test_that("compute_outside_label_layout creates separated outside columns", {
  prepared <- make_layout_ready(mini_volcano_data())
  label_data <- prepared[prepared$gene %in% c("Vip", "Sst", "Gad1", "Slc17a7"), , drop = FALSE]

  layout <- compute_outside_label_layout(label_data, plot_xmax = 6, ymax = 20)

  expect_true(all(c(".label_anchor_x", ".label_text_x", ".label_y", ".label_hjust", ".label_side") %in% names(layout)))
  expect_equal(layout$.label_side[layout$gene == "Sst"], "left")
  expect_equal(layout$.label_side[layout$gene == "Vip"], "right")
  expect_lte(length(unique(layout$.label_hjust[layout$.label_side == "right"])), 1)
  expect_lte(length(unique(layout$.label_hjust[layout$.label_side == "left"])), 1)
  expect_equal(anyDuplicated(paste(layout$.label_side, round(layout$.label_y, 6))), 0)
  expect_lte(max(abs(layout$.label_anchor_x)), 6)
})

test_that("default left labels use one text direction and clear the axis", {
  prepared <- make_layout_ready(mini_volcano_data(), labels = c("Sst", "Pvalb"))
  label_data <- prepared[prepared$gene %in% c("Sst", "Pvalb"), , drop = FALSE]

  layout <- compute_outside_label_layout(label_data, plot_xmax = 6, ymax = 20)
  left_layout <- layout[layout$.label_side == "left", , drop = FALSE]

  expect_gt(nrow(left_layout), 0)
  expect_lte(length(unique(left_layout$.label_hjust)), 1)
  expect_true(all(sign(left_layout$.label_text_x - left_layout$.label_anchor_x) == sign(left_layout$.label_text_x[1] - left_layout$.label_anchor_x[1])))
  expect_true(all(left_layout$.label_text_x < -6 * 0.18))
  expect_true(all(left_layout$.label_text_x > -6 * 0.86))
})

test_that("default right labels use one text direction and stay inside the device", {
  prepared <- make_layout_ready(mini_volcano_data(), labels = c("Vip", "Slc17a7"))
  label_data <- prepared[prepared$gene %in% c("Vip", "Slc17a7"), , drop = FALSE]

  layout <- compute_outside_label_layout(label_data, plot_xmax = 6, ymax = 20)
  right_layout <- layout[layout$.label_side == "right", , drop = FALSE]

  expect_gt(nrow(right_layout), 0)
  expect_lte(length(unique(right_layout$.label_hjust)), 1)
  expect_true(all(sign(right_layout$.label_text_x - right_layout$.label_anchor_x) == sign(right_layout$.label_text_x[1] - right_layout$.label_anchor_x[1])))
  expect_true(all(right_layout$.label_text_x > 6 * 0.18))
  expect_true(all(right_layout$.label_text_x < 6 * 0.86))
})

test_that("adaptive layout permits more detail for a few long labels", {
  prepared <- make_layout_ready(mini_volcano_data(), labels = c("Vip", "Sst"))
  label_data <- prepared[prepared$gene %in% c("Vip", "Sst"), , drop = FALSE]
  label_data$.volcano_label <- paste0(
    label_data$gene,
    "-LongSegment/NeuronMarker;AlphaBetaGammaDelta"
  )

  layout <- compute_outside_label_layout(
    label_data,
    plot_xmax = 7,
    ymax = 10,
    wrap_width = "auto",
    max_lines = "auto"
  )

  expect_equal(unique(layout$.label_max_lines), 2)
  expect_gte(unique(layout$.label_wrap_width), 18)
  expect_true(any(grepl("\n", layout$.volcano_label, fixed = TRUE)))
})

test_that("adaptive layout keeps dense long labels compact and distinct", {
  dense <- data.frame(
    gene = paste0("VeryLongGeneName-With/Many;Separators_", seq_len(40), "_AlphaBetaGammaDelta"),
    log2FC = rep(c(-1.6, 1.6), length.out = 40),
    pvalue = 10^seq(-5, -2, length.out = 40),
    stringsAsFactors = FALSE
  )
  selected <- c(head(dense$gene, 4), dense$gene[seq(9, 24)])
  prepared <- make_layout_ready(dense, labels = selected)
  label_data <- prepared[!is.na(prepared$.volcano_label), , drop = FALSE]

  layout <- compute_outside_label_layout(
    label_data,
    plot_xmax = 6,
    ymax = 8,
    wrap_width = "auto",
    max_lines = "auto"
  )

  expect_equal(unique(layout$.label_max_lines), 1)
  expect_lte(unique(layout$.label_wrap_width), 24)
  expect_equal(anyDuplicated(layout$.volcano_label), 0)
})

test_that("adaptive anchors move outside edge-hugging labeled points", {
  edge <- data.frame(
    gene = c("RightEdge", "LeftEdge", "Center"),
    log2FC = c(8.5, -8.4, 0),
    pvalue = c(1e-5, 2e-5, 0.5),
    stringsAsFactors = FALSE
  )
  prepared <- make_layout_ready(edge, labels = edge$gene)

  layout <- compute_outside_label_layout(
    prepared[prepared$gene != "Center", , drop = FALSE],
    plot_xmax = 10,
    ymax = 8,
    wrap_width = "auto",
    max_lines = "auto"
  )
  right_layout <- layout[layout$.label_side == "right", , drop = FALSE]
  left_layout <- layout[layout$.label_side == "left", , drop = FALSE]

  expect_gt(right_layout$.label_anchor_x, right_layout$.volcano_x)
  expect_lt(left_layout$.label_anchor_x, left_layout$.volcano_x)
  expect_lt(abs(right_layout$.label_anchor_x), 10)
  expect_lt(abs(left_layout$.label_anchor_x), 10)
})

test_that("auto side layout can put an edge-hugging side inside to reduce wasted width", {
  edge <- data.frame(
    gene = paste0("RightEdge", seq_len(4)),
    log2FC = c(8.7, 8.8, 8.9, 9.0),
    pvalue = c(1e-7, 2e-7, 3e-7, 4e-7),
    stringsAsFactors = FALSE
  )
  prepared <- make_layout_ready(edge, labels = edge$gene)

  layout <- compute_outside_label_layout(
    prepared,
    plot_xmax = 10,
    ymax = 8,
    point_data = prepared,
    text_side_right = "auto"
  )
  right_layout <- layout[layout$.label_side == "right", , drop = FALSE]

  expect_equal(unique(right_layout$.label_hjust), 1)
  expect_true(all(right_layout$.label_text_x < right_layout$.label_anchor_x))
  expect_true(all(right_layout$.label_anchor_x > right_layout$.volcano_x))
})

test_that("auto side layout keeps a side outside when inside text would cross a point cloud", {
  labeled <- data.frame(
    gene = paste0("RightLabel", seq_len(4)),
    log2FC = c(0.75, 0.80, 0.85, 0.90),
    pvalue = c(1e-7, 2e-7, 3e-7, 4e-7),
    stringsAsFactors = FALSE
  )
  cloud <- data.frame(
    gene = paste0("cloud", seq_len(80)),
    log2FC = rep(seq(1.55, 2.25, length.out = 20), 4),
    pvalue = rep(10^-seq(6.7, 7.5, length.out = 4), each = 20),
    stringsAsFactors = FALSE
  )
  prepared <- make_layout_ready(rbind(labeled, cloud), labels = labeled$gene)
  label_data <- prepared[prepared$gene %in% labeled$gene, , drop = FALSE]

  layout <- compute_outside_label_layout(
    label_data,
    plot_xmax = 5,
    ymax = 8,
    point_data = prepared,
    text_side_right = "auto"
  )
  right_layout <- layout[layout$.label_side == "right", , drop = FALSE]

  expect_equal(unique(right_layout$.label_hjust), 0)
  expect_true(all(right_layout$.label_text_x > right_layout$.label_anchor_x))
})

test_that("demo-like left labels clear the y-axis safety zone", {
  data <- utils::read.table(
    system.file("extdata", "gene_expression.txt", package = "volcanolabel"),
    header = TRUE,
    sep = "\t",
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  selected <- c("Vip", "Sst", "Gad1", "Gad2", "Slc17a7", "Pvalb", "Rps6ka2")
  prepared <- make_layout_ready(data, labels = selected)
  label_data <- prepared[!is.na(prepared$.volcano_label), , drop = FALSE]
  plot_xmax <- 7.8
  ymax <- max(prepared$.volcano_y) * 1.06

  layout <- compute_outside_label_layout(label_data, plot_xmax = plot_xmax, ymax = ymax)
  rps6ka2 <- layout[layout$gene == "Rps6ka2", , drop = FALSE]

  expect_equal(nrow(rps6ka2), 1)
  expect_gt(rps6ka2$.label_text_x, -plot_xmax * 0.76)
  expect_lt(rps6ka2$.label_text_x, -plot_xmax * 0.18)
  expect_equal(rps6ka2$.label_hjust, 1)
  expect_lt(rps6ka2$.label_text_x, rps6ka2$.label_anchor_x)
})
