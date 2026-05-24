test_that("package exports only plotting-oriented functions", {
  expect_setequal(
    getNamespaceExports("volcanolabel"),
    c(
      "compute_outside_label_layout",
      "save_volcano",
      "theme_volcano",
      "volcano_palette",
      "volcano_plot",
      "wrap_volcano_labels"
    )
  )
  expect_false(exists("read_diff_table", asNamespace("volcanolabel"), inherits = FALSE))
  expect_false(exists("prepare_volcano_data", asNamespace("volcanolabel"), inherits = FALSE))
  expect_false(exists("select_volcano_labels", asNamespace("volcanolabel"), inherits = FALSE))
  expect_false(exists("volcano_cli", asNamespace("volcanolabel"), inherits = FALSE))
})

test_that("volcano_plot consumes already prepared plotting columns", {
  data <- mini_volcano_data()

  plot <- volcano_plot(
    data,
    x = "log2FC",
    y = "neg_log10_p",
    label = "plot_label",
    color = "direction",
    x_cutoff = log2(1.2),
    y_cutoff = -log10(0.05)
  )

  expect_s3_class(plot, "ggplot")
  expect_equal(plot$labels$x, "log2FC")
  expect_equal(plot$labels$y, "neg_log10_p")
  expect_equal(plot$coordinates$clip, "off")
})

test_that("volcano_plot uses y and color columns without deriving analysis fields", {
  data <- data.frame(
    gene = c("A", "B", "C"),
    log2FC = c(-2, 0.2, 2),
    score = c(4, 1, 3),
    display = c("left gene", NA, "right gene"),
    cohort = c("case", "control", "case"),
    stringsAsFactors = FALSE
  )

  plot <- volcano_plot(
    data,
    x = "log2FC",
    y = "score",
    label = "display",
    color = "cohort",
    x_cutoff = 1,
    y_cutoff = 2
  )
  built <- ggplot2::ggplot_build(plot)

  expect_equal(built$data[[1]]$y, data$score)
  expect_equal(as.character(plot$data$regulation), data$cohort)
  expect_false(".volcano_p" %in% names(plot$data))
})

test_that("volcano_plot no longer accepts analysis-style p-value arguments", {
  expect_error(
    volcano_plot(
      mini_volcano_data(),
      x = "log2FC",
      p = "pvalue",
      gene = "gene",
      top_n = 1
    ),
    "unused argument"
  )
})
