stress_case <- function(kind, n = 800, seed = 42) {
  set.seed(seed)
  gene <- sprintf("Gene_%s_%04d", kind, seq_len(n))
  log2fc <- stats::rnorm(n, sd = 0.35)
  pvalue <- stats::runif(n, min = 0.08, max = 1)

  if (kind == "balanced_dense") {
    up <- seq_len(80)
    down <- 81:160
    log2fc[up] <- stats::rnorm(length(up), mean = 1.5, sd = 0.35)
    log2fc[down] <- stats::rnorm(length(down), mean = -1.4, sd = 0.35)
    pvalue[c(up, down)] <- 10^stats::runif(length(c(up, down)), min = -7, max = -2)
  } else if (kind == "no_significant") {
    log2fc <- stats::rnorm(n, sd = 0.05)
    pvalue <- stats::runif(n, min = 0.25, max = 1)
  } else if (kind == "up_only") {
    up <- seq_len(90)
    log2fc[up] <- stats::rnorm(length(up), mean = 1.8, sd = 0.25)
    pvalue[up] <- 10^stats::runif(length(up), min = -6, max = -2)
  } else if (kind == "down_only") {
    down <- seq_len(90)
    log2fc[down] <- stats::rnorm(length(down), mean = -1.8, sd = 0.25)
    pvalue[down] <- 10^stats::runif(length(down), min = -6, max = -2)
  } else if (kind == "long_labels") {
    hit <- seq_len(28)
    gene[hit] <- paste0("VeryLongGeneName-With/Many;Separators_", hit, "_AlphaBetaGammaDelta")
    log2fc[hit] <- rep(c(-1.5, 1.5), length.out = length(hit))
    pvalue[hit] <- 10^stats::runif(length(hit), min = -5, max = -2)
  } else if (kind == "extreme_p") {
    hit <- seq_len(40)
    log2fc[hit] <- rep(c(-2, 2), length.out = length(hit))
    pvalue[hit] <- c(0, rep(.Machine$double.xmin, 3), 10^stats::runif(length(hit) - 4, min = -20, max = -8))
  } else if (kind == "few_long_labels") {
    hit <- seq_len(6)
    gene[hit] <- paste0("SparseVeryLongGeneName-With/Useful;Detail_", hit, "_AlphaBetaGammaDelta")
    log2fc[hit] <- rep(c(-1.6, 1.6), length.out = length(hit))
    pvalue[hit] <- 10^stats::runif(length(hit), min = -5, max = -3)
  } else if (kind == "edge_hugging") {
    hit <- seq_len(12)
    gene[hit] <- paste0("EdgeGene_", hit)
    log2fc[hit] <- rep(c(-8.4, 8.5), length.out = length(hit)) + stats::rnorm(length(hit), sd = 0.08)
    pvalue[hit] <- 10^stats::runif(length(hit), min = -6, max = -3)
  } else {
    stop("unknown stress case: ", kind)
  }

  make_plot_ready(
    data.frame(gene = gene, log2FC = log2fc, pvalue = pvalue, padjs = pvalue, stringsAsFactors = FALSE),
    labels = head(gene, 3)
  )
}

test_that("volcano_plot handles common synthetic edge-case datasets", {
  for (kind in c("balanced_dense", "no_significant", "up_only", "down_only", "long_labels", "extreme_p", "few_long_labels", "edge_hugging")) {
    data <- stress_case(kind, n = 900)
    data$plot_label <- ifelse(data$gene %in% c(head(data$gene, 3), data$gene[seq_len(min(20, nrow(data)))]), data$gene, NA_character_)
    plot <- volcano_plot(
      data,
      x = "log2FC",
      y = "neg_log10_p",
      label = "plot_label",
      color = "direction",
      x_cutoff = log2(1.2),
      y_cutoff = -log10(0.05),
      title = kind
    )
    path <- tempfile(fileext = ".png")
    save_volcano(plot, path, width = 6.4, height = 4.8, dpi = 90)

    expect_s3_class(plot, "ggplot")
    expect_true(file.exists(path))
    expect_gt(file.info(path)$size, 1000)
  }
})

test_that("long-label stress case is capped to compact label text", {
  data <- stress_case("long_labels", n = 300)
  selected <- c(head(data$gene, 3), data$gene[seq_len(12)])
  prepared <- make_layout_ready(data, labels = selected)
  label_data <- prepared[!is.na(prepared$.volcano_label), , drop = FALSE]

  layout <- compute_outside_label_layout(
    label_data,
    plot_xmax = 5,
    ymax = 8,
    wrap_width = 18,
    max_lines = 2
  )

  label_lines <- lengths(strsplit(layout$.volcano_label, "\n", fixed = TRUE))
  expect_lte(max(label_lines), 2)
  expect_true(any(grepl("...", layout$.volcano_label, fixed = TRUE)))
  expect_equal(anyDuplicated(layout$.volcano_label), 0)
})
