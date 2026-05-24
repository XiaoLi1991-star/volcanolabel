mini_volcano_data <- function() {
  out <- data.frame(
    gene = c("Viren", "Sorin", "Pavon", "Gavon1", "Gavon2", "Laxen7", "Quiet"),
    log2FC = c(2.1, -1.8, 0.24, 0.05, -0.40, 1.25, 0),
    pvalue = c(1e-6, 2e-5, 0, 0.3, 0.001, 0.02, 0.8),
    padjs = c(1e-6, 2e-5, 0, 0.3, 0.001, 0.02, 0.8),
    stringsAsFactors = FALSE
  )
  out$neg_log10_p <- -log10(pmax(out$pvalue, .Machine$double.xmin))
  out$direction <- c("Up", "Down", "Normal", "Normal", "Down", "Up", "Normal")
  out$plot_label <- ifelse(out$gene %in% c("Viren", "Sorin", "Gavon1", "Laxen7"), out$gene, NA_character_)
  out
}

demo_label_genes <- function() {
  c(
    "Mefra", "Ashor", "Cavrel", "Halen3", "Kavo2", "Cirob2", "Ralon6", "Esvam", "Pavon",
    "Ervon", "Hivra", "Ilven1", "Tavrel", "Auvon", "Dorin", "Pavri4", "Corin",
    "Viren", "Gavon2", "Sorin", "Laxen7", "Gavon1"
  )
}

make_plot_ready <- function(data,
                            x_cutoff = log2(1.2),
                            y_cutoff = -log10(0.05),
                            labels = head(data$gene, 6)) {
  data$neg_log10_p <- -log10(pmax(data$pvalue, .Machine$double.xmin))
  data$direction <- ifelse(
    data$log2FC >= x_cutoff & data$neg_log10_p >= y_cutoff,
    "Up",
    ifelse(data$log2FC <= -x_cutoff & data$neg_log10_p >= y_cutoff, "Down", "Normal")
  )
  data$plot_label <- ifelse(data$gene %in% labels, data$gene, NA_character_)
  data
}

make_layout_ready <- function(data,
                              x_cutoff = log2(1.2),
                              y_cutoff = -log10(0.05),
                              labels = head(data$gene, 6)) {
  data <- make_plot_ready(data, x_cutoff = x_cutoff, y_cutoff = y_cutoff, labels = labels)
  data$.volcano_x <- data$log2FC
  data$.volcano_y <- data$neg_log10_p
  data$regulation <- data$direction
  data$.volcano_label <- data$plot_label
  data
}
