library(ggplot2)
library(grid)

devtools::load_all(".", quiet = TRUE)

asset_dir <- file.path("docs", "wechat", "assets")
dir.create(asset_dir, recursive = TRUE, showWarnings = FALSE)
unlink(list.files(asset_dir, pattern = "\\.(png|jpg)$|IMAGE_MANIFEST\\.tsv$", full.names = TRUE))

font_family <- ""
main_plot_path <- file.path("man", "figures", "README-auto-ring-redraw-v2.png")
article_main_path <- file.path(asset_dir, "article-01-main-volcano-1080x788.png")
article_before_after_path <- file.path(asset_dir, "article-02-before-after-v2-1080x720.png")
cover_png_path <- file.path(asset_dir, "wechat-cover-v2-900x383.png")
cover_jpg_path <- file.path(asset_dir, "wechat-cover-v2-900x383.jpg")
cover_retina_path <- file.path(asset_dir, "wechat-cover-v2-retina-1800x766.png")
square_png_path <- file.path(asset_dir, "wechat-square-v2-1080x1080.png")
square_jpg_path <- file.path(asset_dir, "wechat-square-v2-1080x1080.jpg")

x_cutoff <- log2(1.2)
y_cutoff <- -log10(0.05)
labels <- c(
  "Mefra", "Ashor", "Cavrel", "Halen3", "Kavo2", "Cirob2", "Ralon6", "Esvam", "Pavon",
  "Ervon", "Hivra", "Ilven1", "Tavrel", "Auvon", "Dorin", "Pavri4", "Corin",
  "Viren", "Gavon2", "Sorin", "Laxen7", "Gavon1"
)
custom_palette <- c(Up = "#B42318", Normal = "#B8C2CC", Down = "#0077B6")

expr <- read.table(
  file.path("inst", "extdata", "gene_expression.txt"),
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE,
  check.names = FALSE
)
expr$neg_log10_p <- -log10(pmax(expr$pvalue, .Machine$double.xmin))
expr$Regulation <- ifelse(
  expr$log2FC >= x_cutoff & expr$neg_log10_p >= y_cutoff,
  "Up",
  ifelse(expr$log2FC <= -x_cutoff & expr$neg_log10_p >= y_cutoff, "Down", "Normal")
)
expr$plot_label <- ifelse(expr$gene %in% labels, expr$gene, NA_character_)

cover_labels <- c("Mefra", "Ashor", "Cirob2", "Ervon", "Dorin", "Corin", "Viren", "Gavon1")
cover_expr <- expr
cover_expr$cover_label <- ifelse(cover_expr$gene %in% cover_labels, cover_expr$gene, NA_character_)

cover_plot <- volcano_plot(
  cover_expr,
  x = "log2FC",
  y = "neg_log10_p",
  label = "cover_label",
  color = "Regulation",
  palette = custom_palette,
  x_cutoff = x_cutoff,
  y_cutoff = y_cutoff,
  point_size = 1.15,
  point_alpha = 0.76,
  label_size = 2.8,
  label_segment_size = 0.28,
  label_segment_alpha = 0.82,
  label_point_ring_color = "auto",
  label_point_ring_size = 2.35,
  label_point_ring_halo_size = 3.0,
  label_anchor_x_left = -3.7,
  label_anchor_x_right = 3.8,
  plot_xmax = 4.85,
  ymax = 6.25,
  xlab = NULL,
  ylab = NULL,
  base_size = 8.5,
  legend_position = "none"
) +
  ggplot2::theme(
    axis.title = ggplot2::element_blank(),
    axis.text = ggplot2::element_blank(),
    axis.ticks = ggplot2::element_blank(),
    axis.line = ggplot2::element_blank(),
    plot.title = ggplot2::element_blank(),
    plot.subtitle = ggplot2::element_blank(),
    plot.margin = ggplot2::margin(0, 0, 0, 0)
  )

save_png <- function(path, width, height, draw, res = 144) {
  ragg::agg_png(path, width = width, height = height, units = "px", res = res, background = "white")
  on.exit(grDevices::dev.off(), add = TRUE)
  grid.newpage()
  draw()
}

save_jpg <- function(path, width, height, draw, res = 144) {
  grDevices::jpeg(path, width = width, height = height, units = "px", res = res, quality = 94)
  on.exit(grDevices::dev.off(), add = TRUE)
  grid.newpage()
  draw()
}

save_png(article_main_path, 1080, 788, function() {
  img <- png::readPNG(main_plot_path)
  grid.raster(img, x = 0.5, y = 0.5, width = 1, height = 1, interpolate = TRUE)
})

draw_round_rect <- function(x, y, width, height, fill, col = NA, r = unit(8, "pt")) {
  grid.roundrect(x, y, width, height, r = r, gp = gpar(fill = fill, col = col))
}

draw_ring_point <- function(x, y, fill, ring, size = 0.035) {
  grid.circle(x, y, r = unit(size * 1.55, "npc"), gp = gpar(fill = NA, col = "#E7ECF2", lwd = 6))
  grid.circle(x, y, r = unit(size * 1.20, "npc"), gp = gpar(fill = NA, col = ring, lwd = 3))
  grid.circle(x, y, r = unit(size * 0.54, "npc"), gp = gpar(fill = fill, col = fill, lwd = 1))
}

draw_cover_plot_panel <- function(x = 0.805, y = 0.52, width = 0.35, height = 0.74) {
  pushViewport(viewport(x = x, y = y, width = width, height = height, clip = "on"))
  grid.draw(ggplotGrob(cover_plot))
  popViewport()
}

draw_cover <- function(scale = 1) {
  grid.rect(gp = gpar(fill = "#F7FAFC", col = NA))
  grid.rect(x = 0.805, y = 0.5, width = 0.39, height = 1, gp = gpar(fill = "#FFFFFF", col = NA))
  grid.text(
    "volcanolabel",
    x = 0.07,
    y = 0.74,
    just = c("left", "center"),
    gp = gpar(col = "#17202B", fontsize = 34 * scale, fontface = "bold", fontfamily = font_family)
  )
  grid.text(
    "把火山图标签稳稳放好",
    x = 0.07,
    y = 0.57,
    just = c("left", "center"),
    gp = gpar(col = "#475467", fontsize = 19 * scale, fontfamily = font_family)
  )
  grid.text(
    "自动外侧标签 · 清晰连接线 · 自适应高亮圈",
    x = 0.07,
    y = 0.44,
    just = c("left", "center"),
    gp = gpar(col = "#667085", fontsize = 11.5 * scale, fontfamily = font_family)
  )
  draw_cover_plot_panel()
  grid.text(
    "开源 R 包",
    x = 0.07,
    y = 0.22,
    just = c("left", "center"),
    gp = gpar(col = "#B42318", fontsize = 13 * scale, fontface = "bold", fontfamily = font_family)
  )
}

save_png(cover_png_path, 900, 383, function() draw_cover(1))
save_jpg(cover_jpg_path, 900, 383, function() draw_cover(1))
save_png(cover_retina_path, 1800, 766, function() draw_cover(2))
save_png(square_png_path, 1080, 1080, function() {
  grid.rect(gp = gpar(fill = "#F7FAFC", col = NA))
  pushViewport(viewport(x = 0.5, y = 0.72, width = 0.90, height = 0.34))
  draw_cover(1.05)
  popViewport()
  img <- png::readPNG(main_plot_path)
  grid.raster(img, x = 0.5, y = 0.31, width = 0.84, height = 0.46, interpolate = TRUE)
})
save_jpg(square_jpg_path, 1080, 1080, function() {
  grid.rect(gp = gpar(fill = "#F7FAFC", col = NA))
  pushViewport(viewport(x = 0.5, y = 0.72, width = 0.90, height = 0.34))
  draw_cover(1.05)
  popViewport()
  img <- png::readPNG(main_plot_path)
  grid.raster(img, x = 0.5, y = 0.31, width = 0.84, height = 0.46, interpolate = TRUE)
})
save_png(file.path(asset_dir, "wechat-share-500x400.png"), 500, 400, function() {
  grid.rect(gp = gpar(fill = "#F7FAFC", col = NA))
  grid.text("volcanolabel", 0.08, 0.83, just = c("left", "center"), gp = gpar(fontsize = 28, fontface = "bold", col = "#17202B", fontfamily = font_family))
  grid.text("火山图标签，不再靠运气", 0.08, 0.72, just = c("left", "center"), gp = gpar(fontsize = 15, col = "#475467", fontfamily = font_family))
  img <- png::readPNG(main_plot_path)
  grid.raster(img, x = 0.5, y = 0.34, width = 0.82, height = 0.50, interpolate = TRUE)
})
save_jpg(file.path(asset_dir, "wechat-share-500x400.jpg"), 500, 400, function() {
  grid.rect(gp = gpar(fill = "#F7FAFC", col = NA))
  grid.text("volcanolabel", 0.08, 0.83, just = c("left", "center"), gp = gpar(fontsize = 28, fontface = "bold", col = "#17202B", fontfamily = font_family))
  grid.text("火山图标签，不再靠运气", 0.08, 0.72, just = c("left", "center"), gp = gpar(fontsize = 15, col = "#475467", fontfamily = font_family))
  img <- png::readPNG(main_plot_path)
  grid.raster(img, x = 0.5, y = 0.34, width = 0.82, height = 0.50, interpolate = TRUE)
})

naive <- ggplot(expr, aes(log2FC, neg_log10_p)) +
  geom_point(aes(color = Regulation, fill = Regulation), shape = 21, size = 1.4, alpha = 0.70, stroke = 0.12) +
  geom_label(
    data = expr[!is.na(expr$plot_label), , drop = FALSE],
    aes(label = plot_label),
    fill = "white",
    color = "#1F2937",
    alpha = 0.92,
    linewidth = 0.22,
    size = 3.2
  ) +
  scale_color_manual(values = custom_palette) +
  scale_fill_manual(values = custom_palette) +
  coord_cartesian(xlim = c(-2.5, 2.5), ylim = c(0, 6.2)) +
  labs(title = "常见问题", subtitle = "标签遮挡、指向不清、点被盖住", x = "log2FC", y = "-log10(pvalue)") +
  theme_volcano(base_size = 10, legend_position = "none")

smart <- volcano_plot(
  expr,
  x = "log2FC",
  y = "neg_log10_p",
  label = "plot_label",
  color = "Regulation",
  palette = custom_palette,
  x_cutoff = x_cutoff,
  y_cutoff = y_cutoff,
  label_point_ring_color = "auto",
  label_anchor_x_left = -3.3,
  label_anchor_x_right = 4.55,
  plot_xmax = 5.65,
  ymax = 6.25,
  xlab = "log2FC",
  ylab = "-log10(pvalue)",
  title = "volcanolabel",
  subtitle = "外侧列、肘形连接线、自适应高亮圈",
  base_size = 10,
  legend_position = "none"
)

save_png(article_before_after_path, 1080, 720, function() {
  grid.rect(gp = gpar(fill = "white", col = NA))
  pushViewport(viewport(x = 0.25, y = 0.50, width = 0.47, height = 0.88))
  grid.draw(ggplotGrob(naive))
  popViewport()
  pushViewport(viewport(x = 0.75, y = 0.50, width = 0.47, height = 0.88))
  grid.draw(ggplotGrob(smart))
  popViewport()
})

save_png(file.path(asset_dir, "article-03-ring-rule-1080x720.png"), 1080, 720, function() {
  grid.rect(gp = gpar(fill = "#F7FAFC", col = NA))
  grid.text("label_point_ring_color = \"auto\"", 0.06, 0.86, just = c("left", "center"), gp = gpar(fontsize = 28, fontface = "bold", col = "#17202B", fontfamily = font_family))
  grid.text("圈色从当前点色自动推导：同色系，但更容易从点云里看出来", 0.06, 0.78, just = c("left", "center"), gp = gpar(fontsize = 16, col = "#667085", fontfamily = font_family))
  ring_cols <- vapply(custom_palette, volcanolabel:::adaptive_ring_color, character(1))
  xs <- c(0.22, 0.50, 0.78)
  for (i in seq_along(xs)) {
    draw_round_rect(xs[i], 0.43, 0.22, 0.42, "#FFFFFF", "#E5E7EB")
    draw_ring_point(xs[i], 0.51, unname(custom_palette[i]), unname(ring_cols[i]), size = 0.045)
    grid.text(names(custom_palette)[i], xs[i], 0.31, gp = gpar(fontsize = 22, fontface = "bold", col = "#17202B", fontfamily = font_family))
    grid.text(paste0(unname(custom_palette[i]), " -> ", unname(ring_cols[i])), xs[i], 0.24, gp = gpar(fontsize = 12, col = "#667085", fontfamily = font_family))
  }
})

save_png(file.path(asset_dir, "article-04-api-flow-1080x720.png"), 1080, 720, function() {
  grid.rect(gp = gpar(fill = "#FFFFFF", col = NA))
  grid.text("volcanolabel 只做绘图层", 0.08, 0.86, just = c("left", "center"), gp = gpar(fontsize = 34, fontface = "bold", col = "#17202B", fontfamily = font_family))
  grid.text("像 ggplot2 一样接收已经准备好的列\n不把差异分析打包进绘图函数", 0.08, 0.77, just = c("left", "center"), gp = gpar(fontsize = 15, lineheight = 1.15, col = "#667085", fontfamily = font_family))
  boxes <- data.frame(
    x = c(0.20, 0.50, 0.80),
    title = c("准备绘图列", "volcano_plot()", "保存 / 修改"),
    body = c("log2FC\n-log10(p)\nRegulation\nplot_label", "自动 label 布局\n自适应 ring\n返回 ggplot", "save_volcano()\n或继续 + theme\n出图投稿"),
    stringsAsFactors = FALSE
  )
  for (i in seq_len(nrow(boxes))) {
    draw_round_rect(boxes$x[i], 0.45, 0.22, 0.42, "#F7FAFC", "#D9E0EA")
    grid.text(boxes$title[i], boxes$x[i], 0.57, gp = gpar(fontsize = 21, fontface = "bold", col = "#17202B", fontfamily = font_family))
    grid.text(boxes$body[i], boxes$x[i], 0.42, gp = gpar(fontsize = 15, col = "#475467", lineheight = 1.2, fontfamily = font_family))
  }
  grid.lines(c(0.33, 0.40), c(0.45, 0.45), arrow = arrow(length = unit(0.18, "inches")), gp = gpar(col = "#9AA4B2", lwd = 2))
  grid.lines(c(0.63, 0.70), c(0.45, 0.45), arrow = arrow(length = unit(0.18, "inches")), gp = gpar(col = "#9AA4B2", lwd = 2))
})

files <- list.files(asset_dir, pattern = "\\.(png|jpg)$", full.names = TRUE)
dimensions <- vapply(files, function(path) {
  if (grepl("\\.png$", path, ignore.case = TRUE)) {
    image <- png::readPNG(path)
    paste0(dim(image)[2], "x", dim(image)[1])
  } else {
    info <- file.info(path)
    stem <- basename(path)
    if (grepl("900x383", stem)) "900x383" else if (grepl("500x400", stem)) "500x400" else if (grepl("1080x1080", stem)) "1080x1080" else "unknown"
  }
}, character(1))
use_map <- c(
  "article-01-main-volcano-1080x788.png" = "article inline main figure",
  "article-02-before-after-v2-1080x720.png" = "before/after comparison",
  "article-03-ring-rule-1080x720.png" = "ring color rule explainer",
  "article-04-api-flow-1080x720.png" = "API workflow explainer",
  "wechat-cover-v2-900x383.jpg" = "main WeChat cover JPG",
  "wechat-cover-v2-900x383.png" = "main WeChat cover",
  "wechat-cover-v2-retina-1800x766.png" = "retina WeChat cover",
  "wechat-share-500x400.jpg" = "share card JPG",
  "wechat-share-500x400.png" = "share card",
  "wechat-square-v2-1080x1080.jpg" = "square secondary cover JPG",
  "wechat-square-v2-1080x1080.png" = "square secondary cover"
)
manifest <- data.frame(
  file = basename(files),
  dimensions = unname(dimensions),
  bytes = file.info(files)$size,
  use = unname(use_map[basename(files)]),
  stringsAsFactors = FALSE
)
manifest <- manifest[order(manifest$file), , drop = FALSE]
write.table(manifest, file.path(asset_dir, "IMAGE_MANIFEST.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
print(manifest)
