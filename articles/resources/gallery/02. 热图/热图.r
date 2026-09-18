input_tsv <- "示例数据.tsv"
output_svg <- "heatmap.svg"
x_label_angle <- 45
y_label_angle <- 0
x_axis_tick_cex <- 2
y_axis_tick_cex <- 2
colorbar_tick_cex <- 1.5

## Continuous palette. The heatmap itself is drawn as vector cells so adjacent
## cells are not smoothed/interpolated in SVG viewers.
heatmap_color_anchors <- c("#2C7BB6", "#ABD9E9", "#FFFFBF", "#FDAE61", "#D7191C")
heatmap_colors <- colorRampPalette(heatmap_color_anchors)(100)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) >= 1) input_tsv <- args[1]
if (length(args) >= 2) output_svg <- args[2]
if (length(args) >= 3) x_axis_tick_cex <- as.numeric(args[3])
if (length(args) >= 4) y_axis_tick_cex <- as.numeric(args[4])
if (length(args) >= 5) colorbar_tick_cex <- as.numeric(args[5])

if (anyNA(c(x_axis_tick_cex, y_axis_tick_cex, colorbar_tick_cex))) {
  stop("Axis and colorbar tick font size parameters must be numeric.")
}

mat_df <- read.delim(
  input_tsv,
  header = TRUE,
  row.names = 1,
  check.names = FALSE,
  na.strings = c("", "NA", "NaN"),
  sep = "\t"
)

mat <- as.matrix(mat_df)
storage.mode(mat) <- "numeric"

if (anyNA(dimnames(mat)) || is.null(rownames(mat)) || is.null(colnames(mat))) {
  stop("Input TSV must contain column names in the first row and row names in the first column.")
}

if (nrow(mat) == 0 || ncol(mat) == 0) {
  stop("Input matrix is empty. Please check the TSV file.")
}

finite_values <- mat[is.finite(mat)]
if (length(finite_values) == 0) {
  stop("No numeric values found in the matrix.")
}

zlim <- range(finite_values, na.rm = TRUE)
if (zlim[1] == zlim[2]) {
  zlim <- zlim + c(-0.5, 0.5)
}

color_breaks <- seq(zlim[1], zlim[2], length.out = length(heatmap_colors) + 1)
na_color <- "#F2F2F2"

draw_axis_labels <- function(side, at, labels, angle, cex = 0.8) {
  usr <- par("usr")
  x_offset <- strheight("M", cex = cex) * 1.4
  y_offset <- strwidth("M", cex = cex) * 1.4

  if (side == 1) {
    axis(side = 1, at = at, labels = FALSE, tick = FALSE)
    adj <- if (angle == 0) c(0.5, 1) else c(1, 1)
    text(
      x = at,
      y = usr[3] - x_offset,
      labels = labels,
      srt = angle,
      adj = adj,
      xpd = NA,
      cex = cex
    )
  }

  if (side == 2) {
    axis(side = 2, at = at, labels = FALSE, tick = FALSE)
    adj <- if (angle == 0) c(1, 0.5) else c(1, 1)
    text(
      x = usr[1] - y_offset,
      y = at,
      labels = labels,
      srt = angle,
      adj = adj,
      xpd = NA,
      cex = cex
    )
  }
}

svg(output_svg, width = 8, height = 7, pointsize = 10)

layout(
  matrix(c(1, 2), nrow = 1),
  widths = c(5.5, 0.7)
)

par(mar = c(7, 8, 3, 1))

plot_mat <- t(mat[nrow(mat):1, , drop = FALSE])
image(
  x = seq_len(ncol(mat)),
  y = seq_len(nrow(mat)),
  z = plot_mat,
  col = heatmap_colors,
  breaks = color_breaks,
  axes = FALSE,
  xlab = "",
  ylab = "",
  useRaster = FALSE
)

na_pos <- which(is.na(mat), arr.ind = TRUE)
if (nrow(na_pos) > 0) {
  rect(
    xleft = na_pos[, "col"] - 0.5,
    ybottom = nrow(mat) - na_pos[, "row"] + 0.5,
    xright = na_pos[, "col"] + 0.5,
    ytop = nrow(mat) - na_pos[, "row"] + 1.5,
    col = na_color,
    border = NA
  )
}

draw_axis_labels(1, seq_len(ncol(mat)), colnames(mat), x_label_angle, x_axis_tick_cex)
draw_axis_labels(2, seq_len(nrow(mat)), rev(rownames(mat)), y_label_angle, y_axis_tick_cex)

box()
title(main = tools::file_path_sans_ext(basename(input_tsv)), line = 1)

par(mar = c(7, 1, 3, 4))
legend_midpoints <- (color_breaks[-1] + color_breaks[-length(color_breaks)]) / 2
image(
  x = 1,
  y = legend_midpoints,
  z = matrix(seq_along(heatmap_colors), nrow = 1),
  col = heatmap_colors,
  axes = FALSE,
  xlab = "",
  ylab = "",
  useRaster = FALSE
)
legend_ticks <- pretty(zlim, n = 5)
legend_ticks <- legend_ticks[legend_ticks >= zlim[1] & legend_ticks <= zlim[2]]
axis(side = 4, at = legend_ticks, las = 2, cex.axis = colorbar_tick_cex)
mtext("Value", side = 4, line = 2.5)
box()

dev.off()

message("Heatmap saved to: ", normalizePath(output_svg, mustWork = FALSE))
