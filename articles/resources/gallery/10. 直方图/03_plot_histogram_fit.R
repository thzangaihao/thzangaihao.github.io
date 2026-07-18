#!/usr/bin/env Rscript

# 基于精简分箱 CSV 和 Python 线性拟合 CSV 绘制发表级直方图与拟合曲线。
# 默认使用 svglite 输出 SVG。svglite 会把坐标轴、标题和图内标注保留为
# SVG <text> 元素，导入 Adobe Illustrator 后可以用文字工具继续编辑。

config <- list(
  bins_csv = "histogram_bins.csv",
  fit_curve_csv = "histogram_fit_curve.csv",
  output = "histogram_fit.svg",

  # 坐标系："normal"、"log_y"、"log_log"。
  coord_system = "normal",

  width = 7.2,
  aspect_ratio = 3 / 2,
  dpi = 600,

  x_min = NA_real_,
  x_max = NA_real_,
  y_min = NA_real_,
  y_max = NA_real_,

  # 最左侧柱子与 Y 轴之间的距离，单位是“一个分箱宽度”。
  # 0   = 最左侧柱子贴近 Y 轴；
  # 0.5 = 留出半个分箱宽度；
  # 1   = 留出一个完整分箱宽度；也可以写 1.5、2 等。
  left_padding_bins = 1,

  x_tick_step = NA_real_,
  y_tick_step = NA_real_,

  main = "",
  xlab = "Value",
  ylab = "Frequency",

  fill = "#006AA7",
  border = "white",
  fit_color = "#FECC02",
  bar_line_width = 0.6,
  fit_line_width = 1.5,

  font_family = "sans",
  base_font_size = 10,
  show_equation = TRUE,

  # SVG 文字输出模式：
  # "editable" 使用 svglite，文字在 Illustrator 中可作为文字对象编辑；
  # "base" 使用 grDevices::svg，兼容性好，但在部分软件中可能被拆成路径。
  svg_text_mode = "editable"
)

fmt <- function(x, digits = 4) formatC(x, digits = digits, format = "g")
ticks <- function(lo, hi, step) {
  if (is.na(step)) return(NULL)
  seq(ceiling(lo / step) * step, hi, by = step)
}
get_script_dir <- function() {
  ofile <- tryCatch(sys.frame(1)$ofile, error = function(e) NULL)
  if (!is.null(ofile) && length(ofile) == 1L && nzchar(ofile)) {
    return(dirname(normalizePath(ofile, winslash = "/", mustWork = FALSE)))
  }
  getwd()
}
is_absolute_path <- function(path) grepl("^(~|[A-Za-z]:|/|\\\\\\\\)", path)
script_relative <- function(path, base_dir) {
  if (is.null(path) || !nzchar(path) || is_absolute_path(path)) return(path)
  file.path(base_dir, path)
}
parse_interval_bounds <- function(interval) {
  # interval may look like "[0, 50)", "(50, 100]", or use a Chinese comma.
  # Extracting numbers is more robust than splitting on a comma.
  pattern <- "[-+]?(?:\\d+\\.?\\d*|\\.\\d+)(?:[eE][-+]?\\d+)?"
  matches <- regmatches(interval, gregexpr(pattern, interval, perl = TRUE))
  bounds <- do.call(rbind, lapply(seq_along(matches), function(i) {
    nums <- suppressWarnings(as.numeric(matches[[i]]))
    if (length(nums) < 2L || any(!is.finite(nums[1:2]))) {
      stop(sprintf("无法从 interval 第 %d 行解析左右边界：%s", i, interval[i]), call. = FALSE)
    }
    nums[1:2]
  }))
  colnames(bounds) <- c("bin_start", "bin_end")
  bounds
}

cfg <- config
script_dir <- get_script_dir()
for (field in c("bins_csv", "fit_curve_csv", "output")) {
  cfg[[field]] <- script_relative(cfg[[field]], script_dir)
}
cfg$coord_system <- tolower(cfg$coord_system)
cfg$svg_text_mode <- tolower(cfg$svg_text_mode)
if (!cfg$coord_system %in% c("normal", "log_y", "log_log")) stop("coord_system 只能为 normal、log_y 或 log_log。", call. = FALSE)
if (!cfg$svg_text_mode %in% c("editable", "base")) stop("svg_text_mode 只能为 editable 或 base。", call. = FALSE)
if (!file.exists(cfg$bins_csv)) stop("分箱 CSV 不存在：", cfg$bins_csv, call. = FALSE)
if (!file.exists(cfg$fit_curve_csv)) stop("拟合曲线 CSV 不存在：", cfg$fit_curve_csv, call. = FALSE)

bins <- utils::read.csv(cfg$bins_csv, stringsAsFactors = FALSE, check.names = FALSE, fileEncoding = "UTF-8-BOM")
curve <- utils::read.csv(cfg$fit_curve_csv, stringsAsFactors = FALSE, check.names = FALSE, fileEncoding = "UTF-8-BOM")
required_bins <- c("interval", "count")
required_curve <- c("x", "predicted_count", "equation", "R^2", "p_value")
missing_bins <- setdiff(required_bins, names(bins))
missing_curve <- setdiff(required_curve, names(curve))
if (length(missing_bins)) stop("分箱 CSV 缺少列：", paste(missing_bins, collapse = ", "), call. = FALSE)
if (length(missing_curve)) stop("拟合 CSV 缺少列：", paste(missing_curve, collapse = ", "), call. = FALSE)

interval_bounds <- parse_interval_bounds(bins$interval)
bin_start <- interval_bounds[, "bin_start"]
bin_end <- interval_bounds[, "bin_end"]
count <- as.numeric(bins$count)
bin_width <- stats::median(bin_end - bin_start, na.rm = TRUE)
if (!any(is.finite(bin_start)) || !any(is.finite(bin_end))) {
  stop("无法从 interval 列解析出有效的分箱边界。请确认格式类似 [0, 50)。", call. = FALSE)
}
if (!is.finite(bin_width) || bin_width <= 0) {
  stop("分箱宽度无法从 interval 列推断，请检查 interval 是否按从小到大排列。", call. = FALSE)
}

bar_xlo <- if (is.na(cfg$x_min)) min(bin_start, na.rm = TRUE) else cfg$x_min
bar_xhi <- if (is.na(cfg$x_max)) max(bin_end, na.rm = TRUE) else cfg$x_max
plot_xlo <- bar_xlo - cfg$left_padding_bins * bin_width
if (cfg$coord_system == "log_log" && plot_xlo <= 0) {
  stop("双对数坐标下 X 轴下限必须大于 0；请调小 left_padding_bins 或设置更大的 x_min。", call. = FALSE)
}

log_y <- cfg$coord_system %in% c("log_y", "log_log")
positive_count <- count[count > 0]
curve_y <- as.numeric(curve$predicted_count)
ylo <- if (!is.na(cfg$y_min)) cfg$y_min else if (log_y) max(min(positive_count) / 2, .Machine$double.eps) else 0
valid_curve <- is.finite(curve$x) & is.finite(curve_y) & (!log_y | curve_y > 0)
auto_yhi <- max(c(positive_count, curve_y[valid_curve]), na.rm = TRUE) * if (log_y) 1.35 else 1.08
yhi <- if (is.na(cfg$y_max)) auto_yhi else cfg$y_max
if (log_y && ylo <= 0) stop("对数 Y 轴要求 y_min 大于 0。", call. = FALSE)
if (ylo >= yhi) stop("y_min 必须小于 y_max。", call. = FALSE)

height <- cfg$width / cfg$aspect_ratio
dir.create(dirname(cfg$output), recursive = TRUE, showWarnings = FALSE)
ext <- tolower(tools::file_ext(cfg$output))
if (ext == "svg" && cfg$svg_text_mode == "editable") {
  if (!requireNamespace("svglite", quietly = TRUE)) {
    stop("svg_text_mode='editable' 需要 R 包 svglite。请安装 svglite，或将 svg_text_mode 改为 'base'。", call. = FALSE)
  }
  svglite::svglite(cfg$output, width = cfg$width, height = height,
                   pointsize = cfg$base_font_size, bg = "white")
} else
if (ext == "svg") grDevices::svg(cfg$output, cfg$width, height, pointsize = cfg$base_font_size, family = cfg$font_family, bg = "white") else
if (ext == "pdf") grDevices::pdf(cfg$output, cfg$width, height, pointsize = cfg$base_font_size, family = cfg$font_family, useDingbats = FALSE) else
if (ext == "png") grDevices::png(cfg$output, cfg$width, height, units = "in", res = cfg$dpi, pointsize = cfg$base_font_size, type = "cairo") else
if (ext %in% c("tif", "tiff")) grDevices::tiff(cfg$output, cfg$width, height, units = "in", res = cfg$dpi, compression = "lzw", pointsize = cfg$base_font_size) else
  stop("输出图片仅支持 svg、pdf、png、tif/tiff。", call. = FALSE)

tryCatch({
  graphics::par(family = cfg$font_family, mar = c(4.3, 4.6, if (nzchar(cfg$main)) 2.2 else 0.8, 0.8),
                mgp = c(2.65, 0.75, 0), tcl = -0.25, las = 1, lend = "round",
                cex.axis = 0.92, cex.lab = 1.05, cex.main = 1.08, font.lab = 2, bty = "l")
  graphics::plot(NA,
                 xlim = c(plot_xlo, bar_xhi), ylim = c(ylo, yhi),
                 log = if (cfg$coord_system == "log_y") "y" else if (cfg$coord_system == "log_log") "xy" else "",
                 xaxs = "i", yaxs = "i",
                 xaxt = if (is.na(cfg$x_tick_step)) "s" else "n",
                 yaxt = if (is.na(cfg$y_tick_step) || log_y) "s" else "n",
                 xlab = cfg$xlab, ylab = cfg$ylab, main = cfg$main)

  drawable <- count > 0 & bin_end > bar_xlo & bin_start < bar_xhi
  bottoms <- if (log_y) rep(ylo, length(count)) else rep(0, length(count))
  graphics::rect(pmax(bin_start[drawable], bar_xlo), bottoms[drawable],
                 pmin(bin_end[drawable], bar_xhi), count[drawable],
                 col = cfg$fill, border = cfg$border, lwd = cfg$bar_line_width)
  graphics::lines(curve$x[valid_curve], curve_y[valid_curve], col = cfg$fit_color, lwd = cfg$fit_line_width)

  if (!is.na(cfg$x_tick_step)) graphics::axis(1, at = ticks(plot_xlo, bar_xhi, cfg$x_tick_step))
  if (!is.na(cfg$y_tick_step) && !log_y) graphics::axis(2, at = ticks(ylo, yhi, cfg$y_tick_step), las = 1)

  if (cfg$show_equation) {
    label <- sprintf("%s\n%s = %s, p = %s",
                     curve$equation[1],
                     "R²",
                     fmt(as.numeric(curve[["R^2"]][1])),
                     format.pval(as.numeric(curve$p_value[1]), digits = 3, eps = 0.001))
    graphics::legend("topright", legend = label, col = cfg$fit_color, lwd = cfg$fit_line_width,
                     bty = "n", cex = 0.86, text.col = "black", inset = 0.015)
  }
  graphics::box(bty = "l", lwd = 0.8)
}, finally = grDevices::dev.off())

message("图片已输出：", normalizePath(cfg$output, winslash = "/", mustWork = FALSE))
