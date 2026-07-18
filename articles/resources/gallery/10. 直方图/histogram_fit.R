#!/usr/bin/env Rscript

# 带拟合曲线的发表级直方图。
# 输入为无表头的单列数值文件。运行前只需修改下面的配置区，然后在 RStudio
# 中点击 Source，或直接运行整个脚本；本脚本不接收命令行参数。
config <- list(
  # 输入文件路径。文件必须是无表头的单列数值 CSV、TSV 或 TXT；相对路径以
  # 当前 R 工作目录为基准。例如："data/example.txt"。
  input = "input_cis.tsv",

  # 图片输出路径。扩展名决定格式；支持 .svg、.pdf、.png、.tif 和 .tiff。
  # 投稿优先推荐 PDF/SVG（矢量图）或 600 dpi TIFF（位图）。
  output = "cis50分箱单对数.svg",

  # 表格输出路径。设为 NULL 时，自动使用图片文件名并添加 "_bins.csv"；
  # 例如 histogram_fit.svg 对应 histogram_fit_bins.csv。
  output_csv = NULL,

  # 拟合模型："exponential" 表示对所选拟合数值建立指数模型；
  # "linear" 表示对所选拟合数值建立线性模型。选择 linear 后，图中和 CSV
  # 始终显示线性方程；其中 y 是各直方图区间的频数。
  fit_model = "exponential",

  # 拟合前如何变换直方图的“箱中心 x”和“频数 y”。该设置只影响拟合，
  # 不改变图片坐标轴："original" 使用原始 x、y；"log_y" 使用 x、ln(y)；
  # "log_log" 使用 ln(x)、ln(y)。对数模式会排除 y=0 的箱，双对数还要求
  # x>0。若要复现 Excel 的“指数趋势线”，请使用 original + exponential。
  fit_value_mode = "original",

  # 坐标系："normal" 为普通坐标；"log_y" 为仅 Y 轴取对数（单对数）；
  # "log_log" 为 X、Y 轴均取对数（双对数）。双对数要求显示范围内 X>0。
  coord_system = "normal",

  # 图片宽度，单位为英寸。常见单栏宽约 3.5 英寸，双栏宽约 7.0 英寸。
  width = 7.2,

  # 图片横纵比，即“宽度 ÷ 高度”。例如 4/3、3/2、16/9；脚本会自动按
  # height = width / aspect_ratio 计算高度，无需另设 height。
  aspect_ratio = 3 / 2,

  # PNG/TIFF 的分辨率（dpi）；SVG/PDF 是矢量图，不受该参数影响。
  dpi = 600,

  # 直方图每个箱的宽度，必须大于 0。数值越小，分箱越细。
  bin_width = 50,

  # X 轴显示范围。NA_real_ 表示按数据自动确定；手动指定时 x_min<x_max。
  x_min = NA_real_,
  x_max = NA_real_,

  # 独立控制参与拟合的箱中心范围。NA_real_ 表示不限制。该范围不会改变
  # 直方图显示范围；若要与 Excel 核对，请确保两边选择了完全相同的中点行。
  fit_x_min = NA_real_,
  fit_x_max = NA_real_,

  # 最左侧柱子与 Y 轴之间的留白，以“箱宽”为单位。1 表示留出一个完整
  # bin_width，0 表示柱子紧贴 Y 轴；可使用 0.5 等小数。
  left_padding_bins = 1,

  # Y 轴显示范围。NA_real_ 表示自动确定。普通坐标的自动下限为 0；
  # 对数 Y 轴的下限必须大于 0，自动值为最小正频数的一半。
  y_min = NA_real_,
  y_max = NA_real_,

  # 主刻度间隔。NA_real_ 表示由 R 自动选择；必须大于 0。
  # 对数 Y 轴始终使用自动的对数刻度，因此 y_tick_step 在该情况下不生效。
  x_tick_step = NA_real_,
  y_tick_step = NA_real_,

  # 图标题及坐标轴标题。main="" 表示不显示主标题，通常更适合论文排版。
  main = "",
  xlab = "Value",
  ylab = "Frequency",

  # 颜色，可使用十六进制色值或 R 颜色名称。fill 为柱体填充色，border 为
  # 柱体边框色，fit_color 为拟合曲线颜色。
  fill = "#006AA7",
  border = "white",
  fit_color = "#FECC02",

  # 柱体边框和拟合曲线的线宽，单位为 R 图形设备的相对线宽。
  bar_line_width = 0.6,
  fit_line_width = 1.5,

  # 字体族。"sans" 通常映射为 Arial/Helvetica；也可填写系统已安装字体名。
  font_family = "sans",

  # 基础字号，单位为 point。坐标标题和刻度字号会在此基础上按比例调整。
  base_font_size = 10,

  # 是否在图内显示拟合方程、R² 和模型整体显著性 p 值。
  show_equation = TRUE,

  # 图中显示哪一种 R²："excel" 表示按 Excel 趋势线口径显示，即原始频数
  # count 与拟合预测频数 predicted_count 的相关系数平方；"model" 表示显示
  # lm() 在实际回归变量上的 R²。指数趋势线下二者可能明显不同。
  r_squared_to_show = "excel",

  # TRUE：忽略空值、非数值和无穷值；FALSE：遇到这些值时停止并报错。
  na_rm = TRUE
)
fmt <- function(x, digits = 4) formatC(x, digits = digits, format = "g")
ticks <- function(lo, hi, step) {
  if (is.na(step)) return(NULL)
  if (step <= 0) stop("刻度步长必须大于 0。", call. = FALSE)
  seq(ceiling(lo / step) * step, hi, by = step)
}
get_script_dir <- function() {
  ofile <- tryCatch(sys.frame(1)$ofile, error = function(e) NULL)
  if (!is.null(ofile) && length(ofile) == 1L && nzchar(ofile)) {
    return(dirname(normalizePath(ofile, winslash = "/", mustWork = FALSE)))
  }
  getwd()
}
is_absolute_path <- function(path) {
  grepl("^(~|[A-Za-z]:|/|\\\\\\\\)", path)
}
script_relative <- function(path, base_dir) {
  if (is.null(path) || !nzchar(path) || is_absolute_path(path)) return(path)
  file.path(base_dir, path)
}

cfg <- config
cfg$fit_model <- tolower(cfg$fit_model)
cfg$fit_value_mode <- tolower(cfg$fit_value_mode)
cfg$coord_system <- tolower(cfg$coord_system)
cfg$r_squared_to_show <- tolower(cfg$r_squared_to_show)
if (!cfg$fit_model %in% c("exponential", "linear")) stop("fit_model 只能为 exponential 或 linear。", call. = FALSE)
if (!cfg$fit_value_mode %in% c("original", "log_y", "log_log")) stop("fit_value_mode 只能为 original、log_y 或 log_log。", call. = FALSE)
if (!cfg$coord_system %in% c("normal", "log_y", "log_log")) stop("coord_system 只能为 normal、log_y 或 log_log。", call. = FALSE)
if (!cfg$r_squared_to_show %in% c("excel", "model")) stop("r_squared_to_show 只能为 excel 或 model。", call. = FALSE)
if (is.null(cfg$input) || !nzchar(cfg$input)) stop("请在 config$input 中设置输入文件路径。", call. = FALSE)
script_dir <- get_script_dir()
cfg$input <- script_relative(cfg$input, script_dir)
cfg$output <- script_relative(cfg$output, script_dir)
if (!is.null(cfg$output_csv) && nzchar(cfg$output_csv)) cfg$output_csv <- script_relative(cfg$output_csv, script_dir)
if (!file.exists(cfg$input)) stop("输入文件不存在：", cfg$input, call. = FALSE)
if (length(cfg$aspect_ratio) != 1L || !is.finite(cfg$aspect_ratio) || cfg$aspect_ratio <= 0) stop("aspect_ratio 必须是大于 0 的有限数值。", call. = FALSE)
if (cfg$width <= 0 || cfg$dpi <= 0 || cfg$bin_width <= 0) stop("width、dpi 和 bin_width 必须大于 0。", call. = FALSE)
if (length(cfg$left_padding_bins) != 1L || !is.finite(cfg$left_padding_bins) || cfg$left_padding_bins < 0) stop("left_padding_bins 必须是大于或等于 0 的有限数值。", call. = FALSE)
if (!is.na(cfg$fit_x_min) && !is.na(cfg$fit_x_max) && cfg$fit_x_min >= cfg$fit_x_max) stop("fit_x_min 必须小于 fit_x_max。", call. = FALSE)
cfg$height <- cfg$width / cfg$aspect_ratio

raw <- scan(cfg$input, what = character(), sep = "", quiet = TRUE)
values <- suppressWarnings(as.numeric(raw)); bad <- !is.finite(values)
if (any(bad) && !cfg$na_rm) stop(sprintf("输入中有 %d 个无效值。", sum(bad)), call. = FALSE)
values <- values[!bad]
if (!length(values)) stop("输入中没有有效数值。", call. = FALSE)

xlo <- if (is.na(cfg$x_min)) min(values) else cfg$x_min
xhi <- if (is.na(cfg$x_max)) max(values) else cfg$x_max
if (xlo >= xhi) { xlo <- xlo - cfg$bin_width / 2; xhi <- xhi + cfg$bin_width / 2 }
if (cfg$coord_system == "log_log" && xlo <= 0) stop("双对数坐标要求 x_min 和所有显示的 X 值大于 0。", call. = FALSE)
br <- seq(floor(xlo / cfg$bin_width) * cfg$bin_width,
          ceiling(xhi / cfg$bin_width) * cfg$bin_width, by = cfg$bin_width)
if (length(br) < 2L) br <- c(xlo, xlo + cfg$bin_width)
h <- hist(values[values >= min(br) & values <= max(br)], breaks = br,
          plot = FALSE, right = FALSE, include.lowest = TRUE)
x <- h$mids; y <- h$counts

# 拟合数据变换与反变换。拟合始终针对直方图箱中心和频数，而不是原始单列样本。
base_keep <- is.finite(x) & is.finite(y)
if (!is.na(cfg$fit_x_min)) base_keep <- base_keep & x >= cfg$fit_x_min
if (!is.na(cfg$fit_x_max)) base_keep <- base_keep & x <= cfg$fit_x_max
if (cfg$fit_value_mode %in% c("log_y", "log_log")) base_keep <- base_keep & y > 0
if (cfg$fit_value_mode == "log_log") base_keep <- base_keep & x > 0
fit_x_all <- if (cfg$fit_value_mode == "log_log") suppressWarnings(log(x)) else x
fit_y_all <- if (cfg$fit_value_mode %in% c("log_y", "log_log")) suppressWarnings(log(y)) else y
inverse_y <- if (cfg$fit_value_mode %in% c("log_y", "log_log")) exp else identity
transform_new_x <- if (cfg$fit_value_mode == "log_log") function(z) suppressWarnings(log(z)) else identity

if (cfg$fit_model == "exponential") {
  # 与 Excel 指数趋势线相同，指数回归通过 ln(拟合Y)=a+b*拟合X 求解；
  # 因此拟合后的 Y 必须大于 0。original 模式会自然排除零频数箱。
  keep <- base_keep & is.finite(fit_x_all) & is.finite(fit_y_all) & fit_y_all > 0
  if (sum(keep) < 3L) stop("当前变换下，指数拟合至少需要 3 个变换后 Y>0 的箱。", call. = FALSE)
  model <- lm(log(fit_y_all[keep]) ~ fit_x_all[keep]); co <- coef(model)
  transformed_predict <- function(z) exp(co[1L] + co[2L] * z)
} else {
  keep <- base_keep & is.finite(fit_x_all) & is.finite(fit_y_all)
  if (sum(keep) < 3L || length(unique(fit_x_all[keep])) < 2L) stop("线性拟合至少需要 3 个有效箱。", call. = FALSE)
  model <- lm(fit_y_all[keep] ~ fit_x_all[keep]); co <- coef(model)
  transformed_predict <- function(z) co[1L] + co[2L] * z
}
predict_fit <- function(z) inverse_y(transformed_predict(transform_new_x(z)))

# 方程按实际回归变量显示：选择 linear 时始终打印线性形式；选择
# exponential 时始终打印指数形式。绘图预测值仍会反变换到原始频数尺度。
sgn <- function(v) if (v < 0) "-" else "+"
if (cfg$fit_model == "linear" && cfg$fit_value_mode == "original")
  equation <- sprintf("y = %s %s %s x", fmt(co[1L], 7), sgn(co[2L]), fmt(abs(co[2L]), 7)) else
if (cfg$fit_model == "linear" && cfg$fit_value_mode == "log_y")
  equation <- sprintf("ln(y) = %s %s %s x", fmt(co[1L], 7), sgn(co[2L]), fmt(abs(co[2L]), 7)) else
if (cfg$fit_model == "linear" && cfg$fit_value_mode == "log_log")
  equation <- sprintf("ln(y) = %s %s %s ln(x)", fmt(co[1L], 7), sgn(co[2L]), fmt(abs(co[2L]), 7)) else
if (cfg$fit_model == "exponential" && cfg$fit_value_mode == "original")
  equation <- sprintf("y = %s exp(%s x)", fmt(exp(co[1L]), 7), fmt(co[2L], 7)) else
if (cfg$fit_model == "exponential" && cfg$fit_value_mode == "log_y")
  equation <- sprintf("ln(y) = %s exp(%s x)", fmt(exp(co[1L]), 7), fmt(co[2L], 7)) else
  equation <- sprintf("ln(y) = %s exp(%s ln(x))", fmt(exp(co[1L]), 7), fmt(co[2L], 7))

sm <- summary(model); model_r2 <- unname(sm$r.squared)
f <- sm$fstatistic
p_value <- if (is.null(f)) NA_real_ else unname(pf(f[1L], f[2L], f[3L], lower.tail = FALSE))

# 输出逐点诊断量。regression_y 是 R² 实际使用的因变量，可直接复制到 Excel
# 与 regression_x 重新回归；未参与拟合的箱保留为 NA。
regression_y <- if (cfg$fit_model == "exponential") suppressWarnings(log(fit_y_all)) else fit_y_all
fitted_regression_y <- residual_regression_y <- rep(NA_real_, length(y))
fitted_regression_y[keep] <- unname(fitted(model))
residual_regression_y[keep] <- unname(residuals(model))
predicted_count <- suppressWarnings(predict_fit(x))
predicted_count[!is.finite(predicted_count)] <- NA_real_
excel_r2 <- suppressWarnings(stats::cor(y[keep], predicted_count[keep], use = "complete.obs")^2)
r2 <- if (cfg$r_squared_to_show == "excel") excel_r2 else model_r2
r2_label <- if (cfg$r_squared_to_show == "excel") "R²(excel)" else "R²(model)"

bin_table <- data.frame(
  bin_start = h$breaks[-length(h$breaks)], bin_end = h$breaks[-1L],
  bin_midpoint = x, count = y, used_for_fit = keep,
  transformed_x = fit_x_all, transformed_y = fit_y_all,
  regression_x = fit_x_all, regression_y = regression_y,
  fitted_regression_y = fitted_regression_y,
  residual_regression_y = residual_regression_y,
  predicted_count = predicted_count,
  fit_model = cfg$fit_model, fit_value_mode = cfg$fit_value_mode,
  coordinate_system = cfg$coord_system,
  r_squared_scale = paste0(if (cfg$fit_model == "exponential") "log(" else "",
                           if (cfg$fit_value_mode == "original") "count" else "log(count)",
                           if (cfg$fit_model == "exponential") ")" else ""),
  # 保留完整精度的回归系数，便于与 Excel LINEST/趋势线结果逐项核对。
  regression_intercept = unname(co[1L]), regression_slope = unname(co[2L]),
  fit_equation = equation,
  displayed_r_squared = r2,
  excel_r_squared = excel_r2,
  model_r_squared = model_r2,
  model_p_value = p_value,
  n_bins_used_for_fit = sum(keep), stringsAsFactors = FALSE, check.names = FALSE
)
if (is.null(cfg$output_csv) || !nzchar(cfg$output_csv)) cfg$output_csv <- paste0(tools::file_path_sans_ext(cfg$output), "_bins.csv")
dir.create(dirname(cfg$output_csv), recursive = TRUE, showWarnings = FALSE)
write.csv(bin_table, cfg$output_csv, row.names = FALSE, fileEncoding = "UTF-8")

# 自动范围显示完整分箱；左侧额外保留 left_padding_bins 个箱宽。
bar_xlo <- if (is.na(cfg$x_min)) min(h$breaks) else xlo
bar_xhi <- if (is.na(cfg$x_max)) max(h$breaks) else xhi
plot_xlo <- bar_xlo - cfg$left_padding_bins * cfg$bin_width
if (cfg$coord_system == "log_log" && plot_xlo <= 0) {
  stop("双对数坐标下，左侧留白后 X 轴下限必须大于 0；请减小 left_padding_bins 或增大 x_min。", call. = FALSE)
}

# 对数 Y 轴不能显示 0；绘图下限取最小正频数的一半。
log_y <- cfg$coord_system %in% c("log_y", "log_log")
positive_y <- y[y > 0]
ylo <- if (!is.na(cfg$y_min)) cfg$y_min else if (log_y) max(min(positive_y) / 2, .Machine$double.eps) else 0
grid_x <- seq(bar_xlo, bar_xhi, length.out = 500L); fit_y <- predict_fit(grid_x)
valid_curve <- is.finite(fit_y) & (!log_y | fit_y > 0) &
  (cfg$fit_value_mode != "log_log" | grid_x > 0)
auto_yhi <- max(c(positive_y, fit_y[valid_curve]), na.rm = TRUE) * if (log_y) 1.35 else 1.08
yhi <- if (is.na(cfg$y_max)) auto_yhi else cfg$y_max
if (log_y && ylo <= 0) stop("对数 Y 轴要求 y_min 大于 0。", call. = FALSE)
if (ylo >= yhi) stop("y_min 必须小于 y_max。", call. = FALSE)

dir.create(dirname(cfg$output), recursive = TRUE, showWarnings = FALSE)
ext <- tolower(tools::file_ext(cfg$output))
if (ext == "svg") svg(cfg$output, cfg$width, cfg$height, pointsize = cfg$base_font_size, family = cfg$font_family, bg = "white") else
if (ext == "pdf") pdf(cfg$output, cfg$width, cfg$height, pointsize = cfg$base_font_size, family = cfg$font_family, useDingbats = FALSE) else
if (ext == "png") png(cfg$output, cfg$width, cfg$height, units = "in", res = cfg$dpi, pointsize = cfg$base_font_size, type = "cairo") else
if (ext %in% c("tif", "tiff")) tiff(cfg$output, cfg$width, cfg$height, units = "in", res = cfg$dpi, compression = "lzw", pointsize = cfg$base_font_size) else
  stop("输出图片仅支持 svg、pdf、png、tif/tiff。", call. = FALSE)

tryCatch({
  par(family = cfg$font_family, mar = c(4.3, 4.6, if (nzchar(cfg$main)) 2.2 else 0.8, 0.8),
      mgp = c(2.65, 0.75, 0), tcl = -0.25, las = 1, lend = "round",
      cex.axis = 0.92, cex.lab = 1.05, cex.main = 1.08, font.lab = 2, bty = "l")
  plot(NA, xlim = c(plot_xlo, bar_xhi), ylim = c(ylo, yhi), log = if (cfg$coord_system == "log_y") "y" else if (cfg$coord_system == "log_log") "xy" else "",
       xaxs = "i", yaxs = "i", xaxt = if (is.na(cfg$x_tick_step)) "s" else "n",
       yaxt = if (is.na(cfg$y_tick_step) || log_y) "s" else "n",
       xlab = cfg$xlab, ylab = cfg$ylab, main = cfg$main)
  bottoms <- if (log_y) rep(ylo, length(y)) else rep(0, length(y))
  drawable <- y > 0 & h$breaks[-1L] > bar_xlo & h$breaks[-length(h$breaks)] < bar_xhi
  rect(pmax(h$breaks[-length(h$breaks)][drawable], bar_xlo), bottoms[drawable],
       pmin(h$breaks[-1L][drawable], bar_xhi), y[drawable],
       col = cfg$fill, border = cfg$border, lwd = cfg$bar_line_width)
  lines(grid_x[valid_curve], fit_y[valid_curve], col = cfg$fit_color, lwd = cfg$fit_line_width)
  if (!is.na(cfg$x_tick_step)) axis(1, at = ticks(plot_xlo, bar_xhi, cfg$x_tick_step))
  if (!is.na(cfg$y_tick_step) && !log_y) axis(2, at = ticks(ylo, yhi, cfg$y_tick_step), las = 1)
  if (cfg$show_equation) {
    label <- sprintf("%s\n%s = %s, p = %s", equation, r2_label, fmt(r2), format.pval(p_value, digits = 3, eps = 0.001))
    legend("topright", legend = label, col = cfg$fit_color, lwd = cfg$fit_line_width,
           bty = "n", cex = 0.86, text.col = "black", inset = 0.015)
  }
  box(bty = "l", lwd = 0.8)
}, finally = dev.off())

message("图片已输出：", normalizePath(cfg$output, winslash = "/", mustWork = FALSE))
message("CSV 已输出：", normalizePath(cfg$output_csv, winslash = "/", mustWork = FALSE))
message(sprintf("拟合使用 %d 个箱子的中点（共 %d 个箱）；中点范围：%s 至 %s。",
                sum(keep), length(x), fmt(min(x[keep]), 7), fmt(max(x[keep]), 7)))
message(sprintf("拟合：%s；%s = %s；R²(model) = %s；p = %s",
                equation, r2_label, fmt(r2), fmt(model_r2),
                format.pval(p_value, digits = 3, eps = 0.001)))
