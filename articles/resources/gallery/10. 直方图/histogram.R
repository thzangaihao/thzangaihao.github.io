#!/usr/bin/env Rscript

# 单列数值直方图（输入文件无表头）
# 日常使用时可直接修改本配置区，也可用同名命令行参数临时覆盖。
config <- list(
  input       = NULL,             # 输入 CSV/TSV/TXT 路径（必填）
  output      = "histogram.svg",  # 输出 SVG 路径
  output_csv  = NULL,             # 分箱统计 CSV；NULL 时按 SVG 文件名自动生成
  width       = 8,                # 图片宽度（英寸）
  height      = 5,                # 图片高度（英寸）
  bin_width   = 1,                # 直方图分箱步长
  x_min       = NA_real_,         # X 轴最小值；NA 表示根据数据自动确定
  x_max       = NA_real_,         # X 轴最大值；NA 表示根据数据自动确定
  y_min       = 0,                # Y 轴最小值
  y_max       = NA_real_,         # Y 轴最大值；NA 表示自动确定
  x_tick_step = NA_real_,         # X 轴刻度步长；NA 表示自动确定
  y_tick_step = NA_real_,         # Y 轴刻度步长；NA 表示自动确定
  main        = "Histogram",      # 主标题；设为 "" 可隐藏
  xlab        = "Value",          # X 轴标题
  ylab        = "Frequency",      # Y 轴标题
  fill        = "#4C78A8",        # 柱体填充色
  border      = "white",          # 柱体边框色
  line_width  = 1,                # 柱体边框粗细
  font_size   = 1,                # 全局文字缩放倍数
  font_family = "sans",           # 字体族
  dpi         = 96,               # SVG 名义分辨率
  na_rm       = TRUE               # 是否忽略空值/非数值；FALSE 时直接报错
)

usage <- function() {
  cat(paste0(
    "用法:\n",
    "  Rscript histogram.R --input=数据.tsv [--output=直方图.svg] [--output-csv=分箱统计.csv] [选项]\n\n",
    "常用选项:\n",
    "  --width=8 --height=5 --bin-width=1\n",
    "  --x-min=0 --x-max=100 --y-min=0 --y-max=50\n",
    "  --x-tick-step=10 --y-tick-step=5\n",
    "  --main=标题 --xlab=数值 --ylab=频数\n",
    "  --fill=#4C78A8 --border=white --line-width=1\n",
    "  --font-size=1 --font-family=sans --dpi=96 --na-rm=true\n"
  ))
}

parse_bool <- function(x, name) {
  value <- tolower(x)
  if (value %in% c("true", "t", "1", "yes", "y")) return(TRUE)
  if (value %in% c("false", "f", "0", "no", "n")) return(FALSE)
  stop(sprintf("参数 --%s 必须是 true 或 false。", name), call. = FALSE)
}

parse_number <- function(x, name) {
  value <- suppressWarnings(as.numeric(x))
  if (length(value) != 1L || is.na(value)) {
    stop(sprintf("参数 --%s 必须是数值。", name), call. = FALSE)
  }
  value
}

parse_args <- function(args, defaults) {
  if (any(args %in% c("-h", "--help"))) {
    usage()
    quit(status = 0)
  }

  aliases <- c(
    input = "input", output = "output", "output-csv" = "output_csv",
    width = "width", height = "height",
    "bin-width" = "bin_width", "x-min" = "x_min", "x-max" = "x_max",
    "y-min" = "y_min", "y-max" = "y_max",
    "x-tick-step" = "x_tick_step", "y-tick-step" = "y_tick_step",
    main = "main", xlab = "xlab", ylab = "ylab", fill = "fill",
    border = "border", "line-width" = "line_width",
    "font-size" = "font_size", "font-family" = "font_family",
    dpi = "dpi", "na-rm" = "na_rm"
  )
  numeric_fields <- c(
    "width", "height", "bin_width", "x_min", "x_max", "y_min", "y_max",
    "x_tick_step", "y_tick_step", "line_width", "font_size", "dpi"
  )

  for (arg in args) {
    matched <- regexec("^--([^=]+)=(.*)$", arg)
    parts <- regmatches(arg, matched)[[1L]]
    if (length(parts) != 3L) {
      stop(sprintf("无法识别参数 '%s'；参数应写成 --名称=值。", arg), call. = FALSE)
    }
    cli_name <- parts[2L]
    if (!cli_name %in% names(aliases)) {
      stop(sprintf("未知参数 --%s。使用 --help 查看帮助。", cli_name), call. = FALSE)
    }
    field <- unname(aliases[[cli_name]])
    raw_value <- parts[3L]
    if (field %in% numeric_fields) {
      defaults[[field]] <- parse_number(raw_value, cli_name)
    } else if (field == "na_rm") {
      defaults[[field]] <- parse_bool(raw_value, cli_name)
    } else {
      defaults[[field]] <- raw_value
    }
  }
  defaults
}

make_ticks <- function(lower, upper, step, axis_name) {
  if (is.na(step)) return(NULL)
  if (step <= 0) stop(sprintf("%s 轴刻度步长必须大于 0。", axis_name), call. = FALSE)
  first <- ceiling(lower / step) * step
  if (first > upper) return(numeric(0))
  seq(first, upper, by = step)
}

cfg <- parse_args(commandArgs(trailingOnly = TRUE), config)

if (is.null(cfg$input) || !nzchar(cfg$input)) {
  usage()
  stop("缺少输入文件，请设置 config$input 或使用 --input=文件路径。", call. = FALSE)
}
if (!file.exists(cfg$input)) stop("输入文件不存在：", cfg$input, call. = FALSE)
if (cfg$width <= 0 || cfg$height <= 0) stop("图片宽度和高度必须大于 0。", call. = FALSE)
if (cfg$bin_width <= 0) stop("分箱步长 bin_width 必须大于 0。", call. = FALSE)
if (cfg$dpi <= 0) stop("dpi 必须大于 0。", call. = FALSE)

# sep="" 会将连续空白、逗号或制表符均视作分隔符，适配单列 CSV/TSV/TXT。
raw <- scan(cfg$input, what = character(), sep = "", quiet = TRUE)
values <- suppressWarnings(as.numeric(raw))
bad <- is.na(values) | !is.finite(values)
if (any(bad) && !cfg$na_rm) {
  stop(sprintf("输入中有 %d 个空值、非数值或无穷值。", sum(bad)), call. = FALSE)
}
values <- values[!bad]
if (!length(values)) stop("输入文件中没有可绘制的有限数值。", call. = FALSE)

x_lower <- if (is.na(cfg$x_min)) min(values) else cfg$x_min
x_upper <- if (is.na(cfg$x_max)) max(values) else cfg$x_max
if (x_lower >= x_upper) {
  x_lower <- x_lower - cfg$bin_width / 2
  x_upper <- x_upper + cfg$bin_width / 2
}
if (!is.na(cfg$x_min) && !is.na(cfg$x_max) && cfg$x_min >= cfg$x_max) {
  stop("x_min 必须小于 x_max。", call. = FALSE)
}

# 分箱边界对齐到 bin_width 的整数倍，并完整覆盖绘图范围。
break_start <- floor(x_lower / cfg$bin_width) * cfg$bin_width
break_end <- ceiling(x_upper / cfg$bin_width) * cfg$bin_width
if (break_end <= break_start) break_end <- break_start + cfg$bin_width
breaks <- seq(break_start, break_end, by = cfg$bin_width)
plot_values <- values[values >= break_start & values <= break_end]
if (!length(plot_values)) stop("指定的 X 轴范围内没有数据。", call. = FALSE)
hist_info <- hist(
  plot_values, breaks = breaks, plot = FALSE, right = FALSE,
  include.lowest = TRUE
)

# 输出与直方图完全一致的分箱统计。除最后一个区间外，区间均为 [左边界, 右边界)。
if (is.null(cfg$output_csv) || !nzchar(cfg$output_csv)) {
  cfg$output_csv <- paste0(tools::file_path_sans_ext(cfg$output), "_bins.csv")
}
csv_dir <- dirname(cfg$output_csv)
if (!dir.exists(csv_dir)) dir.create(csv_dir, recursive = TRUE, showWarnings = FALSE)

bin_count <- length(hist_info$counts)
bin_table <- data.frame(
  bin_start = hist_info$breaks[seq_len(bin_count)],
  bin_end = hist_info$breaks[seq_len(bin_count) + 1L],
  interval = sprintf(
    ifelse(seq_len(bin_count) == bin_count, "[%s, %s]", "[%s, %s)"),
    format(hist_info$breaks[seq_len(bin_count)], trim = TRUE, scientific = FALSE),
    format(hist_info$breaks[seq_len(bin_count) + 1L], trim = TRUE, scientific = FALSE)
  ),
  count = hist_info$counts,
  check.names = FALSE
)
utils::write.csv(bin_table, cfg$output_csv, row.names = FALSE, fileEncoding = "UTF-8")

y_upper <- if (is.na(cfg$y_max)) max(hist_info$counts) else cfg$y_max
if (cfg$y_min >= y_upper) {
  if (is.na(cfg$y_max) && cfg$y_min == 0 && y_upper == 0) y_upper <- 1 else
    stop("y_min 必须小于 y_max（或自动计算出的最大频数）。", call. = FALSE)
}

out_dir <- dirname(cfg$output)
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

grDevices::svg(
  filename = cfg$output, width = cfg$width, height = cfg$height,
  pointsize = 12 * cfg$font_size, family = cfg$font_family,
  bg = "white", antialias = "default"
)
tryCatch({
  par(family = cfg$font_family, cex = cfg$font_size)
  hist(
    plot_values, breaks = breaks, freq = TRUE, col = cfg$fill, border = cfg$border,
    lwd = cfg$line_width, main = cfg$main, xlab = cfg$xlab, ylab = cfg$ylab,
    xlim = c(x_lower, x_upper), ylim = c(cfg$y_min, y_upper),
    xaxt = if (is.na(cfg$x_tick_step)) "s" else "n",
    yaxt = if (is.na(cfg$y_tick_step)) "s" else "n",
    include.lowest = TRUE, right = FALSE
  )
  x_ticks <- make_ticks(x_lower, x_upper, cfg$x_tick_step, "X")
  y_ticks <- make_ticks(cfg$y_min, y_upper, cfg$y_tick_step, "Y")
  if (!is.null(x_ticks)) axis(1, at = x_ticks)
  if (!is.null(y_ticks)) axis(2, at = y_ticks, las = 1)
}, finally = {
  grDevices::dev.off()
})

message("SVG 已输出：", normalizePath(cfg$output, winslash = "/", mustWork = FALSE))
message("分箱统计 CSV 已输出：", normalizePath(cfg$output_csv, winslash = "/", mustWork = FALSE))
