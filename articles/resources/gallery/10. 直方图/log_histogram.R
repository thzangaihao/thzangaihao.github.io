 #!/usr/bin/env Rscript

# 单列正数数据的对数分箱直方图（输入文件无表头）
# 分箱在对数尺度上等宽，原始尺度上的箱宽按等比数列增长。
config <- list(
  input           = NULL,                    # 输入 CSV/TSV/TXT 路径（必填）
  output          = "log_histogram.svg",     # 输出 SVG 路径
  output_csv      = NULL,                    # NULL 时自动生成 *_bins.csv
  width           = 8,
  height          = 5,
  log_base        = 10,                      # 对数底数，必须大于 1
  bins_per_decade = 10,                      # 每个对数单位（数量级）的箱数
  x_min           = NA_real_,                # 必须大于 0；NA 表示自动确定
  x_max           = NA_real_,                # 必须大于 0；NA 表示自动确定
  y_min           = 0,
  y_max           = NA_real_,
  y_tick_step     = NA_real_,
  main            = "Log-binned Histogram",
  xlab            = "Value (log scale)",
  ylab            = "Frequency",
  fill            = "#4C78A8",
  border          = "white",
  line_width      = 1,
  font_size       = 1,
  font_family     = "sans",
  na_rm           = TRUE                     # 是否忽略空值、非数值、无穷值及非正数
)

usage <- function() {
  cat(paste0(
    "用法:\n",
    "  Rscript log_histogram.R --input=数据.tsv [--output=直方图.svg] ",
    "[--output-csv=分箱统计.csv] [选项]\n\n",
    "常用选项:\n",
    "  --log-base=10 --bins-per-decade=10\n",
    "  --x-min=0.01 --x-max=1000 --y-min=0 --y-max=50\n",
    "  --y-tick-step=5 --main=标题 --xlab=数值 --ylab=频数\n",
    "  --fill=#4C78A8 --border=white --line-width=1\n",
    "  --width=8 --height=5 --font-size=1 --font-family=sans --na-rm=true\n"
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
  if (length(value) != 1L || is.na(value) || !is.finite(value)) {
    stop(sprintf("参数 --%s 必须是有限数值。", name), call. = FALSE)
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
    width = "width", height = "height", "log-base" = "log_base",
    "bins-per-decade" = "bins_per_decade", "x-min" = "x_min",
    "x-max" = "x_max", "y-min" = "y_min", "y-max" = "y_max",
    "y-tick-step" = "y_tick_step", main = "main", xlab = "xlab",
    ylab = "ylab", fill = "fill", border = "border",
    "line-width" = "line_width", "font-size" = "font_size",
    "font-family" = "font_family", "na-rm" = "na_rm"
  )
  numeric_fields <- c(
    "width", "height", "log_base", "bins_per_decade", "x_min", "x_max",
    "y_min", "y_max", "y_tick_step", "line_width", "font_size"
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

make_y_ticks <- function(lower, upper, step) {
  if (is.na(step)) return(NULL)
  if (step <= 0) stop("Y 轴刻度步长必须大于 0。", call. = FALSE)
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
if (cfg$log_base <= 1) stop("log_base 必须大于 1。", call. = FALSE)
if (cfg$bins_per_decade <= 0) stop("bins_per_decade 必须大于 0。", call. = FALSE)
if (cfg$line_width < 0 || cfg$font_size <= 0) stop("线宽不能为负，字号必须大于 0。", call. = FALSE)
if ((!is.na(cfg$x_min) && cfg$x_min <= 0) || (!is.na(cfg$x_max) && cfg$x_max <= 0)) {
  stop("对数坐标下 x_min 和 x_max 必须大于 0。", call. = FALSE)
}
if (!is.na(cfg$x_min) && !is.na(cfg$x_max) && cfg$x_min >= cfg$x_max) {
  stop("x_min 必须小于 x_max。", call. = FALSE)
}

raw <- scan(cfg$input, what = character(), sep = "", quiet = TRUE)
values <- suppressWarnings(as.numeric(raw))
bad <- is.na(values) | !is.finite(values) | values <= 0
if (any(bad) && !cfg$na_rm) {
  stop(sprintf("输入中有 %d 个空值、非数值、无穷值或非正数。", sum(bad)), call. = FALSE)
}
if (any(bad) && cfg$na_rm) message(sprintf("已忽略 %d 个无法进行对数分箱的值。", sum(bad)))
values <- values[!bad]
if (!length(values)) stop("输入文件中没有可用于对数分箱的正有限数值。", call. = FALSE)

x_lower <- if (is.na(cfg$x_min)) min(values) else cfg$x_min
x_upper <- if (is.na(cfg$x_max)) max(values) else cfg$x_max
if (x_lower == x_upper) {
  factor <- cfg$log_base ^ (0.5 / cfg$bins_per_decade)
  x_lower <- x_lower / factor
  x_upper <- x_upper * factor
}

# 在 log_base 对数尺度上对齐边界并等距切分。
log_lower <- log(x_lower, base = cfg$log_base)
log_upper <- log(x_upper, base = cfg$log_base)
log_step <- 1 / cfg$bins_per_decade
break_start <- floor(log_lower / log_step) * log_step
break_end <- ceiling(log_upper / log_step) * log_step
if (break_end <= break_start) break_end <- break_start + log_step
log_breaks <- seq(break_start, break_end, by = log_step)
breaks <- cfg$log_base ^ log_breaks

plot_values <- values[values >= x_lower & values <= x_upper]
if (!length(plot_values)) stop("指定的 X 轴范围内没有数据。", call. = FALSE)
hist_info <- hist(
  plot_values, breaks = breaks, plot = FALSE, right = FALSE,
  include.lowest = TRUE
)

if (is.null(cfg$output_csv) || !nzchar(cfg$output_csv)) {
  cfg$output_csv <- paste0(tools::file_path_sans_ext(cfg$output), "_bins.csv")
}
csv_dir <- dirname(cfg$output_csv)
if (!dir.exists(csv_dir)) dir.create(csv_dir, recursive = TRUE, showWarnings = FALSE)

bin_count <- length(hist_info$counts)
left <- hist_info$breaks[seq_len(bin_count)]
right <- hist_info$breaks[seq_len(bin_count) + 1L]
bin_table <- data.frame(
  bin_start = left,
  bin_end = right,
  interval = sprintf(
    ifelse(seq_len(bin_count) == bin_count, "[%s, %s]", "[%s, %s)"),
    format(left, trim = TRUE, scientific = TRUE, digits = 10),
    format(right, trim = TRUE, scientific = TRUE, digits = 10)
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
    xlim = c(x_lower, x_upper), ylim = c(cfg$y_min, y_upper), log = "x",
    yaxt = if (is.na(cfg$y_tick_step)) "s" else "n",
    include.lowest = TRUE, right = FALSE
  )
  y_ticks <- make_y_ticks(cfg$y_min, y_upper, cfg$y_tick_step)
  if (!is.null(y_ticks)) axis(2, at = y_ticks, las = 1)
}, finally = {
  grDevices::dev.off()
})

message("SVG 已输出：", normalizePath(cfg$output, winslash = "/", mustWork = FALSE))
message("分箱统计 CSV 已输出：", normalizePath(cfg$output_csv, winslash = "/", mustWork = FALSE))
