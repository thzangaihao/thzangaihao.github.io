#!/usr/bin/env Rscript

# =============================================================================
# 参数设置区：只需修改本区内容，脚本其余部分无需改动
# =============================================================================
config <- list(
  # 输入一个 BigWig 文件；相对路径以本脚本所在目录为基准。
  input_file = "sample.bw",

  # 输出 SVG 文件；相对路径以本脚本所在目录为基准。
  output_file = "chromosome_coverage.svg",

  # 统计窗口大小（bp）。每个窗口计算 BigWig 信号的加权平均值。
  bin_size = 100000L,

  # 需要绘制的染色体及其排列顺序。
  # 设置为 NULL 时，绘制 BigWig 文件中的全部染色体。
  chromosomes = NULL,
  # chromosomes = c("chr1", "chr2", "chr3"),

  # 多面板布局：每行的染色体面板数。
  panel_columns = 2L,

  # SVG 尺寸（英寸）。总高度 = panel_height * 面板行数。
  figure_width = 12,
  panel_height = 3.2,

  # coverage 线条颜色和宽度。
  line_color = "#0072B2",
  line_width = 1.2,

  # 图形文字与样式。建议使用 Illustrator 能识别的本机字体。
  font_family = "Arial",
  base_font_size = 10,
  title = "Chromosome coverage",
  x_label = "Chromosome position (Mb)",
  y_label = "Mean BigWig signal",
  grid_color = "#E6E6E6",
  background_color = "white"
)
# =============================================================================
# 参数设置区结束
# =============================================================================

get_script_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg)) {
    return(dirname(normalizePath(sub("^--file=", "", file_arg[[1L]]), mustWork = FALSE)))
  }
  getwd()
}

is_absolute_path <- function(path) {
  grepl("^(~|[A-Za-z]:|/|\\\\)", path)
}

resolve_path <- function(path, base_dir) {
  if (is_absolute_path(path)) path else file.path(base_dir, path)
}

mean_signal_in_bins <- function(bigwig, chromosome, chromosome_length, bin_size) {
  starts <- seq.int(1, chromosome_length, by = bin_size)
  ends <- pmin(starts + bin_size - 1, chromosome_length)
  bins <- GenomicRanges::GRanges(
    seqnames = chromosome,
    ranges = IRanges::IRanges(start = starts, end = ends)
  )

  signal <- rtracklayer::import(
    bigwig,
    format = "BigWig",
    which = GenomicRanges::GRanges(
      seqnames = chromosome,
      ranges = IRanges::IRanges(start = 1, end = chromosome_length)
    )
  )

  weighted_sum <- numeric(length(bins))
  if (length(signal)) {
    hits <- GenomicRanges::findOverlaps(bins, signal, ignore.strand = TRUE)
    if (length(hits)) {
      query_id <- S4Vectors::queryHits(hits)
      subject_id <- S4Vectors::subjectHits(hits)
      overlap_width <- pmax(
        0,
        pmin(GenomicRanges::end(bins)[query_id], GenomicRanges::end(signal)[subject_id]) -
          pmax(GenomicRanges::start(bins)[query_id], GenomicRanges::start(signal)[subject_id]) + 1
      )
      score <- as.numeric(S4Vectors::mcols(signal)$score[subject_id])
      score[!is.finite(score)] <- 0
      summed <- rowsum(score * overlap_width, group = query_id, reorder = FALSE)
      weighted_sum[as.integer(rownames(summed))] <- summed[, 1L]
    }
  }

  data.frame(
    midpoint = (starts + ends) / 2,
    coverage = weighted_sum / (ends - starts + 1)
  )
}

cfg <- config
script_dir <- get_script_dir()
cfg$input_file <- resolve_path(cfg$input_file, script_dir)
cfg$output_file <- resolve_path(cfg$output_file, script_dir)

required_packages <- c("rtracklayer", "svglite")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1L), quietly = TRUE)
]
if (length(missing_packages)) {
  stop(
    "缺少 R 包：", paste(missing_packages, collapse = ", "),
    "。可使用 BiocManager::install('rtracklayer') 和 install.packages('svglite') 安装。",
    call. = FALSE
  )
}
if (length(cfg$input_file) != 1L || !nzchar(cfg$input_file)) {
  stop("input_file 只能指定一个 BigWig 文件。", call. = FALSE)
}
if (!file.exists(cfg$input_file)) stop("输入文件不存在：", cfg$input_file, call. = FALSE)
if (tolower(tools::file_ext(cfg$output_file)) != "svg") {
  stop("output_file 必须使用 .svg 扩展名。", call. = FALSE)
}
if (!is.finite(cfg$bin_size) || cfg$bin_size < 1) {
  stop("bin_size 必须是正整数。", call. = FALSE)
}
if (!is.finite(cfg$panel_columns) || cfg$panel_columns < 1) {
  stop("panel_columns 必须是正整数。", call. = FALSE)
}
if (!is.finite(cfg$figure_width) || cfg$figure_width <= 0 ||
    !is.finite(cfg$panel_height) || cfg$panel_height <= 0) {
  stop("figure_width 和 panel_height 必须是正数。", call. = FALSE)
}

bigwig <- rtracklayer::BigWigFile(cfg$input_file)
sequence_lengths <- GenomeInfoDb::seqlengths(GenomeInfoDb::seqinfo(bigwig))
available_chromosomes <- names(sequence_lengths)[is.finite(sequence_lengths) & sequence_lengths > 0]

if (is.null(cfg$chromosomes)) {
  chromosomes <- available_chromosomes
} else {
  absent <- setdiff(cfg$chromosomes, available_chromosomes)
  if (length(absent)) {
    warning("以下染色体不存在或缺少有效长度，已跳过：", paste(absent, collapse = ", "))
  }
  chromosomes <- cfg$chromosomes[cfg$chromosomes %in% available_chromosomes]
}
if (!length(chromosomes)) stop("BigWig 文件中没有可绘制的染色体。", call. = FALSE)

panel_columns <- min(as.integer(cfg$panel_columns), length(chromosomes))
panel_rows <- ceiling(length(chromosomes) / panel_columns)
figure_height <- cfg$panel_height * panel_rows
dir.create(dirname(cfg$output_file), recursive = TRUE, showWarnings = FALSE)

# svglite 会将文字写为 SVG <text> 元素，导入 Adobe Illustrator 后仍可编辑。
svglite::svglite(
  file = cfg$output_file,
  width = cfg$figure_width,
  height = figure_height,
  bg = cfg$background_color,
  pointsize = cfg$base_font_size
)

tryCatch({
  graphics::par(
    mfrow = c(panel_rows, panel_columns),
    mar = c(4.1, 4.5, 2.3, 0.8),
    oma = c(0, 0, if (nzchar(cfg$title)) 2 else 0, 0),
    family = cfg$font_family,
    las = 1,
    bty = "l"
  )

  for (chromosome in chromosomes) {
    chromosome_length <- as.numeric(sequence_lengths[[chromosome]])
    profile <- mean_signal_in_bins(bigwig, chromosome, chromosome_length, cfg$bin_size)
    coverage_max <- suppressWarnings(max(profile$coverage, na.rm = TRUE))
    if (!is.finite(coverage_max) || coverage_max <= 0) coverage_max <- 1

    graphics::plot(
      NA,
      xlim = c(0, chromosome_length / 1e6),
      ylim = c(0, coverage_max * 1.05),
      xlab = cfg$x_label,
      ylab = cfg$y_label,
      main = chromosome,
      xaxs = "i",
      yaxs = "i"
    )
    graphics::grid(col = cfg$grid_color, lty = 1)
    graphics::lines(
      profile$midpoint / 1e6,
      profile$coverage,
      col = cfg$line_color,
      lwd = cfg$line_width
    )
  }

  empty_panels <- panel_rows * panel_columns - length(chromosomes)
  if (empty_panels > 0) {
    for (i in seq_len(empty_panels)) graphics::plot.new()
  }
  if (nzchar(cfg$title)) graphics::mtext(cfg$title, outer = TRUE, side = 3, font = 2)
}, finally = grDevices::dev.off())

message("SVG 已输出：", normalizePath(cfg$output_file, winslash = "/", mustWork = FALSE))
