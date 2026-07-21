#!/usr/bin/env Rscript

# 将一个 BigWig 文件中的所有染色体首尾拼接，绘制一条全基因组 coverage 轨道。
# 使用 svglite 输出 SVG，文字可在 Adobe Illustrator 中继续编辑。

# =============================================================================
# 参数设置区：只需修改本区内容
# =============================================================================
config <- list(
  # 输入一个 BigWig 文件；相对路径以本脚本所在目录为基准。
  input_file = "sample.bw",

  # 输出 SVG 文件；相对路径以本脚本所在目录为基准。
  output_file = "genome_coverage_track.svg",

  # 染色体及排列顺序。NULL 表示使用 BigWig 文件中的全部染色体及原始顺序。
  chromosomes = NULL,
  # chromosomes = c("chr1", "chr2", "chr3", "chr4", "chr5"),

  # 统计窗口大小（bp）。每个窗口计算 BigWig 信号的加权平均值。
  bin_size = 250000L,

  # 相邻染色体之间的空隙（bp）；设置为 0 表示首尾相接。
  chromosome_gap = 0,

  # Y 轴范围。y_max = Inf 表示不限制，自动使用数据最大值。
  # 设置为有限数值（如 50）时，超过上限的信号会被裁切。
  y_min = 0,
  y_max = Inf,

  # SVG 尺寸（英寸）。
  figure_width = 12,
  figure_height = 2.4,

  # coverage 填充和轮廓线样式。
  fill_color = "#4C78A8",
  fill_alpha = 0.75,
  line_color = "#1F4E79",
  line_width = 0.6,

  # 染色体边界和标签样式。
  boundary_color = "#808080",
  boundary_line_type = "dotted",
  boundary_line_width = 0.5,
  show_chromosome_labels = TRUE,
  chromosome_label_size = 0.85,

  # 坐标轴与字体设置。
  y_label = "BigWig signal",
  show_y_axis = TRUE,
  font_family = "Arial",
  base_font_size = 10,
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
    chromosome = chromosome,
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
if (!is.finite(cfg$bin_size) || cfg$bin_size < 1) stop("bin_size 必须是正整数。", call. = FALSE)
if (!is.finite(cfg$chromosome_gap) || cfg$chromosome_gap < 0) {
  stop("chromosome_gap 必须大于或等于 0。", call. = FALSE)
}
if (!is.finite(cfg$y_min)) stop("y_min 必须是有限数值。", call. = FALSE)
if (length(cfg$y_max) != 1L || is.na(cfg$y_max) || cfg$y_max <= cfg$y_min) {
  stop("y_max 必须为 Inf，或为大于 y_min 的有限数值。", call. = FALSE)
}
if (!is.finite(cfg$figure_width) || cfg$figure_width <= 0 ||
    !is.finite(cfg$figure_height) || cfg$figure_height <= 0) {
  stop("figure_width 和 figure_height 必须是正数。", call. = FALSE)
}
if (!is.finite(cfg$fill_alpha) || cfg$fill_alpha < 0 || cfg$fill_alpha > 1) {
  stop("fill_alpha 必须在 0 到 1 之间。", call. = FALSE)
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

chromosome_lengths <- as.numeric(sequence_lengths[chromosomes])
names(chromosome_lengths) <- chromosomes
chromosome_starts <- c(0, head(cumsum(chromosome_lengths + cfg$chromosome_gap), -1L))
chromosome_ends <- chromosome_starts + chromosome_lengths
chromosome_centers <- (chromosome_starts + chromosome_ends) / 2
genome_end <- tail(chromosome_ends, 1L)

pieces <- lapply(seq_along(chromosomes), function(i) {
  out <- mean_signal_in_bins(
    bigwig,
    chromosomes[[i]],
    chromosome_lengths[[i]],
    cfg$bin_size
  )
  out$genome_position <- chromosome_starts[[i]] + out$midpoint
  out
})
profile <- do.call(rbind, pieces)

auto_y_max <- suppressWarnings(max(profile$coverage, na.rm = TRUE))
if (!is.finite(auto_y_max) || auto_y_max <= cfg$y_min) auto_y_max <- cfg$y_min + 1
plot_y_max <- if (is.finite(cfg$y_max)) cfg$y_max else auto_y_max

dir.create(dirname(cfg$output_file), recursive = TRUE, showWarnings = FALSE)
svglite::svglite(
  file = cfg$output_file,
  width = cfg$figure_width,
  height = cfg$figure_height,
  bg = cfg$background_color,
  pointsize = cfg$base_font_size
)

tryCatch({
  graphics::par(
    mar = c(if (isTRUE(cfg$show_chromosome_labels)) 2.2 else 0.6, 4.5, 0.4, 0.4),
    family = cfg$font_family,
    las = 1,
    bty = "n",
    xpd = FALSE
  )
  graphics::plot(
    NA,
    xlim = c(0, genome_end),
    ylim = c(cfg$y_min, plot_y_max),
    xaxs = "i",
    yaxs = "i",
    xaxt = "n",
    yaxt = if (isTRUE(cfg$show_y_axis)) "s" else "n",
    xlab = "",
    ylab = cfg$y_label
  )

  for (piece in split(profile, factor(profile$chromosome, levels = chromosomes))) {
    graphics::polygon(
      c(piece$genome_position, rev(piece$genome_position)),
      c(piece$coverage, rep(cfg$y_min, nrow(piece))),
      col = grDevices::adjustcolor(cfg$fill_color, alpha.f = cfg$fill_alpha),
      border = NA
    )
    graphics::lines(
      piece$genome_position,
      piece$coverage,
      col = cfg$line_color,
      lwd = cfg$line_width
    )
  }

  if (length(chromosome_ends) > 1L) {
    graphics::abline(
      v = chromosome_ends[-length(chromosome_ends)],
      col = cfg$boundary_color,
      lty = cfg$boundary_line_type,
      lwd = cfg$boundary_line_width
    )
  }
  graphics::box(bty = "l")

  if (isTRUE(cfg$show_chromosome_labels)) {
    graphics::axis(
      side = 1,
      at = chromosome_centers,
      labels = chromosomes,
      tick = FALSE,
      line = 0.3,
      cex.axis = cfg$chromosome_label_size
    )
  }
}, finally = grDevices::dev.off())

message("全基因组 coverage 轨道已输出：", normalizePath(cfg$output_file, winslash = "/", mustWork = FALSE))
