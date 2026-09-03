#!/usr/bin/env Rscript

# 基于 karyoploteR 绘制 FASTA 染色体线性图
# 运行方式：Rscript 03_plot_karyoploteR_linear_genome.R
# 只需修改下面的参数配置区，无需传入命令行参数。

# =============================================================================
# 参数配置区
# =============================================================================
config <- list(
  # 输入 FASTA；相对路径以本脚本所在目录为基准。
  fasta_file = "genome.fa",

  # 输出 SVG；相对路径以本脚本所在目录为基准。
  output_file = "genome_linear_karyotype.svg",

  # 是否强制重建 FASTA 索引。FALSE 时仅在 .fai 缺失或早于 FASTA 时重建。
  rebuild_index = FALSE,

  # 需要绘制的序列名称及顺序。NULL 表示使用 FASTA 中的全部序列。
  chromosomes = NULL,
  # chromosomes = c("chr1", "chr2", "chr3", "chrX"),

  # 可选名称过滤正则表达式；NULL 表示不过滤。
  chromosome_regex = NULL,
  # chromosome_regex = "^chr([0-9]+|X|Y)$",

  # 只绘制长度不小于该值的序列（bp）。设为 0 保留全部序列。
  min_length = 0,

  # chromosomes = NULL 时的排列方式："fasta"、"natural" 或 "length"。
  chromosome_order = "natural",

  # 图形尺寸（英寸）。染色体较多时可增大 height。
  width = 12,
  height = 8,

  # SVG 字体与字号。
  font_family = "Arial",
  base_font_size = 10,

  # 染色体样式。
  chromosome_color = "#4C78A8",
  chromosome_border = "#263238",

  # 是否显示染色体名称和坐标刻度。
  show_chromosome_names = TRUE,
  show_base_numbers = TRUE,
  tick_distance = 10000000,       # 主刻度间距（bp）
  minor_tick_distance = 2000000,  # 次刻度间距（bp）
  tick_label_cex = 0.75,
  tick_units = "Mb",             # 可选："bp"、"Kb"、"Mb"、"Gb"

  # 页面边距；依次为下、左、上、右。
  margins = c(0.08, 0.12, 0.06, 0.03)
)
# =============================================================================
# 参数配置区结束
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

natural_order <- function(x) {
  # 将连续数字补零，使 chr2 排在 chr10 前面。
  key <- vapply(x, function(value) {
    parts <- regmatches(value, gregexpr("[0-9]+|[^0-9]+", value, perl = TRUE))[[1L]]
    paste(vapply(parts, function(part) {
      if (grepl("^[0-9]+$", part)) sprintf("%020.0f", as.numeric(part)) else tolower(part)
    }, character(1L)), collapse = "")
  }, character(1L))
  order(key, seq_along(x))
}

cfg <- config
script_dir <- get_script_dir()
cfg$fasta_file <- resolve_path(cfg$fasta_file, script_dir)
cfg$output_file <- resolve_path(cfg$output_file, script_dir)

required_packages <- c("Rsamtools", "karyoploteR")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1L), quietly = TRUE)
]
if (length(missing_packages)) {
  stop(
    "缺少 R/Bioconductor 包：", paste(missing_packages, collapse = ", "),
    "\n请先运行：\n",
    "if (!requireNamespace('BiocManager', quietly=TRUE)) install.packages('BiocManager')\n",
    "BiocManager::install(c('Rsamtools', 'karyoploteR'))",
    call. = FALSE
  )
}

if (!file.exists(cfg$fasta_file)) stop("FASTA 文件不存在：", cfg$fasta_file, call. = FALSE)
if (tolower(tools::file_ext(cfg$output_file)) != "svg") {
  stop("output_file 必须使用 .svg 扩展名。", call. = FALSE)
}
if (!cfg$chromosome_order %in% c("fasta", "natural", "length")) {
  stop("chromosome_order 只能是 fasta、natural 或 length。", call. = FALSE)
}
if (!cfg$tick_units %in% c("bp", "Kb", "Mb", "Gb")) {
  stop("tick_units 只能是 bp、Kb、Mb 或 Gb。", call. = FALSE)
}
numeric_positive <- c("width", "height", "base_font_size", "tick_distance", "tick_label_cex")
for (name in numeric_positive) {
  if (length(cfg[[name]]) != 1L || !is.finite(cfg[[name]]) || cfg[[name]] <= 0) {
    stop(name, " 必须是正数。", call. = FALSE)
  }
}
if (!is.finite(cfg$minor_tick_distance) || cfg$minor_tick_distance < 0) {
  stop("minor_tick_distance 必须是非负数。", call. = FALSE)
}
if (!is.finite(cfg$min_length) || cfg$min_length < 0) {
  stop("min_length 必须是非负数。", call. = FALSE)
}
if (length(cfg$margins) != 4L || any(!is.finite(cfg$margins)) || any(cfg$margins < 0)) {
  stop("margins 必须是由 4 个非负数构成的向量。", call. = FALSE)
}

fai_file <- paste0(cfg$fasta_file, ".fai")
fasta_time <- file.info(cfg$fasta_file)$mtime
index_time <- if (file.exists(fai_file)) file.info(fai_file)$mtime else as.POSIXct(NA)
need_index <- isTRUE(cfg$rebuild_index) || !file.exists(fai_file) ||
  is.na(index_time) || index_time < fasta_time

if (need_index) {
  message("正在构建 FASTA 索引：", fai_file)
  tryCatch(
    Rsamtools::indexFa(cfg$fasta_file),
    error = function(e) {
      stop(
        "FASTA 索引构建失败：", conditionMessage(e),
        "\n请确认 FASTA 为未压缩文件、序列名称唯一，并且输入目录可写。",
        call. = FALSE
      )
    }
  )
} else {
  message("使用已有 FASTA 索引：", fai_file)
}

fasta_index <- Rsamtools::scanFaIndex(cfg$fasta_file)
sequence_names <- as.character(GenomeInfoDb::seqnames(fasta_index))
sequence_lengths <- as.numeric(IRanges::width(fasta_index))
names(sequence_lengths) <- sequence_names

keep <- sequence_lengths >= cfg$min_length
if (!is.null(cfg$chromosome_regex)) {
  keep <- keep & grepl(cfg$chromosome_regex, sequence_names, perl = TRUE)
}
available_names <- sequence_names[keep]

if (is.null(cfg$chromosomes)) {
  chromosomes <- available_names
  if (cfg$chromosome_order == "natural") {
    chromosomes <- chromosomes[natural_order(chromosomes)]
  } else if (cfg$chromosome_order == "length") {
    chromosomes <- chromosomes[order(sequence_lengths[chromosomes], decreasing = TRUE)]
  }
} else {
  chromosomes <- unique(as.character(cfg$chromosomes))
  missing_names <- setdiff(chromosomes, sequence_names)
  filtered_names <- setdiff(intersect(chromosomes, sequence_names), available_names)
  if (length(missing_names)) warning("FASTA 中不存在，已跳过：", paste(missing_names, collapse = ", "))
  if (length(filtered_names)) warning("未通过长度或正则过滤，已跳过：", paste(filtered_names, collapse = ", "))
  chromosomes <- chromosomes[chromosomes %in% available_names]
}
if (!length(chromosomes)) stop("没有符合配置条件的序列可供绘制。", call. = FALSE)

genome <- data.frame(
  chr = chromosomes,
  start = rep(1, length(chromosomes)),
  end = unname(sequence_lengths[chromosomes]),
  stringsAsFactors = FALSE
)

dir.create(dirname(cfg$output_file), recursive = TRUE, showWarnings = FALSE)
grDevices::svg(
  filename = cfg$output_file,
  width = cfg$width,
  height = cfg$height,
  family = cfg$font_family,
  pointsize = cfg$base_font_size,
  bg = "white",
  onefile = TRUE
)

tryCatch({
  plot_params <- karyoploteR::getDefaultPlotParams(plot.type = 1)
  plot_params$ideogramheight <- 15
  plot_params$data1height <- 0
  plot_params$bottommargin <- cfg$margins[[1L]]
  plot_params$leftmargin <- cfg$margins[[2L]]
  plot_params$topmargin <- cfg$margins[[3L]]
  plot_params$rightmargin <- cfg$margins[[4L]]

  kp <- karyoploteR::plotKaryotype(
    genome = genome,
    chromosomes = chromosomes,
    plot.type = 1,
    plot.params = plot_params,
    ideogram.plotter = NULL,
    labels.plotter = NULL
  )

  karyoploteR::kpRect(
    kp,
    chr = genome$chr,
    x0 = genome$start,
    x1 = genome$end,
    y0 = 0.22,
    y1 = 0.78,
    col = cfg$chromosome_color,
    border = cfg$chromosome_border,
    data.panel = "ideogram"
  )

  if (isTRUE(cfg$show_chromosome_names)) {
    karyoploteR::kpAddChromosomeNames(kp, cex = 0.9)
  }
  if (isTRUE(cfg$show_base_numbers)) {
    karyoploteR::kpAddBaseNumbers(
      kp,
      tick.dist = cfg$tick_distance,
      minor.tick.dist = cfg$minor_tick_distance,
      tick.len = 5,
      minor.tick.len = 2.5,
      cex = cfg$tick_label_cex,
      add.units = TRUE,
      units = cfg$tick_units
    )
  }
}, finally = grDevices::dev.off())

message("已绘制 ", length(chromosomes), " 条序列。")
message("SVG 已输出：", normalizePath(cfg$output_file, winslash = "/", mustWork = FALSE))
