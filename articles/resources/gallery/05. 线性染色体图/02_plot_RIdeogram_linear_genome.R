#!/usr/bin/env Rscript

# 基于 RIdeogram 绘制 FASTA 染色体线性图
# 运行方式：Rscript 02_plot_RIdeogram_linear_genome.R
# 只需修改下面的参数配置区，无需传入命令行参数。

# =============================================================================
# 参数配置区
# =============================================================================
config <- list(
  # 输入 FASTA；相对路径以本脚本所在目录为基准。
  fasta_file = "genome.fa",

  # 输出文件。RIdeogram 仅直接支持 SVG。
  output_file = "genome_RIdeogram.svg",

  # FALSE：仅在 .fai 缺失或旧于 FASTA 时重建；TRUE：每次强制重建。
  rebuild_index = FALSE,

  # 指定绘制的序列及顺序；NULL 表示自动选择。
  chromosomes = NULL,
  # chromosomes = c("chr1", "chr2", "chr3", "chrX"),

  # 可选的序列名过滤正则；NULL 表示不过滤。
  chromosome_regex = NULL,
  # chromosome_regex = "^chr([0-9]+|X|Y)$",

  # 过滤掉短于此长度的序列（bp）。
  min_length = 10000,

  # chromosomes = NULL 时的排序："fasta"、"natural" 或 "length"。
  chromosome_order = "natural",

  # 最多绘制多少条序列；NULL 表示不限制。组装结果含大量 contig 时建议设置。
  max_chromosomes = NULL,

  # RIdeogram 的绘图区宽度参数，不是英寸。染色体越多通常应设置得越大。
  plot_width = 170,

  # RIdeogram 默认使用白色填充，在白色背景中很难看清；生成后自动替换样式。
  chromosome_fill = "#4C78A8",
  chromosome_border = "#263238",

  # 可选的胶囊内部热图数据。NULL 表示不绘制热图。
  # 文件至少包含 Chr、Start、End、Value 四列，支持 .tsv/.txt/.csv。
  heatmap_file = NULL,
  # heatmap_file = "gene_density.tsv",

  # 如果你的列名不同，在这里填写实际列名。
  heatmap_columns = c(
    chromosome = "Chr",
    start = "Start",
    end = "End",
    value = "Value"
  ),

  # 热图低值、中值和高值颜色。
  heatmap_colors = c("#4575B4", "#FFFFBF", "#D73027"),

  # TRUE：超出染色体范围的区间截断到合法范围；FALSE：遇到越界即报错。
  clip_heatmap_ranges = TRUE,

  # 图例位置。当前仅画核型、没有图例，但保留为 RIdeogram 原生配置项。
  legend_x = 160,
  legend_y = 35
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
  key <- vapply(x, function(value) {
    parts <- regmatches(value, gregexpr("[0-9]+|[^0-9]+", value, perl = TRUE))[[1L]]
    paste(vapply(parts, function(part) {
      if (grepl("^[0-9]+$", part)) sprintf("%020.0f", as.numeric(part)) else tolower(part)
    }, character(1L)), collapse = "")
  }, character(1L))
  order(key, seq_along(x))
}

read_fasta_index <- function(fai_file) {
  # .fai 五列依次为：序列 ID、碱基长度、文件偏移、每行碱基数、每行字节数。
  fai <- tryCatch(
    utils::read.delim(
      fai_file,
      header = FALSE,
      sep = "\t",
      quote = "",
      comment.char = "",
      colClasses = c("character", "numeric", "numeric", "numeric", "numeric"),
      stringsAsFactors = FALSE
    ),
    error = function(e) stop("无法读取 FASTA 索引：", conditionMessage(e), call. = FALSE)
  )

  if (!nrow(fai) || ncol(fai) != 5L) {
    stop("FASTA 索引格式异常：.fai 应包含五列且至少有一条序列。", call. = FALSE)
  }
  names(fai) <- c("name", "length", "offset", "line_bases", "line_width")
  if (anyNA(fai$name) || any(!nzchar(fai$name))) stop("索引中存在空序列 ID。", call. = FALSE)

  duplicate_names <- unique(fai$name[duplicated(fai$name)])
  if (length(duplicate_names)) {
    stop(
      "FASTA 中存在重复序列 ID（标题中第一个空白前的字段必须唯一）：",
      paste(duplicate_names, collapse = ", "),
      call. = FALSE
    )
  }
  numeric_columns <- fai[c("length", "offset", "line_bases", "line_width")]
  if (anyNA(numeric_columns) || any(!is.finite(as.matrix(numeric_columns))) ||
      any(as.matrix(numeric_columns) < 0)) {
    stop("FASTA 索引中存在非法的长度、偏移或行宽。", call. = FALSE)
  }
  if (any(fai$length != floor(fai$length))) {
    stop("FASTA 索引中存在非整数序列长度。", call. = FALSE)
  }
  fai
}

standardize_svg <- function(svg_file, fill_color, border_color) {
  # RIdeogram 0.2.2 生成的根节点可能没有 SVG 命名空间和 viewBox，
  # 这会导致部分浏览器、网页 img 标签和文件预览器显示空白。
  svg_text <- paste(readLines(svg_file, warn = FALSE, encoding = "UTF-8"), collapse = "\n")
  if (!grepl("<svg\\b", svg_text, perl = TRUE)) {
    stop("输出文件中没有找到 <svg> 根节点。", call. = FALSE)
  }

  if (!grepl("<svg\\b[^>]*\\bxmlns=", svg_text, perl = TRUE)) {
    svg_text <- sub(
      "<svg\\b",
      '<svg xmlns="http://www.w3.org/2000/svg"',
      svg_text,
      perl = TRUE
    )
  }

  svg_tag <- regmatches(svg_text, regexpr("<svg\\b[^>]*>", svg_text, perl = TRUE))
  if (!grepl("\\bviewBox=", svg_tag, perl = TRUE)) {
    dimensions <- regexec(
      'width="([0-9.]+)"[^>]*height="([0-9.]+)"',
      svg_tag,
      perl = TRUE
    )
    values <- regmatches(svg_tag, dimensions)[[1L]]
    if (length(values) == 3L) {
      replacement <- sub(
        ">$",
        sprintf(' viewBox="0 0 %s %s">', values[[2L]], values[[3L]]),
        svg_tag
      )
      svg_text <- sub("<svg\\b[^>]*>", replacement, svg_text, perl = TRUE)
    }
  }

  # 无着丝粒三列模式中，这两类样式分别对应染色体填充和外边框。
  svg_text <- gsub(
    "fill:white; stroke:white;",
    paste0("fill:", fill_color, "; stroke:", fill_color, ";"),
    svg_text,
    fixed = TRUE
  )
  svg_text <- gsub("stroke:grey;", paste0("stroke:", border_color, ";"), svg_text, fixed = TRUE)

  writeLines(svg_text, con = svg_file, useBytes = TRUE)
}

read_heatmap_data <- function(file, column_map, chromosome_lengths, clip_ranges) {
  extension <- tolower(tools::file_ext(file))
  separator <- if (extension == "csv") "," else "\t"
  heatmap <- tryCatch(
    utils::read.table(
      file,
      header = TRUE,
      sep = separator,
      quote = '"',
      comment.char = "",
      check.names = FALSE,
      stringsAsFactors = FALSE
    ),
    error = function(e) stop("无法读取热图数据：", conditionMessage(e), call. = FALSE)
  )

  required_keys <- c("chromosome", "start", "end", "value")
  if (is.null(names(column_map)) || !all(required_keys %in% names(column_map))) {
    stop("heatmap_columns 必须包含 chromosome、start、end、value 四个命名项。", call. = FALSE)
  }
  source_columns <- unname(column_map[required_keys])
  missing_columns <- setdiff(source_columns, names(heatmap))
  if (length(missing_columns)) {
    stop("热图文件缺少列：", paste(missing_columns, collapse = ", "), call. = FALSE)
  }

  heatmap <- heatmap[source_columns]
  names(heatmap) <- c("Chr", "Start", "End", "Value")
  heatmap$Chr <- as.character(heatmap$Chr)
  heatmap$Start <- suppressWarnings(as.numeric(heatmap$Start))
  heatmap$End <- suppressWarnings(as.numeric(heatmap$End))
  heatmap$Value <- suppressWarnings(as.numeric(heatmap$Value))

  if (!nrow(heatmap)) stop("热图文件没有数据行。", call. = FALSE)
  if (anyNA(heatmap) || any(!is.finite(as.matrix(heatmap[c("Start", "End", "Value")]))) ) {
    stop("热图数据包含空值、非数字坐标或非数字 Value。", call. = FALSE)
  }
  if (any(heatmap$Start < 0) || any(heatmap$End <= heatmap$Start)) {
    stop("热图区间必须满足 Start >= 0 且 End > Start。", call. = FALSE)
  }

  unknown <- setdiff(unique(heatmap$Chr), names(chromosome_lengths))
  if (length(unknown)) {
    warning("热图中以下序列未被绘制，相关行已跳过：", paste(unknown, collapse = ", "))
    heatmap <- heatmap[heatmap$Chr %in% names(chromosome_lengths), , drop = FALSE]
  }
  if (!nrow(heatmap)) stop("过滤后没有可绘制的热图数据。", call. = FALSE)

  maximum_end <- unname(chromosome_lengths[heatmap$Chr])
  out_of_range <- heatmap$End > maximum_end
  if (any(out_of_range)) {
    if (!isTRUE(clip_ranges)) {
      examples <- unique(heatmap$Chr[out_of_range])
      stop("热图区间超出染色体长度：", paste(head(examples, 10L), collapse = ", "), call. = FALSE)
    }
    warning("有 ", sum(out_of_range), " 个热图区间超出染色体末端，已自动截断。")
    heatmap$End[out_of_range] <- maximum_end[out_of_range]
  }
  heatmap <- heatmap[heatmap$End > heatmap$Start, , drop = FALSE]
  if (!nrow(heatmap)) stop("截断后没有有效的热图区间。", call. = FALSE)
  heatmap
}

cfg <- config
script_dir <- get_script_dir()
cfg$fasta_file <- resolve_path(cfg$fasta_file, script_dir)
cfg$output_file <- resolve_path(cfg$output_file, script_dir)
if (!is.null(cfg$heatmap_file)) cfg$heatmap_file <- resolve_path(cfg$heatmap_file, script_dir)

required_packages <- c("Rsamtools", "RIdeogram")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1L), quietly = TRUE)
]
if (length(missing_packages)) {
  stop(
    "缺少 R 包：", paste(missing_packages, collapse = ", "),
    "\n安装方法：\n",
    "if (!requireNamespace('BiocManager', quietly=TRUE)) install.packages('BiocManager')\n",
    "BiocManager::install('Rsamtools')\n",
    "install.packages('RIdeogram')",
    call. = FALSE
  )
}

if (length(cfg$fasta_file) != 1L || !nzchar(cfg$fasta_file) || !file.exists(cfg$fasta_file)) {
  stop("FASTA 文件不存在：", cfg$fasta_file, call. = FALSE)
}
if (length(cfg$output_file) != 1L || tolower(tools::file_ext(cfg$output_file)) != "svg") {
  stop("output_file 必须是一个 .svg 文件。", call. = FALSE)
}
if (!cfg$chromosome_order %in% c("fasta", "natural", "length")) {
  stop("chromosome_order 只能是 fasta、natural 或 length。", call. = FALSE)
}
if (length(cfg$min_length) != 1L || !is.finite(cfg$min_length) || cfg$min_length < 0) {
  stop("min_length 必须是非负数。", call. = FALSE)
}
if (length(cfg$plot_width) != 1L || !is.finite(cfg$plot_width) || cfg$plot_width <= 0) {
  stop("plot_width 必须是正数。", call. = FALSE)
}
for (name in c("legend_x", "legend_y")) {
  if (length(cfg[[name]]) != 1L || !is.finite(cfg[[name]]) || cfg[[name]] < 0) {
    stop(name, " 必须是非负数。", call. = FALSE)
  }
}
if (!is.null(cfg$max_chromosomes) &&
    (length(cfg$max_chromosomes) != 1L || !is.finite(cfg$max_chromosomes) ||
     cfg$max_chromosomes < 1 || cfg$max_chromosomes != floor(cfg$max_chromosomes))) {
  stop("max_chromosomes 必须是 NULL 或正整数。", call. = FALSE)
}
for (name in c("chromosome_fill", "chromosome_border")) {
  if (length(cfg[[name]]) != 1L || is.na(cfg[[name]]) || !nzchar(cfg[[name]]) ||
      !grepl("^[#A-Za-z0-9(),. %+_-]+$", cfg[[name]])) {
    stop(name, " 不是有效的安全颜色字符串。", call. = FALSE)
  }
}
if (length(cfg$heatmap_colors) < 2L || anyNA(cfg$heatmap_colors) ||
    any(!nzchar(cfg$heatmap_colors))) {
  stop("heatmap_colors 至少需要两个非空颜色。", call. = FALSE)
}
if (!is.null(cfg$heatmap_file) &&
    (length(cfg$heatmap_file) != 1L || !nzchar(cfg$heatmap_file) ||
     !file.exists(cfg$heatmap_file))) {
  stop("热图文件不存在：", cfg$heatmap_file, call. = FALSE)
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
        "\n请确认 FASTA 未使用普通 gzip 压缩、序列 ID 唯一且输入目录可写。",
        call. = FALSE
      )
    }
  )
} else {
  message("使用已有 FASTA 索引：", fai_file)
}

# 只读取 .fai 第 2 列作为序列长度。FASTA 的 > 标题行及标题描述不会计入长度。
fasta_index <- read_fasta_index(fai_file)
sequence_names <- fasta_index$name
sequence_lengths <- fasta_index$length
names(sequence_lengths) <- sequence_names

empty_sequences <- sequence_names[sequence_lengths == 0]
if (length(empty_sequences)) {
  warning("以下空序列无法绘制，已跳过：", paste(empty_sequences, collapse = ", "))
}

keep <- sequence_lengths > 0 & sequence_lengths >= cfg$min_length
if (!is.null(cfg$chromosome_regex)) {
  if (length(cfg$chromosome_regex) != 1L || is.na(cfg$chromosome_regex)) {
    stop("chromosome_regex 必须是 NULL 或单个正则表达式字符串。", call. = FALSE)
  }
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
  if (length(filtered_names)) warning("未通过过滤条件，已跳过：", paste(filtered_names, collapse = ", "))
  chromosomes <- chromosomes[chromosomes %in% available_names]
}

if (!is.null(cfg$max_chromosomes) && length(chromosomes) > cfg$max_chromosomes) {
  warning("序列数量超过 max_chromosomes，仅绘制排序后的前 ", cfg$max_chromosomes, " 条。")
  chromosomes <- head(chromosomes, cfg$max_chromosomes)
}
if (!length(chromosomes)) stop("没有符合配置条件的非空序列可供绘制。", call. = FALSE)

# RIdeogram 无着丝粒信息时使用三列核型。Start=0，End 为真实碱基长度。
karyotype <- data.frame(
  Chr = chromosomes,
  Start = rep(0, length(chromosomes)),
  End = unname(sequence_lengths[chromosomes]),
  stringsAsFactors = FALSE,
  check.names = FALSE
)

heatmap_data <- NULL
if (!is.null(cfg$heatmap_file)) {
  selected_lengths <- sequence_lengths[chromosomes]
  heatmap_data <- read_heatmap_data(
    cfg$heatmap_file,
    cfg$heatmap_columns,
    selected_lengths,
    cfg$clip_heatmap_ranges
  )
  message("已载入 ", nrow(heatmap_data), " 个热图区间。")
}

dir.create(dirname(cfg$output_file), recursive = TRUE, showWarnings = FALSE)
tryCatch(
  RIdeogram::ideogram(
    karyotype = karyotype,
    overlaid = heatmap_data,
    colorset1 = cfg$heatmap_colors,
    width = cfg$plot_width,
    Lx = cfg$legend_x,
    Ly = cfg$legend_y,
    output = cfg$output_file
  ),
  error = function(e) {
    if (file.exists(cfg$output_file) && file.info(cfg$output_file)$size == 0) {
      unlink(cfg$output_file)
    }
    stop("RIdeogram 绘图失败：", conditionMessage(e), call. = FALSE)
  }
)

if (!file.exists(cfg$output_file) || file.info(cfg$output_file)$size <= 0) {
  stop("RIdeogram 未生成有效的 SVG 文件。", call. = FALSE)
}

tryCatch(
  standardize_svg(cfg$output_file, cfg$chromosome_fill, cfg$chromosome_border),
  error = function(e) stop("SVG 标准化失败：", conditionMessage(e), call. = FALSE)
)

message("索引检查完成：长度来自 .fai 第 2 列，仅统计序列碱基。")
message("已绘制 ", length(chromosomes), " 条序列。")
message("SVG 已输出：", normalizePath(cfg$output_file, winslash = "/", mustWork = TRUE))
