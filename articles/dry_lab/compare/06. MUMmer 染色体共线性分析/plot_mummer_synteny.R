#!/usr/bin/env Rscript

# ============================================================================
#                       用户参数配置区（只需修改这里）
# ============================================================================

# 分析脚本生成的 *.plot.tsv。留空时自动搜索当前目录：仅有一个则直接使用，
# 有多个则显示编号供选择。
INPUT_FILE <- ""

# 输出文件前缀。留空时自动使用输入文件名，例如 sample.plot.tsv -> sample
OUTPUT_PREFIX <- ""

# 输出类型
DRAW_SVG <- TRUE                 # ggplot2 静态矢量图
DRAW_HTML <- TRUE                # plotly 可缩放交互图

# 数据量与画布
MAX_ALIGNMENTS <- 100000         # 超过时保留比对长度最长的记录；设为 Inf 表示不限制
PLOT_WIDTH <- 12                 # SVG 宽度（英寸）
PLOT_HEIGHT <- 10                # SVG 高度（英寸）
MAX_AXIS_LABELS <- 40            # 每条坐标轴最多显示多少个序列标签

# 样式
FORWARD_COLOR <- "#2563EB"      # 正向比对：蓝色
REVERSE_COLOR <- "#DC2626"      # 反向比对：红色
LINE_WIDTH <- 0.32
LINE_ALPHA <- 0.72
BASE_FONT_SIZE <- 11

# ============================================================================
#                     程序主体（一般不需要修改）
# ============================================================================

choose_input_file <- function(configured_file) {
  if (nzchar(configured_file)) {
    path <- normalizePath(configured_file, mustWork = FALSE)
    if (!file.exists(path)) stop("找不到 INPUT_FILE：", path)
    return(path)
  }

  files <- list.files(".", pattern = "\\.plot\\.tsv$", recursive = TRUE, full.names = TRUE)
  if (length(files) == 0) stop("当前目录及子目录中没有找到 *.plot.tsv。")
  if (length(files) == 1) return(normalizePath(files[[1]]))

  cat("找到以下绘图数据文件：\n")
  for (i in seq_along(files)) cat(sprintf("  [%d] %s\n", i, files[[i]]))
  repeat {
    answer <- readline("请输入文件编号（q 退出）：")
    if (tolower(answer) == "q") quit(save = "no", status = 0)
    index <- suppressWarnings(as.integer(answer))
    if (!is.na(index) && index >= 1 && index <= length(files)) {
      return(normalizePath(files[[index]]))
    }
    cat("输入无效，请重新输入。\n")
  }
}

read_metadata <- function(input_file) {
  metadata_file <- sub("\\.plot\\.tsv$", ".plot.json", input_file)
  if (!file.exists(metadata_file) || !requireNamespace("jsonlite", quietly = TRUE)) return(list())
  tryCatch(jsonlite::fromJSON(metadata_file), error = function(error) {
    warning("无法读取元数据：", conditionMessage(error))
    list()
  })
}

make_layout <- function(ids, starts, ends) {
  ordered_ids <- unique(ids)
  sizes <- vapply(ordered_ids, function(id) {
    max(c(starts[ids == id], ends[ids == id]), na.rm = TRUE)
  }, numeric(1))
  offsets <- c(0, head(cumsum(sizes), -1))
  names(offsets) <- ordered_ids
  centers <- offsets + sizes / 2
  boundaries <- cumsum(sizes)
  list(ids = ordered_ids, sizes = sizes, offsets = offsets,
       centers = centers, boundaries = boundaries, total = sum(sizes))
}

limit_axis_labels <- function(centers, maximum) {
  if (length(centers) <= maximum) return(seq_along(centers))
  unique(round(seq(1, length(centers), length.out = maximum)))
}

input_file <- choose_input_file(INPUT_FILE)
if (!is.finite(MAX_ALIGNMENTS) && MAX_ALIGNMENTS != Inf) stop("MAX_ALIGNMENTS 设置无效。")
if (is.finite(MAX_ALIGNMENTS) && MAX_ALIGNMENTS < 1) stop("MAX_ALIGNMENTS 必须大于 0。")
if (!DRAW_SVG && !DRAW_HTML) stop("DRAW_SVG 和 DRAW_HTML 不能同时为 FALSE。")

cat("读取绘图数据：", input_file, "\n")
aln <- read.delim(input_file, header = TRUE, stringsAsFactors = FALSE,
                  check.names = FALSE, quote = "", comment.char = "")
required_columns <- c(
  "ref_start", "ref_end", "query_start", "query_end", "ref_align_length",
  "query_align_length", "identity", "ref_id", "query_id"
)
missing_columns <- setdiff(required_columns, names(aln))
if (length(missing_columns) > 0) stop("绘图文件缺少列：", paste(missing_columns, collapse = ", "))
if (nrow(aln) == 0) stop("绘图文件中没有比对记录。")

numeric_columns <- required_columns[seq_len(7)]
aln[numeric_columns] <- lapply(aln[numeric_columns], as.numeric)
aln <- aln[complete.cases(aln[required_columns]), , drop = FALSE]
if (nrow(aln) == 0) stop("移除无效记录后没有可绘制的数据。")

if (is.finite(MAX_ALIGNMENTS) && nrow(aln) > MAX_ALIGNMENTS) {
  aln <- aln[order(aln$ref_align_length, decreasing = TRUE), , drop = FALSE]
  aln <- aln[seq_len(MAX_ALIGNMENTS), , drop = FALSE]
  cat(sprintf("记录较多，已保留最长的 %d 条比对。\n", MAX_ALIGNMENTS))
}

metadata <- read_metadata(input_file)
reference_name <- if (!is.null(metadata$reference_name)) metadata$reference_name else "Reference"
query_name <- if (!is.null(metadata$query_name)) metadata$query_name else "Query"
output_prefix <- if (nzchar(OUTPUT_PREFIX)) OUTPUT_PREFIX else sub("\\.plot\\.tsv$", "", input_file)
output_prefix <- normalizePath(dirname(output_prefix), mustWork = TRUE) |>
  file.path(basename(output_prefix))

ref_layout <- make_layout(aln$ref_id, aln$ref_start, aln$ref_end)
query_layout <- make_layout(aln$query_id, aln$query_start, aln$query_end)
aln$x <- aln$ref_start + ref_layout$offsets[aln$ref_id]
aln$xend <- aln$ref_end + ref_layout$offsets[aln$ref_id]
aln$y <- aln$query_start + query_layout$offsets[aln$query_id]
aln$yend <- aln$query_end + query_layout$offsets[aln$query_id]
aln$orientation <- ifelse(
  (aln$ref_end - aln$ref_start) * (aln$query_end - aln$query_start) >= 0,
  "Forward", "Reverse"
)
aln$tooltip <- sprintf(
  paste0("Reference: %s<br>Ref: %s–%s<br>",
         "Query: %s<br>Query: %s–%s<br>Identity: %.2f%%<br>Length: %s bp"),
  aln$ref_id, aln$ref_start, aln$ref_end, aln$query_id,
  aln$query_start, aln$query_end, aln$identity, aln$ref_align_length
)

if (!requireNamespace("ggplot2", quietly = TRUE)) {
  stop("缺少必需的 R 包 ggplot2。请运行：install.packages('ggplot2')")
}
suppressPackageStartupMessages(library(ggplot2))

ref_labels <- limit_axis_labels(ref_layout$centers, MAX_AXIS_LABELS)
query_labels <- limit_axis_labels(query_layout$centers, MAX_AXIS_LABELS)
plot_object <- ggplot(aln) +
  geom_segment(
    aes(x = x, xend = xend, y = y, yend = yend,
        colour = orientation, text = tooltip),
    linewidth = LINE_WIDTH, alpha = LINE_ALPHA, lineend = "round"
  ) +
  geom_vline(xintercept = head(ref_layout$boundaries, -1),
             colour = "#E5E7EB", linewidth = 0.22) +
  geom_hline(yintercept = head(query_layout$boundaries, -1),
             colour = "#E5E7EB", linewidth = 0.22) +
  scale_colour_manual(
    values = c(Forward = FORWARD_COLOR, Reverse = REVERSE_COLOR),
    breaks = c("Forward", "Reverse")
  ) +
  scale_x_continuous(
    breaks = ref_layout$centers[ref_labels], labels = ref_layout$ids[ref_labels],
    expand = expansion(mult = c(0.005, 0.005))
  ) +
  scale_y_continuous(
    breaks = query_layout$centers[query_labels], labels = query_layout$ids[query_labels],
    expand = expansion(mult = c(0.005, 0.005))
  ) +
  coord_cartesian(clip = "off") +
  labs(
    title = paste(reference_name, "vs", query_name),
    subtitle = paste0(
      "MUMmer whole-genome alignment · ",
      format(nrow(aln), big.mark = ",", scientific = FALSE), " alignments"
    ),
    x = reference_name, y = query_name, colour = "Orientation"
  ) +
  theme_minimal(base_size = BASE_FONT_SIZE) +
  theme(
    panel.grid = element_blank(),
    panel.border = element_rect(colour = "#9CA3AF", fill = NA, linewidth = 0.4),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8, colour = "#374151"),
    axis.text.y = element_text(size = 8, colour = "#374151"),
    axis.title = element_text(face = "bold"),
    plot.title = element_text(face = "bold", size = 15, colour = "#111827"),
    plot.subtitle = element_text(colour = "#6B7280"),
    legend.position = "top",
    legend.title = element_blank(),
    plot.margin = margin(10, 18, 10, 10)
  )

if (DRAW_SVG) {
  svg_file <- paste0(output_prefix, ".synteny.svg")
  tryCatch({
    ggsave(svg_file, plot_object, device = grDevices::svg, width = PLOT_WIDTH,
           height = PLOT_HEIGHT, units = "in", limitsize = FALSE)
    cat("SVG：", svg_file, "\n")
  }, error = function(error) {
    warning("SVG 生成失败，但不会阻止 HTML：", conditionMessage(error))
  })
}

if (DRAW_HTML) {
  if (!requireNamespace("plotly", quietly = TRUE) ||
      !requireNamespace("htmlwidgets", quietly = TRUE)) {
    warning("缺少 plotly 或 htmlwidgets，跳过 HTML；SVG 不受影响。")
  } else {
    html_file <- paste0(output_prefix, ".synteny.html")
    tryCatch({
      interactive_plot <- plotly::ggplotly(plot_object, tooltip = "text", dynamicTicks = TRUE)
      interactive_plot <- plotly::layout(
        interactive_plot,
        hovermode = "closest",
        legend = list(orientation = "h", x = 0.5, xanchor = "center", y = 1.08)
      )
      htmlwidgets::saveWidget(interactive_plot, html_file, selfcontained = TRUE)
      cat("交互图：", html_file, "\n")
    }, error = function(error) {
      warning("HTML 生成失败，但 SVG 和绘图数据不受影响：", conditionMessage(error))
    })
  }
}

cat("绘图完成。\n")
