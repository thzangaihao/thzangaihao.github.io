# ============================================================
# 分组条形图 + 原始散点 + 标准差 + 两两显著性检验
# 输入格式：第一列为组名，其余各列均为重复样本。
# 默认以第一条数据行作为对照组（例如 CK）。
# ============================================================

# ------------------------- 用户设置区 -------------------------
INPUT_FILE  <- "C:/Users/thzan/Desktop/experimental_data.csv"
OUTPUT_FILE <- "experimental_data_barplot.svg"

# 图片尺寸（单位：英寸）及绘图方向。
# ORIENTATION："vertical" 为纵向条形图；"horizontal" 为横向条形图。
ORIENTATION <- "vertical"
FIG_WIDTH   <- 10
FIG_HEIGHT  <- 7

# 是否颠倒分组的显示方向。
# TRUE：纵向条形图从右向左排列；横向条形图从上向下排列。
# FALSE：使用 ggplot2 的默认因子排列方向。
REVERSE_GROUP_ORDER <- FALSE

# 配色：只填写一种颜色时所有组共用该颜色；填写多种颜色时按数据行顺序匹配。
# 也支持命名颜色，例如 c(CK="#4C78A8", 处理组="#F58518")。
BAR_COLORS <- c(
  "#F4C95D", "#3B6EA8", "#3B6EA8", "#3B6EA8", "#3B6EA8", "#3B6EA8", "#3B6EA8", "#3B6EA8", "#3B6EA8"
)
BAR_ALPHA   <- 0.82
POINT_COLOR <- "#202020"
POINT_SIZE  <- 2.3
POINT_ALPHA <- 0.80
POINT_JITTER_WIDTH <- 0.12

# 标题、坐标轴标签和字体。
PLOT_TITLE <- "Experimental results"
X_LABEL    <- NULL
Y_LABEL    <- "Value"
BASE_FONT_FAMILY <- "sans"
BASE_FONT_SIZE   <- 12

# 误差棒及统计检验设置。
ERROR_BAR <- "SD"             # 当前使用标准差：均值 ± SD
CONTROL_ROW <- 1L              # 第几条数据行作为对照组
ALPHA_VAR   <- 0.05            # F 检验 p >= 此值时认为方差齐
ALPHA_SIG   <- 0.05            # 显著性阈值，同时写入 SVG 图注
P_ADJUST_METHOD <- "none"      # 可设为 "none"、"holm"、"bonferroni"、"BH" 等
T_TEST_ALTERNATIVE <- "two.sided"

# 显著性符号；开启多重比较校正时，按校正后的 p 值标注。
P_CUTS   <- c(0, 0.001, 0.01, 0.05, Inf)
P_LABELS <- c("***", "**", "*", "ns")

# 显著性符号与该组误差棒顶端或最高散点之间的距离，占总数据范围的比例。
ANNOTATION_OFFSET <- 0.035

# 散点抖动的随机种子；只影响点的位置，不影响统计结果。
JITTER_SEED <- 20260903
# ----------------------------------------------------------------

if (!requireNamespace("ggplot2", quietly = TRUE)) {
  stop("Package 'ggplot2' is required. Install it with install.packages('ggplot2').")
}

if (!file.exists(INPUT_FILE)) stop("Input file does not exist: ", INPUT_FILE)
if (!ORIENTATION %in% c("vertical", "horizontal")) {
  stop("ORIENTATION must be either 'vertical' or 'horizontal'.")
}
if (length(P_LABELS) != length(P_CUTS) - 1L) {
  stop("length(P_LABELS) must equal length(P_CUTS) - 1.")
}

# check.names=FALSE 用于保留原始样本表头；fileEncoding 用于处理 UTF-8 BOM。
raw <- read.csv(
  INPUT_FILE, header = TRUE, check.names = FALSE,
  stringsAsFactors = FALSE, fileEncoding = "UTF-8-BOM",
  na.strings = c("", "NA", "NaN")
)
if (ncol(raw) < 2L) stop("The CSV needs one group column and at least one replicate column.")
if (nrow(raw) < 2L) stop("At least two groups (control + treatment) are required.")

group_col <- names(raw)[1L]
group_names <- trimws(as.character(raw[[1L]]))
if (anyNA(group_names) || any(group_names == "") || anyDuplicated(group_names)) {
  stop("Group names in the first column must be non-empty and unique.")
}
if (CONTROL_ROW < 1L || CONTROL_ROW > nrow(raw)) stop("CONTROL_ROW is out of range.")

# 将所有非空重复值转换为数值；遇到非数值文本时停止运行并提示。
replicate_block <- raw[-1L]
numeric_block <- lapply(replicate_block, function(x) suppressWarnings(as.numeric(x)))
bad_cells <- Map(function(original, converted) {
  !is.na(original) & trimws(as.character(original)) != "" & is.na(converted)
}, replicate_block, numeric_block)
if (any(unlist(bad_cells, use.names = FALSE))) {
  stop("A non-numeric value was found in the replicate columns.")
}
numeric_block <- as.data.frame(numeric_block, check.names = FALSE)

long_list <- lapply(seq_len(nrow(raw)), function(i) {
  values <- as.numeric(unlist(numeric_block[i, , drop = FALSE], use.names = FALSE))
  keep <- is.finite(values)
  data.frame(
    group = group_names[i],
    replicate = names(numeric_block)[keep],
    value = values[keep],
    stringsAsFactors = FALSE
  )
})
long <- do.call(rbind, long_list)
long$group <- factor(long$group, levels = group_names)

sample_sizes <- vapply(long_list, nrow, integer(1L))
if (any(sample_sizes < 2L)) {
  stop("Every group needs at least two finite observations for SD and t-tests. Problem group(s): ",
       paste(group_names[sample_sizes < 2L], collapse = ", "))
}

summary_df <- do.call(rbind, lapply(seq_along(long_list), function(i) {
  x <- long_list[[i]]$value
  data.frame(
    group = group_names[i], mean = mean(x), sd = stats::sd(x), n = length(x),
    stringsAsFactors = FALSE
  )
}))
summary_df$group <- factor(summary_df$group, levels = group_names)

# 每个非对照组分别与对照组比较。对于两组数据，var.test 执行经典的
# 双侧 F 方差齐性检验。
control_name <- group_names[CONTROL_ROW]
control_values <- long_list[[CONTROL_ROW]]$value
treatment_rows <- setdiff(seq_len(nrow(raw)), CONTROL_ROW)

tests <- do.call(rbind, lapply(treatment_rows, function(i) {
  treatment_values <- long_list[[i]]$value
  variance_test <- stats::var.test(
    treatment_values, control_values, alternative = "two.sided"
  )
  equal_variance <- is.finite(variance_test$p.value) && variance_test$p.value >= ALPHA_VAR
  t_result <- stats::t.test(
    treatment_values, control_values,
    alternative = T_TEST_ALTERNATIVE,
    var.equal = equal_variance
  )
  data.frame(
    group = group_names[i],
    variance_p = variance_test$p.value,
    equal_variance = equal_variance,
    test_method = if (equal_variance) "Student" else "Welch",
    p_raw = t_result$p.value,
    stringsAsFactors = FALSE
  )
}))
tests$p_used <- stats::p.adjust(tests$p_raw, method = P_ADJUST_METHOD)
tests$significance <- as.character(cut(
  tests$p_used, breaks = P_CUTS, labels = P_LABELS,
  include.lowest = TRUE, right = TRUE
))

tests$annotation <- tests$significance

data_min <- min(c(long$value, 0), na.rm = TRUE)
data_max <- max(summary_df$mean + summary_df$sd, na.rm = TRUE)
data_range <- max(data_max - data_min, abs(data_max) * 0.1, 1e-8)

# 每组显著性符号分别紧贴在“误差棒顶端”和“最高散点”中较高者的上方。
group_top <- vapply(seq_along(long_list), function(i) {
  max(max(long_list[[i]]$value), summary_df$mean[i] + summary_df$sd[i])
}, numeric(1L))
names(group_top) <- group_names
tests$y <- unname(group_top[tests$group]) + ANNOTATION_OFFSET * data_range

# 将完整分析方法写入统计 CSV，不在图片下方显示。
adjustment_note <- if (P_ADJUST_METHOD == "none") {
  "不进行多重比较校正，使用原始 p 值"
} else {
  paste0("使用 ", P_ADJUST_METHOD, " 方法校正 p 值")
}
tests$方差检验说明 <- paste0(
  "各处理组分别与 ", control_name,
  " 进行双侧 F 方差齐性检验；p >= ", ALPHA_VAR, " 判定为方差齐"
)
tests$t检验选择规则 <- paste0(
  "方差齐时使用双侧 Student 等方差 t 检验，方差不齐时使用双侧 Welch t 检验；",
  adjustment_note
)
tests$显著性符号规则 <- paste0(
  "***: p<0.001；**: p<0.01；*: p<", ALPHA_SIG,
  "；ns: p>=", ALPHA_SIG
)
tests$作图说明 <- "条形为均值，误差棒为标准差（SD），散点为全部有效观测值"

# 根据数据中的组顺序确定配色。
if (!is.null(names(BAR_COLORS)) && all(group_names %in% names(BAR_COLORS))) {
  palette <- BAR_COLORS[group_names]
} else {
  palette <- rep(BAR_COLORS, length.out = length(group_names))
  names(palette) <- group_names
}

# 根据开关确定最终显示顺序。
display_levels <- if (REVERSE_GROUP_ORDER) rev(group_names) else group_names
long$group <- factor(as.character(long$group), levels = display_levels)
summary_df$group <- factor(as.character(summary_df$group), levels = display_levels)
tests$group <- factor(as.character(tests$group), levels = display_levels)

set.seed(JITTER_SEED)
p <- ggplot2::ggplot(summary_df, ggplot2::aes(x = group, y = mean, fill = group)) +
  ggplot2::geom_col(width = 0.68, alpha = BAR_ALPHA, colour = "black", linewidth = 0.35) +
  ggplot2::geom_errorbar(
    ggplot2::aes(ymin = mean - sd, ymax = mean + sd),
    width = 0.18, linewidth = 0.55
  ) +
  ggplot2::geom_point(
    data = long,
    ggplot2::aes(x = group, y = value),
    inherit.aes = FALSE,
    position = ggplot2::position_jitter(width = POINT_JITTER_WIDTH, height = 0, seed = JITTER_SEED),
    colour = POINT_COLOR, size = POINT_SIZE, alpha = POINT_ALPHA
  ) +
  ggplot2::geom_text(
    data = tests,
    ggplot2::aes(x = group, y = y, label = annotation),
    inherit.aes = FALSE, vjust = 0, lineheight = 0.95, size = 3.5
  ) +
  ggplot2::scale_fill_manual(values = palette, drop = FALSE) +
  ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0.05, 0.08))) +
  ggplot2::coord_cartesian(
    ylim = c(data_min, max(tests$y) + 0.10 * data_range), clip = "off"
  ) +
  ggplot2::labs(
    title = PLOT_TITLE, x = X_LABEL, y = Y_LABEL
  ) +
  ggplot2::theme_classic(base_size = BASE_FONT_SIZE, base_family = BASE_FONT_FAMILY) +
  ggplot2::theme(
    legend.position = "none",
    axis.text.x = ggplot2::element_text(angle = 35, hjust = 1),
    plot.margin = ggplot2::margin(8, 14, 8, 8)
  )

if (ORIENTATION == "horizontal") {
  p <- p + ggplot2::coord_flip(
    ylim = c(data_min, max(tests$y) + 0.10 * data_range), clip = "off"
  ) + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 0, hjust = 0.5))
}

# 使用 R 自带的 Cairo SVG 图形设备。Cairo 会将文字保存为矢量字形路径，
# 而不是栅格图片；这样缩放 SVG 时文字和图形都不会失真。
if (!isTRUE(capabilities("cairo"))) {
  stop("当前 R 环境不支持 Cairo，无法保证 SVG 文字以矢量路径输出。")
}
grDevices::svg(
  filename = OUTPUT_FILE, width = FIG_WIDTH, height = FIG_HEIGHT,
  family = BASE_FONT_FAMILY, onefile = TRUE
)
print(p)
grDevices::dev.off()

# 在图片旁同时保存统计结果表，便于复核每个显著性符号。
stats_file <- sub("\\.svg$", "_statistics.csv", OUTPUT_FILE, ignore.case = TRUE)
if (identical(stats_file, OUTPUT_FILE)) stats_file <- paste0(OUTPUT_FILE, "_statistics.csv")
utils::write.csv(tests, stats_file, row.names = FALSE, fileEncoding = "UTF-8")

message("SVG saved to: ", normalizePath(OUTPUT_FILE, winslash = "/", mustWork = FALSE))
message("Statistics saved to: ", normalizePath(stats_file, winslash = "/", mustWork = FALSE))
