# ==========================================
# 📊 GO 功能分类统计图 (发表级右对齐条形图)
# ==========================================

# 1. 载入工具包
library(ggplot2)
library(dplyr)
library(forcats)
library(this.path)
library(grid)

# 2. 路径处理
tryCatch({
  setwd(this.path::this.dir())
}, error = function(e) {
  message("⚠️ 路径自动设置失败，请手动确认工作目录")
})

# 3. 读取数据 (使用你上传的“小库GO结构.csv”)
file_name <- "你的数据路径"
if (!file.exists(file_name)) {
  stop(paste("❌ 找不到文件:", file_name))
}
df <- read.csv(file_name)

# 4. 数据预处理
# 统计每个 Category 下的条目，并按 Count 排序
plot_data <- df %>%
  group_by(Category) %>%
  arrange(desc(Count)) %>%
  slice_head(n = 10) %>% # 每个分类选前10个展示，若想全显可删掉此行
  ungroup() %>%
  # 关键：按 Category 排序，再按 Count 排序 Description
  mutate(Description = fct_reorder(Description, Count))

# 5. 绘图逻辑
p <- ggplot(plot_data, aes(x = Count, y = Description, fill = Category)) +
  theme_bw() +
  
  # 绘制条形图
  geom_bar(stat = "identity", width = 0.7, alpha = 0.8, color = "black", linewidth = 0.3) +
  
  # 分栏：BP/CC/MF 垂直排列
  facet_grid(Category ~ ., scales = "free_y", space = "free_y") +
  
  # 设置经典科研配色
  scale_fill_manual(values = c("BP" = "#89CFF0", "CC" = "#A1D99B", "MF" = "#FDB462")) +
  
  # 精细主题调节
  theme(
    # 分类标签背景
    strip.background = element_rect(fill = "#EEEEEE", color = "black"),
    strip.text = element_text(face = "bold", size = 11),
    
    # 坐标轴文字
    axis.text.y = element_text(size = 9, color = "black"),
    axis.text.x = element_text(size = 10, color = "black"),
    axis.title = element_text(face = "bold", size = 12),
    
    # 移除网格线
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    
    # 图例位置（放在底部或隐藏，条形图通常不需要额外图例因为已有分栏标签）
    legend.position = "none",
    
    # 边缘紧凑处理
    plot.margin = margin(t = 10, r = 10, b = 10, l = 5, unit = "pt")
  ) +
  
  labs(x = "Number of Genes", y = "GO Term Description", 
       title = "GO Functional Classification")

# 6. ✨ 核心魔法：锁定物理宽度并计算自适应总宽
g <- ggplotGrob(p)

# 找到绘图面板列，强制锁定为 4 英寸宽
panel_cols <- unique(g$layout$l[grepl("panel", g$layout$name)])
g$widths[panel_cols] <- unit(4, "in")

# 自动计算能够容纳左侧文字的最窄总宽度
exact_width <- convertWidth(sum(g$widths), "in", valueOnly = TRUE)

# 7. 保存图像
output_file <- "GO_Classification_Bar.png"
ggsave(output_file, plot = g, 
       width = exact_width, height = 8, units = "in", dpi = 300, bg = "white")

print(sprintf("✅ 统计图绘制完成！图片总宽度为 %.2f 英寸", exact_width))