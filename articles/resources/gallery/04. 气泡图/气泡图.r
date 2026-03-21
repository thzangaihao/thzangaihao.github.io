# ==========================================
# 🎨 发表级 GO 富集气泡图 (完美右对齐版)
# ==========================================

# 1. 载入必须的工具包
# 如果报错，请运行: install.packages(c("ggplot2", "dplyr", "forcats", "this.path", "grid"))
library(ggplot2)
library(dplyr)
library(forcats)
library(this.path)  # 解决路径自动化
library(grid)       # ✨ 替换了 egg，用于底层精细提取宽度，同时保护完美的气泡纵向间距

# 2. 路径自动化：自动识别脚本所在位置
tryCatch({
  setwd(this.path::this.dir())
}, error = function(e) {
  message("⚠️ 无法自动获取路径，请确保在 VS Code/RStudio 中设置了正确的运行目录")
})

# 3. 读取数据
file_name <- "气泡图示例数据.csv"
if (!file.exists(file_name)) {
  stop(paste("❌ 找不到文件:", file_name, "。请检查文件是否在脚本同级目录下。"))
}
df <- read.csv(file_name)

# 4. 数据预处理
plot_data <- df %>%
  group_by(Category) %>%
  arrange(P_value) %>%
  slice_head(n = 10) %>%
  ungroup() %>%
  mutate(Description = fct_reorder(Description, RichFactor))

# 5. 核心绘图逻辑
p <- ggplot(plot_data, aes(x = RichFactor, y = Description)) +
  theme_bw() +
  geom_point(aes(size = Count, color = P_value), alpha = 0.8, stroke = 0.5) +
  
  # facet_grid 保留自由分配高度的特性
  facet_grid(Category ~ ., scales = "free_y", space = "free_y") +
  
  scale_color_viridis_c(direction = -1, limits = c(0, 0.05), oob = scales::squish) +
  
  # 💡 强力建议：如果你要拼贴多张不同的 GO 图，最好把下面的 max(df$Count) 替换成一个全局的最大固定值（比如 100）
  # 这样不仅可以保证气泡大小绝对一致，也能保证右侧图例的物理宽度“一毫米不差”。
  scale_size_continuous(range = c(3, 10), limits = c(0, max(df$Count))) +
  
  theme(
    strip.background = element_rect(fill = "#E0E0E0", color = "black", linewidth = 0.5),
    strip.text = element_text(face = "bold", size = 12, color = "black"),
    axis.text.y = element_text(size = 10, color = "black"),
    axis.text.x = element_text(size = 10, color = "black"),
    axis.title = element_text(face = "bold", size = 12),
    legend.title = element_text(face = "bold", size = 10),
    legend.key.height = unit(0.8, "cm"),
    panel.grid.minor = element_blank(),
    
    # ✨ 关键：消除左右两侧多余的白边，让绘图元素绝对贴着图片边缘
    plot.margin = margin(t = 10, r = 5, b = 10, l = 5, unit = "pt")
  ) +
  labs(x = "Rich Factor", y = NULL, color = "P-value", size = "Gene Count")

# 6. ✨ 核心魔法：仅强行锁定“绘图区”的物理宽度，保护比例高度
# 将 ggplot 转化为底层的图形表格对象 (gtable)
g <- ggplotGrob(p)

# 找到其中对应气泡图面板 (panel) 的列，强制将其物理宽度设定为精确的 4 英寸
panel_cols <- unique(g$layout$l[grepl("panel", g$layout$name)])
g$widths[panel_cols] <- unit(4, "in")

# 7. ✨ 核心魔法：“自适应无缝紧贴”算法
# 动态计算经过宽度修改后，整个图像的【精确绝对总宽度】（左侧长文字宽度 + 4英寸面板 + 右侧图例）
exact_width <- convertWidth(sum(g$widths), "in", valueOnly = TRUE)

# 保存图像：不再使用固定的 width = 12 瞎猜了！
# 使用精确的 exact_width，生成的图片没有任何多余留白，左侧正好贴着最长的字，右侧正好贴着图例！
output_file <- "GO_Bubble_Refined.png"
ggsave(output_file, plot = g, 
       width = exact_width, height = 10, units = "in", dpi = 300, bg = "white")

print(sprintf("✅ 绘图修缮完成！图片精确宽度为 %.2f 英寸，已保存至: %s/%s", exact_width, getwd(), output_file))