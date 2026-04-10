library(ggplot2)
library(dplyr)
library(tidyr)
library(this.path)
library(svglite)
library(ggh4x) 

# 自动设置工作目录为脚本所在位置
setwd(this.path::this.dir())

# ==========================================
# 1. 序列与位置参数输入区
# ==========================================
# 定义局部核苷酸序列的起始位置
start_pos <- 1147

# 手动定义氨基酸序列的起始位置（避免内含子打断导致计算错误）
aa_start_pos <- 342

# ---> 在此统一控制所有方框的边框粗细 <---
border_thickness <- 0.5

# 输入 Wild Type (野生型) 和 Mutant (突变型) 的序列片段
wt_dna  <- c("G", "A", "T", "C", "C", "C", "T", "T", "C")
mut_dna <- c("G", "A", "T", "T", "C", "C", "T", "T", "C") 

wt_pro  <- c("Asp", "Pro", "Phe")
mut_pro <- c("Asp", "Ser", "Phe") 

# ==========================================
# 2. 数据处理与坐标映射
# ==========================================
# 生成真实的碱基坐标
dna_positions <- start_pos:(start_pos + length(wt_dna) - 1)

# DNA 数据框
dna_df <- data.frame(
  Position = dna_positions,
  WT = wt_dna,
  Mutant = mut_dna
) %>%
  pivot_longer(cols = c(WT, Mutant), names_to = "Type", values_to = "Base") %>%
  group_by(Position) %>%
  mutate(Is_Mutated = ifelse(length(unique(Base)) > 1, "Yes", "No")) %>%
  ungroup() %>%
  mutate(Molecule = "DNA", X_pos = Position)

# Protein 数据框 (X坐标对齐到每个密码子 3 个碱基的中心)
pro_positions <- (1:length(wt_pro)) * 3 - 2 + start_pos

# 根据手动输入的 aa_start_pos 生成氨基酸编号
aa_numbers <- aa_start_pos:(aa_start_pos + length(wt_pro) - 1)

pro_df <- data.frame(
  Position = pro_positions,
  WT = wt_pro,
  Mutant = mut_pro
) %>%
  pivot_longer(cols = c(WT, Mutant), names_to = "Type", values_to = "AA") %>%
  group_by(Position) %>%
  mutate(Is_Mutated = ifelse(length(unique(AA)) > 1, "Yes", "No")) %>%
  ungroup() %>%
  mutate(Molecule = "Protein", X_pos = Position) 

# 控制上下排显示的 Y 轴高度 (WT 在上，Mutant 在下)
y_mapping <- c("WT" = 2.2, "Mutant" = 1)
dna_df$Y_pos <- y_mapping[dna_df$Type]
pro_df$Y_pos <- y_mapping[pro_df$Type]

# ==========================================
# 3. 开始绘图
# ==========================================
# 定义极值，用于放置端点标签
min_x <- min(dna_df$X_pos)
max_x <- max(dna_df$X_pos)

p <- ggplot() +
  
  # ---------- DNA 层 ----------
  # 移除了 linewidth 的 aes 映射，直接使用全局 border_thickness
  geom_tile(data = dna_df, 
            aes(x = X_pos, y = Y_pos, fill = Base, color = Is_Mutated), 
            linewidth = border_thickness, width = 0.9, height = 0.45) +
  geom_text(data = dna_df, aes(x = X_pos, y = Y_pos, label = Base), 
            size = 6, fontface = "bold", color = "black") +
  
  # ---------- Protein 层 ----------
  # 移除了 linewidth 的 aes 映射，直接使用全局 border_thickness
  geom_tile(data = pro_df, 
            aes(x = X_pos, y = Y_pos + 0.52, color = Is_Mutated), 
            linewidth = border_thickness, fill = "white", width = 2.8, height = 0.45) +
  geom_text(data = pro_df, aes(x = X_pos, y = Y_pos + 0.52, label = AA), 
            size = 5.5, fontface = "bold") +
  
  # ---------- 5'/3' (DNA端点) 和 N'/C' (蛋白端点) ----------
  annotate("text", x = min_x - 0.8, y = 2.2, label = "5'", fontface = "bold", size = 5) +
  annotate("text", x = max_x + 0.8, y = 2.2, label = "3'", fontface = "bold", size = 5) +
  annotate("text", x = min_x - 0.8, y = 2.72, label = "N'", fontface = "bold", size = 5) +
  annotate("text", x = max_x + 0.8, y = 2.72, label = "C'", fontface = "bold", size = 5) +
  
  annotate("text", x = min_x - 0.8, y = 1, label = "5'", fontface = "bold", size = 5) +
  annotate("text", x = max_x + 0.8, y = 1, label = "3'", fontface = "bold", size = 5) +
  annotate("text", x = min_x - 0.8, y = 1.52, label = "N'", fontface = "bold", size = 5) +
  annotate("text", x = max_x + 0.8, y = 1.52, label = "C'", fontface = "bold", size = 5) +
  
  # ---------- 映射与配色 ----------
  scale_fill_manual(values = c("G" = "#9f3f33", "A" = "#b0c9d5", "T" = "#4a6d97", "C" = "#d5a773")) +
  # 突变位点边框变红，未突变保持黑色
  scale_color_manual(values = c("No" = "black", "Yes" = "#ff0059")) +
  # 注意：由于使用了全局 linewidth，移除了原有的 scale_linewidth_manual
  
  scale_y_continuous(breaks = c(1, 2.2), labels = c("Mutant", "Wild Type"), limits = c(0.5, 3)) +
  
  # ---------- 双 X 轴设置 ----------
  scale_x_continuous(
    breaks = seq(min_x, max_x, by = 1), 
    name = "TGAM01_v205164",
    sec.axis = sec_axis(
      transform = ~ ., 
      breaks = pro_positions, 
      labels = aa_numbers
    )
  ) +
  
  # ---------- 极简发文主题 ----------
  theme_classic(base_size = 14) + 
  theme(
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.y = element_text(color = "black", face = "bold", size = 14, hjust = 1, margin = margin(r = 15)),
    axis.title.y = element_blank(),
    
    axis.line.x.bottom = element_line(linewidth = 1, colour = "black"),
    axis.ticks.x.bottom = element_line(linewidth = 1, colour = "black"),
    axis.text.x.bottom = element_text(color = "black", size = 11, face = "bold"),
    axis.title.x.bottom = element_text(color = "black", size = 14, face = "bold", margin = margin(t = 12)),
    
    axis.line.x.top = element_line(linewidth = 1, colour = "black"),
    axis.ticks.x.top = element_line(linewidth = 1, colour = "black"),
    axis.text.x.top = element_text(color = "black", size = 12, face = "bold"),
    axis.title.x.top = element_text(color = "black", size = 14, face = "bold", margin = margin(b = 12)),
    
    legend.position = "none",
    panel.grid = element_blank(),
    plot.margin = margin(t = 20, r = 20, b = 20, l = 20)
  ) +
  coord_cartesian(xlim = c(min_x - 1.2, max_x + 1.2)) +
  force_panelsizes(rows = unit(6, "cm"), cols = unit(15, "cm"))

# ==========================================
# 4. 保存高质量图片
# ==========================================
ggsave("TGAM01_v205164_UnifiedBorder.svg", plot = p, 
       width = 22, height = 9, units = "cm", 
       dpi = 300, device = "svg")

message("绘图完成！边框粗细已统一控制，保存为 TGAM01_v205164_UnifiedBorder.svg")