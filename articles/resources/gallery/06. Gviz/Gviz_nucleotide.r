# ==========================================
# 测序 Reads 双样本验证图 (极致纯净版)
# - Reads 与 Coverage 轨道彻底拆分，告别大灰块
# - 左侧标题与坐标轴字体大小独立调节
# - Reads 内部去除文字，仅保留纯净色块
# ==========================================

library(Gviz)
library(GenomicRanges)
library(Rsamtools)
library(this.path)

# 关闭 UCSC 命名限制
options(ucscChromosomeNames = FALSE) 

# 自动设置工作目录
setwd(this.path::this.dir())

# ==========================================
# 1. 参数配置区 
# ==========================================
wt_bam    <- "CK-7-9_clean\\CK-7-9_clean.bam"             
mut_bam   <- "UV2-7-9-8_clean\\UV2-7-9-8_clean.bam" 
ref_fasta <- "CK-7-9_cleaned_v1.fasta" 
chr_name  <- "NODE_63_length_216026_cov_23.186387" 

mut_pos <- 35204
window_size <- 20  
plot_start  <- mut_pos - window_size
plot_end    <- mut_pos + window_size

# ==========================================
# 2. 构建绘图轨道
# ==========================================

# 【轨道 1】：顶部的基因组坐标尺
axisTrack <- GenomeAxisTrack(col = "black", fontcolor = "black", fontsize = 14, lwd = 1.2)

# 【轨道 2】：底层的参考序列字母
seqTrack <- SequenceTrack(
  ref_fasta, 
  chromosome = chr_name, 
  name = "Ref",
  noLetters = FALSE,    
  add53 = TRUE,         
  add35 = TRUE,         
  cex = 1.2,            
  fontface = "bold",    
  background.title = "transparent", # 透明背景
  col.title = "black"
)

# ==========================================
# 核心拆分：WT 样本分为 Coverage 和 Reads 两个独立轨道
# ==========================================
# 【轨道 3】：WT 覆盖度柱状图 (Coverage)
wtCovTrack <- AlignmentsTrack(
  wt_bam, isPaired = TRUE, chromosome = chr_name, referenceSequence = seqTrack, 
  type = "coverage",         # 【优化】：仅显示覆盖度
  showMismatches = TRUE,     
  name = "Cov",          # 换行显示更美观
  background.title = "transparent", # 去除大灰块
  col.title = "black",
  col.axis = "black",
  cex.title = 0.9,           # 【独立调节】：左侧标题的字体大小
  cex.axis = 0.7             # 【独立调节】：Y轴数字刻度的字体大小
)

# 【轨道 4】：WT 具体比对序列 (Reads)
wtReadsTrack <- AlignmentsTrack(
  wt_bam, isPaired = TRUE, chromosome = chr_name, referenceSequence = seqTrack, 
  type = "pileup",           # 【优化】：仅显示 Reads
  showMismatches = TRUE,     
  fill.reads = "grey85", col.reads = "grey60", col.mismatch = "black",
  name = "Reads",
  background.title = "transparent", # 去除大灰块
  col.title = "black",
  cex.title = 0.9            # 【独立调节】：左侧标题的字体大小 (此处无Y轴)
)

# ==========================================
# 核心拆分：Mutant 样本分为 Coverage 和 Reads 两个独立轨道
# ==========================================
# 【轨道 5】：Mutant 覆盖度柱状图 (Coverage)
mutCovTrack <- AlignmentsTrack(
  mut_bam, isPaired = TRUE, chromosome = chr_name, referenceSequence = seqTrack, 
  type = "coverage",
  showMismatches = TRUE,     
  name = "Cov",
  background.title = "transparent",
  col.title = "black",
  col.axis = "black",
  cex.title = 0.9,           
  cex.axis = 0.7             
)

# 【轨道 6】：Mutant 具体比对序列 (Reads)
mutReadsTrack <- AlignmentsTrack(
  mut_bam, isPaired = TRUE, chromosome = chr_name, referenceSequence = seqTrack, 
  type = "pileup",
  showMismatches = TRUE,     
  fill.reads = "grey85", col.reads = "grey60", col.mismatch = "black",    
  name = "Reads",
  background.title = "transparent",
  col.title = "black",
  cex.title = 0.9            
)

# 【轨道 7】：中心突变高亮带 (现在需要包裹 4 个轨道)
htTrack <- HighlightTrack(
  trackList = list(wtCovTrack, wtReadsTrack, mutCovTrack, mutReadsTrack), 
  start = mut_pos,
  end = mut_pos,
  chromosome = chr_name,
  col = "#000000",            
  fill = "#000000",           
  alpha = 0.1,                
  inBackground = FALSE        
)

# ==========================================
# 3. 渲染出图
# ==========================================
svg("5164.svg", width = 8, height = 7)

plotTracks(
  trackList = list(axisTrack, htTrack, seqTrack),
  from = plot_start,
  to = plot_end,
  chromosome = chr_name,
  
  # 高度分配比例：坐标尺(1) : WTCov(1.5) : WTReads(3) : MutCov(1.5) : MutReads(3) : 底部序列(1)
  sizes = c(1, 1.5, 3, 1.5, 3, 1), 
  
  background.panel = "white",
  col.grid = "transparent",
  
  # 调节错配方块内部的字母显示
  cex.mismatch = 1.0
)

dev.off() 

message("极简拆分版绘图完成！已保存为 IGV_Style_Validation_Pure.svg")