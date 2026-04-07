# ==========================================
# 单样本基因级别突变分布图 (自动定宽 & AI排版优化版)
# - 修复底层版本警告信息，控制台极致纯净
# - 深度解析 GFF3，精美呈现【外显子-内含子】箭头结构
# ==========================================

library(Gviz)
library(GenomicRanges)
library(Rsamtools)
library(rtracklayer) 
library(this.path)

options(ucscChromosomeNames = FALSE) 
setwd(this.path::this.dir())

# ==========================================
# 1. 参数配置区 
# ==========================================
sample_bam <- "CK-7-9_clean\\CK-7-9_clean.bam"             
ref_fasta  <- "CK-7-9_cleaned_v1.fasta" 
anno_file  <- "CK-7-9_cleaned_v1_annotated.gff3" 

target_gene <- "TGAM01_v208328"  
flank_size  <- 50              
inch_per_bp <- 0.02 

# ==========================================
# 2. 解析注释文件 (升级版：提取外显子与内含子)
# ==========================================
anno <- import(anno_file)

# 安全提取函数
get_attr <- function(x) {
  if (is.null(x)) return(rep("", length(anno)))
  vec <- as.character(x)
  vec[is.na(vec)] <- ""
  return(vec)
}

anno_type   <- get_attr(anno$type)
anno_name   <- get_attr(anno$Name)
anno_id     <- get_attr(anno$ID)
anno_parent <- get_attr(anno$Parent)

# 【核心优化】：不仅找基因本身，还要把 Parent 属于该基因的外显子 (exon/CDS) 全部抓出来
target_idx <- which(
  anno_name == target_gene | 
  anno_id == target_gene | 
  anno_id == paste0("gene-", target_gene) |
  grepl(target_gene, anno_parent)
)

gene_region_all <- anno[target_idx]

if (length(gene_region_all) == 0) {
  available_names <- unique(anno_name[anno_name != ""])
  message("当前注释文件中可用的 Name 示例：", paste(head(available_names, 5), collapse = ", "))
  stop("找不到该基因，请核对是否拼写错误！")
}

# 自动提取该基因的整体起止点（用于确定画布）
chr_name   <- as.character(seqnames(gene_region_all)[1])
gene_start <- min(start(gene_region_all))
gene_end   <- max(end(gene_region_all))

plot_start <- gene_start - flank_size
plot_end   <- gene_end + flank_size
plot_width_bases <- plot_end - plot_start + 1
canvas_width <- max(4, plot_width_bases * inch_per_bp) 

message(sprintf("🎯 目标基因: %s | 所在 Contig: %s\n📏 绘图范围: %d - %d | 序列总长: %d bp\n🎨 自动计算画布宽度为: %.2f 英寸", 
                target_gene, chr_name, plot_start, plot_end, plot_width_bases, canvas_width))

# ==========================================
# 3. 构建绘图轨道
# ==========================================

axisTrack <- GenomeAxisTrack(col = "black", fontcolor = "black", fontsize = 14, lwd = 1.2)

# 【核心优化】：分离外显子用于精细绘图
exons <- gene_region_all[get_attr(gene_region_all$type) %in% c("exon", "CDS")]

if (length(exons) > 0) {
  plot_gr <- exons
  track_group <- get_attr(exons$Parent)          # 按照转录本分组，Gviz 会自动画出内含子连线
  track_id    <- rep(target_gene, length(exons)) # 强制统一显示基因名
} else {
  # 兜底：如果此基因没有注释外显子，降级显示为大方块
  plot_gr <- gene_region_all[1]
  track_group <- target_gene
  track_id    <- target_gene
}

# 基因结构注释轨道 (优化为 IGV 风格：短箭头 + 粗内含子)
geneTrack <- AnnotationTrack(
  plot_gr, 
  name = "Gene",
  group = track_group,           # 激活分组连线特性 (内含子)
  id = track_id,
  groupAnnotation = "id",        # 将基因名显示在结构上
  just.group = "above",          # 名字悬浮在上方
  showId = TRUE,
  
  # === 形状与线条优化 (IGV Style) ===
  shape = "arrow",               
  arrowHeadWidth = 10,           # 【新增】控制箭头的相对长度 (默认通常是 30，改小让箭头变短)
  arrowHeadMaxWidth = 10,        # 【新增】强制限制箭头的最大绝对像素长度，防止长外显子箭头过长
  lwd = 1.5,                     # 【新增】线宽参数：加粗内含子连线（同时也会让外显子的黑色边框更清晰）
  
  chromosome = chr_name,
  fill = "#82b2d2",              # 高级灰蓝色
  col = "black",                 # 边框和内含子连线的颜色
  fontcolor.group = "black",     
  fontsize.group = 16,           
  background.title = "transparent",
  col.title = "black",
  cex.title = 0.9
)

seqTrack <- SequenceTrack(
  ref_fasta, chromosome = chr_name, name = "Ref",
  noLetters = FALSE, add53 = TRUE, add35 = TRUE,         
  cex = 1.2, fontface = "bold",    
  background.title = "transparent", col.title = "black"
)

sampleCovTrack <- AlignmentsTrack(
  sample_bam, isPaired = TRUE, chromosome = chr_name, referenceSequence = seqTrack, 
  type = "coverage", showMismatches = TRUE,     
  name = "Cov",
  background.title = "transparent", col.title = "black", col.axis = "black",
  cex.title = 0.9, cex.axis = 0.7             
)

sampleReadsTrack <- AlignmentsTrack(
  sample_bam, isPaired = TRUE, chromosome = chr_name, referenceSequence = seqTrack, 
  type = "pileup", showMismatches = TRUE,     
  fill.reads = "grey85", col.reads = "grey60", col.mismatch = "black",    
  name = "Reads",
  background.title = "transparent", col.title = "black",
  cex.title = 0.9            
)

# ==========================================
# 4. 渲染出图 (加入 suppressWarnings 屏蔽底层提示)
# ==========================================
output_name <- paste0(target_gene, "_Validation.svg")
svg(output_name, width = canvas_width, height = 5) 

# 强行屏蔽底层版本更迭带来的无关警告
suppressWarnings({
  plotTracks(
    trackList = list(axisTrack, geneTrack, sampleCovTrack, sampleReadsTrack, seqTrack),
    from = plot_start,
    to = plot_end,
    chromosome = chr_name,
    sizes = c(1, 0.8, 1.5, 3.5, 1), # 稍微加宽一点基因轨道的高度，让模型更饱满
    background.panel = "white",
    col.grid = "transparent",
    cex.mismatch = 1.0
  )
})

dev.off() 

message("✅ 极简拆分版绘图完成！已保存为: ", output_name)