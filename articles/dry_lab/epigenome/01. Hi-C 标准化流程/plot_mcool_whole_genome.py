#!/usr/bin/env python3
"""Plot all chromosomes/scaffolds from an .mcool file in one heatmap."""

from pathlib import Path

import cooler
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LinearSegmentedColormap


# =============================================================================
# 用户配置区：修改后直接运行 `python plot_mcool_whole_genome.py`
# =============================================================================

MCOOL = "./contacts.mcool"
RESOLUTION = 100_000

# 根据扩展名选择格式，支持 .svg、.png 和 .pdf。
OUTPUT_FILE = "contact_map.whole_genome.100kb.svg"

# None 表示包含 mcool 中的全部染色体/scaffold。
# 如果只想画指定染色体，可写成 ["chr1", "chr2", "chr3"]。
CHROMS = None

# TRUE 使用平衡后的矩阵；若文件没有 weight，可改为 False。
BALANCED = True

# 颜色使用 log1p 变换，并将高于此分位数的信号压到颜色上限。
CLIP_QUANTILE = 0.995

# Matplotlib 内置色板，例如 Reds、Blues、Purples、viridis、magma。
# 仅在 CUSTOM_COLORS = None 时使用。
COLORMAP = "viridis"

# 自定义由低到高的渐变色。设为 None 使用上面的 COLORMAP。
# 示例：白 → 浅黄 → 橙 → 深红
CUSTOM_COLORS = None
# CUSTOM_COLORS = ["#FFFFFF", "#FFF3BF", "#FDAE6B", "#E6550D", "#7F0000"]

# 所有染色体都会进入矩阵，但只有达到该长度的染色体显示文字标签，
# 避免大量短 scaffold 的名字重叠。设为 0 可显示全部名称。
LABEL_MIN_LENGTH = 1_000_000

FIGURE_SIZE = (14, 14)
DPI = 220

# 稠密矩阵保护。100 kb 的水稻基因组通常远低于该值。
MAX_BINS = 15_000

# =============================================================================
# 绘图区域：一般不需要修改
# =============================================================================


def fail(message: str) -> None:
    raise SystemExit(f"错误：{message}")


mcool_path = Path(MCOOL).expanduser().resolve()
if not mcool_path.is_file():
    fail(f"找不到 mcool 文件：{mcool_path}")
if RESOLUTION <= 0:
    fail("RESOLUTION 必须大于 0")
if not 0 < CLIP_QUANTILE <= 1:
    fail("CLIP_QUANTILE 必须位于 (0, 1] 区间")

uri = f"{mcool_path}::/resolutions/{RESOLUTION}"
try:
    contact = cooler.Cooler(uri)
except Exception as error:
    fail(
        f"无法打开 {uri}\n{error}\n"
        f"请先运行 `cooler ls {mcool_path}` 检查可用分辨率。"
    )

chromsizes = contact.chromsizes
available_chroms = list(chromsizes.index)
if CHROMS is None:
    selected_chroms = available_chroms
else:
    missing = [chrom for chrom in CHROMS if chrom not in chromsizes.index]
    if missing:
        fail("以下染色体不在 mcool 中：" + ", ".join(missing))
    selected_chroms = list(CHROMS)

bins = contact.bins()[:]
selected_mask = bins["chrom"].isin(selected_chroms).to_numpy()
selected_indices = np.flatnonzero(selected_mask)
n_bins = len(selected_indices)
if n_bins == 0:
    fail("选择结果中没有任何 bins")
if n_bins > MAX_BINS:
    estimated_gib = n_bins * n_bins * 8 / 1024**3
    suggested = int(np.ceil(RESOLUTION * n_bins / MAX_BINS / 10_000) * 10_000)
    fail(
        f"矩阵包含 {n_bins:,} 个 bins，单个矩阵约需 {estimated_gib:.1f} GiB。\n"
        f"请改用约 {suggested:,} bp 或更粗的已有分辨率，"
        "或确认内存充足后提高 MAX_BINS。"
    )

print(f"读取矩阵：{uri}")
print(f"染色体/scaffold 数：{len(selected_chroms):,}")
print(f"矩阵 bins：{n_bins:,} × {n_bins:,}")

try:
    full_matrix = contact.matrix(balance=BALANCED)[:]
except Exception as error:
    if BALANCED:
        fail(f"读取 balanced 矩阵失败：{error}\n可将 BALANCED 改为 False。")
    raise

matrix = full_matrix[np.ix_(selected_indices, selected_indices)]
del full_matrix
matrix = np.asarray(matrix, dtype=np.float64)
matrix[~np.isfinite(matrix)] = 0
matrix[matrix < 0] = 0
matrix = np.log1p(matrix)

positive = matrix[matrix > 0]
if positive.size == 0:
    fail("矩阵没有非零互作信号")
color_max = float(np.quantile(positive, CLIP_QUANTILE))
if not np.isfinite(color_max) or color_max <= 0:
    color_max = float(positive.max())

selected_bins = bins.iloc[selected_indices].reset_index(drop=True)
chrom_array = selected_bins["chrom"].astype(str).to_numpy()
change_points = np.flatnonzero(chrom_array[1:] != chrom_array[:-1]) + 1
starts = np.r_[0, change_points]
ends = np.r_[change_points, n_bins]
names = chrom_array[starts]
midpoints = (starts + ends) / 2
lengths = chromsizes.loc[names].to_numpy()
label_mask = lengths >= LABEL_MIN_LENGTH

if CUSTOM_COLORS is None:
    color_map = COLORMAP
else:
    if len(CUSTOM_COLORS) < 2:
        fail("CUSTOM_COLORS 至少需要两个颜色")
    color_map = LinearSegmentedColormap.from_list(
        "hic_custom", CUSTOM_COLORS, N=256
    )

figure, axis = plt.subplots(figsize=FIGURE_SIZE)
image = axis.imshow(
    matrix,
    cmap=color_map,
    vmin=0,
    vmax=color_max,
    origin="upper",
    interpolation="none",
    rasterized=True,
)

for boundary in ends[:-1]:
    axis.axvline(boundary - 0.5, color="#444444", linewidth=0.35, alpha=0.55)
    axis.axhline(boundary - 0.5, color="#444444", linewidth=0.35, alpha=0.55)

axis.set_xticks(midpoints[label_mask])
axis.set_xticklabels(names[label_mask], rotation=90, fontsize=7)
axis.set_yticks(midpoints[label_mask])
axis.set_yticklabels(names[label_mask], fontsize=7)
axis.tick_params(length=0)
axis.set_xlabel("Chromosome / scaffold")
axis.set_ylabel("Chromosome / scaffold")
axis.set_title(f"Whole-genome Hi-C contact map ({RESOLUTION:,} bp)")

colorbar = figure.colorbar(image, ax=axis, fraction=0.046, pad=0.03)
value_name = "balanced contact" if BALANCED else "raw contact count"
colorbar.set_label(f"log1p({value_name}); upper {CLIP_QUANTILE * 100:g}% clipped")

figure.tight_layout()
output = Path(OUTPUT_FILE).expanduser().resolve()
if output.suffix.lower() not in {".svg", ".png", ".pdf"}:
    fail("OUTPUT_FILE 扩展名必须是 .svg、.png 或 .pdf")
output.parent.mkdir(parents=True, exist_ok=True)
figure.savefig(output, dpi=DPI, bbox_inches="tight", facecolor="white")
plt.close(figure)

print(f"完成：{output}")
