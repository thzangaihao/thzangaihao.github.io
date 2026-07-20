#!/usr/bin/env python3
"""交互式绘制 .mcool 文件中的全基因组 Hi-C 热图。"""

# 此文件与文章目录中的脚本保持一致。

from datetime import datetime
from pathlib import Path

import cooler
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LinearSegmentedColormap


# =============================================================================
# 用户配置区
# =============================================================================
# MCOOL
#   在脚本所在目录及全部子目录中递归搜索的文件名模式。
#   找到多个文件时会列出文件并要求输入序号；只有一个时自动选择。
MCOOL = "*.mcool"

# RESOLUTION
#   从 .mcool 文件读取的分辨率，单位为 bp。
#   该分辨率必须已经存在于 .mcool 文件中。
RESOLUTION = 100_000

# OUTPUT_FILE
#   默认输出格式，支持 svg、png 和 pdf。
#   运行时可以交互输入输出文件名；直接回车则以当前时间命名。
#   输入文件名时若省略扩展名，会自动添加这里设置的默认扩展名。
OUTPUT_FILE = "svg"

# CHROMS
#   None：运行时列出文件中的染色体 / scaffold，并交互选择。
#   交互输入支持 1,2,3、1-3、1,2,4-6 和 all。
#   也可预设为 "all"、"1,2,4-6"，或染色体名列表
#   （例如 ["chr1", "chr2"]），从而跳过染色体选择步骤。
CHROMS = None

# BALANCED
#   True 使用平衡后的矩阵；如果文件中没有 weight，可改为 False。
BALANCED = True

# CLIP_QUANTILE
#   对矩阵执行 log1p 转换后，将高于该分位数的信号压到颜色上限，
#   避免少量极高值降低主体区域的颜色对比度。
CLIP_QUANTILE = 0.995

# COLORMAP
#   Matplotlib 内置色板，仅在 CUSTOM_COLORS = None 时生效。
COLORMAP = "Reds"

# CUSTOM_COLORS
#   自定义从低值到高值的渐变色。例如：
#   ["#FFFFFF", "#FFF3BF", "#FDAE6B", "#E6550D", "#7F0000"]
CUSTOM_COLORS = None

# LABEL_MIN_LENGTH
#   只有长度达到该值的染色体 / scaffold 才显示底部文字标签；
#   所有被选择的染色体仍会进入矩阵。设置为 0 可显示全部名称。
LABEL_MIN_LENGTH = 1_000_000

# DPI
#   PNG 等栅格输出的分辨率。SVG/PDF 主体仍为矢量容器。
DPI = 220

# MAX_BINS
#   稠密矩阵的安全上限。提高该值前请确认有足够内存。
MAX_BINS = 15_000

# 三角热图自身的画布尺寸，单位为英寸。三角形高度约为底边的一半，
# 因而使用比正方形热图更扁的画布。
FIGURE_SIZE = (14, 8)


def fail(message: str) -> None:
    raise SystemExit(f"错误：{message}")


def choose_mcool(script_dir: Path) -> Path:
    matches = sorted(
        (path for path in script_dir.rglob(MCOOL) if path.is_file()),
        key=lambda path: str(path.relative_to(script_dir)).lower(),
    )
    if not matches:
        fail(f"在 {script_dir} 及其子目录中找不到匹配 {MCOOL!r} 的文件")
    print("\n可用的 mcool 文件：")
    for index, path in enumerate(matches, 1):
        print(f"  [{index}] {path.relative_to(script_dir)}")
    if len(matches) == 1:
        print("仅找到一个文件，已自动选择。")
        return matches[0].resolve()
    while True:
        answer = input(f"请选择文件 [1-{len(matches)}，默认 1]：").strip() or "1"
        if answer.isdigit() and 1 <= int(answer) <= len(matches):
            return matches[int(answer) - 1].resolve()
        print("输入无效，请输入列表中的序号。")


def parse_chrom_selection(selection: str, available_chroms: list[str]) -> list[str]:
    selection = selection.strip().lower()
    if selection == "all":
        return list(available_chroms)
    if not selection:
        raise ValueError("选择不能为空")
    selected_indices: set[int] = set()
    for item in selection.replace("，", ",").split(","):
        item = item.strip()
        if not item:
            raise ValueError("存在空的选择项")
        if "-" in item:
            parts = item.split("-")
            if len(parts) != 2 or not all(part.strip().isdigit() for part in parts):
                raise ValueError(f"无法识别范围 {item!r}")
            start, end = (int(part.strip()) for part in parts)
            if start > end:
                raise ValueError(f"范围起点不能大于终点：{item}")
            selected_indices.update(range(start, end + 1))
        elif item.isdigit():
            selected_indices.add(int(item))
        else:
            raise ValueError(f"无法识别选择项 {item!r}")
    invalid = sorted(i for i in selected_indices if not 1 <= i <= len(available_chroms))
    if invalid:
        raise ValueError("序号超出范围：" + ", ".join(map(str, invalid)))
    return [chrom for i, chrom in enumerate(available_chroms, 1) if i in selected_indices]


def choose_chroms(available_chroms: list[str]) -> list[str]:
    print("\n可绘制的染色体 / scaffold：")
    for index, chrom in enumerate(available_chroms, 1):
        print(f"  [{index}] {chrom}")
    while True:
        answer = input("请选择（如 1,2,3 / 1-3 / 1,2,4-6 / all）：").strip()
        try:
            return parse_chrom_selection(answer, available_chroms)
        except ValueError as error:
            print(f"输入无效：{error}")


def choose_output(script_dir: Path) -> Path:
    default_extension = OUTPUT_FILE.lower().lstrip(".")
    if default_extension not in {"svg", "png", "pdf"}:
        fail("OUTPUT_FILE 必须是 svg、png 或 pdf")
    default_name = datetime.now().strftime(f"%Y%m%d_%H%M%S.{default_extension}")
    answer = input(f"\n输出文件名 [默认 {default_name}]：").strip()
    output = Path(answer or default_name).expanduser()
    if not output.suffix:
        output = output.with_suffix(f".{default_extension}")
    if output.suffix.lower() not in {".svg", ".png", ".pdf"}:
        fail("输出文件扩展名必须是 .svg、.png 或 .pdf")
    if not output.is_absolute():
        output = script_dir / output
    return output.resolve()


def main() -> None:
    script_dir = Path(__file__).resolve().parent
    mcool_path = choose_mcool(script_dir)
    if RESOLUTION <= 0:
        fail("RESOLUTION 必须大于 0")
    if not 0 < CLIP_QUANTILE <= 1:
        fail("CLIP_QUANTILE 必须位于 (0, 1] 区间")
    uri = f"{mcool_path}::/resolutions/{RESOLUTION}"
    try:
        contact = cooler.Cooler(uri)
    except Exception as error:
        fail(f"无法打开 {uri}\n{error}\n请运行 `cooler ls {mcool_path}` 检查可用分辨率。")
    chromsizes = contact.chromsizes
    available_chroms = list(chromsizes.index)
    if CHROMS is None:
        selected_chroms = choose_chroms(available_chroms)
    elif isinstance(CHROMS, str):
        try:
            selected_chroms = parse_chrom_selection(CHROMS, available_chroms)
        except ValueError as error:
            fail(f"CHROMS 配置无效：{error}")
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
        fail(f"矩阵包含 {n_bins:,} 个 bins，单个矩阵约需 {estimated_gib:.1f} GiB。\n请改用约 {suggested:,} bp 或更粗的已有分辨率，或确认内存充足后提高 MAX_BINS。")
    output = choose_output(script_dir)
    print(f"\n读取矩阵：{uri}")
    print(f"染色体 / scaffold 数：{len(selected_chroms):,}")
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
            fail("CUSTOM_COLORS 至少需要两种颜色")
        color_map = LinearSegmentedColormap.from_list("hic_custom", CUSTOM_COLORS, N=256)
    figure, axis = plt.subplots(figsize=FIGURE_SIZE)
    image = axis.imshow(matrix, cmap=color_map, vmin=0, vmax=color_max, origin="upper", interpolation="none", rasterized=True)
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
    axis.set_title(f"Hi-C contact map ({RESOLUTION:,} bp)")
    colorbar = figure.colorbar(image, ax=axis, fraction=0.046, pad=0.03)
    value_name = "balanced contact" if BALANCED else "raw contact count"
    colorbar.set_label(f"log1p({value_name}); upper {CLIP_QUANTILE * 100:g}% clipped")
    figure.tight_layout()
    output.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(output, dpi=DPI, bbox_inches="tight", facecolor="white")
    plt.close(figure)
    print(f"完成：{output}")


if __name__ == "__main__":
    main()
