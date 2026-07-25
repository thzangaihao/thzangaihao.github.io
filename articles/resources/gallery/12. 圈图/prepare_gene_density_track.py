#!/usr/bin/env python3
"""从 GFF/GFF3/GTF 注释生成 Circos 基因密度数据和独立轨道配置。

本脚本不会修改任何已有 Circos 文件，只生成：

1. <prefix>_gene_density.txt：Circos 四列数值轨道；
2. <prefix>_gene_density_plot.conf：只含一个 <plot> 块的轨道配置。

GFF/GTF 坐标是 1-based 闭区间，本脚本会转换为 Circos 常用的
0-based 半开区间。脚本只使用 Python 标准库。
"""

from __future__ import annotations

import argparse
import gzip
import math
import re
import sys
from collections import defaultdict
from pathlib import Path
from typing import Iterator, TextIO


def positive_integer(text: str) -> int:
    """确保窗口和步长为大于 0 的整数。"""
    value = int(text)
    if value <= 0:
        raise argparse.ArgumentTypeError("必须是大于 0 的整数")
    return value


def open_text(path: Path) -> TextIO:
    """打开普通文本文件或扩展名为 .gz 的 gzip 压缩文件。"""
    if path.suffix.lower() == ".gz":
        return gzip.open(path, "rt", encoding="utf-8")
    return path.open("rt", encoding="utf-8")


def percentile(values: list[float], percentage: float) -> float:
    """使用线性插值计算百分位数，不依赖 NumPy。"""
    ordered = sorted(values)
    if len(ordered) == 1:
        return ordered[0]
    position = (len(ordered) - 1) * percentage / 100.0
    lower = math.floor(position)
    upper = math.ceil(position)
    if lower == upper:
        return ordered[lower]
    fraction = position - lower
    return ordered[lower] * (1.0 - fraction) + ordered[upper] * fraction


def safe_prefix(path: Path) -> str:
    """从注释文件名生成安全的输出前缀。"""
    name = path.name
    if name.lower().endswith(".gz"):
        name = name[:-3]
    for suffix in (".gff3", ".gff", ".gtf"):
        if name.lower().endswith(suffix):
            name = name[: -len(suffix)]
            break
    cleaned = re.sub(r"[^A-Za-z0-9_.-]+", "_", name).strip("._")
    return cleaned or "annotation"


def parse_feature_types(text: str) -> set[str]:
    """解析以逗号分隔的 feature type，比较时忽略大小写。"""
    feature_types = {item.strip().lower() for item in text.split(",") if item.strip()}
    if not feature_types:
        raise argparse.ArgumentTypeError("至少需要一个 feature type")
    return feature_types


def read_annotation(
    path: Path,
    feature_types: set[str],
) -> tuple[dict[str, list[tuple[int, int]]], dict[str, int], list[str]]:
    """读取 GFF/GFF3/GTF，返回区间、观察到的最大坐标和染色体顺序。"""
    intervals: dict[str, list[tuple[int, int]]] = defaultdict(list)
    maximum_ends: dict[str, int] = {}
    chromosome_order: list[str] = []
    chromosomes_seen: set[str] = set()

    with open_text(path) as handle:
        for line_number, raw_line in enumerate(handle, start=1):
            line = raw_line.rstrip("\r\n")
            if not line:
                continue

            # GFF3 的 ##sequence-region 可提供完整染色体长度。
            if line.startswith("##sequence-region"):
                fields = line.split()
                if len(fields) >= 4:
                    chromosome = fields[1]
                    try:
                        region_end = int(fields[3])
                    except ValueError:
                        raise ValueError(
                            f"第 {line_number} 行的 ##sequence-region 坐标无效"
                        )
                    maximum_ends[chromosome] = max(
                        maximum_ends.get(chromosome, 0), region_end
                    )
                    if chromosome not in chromosomes_seen:
                        chromosome_order.append(chromosome)
                        chromosomes_seen.add(chromosome)
                continue
            if line.startswith("#"):
                continue

            fields = line.split("\t")
            if len(fields) != 9:
                raise ValueError(
                    f"第 {line_number} 行不是标准九列 GFF/GTF：检测到 {len(fields)} 列"
                )
            chromosome, feature_type = fields[0], fields[2].lower()
            if feature_type not in feature_types:
                continue
            try:
                # GFF/GTF：start/end 均为 1-based 且包含端点。
                start_1based = int(fields[3])
                end_1based = int(fields[4])
            except ValueError:
                raise ValueError(f"第 {line_number} 行的起止坐标不是整数")
            if start_1based < 1 or end_1based < start_1based:
                raise ValueError(f"第 {line_number} 行的起止坐标无效")

            start_0based = start_1based - 1
            end_0based = end_1based
            intervals[chromosome].append((start_0based, end_0based))
            maximum_ends[chromosome] = max(
                maximum_ends.get(chromosome, 0), end_0based
            )
            if chromosome not in chromosomes_seen:
                chromosome_order.append(chromosome)
                chromosomes_seen.add(chromosome)

    return intervals, maximum_ends, chromosome_order


def read_chromosome_sizes(path: Path) -> tuple[dict[str, int], list[str]]:
    """读取两列 chrom.sizes 或 Circos karyotype.txt。

    支持格式：
    chromosome  length
    chr - ID LABEL START END COLOR
    """
    sizes: dict[str, int] = {}
    order: list[str] = []
    with open_text(path) as handle:
        for line_number, raw_line in enumerate(handle, start=1):
            line = raw_line.strip()
            if not line or line.startswith("#"):
                continue
            fields = line.split()
            try:
                if len(fields) >= 6 and fields[0] == "chr" and fields[1] == "-":
                    chromosome = fields[2]
                    length = int(fields[5]) - int(fields[4])
                elif len(fields) >= 2:
                    chromosome = fields[0]
                    length = int(fields[1])
                else:
                    raise ValueError
            except ValueError:
                raise ValueError(
                    f"染色体长度文件第 {line_number} 行格式无效：{line}"
                )
            if length <= 0:
                raise ValueError(
                    f"染色体长度文件第 {line_number} 行的长度必须大于 0"
                )
            if chromosome in sizes:
                raise ValueError(f"染色体长度文件存在重复名称：{chromosome}")
            sizes[chromosome] = length
            order.append(chromosome)
    return sizes, order


def merge_intervals(intervals: list[tuple[int, int]]) -> list[tuple[int, int]]:
    """合并重叠区间，防止重叠基因或转录本被重复计算覆盖碱基数。"""
    if not intervals:
        return []
    ordered = sorted(intervals)
    merged: list[tuple[int, int]] = [ordered[0]]
    for start, end in ordered[1:]:
        previous_start, previous_end = merged[-1]
        if start <= previous_end:
            merged[-1] = (previous_start, max(previous_end, end))
        else:
            merged.append((start, end))
    return merged


def calculate_count_density(
    intervals: list[tuple[int, int]],
    chromosome_length: int,
    window_size: int,
    step_size: int,
    count_method: str,
) -> list[tuple[int, int, float]]:
    """计算窗口基因数；midpoint 保证每个基因只计入一个窗口。"""
    windows = [
        [start, min(start + window_size, chromosome_length), 0.0]
        for start in range(0, chromosome_length, step_size)
    ]
    if count_method == "midpoint":
        for feature_start, feature_end in intervals:
            midpoint = (feature_start + feature_end - 1) // 2
            # 对非重叠窗口，这是 midpoint 所属窗口的直接索引。
            if step_size == window_size:
                index = midpoint // step_size
                if 0 <= index < len(windows):
                    windows[index][2] += 1.0
            else:
                # 重叠或有间隔时，明确测试每个候选窗口。
                first = max(0, (midpoint - window_size + 1) // step_size)
                last = min(len(windows) - 1, midpoint // step_size)
                for index in range(first, last + 1):
                    if windows[index][0] <= midpoint < windows[index][1]:
                        windows[index][2] += 1.0
    else:
        # overlap 会让跨越窗口边界的同一基因在多个窗口中各计一次。
        for feature_start, feature_end in intervals:
            first = max(0, feature_start // step_size)
            last = min(len(windows) - 1, (feature_end - 1) // step_size)
            # 当 window_size != step_size 时扩大候选范围，再检查真实重叠。
            if window_size > step_size:
                first = max(0, (feature_start - window_size + 1) // step_size)
            for index in range(first, last + 1):
                window_start, window_end = windows[index][0], windows[index][1]
                if feature_start < window_end and feature_end > window_start:
                    windows[index][2] += 1.0
    return [(int(start), int(end), value) for start, end, value in windows]


def calculate_coverage_density(
    intervals: list[tuple[int, int]],
    chromosome_length: int,
    window_size: int,
    step_size: int,
    metric: str,
) -> list[tuple[int, int, float]]:
    """计算去重后的基因覆盖碱基数或窗口覆盖比例。"""
    merged = merge_intervals(intervals)
    result: list[tuple[int, int, float]] = []
    interval_index = 0

    for window_start in range(0, chromosome_length, step_size):
        window_end = min(window_start + window_size, chromosome_length)
        while interval_index < len(merged) and merged[interval_index][1] <= window_start:
            interval_index += 1
        covered_bases = 0
        scan_index = interval_index
        while scan_index < len(merged) and merged[scan_index][0] < window_end:
            feature_start, feature_end = merged[scan_index]
            covered_bases += max(
                0, min(feature_end, window_end) - max(feature_start, window_start)
            )
            scan_index += 1
        value = float(covered_bases)
        if metric == "fraction":
            value /= window_end - window_start
        result.append((window_start, window_end, value))
    return result


def make_track_config(
    data_reference: str,
    metric: str,
    feature_types: set[str],
    plot_type: str,
    r0: str,
    r1: str,
    color: str,
    display_min: float,
    display_max: float,
    window_size: int,
    step_size: int,
) -> str:
    """创建带详细注释的独立 Circos 轨道配置。"""
    fill_parameters = ""
    if plot_type == "histogram":
        fill_parameters = f"""# 柱形填充色，名称需要在 colors.conf 中定义。
fill_color = {color}
# 关闭每个柱子的边框，避免产生密集黑线。
thickness = 0p
extend_bin = no
"""
    else:
        fill_parameters = """# 折线宽度，p 表示像素。
thickness = 2p
"""

    feature_description = ",".join(sorted(feature_types))
    metric_description = {
        "count": "每个窗口中的 feature 数量",
        "bases": "每个窗口被 feature 覆盖的去重碱基数",
        "fraction": "每个窗口被 feature 覆盖的比例（0–1）",
    }[metric]
    legend_symbol = "line" if plot_type == "line" else "box"
    return f"""# GFF/GFF3/GTF 基因密度轨道配置
# 统计 feature type：{feature_description}
# 密度含义：{metric_description}
# 窗口大小：{window_size} bp；步长：{step_size} bp
# 本文件只包含一个 <plot> 块，应放进 config/tracks.conf 的 <plots> 内。
# legend_label = Gene density
# legend_symbol = {legend_symbol}
<plot>
# histogram 为柱形轨道；line 为折线轨道。
type = {plot_type}
# 路径相对于运行 circos 命令时的当前项目目录。
file = {data_reference}
# r0/r1 分别控制轨道内半径和外半径。
r0 = {r0}
r1 = {r1}
# 数据显示范围；自动上限默认使用数据的第 99 百分位数。
min = {display_min:.8g}
max = {display_max:.8g}
# 颜色名称需要在 config/colors.conf 中存在。
color = {color}
{fill_parameters}
<axes>
<axis>
# 每隔轨道高度的 20% 绘制一条浅灰色辅助线。
spacing = 0.2r
color = lgrey
thickness = 1p
</axis>
</axes>
</plot>
"""


def make_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="从 GFF/GFF3/GTF 生成 Circos 基因密度数据与轨道配置。",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("annotation", type=Path, help="输入 GFF、GFF3、GTF 或 .gz 文件")
    parser.add_argument(
        "-o", "--output-dir", type=Path, default=Path("circos_gene_density"),
        help="数据文件和配置文件的输出目录",
    )
    parser.add_argument("--prefix", help="输出文件前缀；默认从注释文件名生成")
    parser.add_argument(
        "-w", "--window-size", type=positive_integer, default=100_000,
        help="密度统计窗口长度，单位 bp",
    )
    parser.add_argument(
        "-s", "--step-size", type=positive_integer,
        help="相邻窗口起点距离；默认等于窗口长度",
    )
    parser.add_argument(
        "--feature-type", type=parse_feature_types, default={"gene"},
        help="要统计的第三列类型，可用逗号指定多个，例如 gene,pseudogene",
    )
    parser.add_argument(
        "--metric", choices=("count", "bases", "fraction"), default="count",
        help="输出基因数、去重覆盖碱基数或窗口覆盖比例",
    )
    parser.add_argument(
        "--count-method", choices=("midpoint", "overlap"), default="midpoint",
        help="count 模式按基因中点归入窗口，或每个重叠窗口均计数",
    )
    parser.add_argument(
        "--chrom-sizes", type=Path,
        help="可选的两列 chrom.sizes 或 karyotype.txt，用于获得完整染色体长度",
    )
    parser.add_argument("--include", metavar="REGEX", help="只处理名称匹配的染色体")
    parser.add_argument("--exclude", metavar="REGEX", help="排除名称匹配的染色体")
    parser.add_argument(
        "--missing-chromosome", choices=("error", "skip"), default="error",
        help="注释染色体不在 --chrom-sizes 中时终止或跳过",
    )
    parser.add_argument(
        "--plot-type", choices=("histogram", "line"), default="histogram",
        help="配置文件中的 Circos 轨道类型",
    )
    parser.add_argument("--r0", default="0.66r", help="轨道内半径")
    parser.add_argument("--r1", default="0.78r", help="轨道外半径")
    parser.add_argument(
        "--color", default="circle_color_04",
        help="轨道颜色名称，需要在 colors.conf 中定义",
    )
    parser.add_argument("--display-min", type=float, help="手动指定配置显示下限")
    parser.add_argument("--display-max", type=float, help="手动指定配置显示上限")
    parser.add_argument(
        "--max-percentile", type=float, default=99.0,
        help="自动显示上限使用的百分位数",
    )
    parser.add_argument(
        "--config-data-path",
        help="配置 file 参数使用的路径；默认 data/<数据文件名>",
    )
    parser.add_argument(
        "--omit-zero-windows", action="store_true",
        help="不输出密度为 0 的窗口",
    )
    parser.add_argument("--overwrite", action="store_true", help="覆盖已有同名输出")
    return parser


def main() -> int:
    args = make_parser().parse_args()
    annotation_path = args.annotation.resolve()
    output_dir = args.output_dir.resolve()
    step_size = args.step_size or args.window_size

    if not annotation_path.is_file():
        print(f"错误：找不到注释文件：{annotation_path}", file=sys.stderr)
        return 2
    if args.chrom_sizes and not args.chrom_sizes.is_file():
        print(f"错误：找不到染色体长度文件：{args.chrom_sizes}", file=sys.stderr)
        return 2
    if not 0 < args.max_percentile <= 100:
        print("错误：--max-percentile 必须大于 0 且不超过 100。", file=sys.stderr)
        return 2
    if (
        args.display_min is not None
        and args.display_max is not None
        and args.display_min >= args.display_max
    ):
        print("错误：--display-min 必须小于 --display-max。", file=sys.stderr)
        return 2

    try:
        include_pattern = re.compile(args.include) if args.include else None
        exclude_pattern = re.compile(args.exclude) if args.exclude else None
    except re.error as error:
        print(f"错误：正则表达式无效：{error}", file=sys.stderr)
        return 2

    prefix = args.prefix or safe_prefix(annotation_path)
    if not re.fullmatch(r"[A-Za-z0-9_.-]+", prefix):
        print("错误：--prefix 只能包含字母、数字、点、下划线和连字符。", file=sys.stderr)
        return 2
    data_path = output_dir / f"{prefix}_gene_density.txt"
    config_path = output_dir / f"{prefix}_gene_density_plot.conf"
    existing = [path for path in (data_path, config_path) if path.exists()]
    if existing and not args.overwrite:
        names = "、".join(path.name for path in existing)
        print(f"错误：{names} 已存在；如需覆盖请添加 --overwrite。", file=sys.stderr)
        return 2

    try:
        intervals, observed_sizes, annotation_order = read_annotation(
            annotation_path, args.feature_type
        )
        if args.chrom_sizes:
            chromosome_sizes, chromosome_order = read_chromosome_sizes(
                args.chrom_sizes.resolve()
            )
            absent = sorted(set(intervals) - set(chromosome_sizes))
            if absent and args.missing_chromosome == "error":
                preview = ", ".join(absent[:10])
                raise ValueError(
                    f"{len(absent)} 条注释染色体不在长度文件中：{preview}"
                )
        else:
            chromosome_sizes = observed_sizes
            chromosome_order = annotation_order
    except (OSError, UnicodeError, ValueError) as error:
        print(f"错误：读取输入失败：{error}", file=sys.stderr)
        return 1

    feature_count = sum(len(items) for items in intervals.values())
    if feature_count == 0:
        requested = ",".join(sorted(args.feature_type))
        print(f"错误：没有找到 feature type：{requested}", file=sys.stderr)
        return 1

    output_dir.mkdir(parents=True, exist_ok=True)
    values_written: list[float] = []
    windows_written = 0
    chromosomes_written = 0

    try:
        with data_path.open("w", encoding="utf-8", newline="\n") as output:
            for chromosome in chromosome_order:
                if chromosome not in chromosome_sizes:
                    continue
                if include_pattern and not include_pattern.search(chromosome):
                    continue
                if exclude_pattern and exclude_pattern.search(chromosome):
                    continue
                chromosome_length = chromosome_sizes[chromosome]
                chromosome_intervals = intervals.get(chromosome, [])

                # 防止错误注释坐标超过用户提供的染色体长度。
                if chromosome_intervals:
                    largest_end = max(end for _, end in chromosome_intervals)
                    if largest_end > chromosome_length:
                        raise ValueError(
                            f"{chromosome} 的注释终点 {largest_end} "
                            f"超过染色体长度 {chromosome_length}"
                        )

                if args.metric == "count":
                    windows = calculate_count_density(
                        chromosome_intervals,
                        chromosome_length,
                        args.window_size,
                        step_size,
                        args.count_method,
                    )
                else:
                    windows = calculate_coverage_density(
                        chromosome_intervals,
                        chromosome_length,
                        args.window_size,
                        step_size,
                        args.metric,
                    )

                chromosomes_written += 1
                for start, end, value in windows:
                    if value == 0 and args.omit_zero_windows:
                        continue
                    value_text = (
                        str(int(value)) if args.metric in {"count", "bases"}
                        else f"{value:.8g}"
                    )
                    output.write(f"{chromosome}\t{start}\t{end}\t{value_text}\n")
                    values_written.append(value)
                    windows_written += 1
    except (OSError, ValueError) as error:
        data_path.unlink(missing_ok=True)
        config_path.unlink(missing_ok=True)
        print(f"错误：生成密度失败：{error}", file=sys.stderr)
        return 1

    if not values_written:
        data_path.unlink(missing_ok=True)
        print("错误：筛选后没有可输出的窗口。", file=sys.stderr)
        return 1

    display_min = args.display_min if args.display_min is not None else 0.0
    display_max = (
        args.display_max
        if args.display_max is not None
        else percentile(values_written, args.max_percentile)
    )
    if display_max <= display_min:
        display_max = display_min + 1.0

    data_reference = (
        args.config_data_path
        if args.config_data_path
        else f"data/{data_path.name}"
    )
    config_text = make_track_config(
        data_reference=data_reference,
        metric=args.metric,
        feature_types=args.feature_type,
        plot_type=args.plot_type,
        r0=args.r0,
        r1=args.r1,
        color=args.color,
        display_min=display_min,
        display_max=display_max,
        window_size=args.window_size,
        step_size=step_size,
    )
    try:
        config_path.write_text(config_text, encoding="utf-8")
    except OSError as error:
        data_path.unlink(missing_ok=True)
        print(f"错误：无法写入配置文件：{error}", file=sys.stderr)
        return 1

    print(f"完成：读取 {feature_count} 个 feature")
    print(f"输出：{chromosomes_written} 条染色体、{windows_written} 个窗口")
    print(f"数据文件：{data_path}")
    print(f"配置文件：{config_path}")
    print(f"配置显示范围：{display_min:.8g} 至 {display_max:.8g}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
