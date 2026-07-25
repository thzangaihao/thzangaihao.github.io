#!/usr/bin/env python3
"""将 bigWig 信号按窗口汇总为 Circos 数据文件和独立轨道配置。

本脚本不会读取或修改 circos.conf、karyotype.txt、colors.conf，也不会修改
prepare_circos_gc.py 和 create_circos_project.py。它只生成：

1. <prefix>_coverage.txt：Circos 四列数值轨道；
2. <prefix>_coverage_plot.conf：只含一个 <plot> 块的轨道配置。

依赖：pyBigWig（安装命令：conda install -c bioconda pybigwig）。
"""

from __future__ import annotations

import argparse
import math
import re
import sys
from pathlib import Path
from typing import Any


def positive_integer(text: str) -> int:
    """确保窗口和步长为大于 0 的整数。"""
    value = int(text)
    if value <= 0:
        raise argparse.ArgumentTypeError("必须是大于 0 的整数")
    return value


def percentile(values: list[float], percentage: float) -> float:
    """使用线性插值计算百分位数，不依赖 NumPy。"""
    ordered = sorted(values)
    if len(ordered) == 1:
        return ordered[0]
    position = (len(ordered) - 1) * percentage / 100.0
    lower_index = math.floor(position)
    upper_index = math.ceil(position)
    if lower_index == upper_index:
        return ordered[lower_index]
    fraction = position - lower_index
    return (
        ordered[lower_index] * (1.0 - fraction)
        + ordered[upper_index] * fraction
    )


def safe_prefix(path: Path) -> str:
    """从文件名产生适合作为输出文件前缀的字符串。"""
    name = path.name
    for suffix in (".bigwig", ".bw2", ".bw"):
        if name.lower().endswith(suffix):
            name = name[: -len(suffix)]
            break
    cleaned = re.sub(r"[^A-Za-z0-9_.-]+", "_", name).strip("._")
    return cleaned or "bigwig"


def make_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="从 bigWig/.bw/.bw2 生成 Circos 覆盖度数据和轨道配置。",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("bigwig", type=Path, help="输入 bigWig、.bw 或 .bw2 文件")
    parser.add_argument(
        "-o", "--output-dir", type=Path, default=Path("circos_bigwig"),
        help="数据文件和配置文件的输出目录",
    )
    parser.add_argument(
        "--prefix",
        help="输出文件前缀；默认由 bigWig 文件名自动生成",
    )
    parser.add_argument(
        "-w", "--window-size", type=positive_integer, default=100_000,
        help="每个统计窗口的长度，单位 bp",
    )
    parser.add_argument(
        "-s", "--step-size", type=positive_integer,
        help="相邻窗口起点的距离；默认等于窗口长度",
    )
    parser.add_argument(
        "--statistic", choices=("mean", "max", "min", "sum", "std"), default="mean",
        help="每个窗口采用的 bigWig 汇总方法",
    )
    parser.add_argument(
        "--missing", choices=("zero", "skip"), default="zero",
        help="无信号窗口写为 0，或者完全跳过",
    )
    parser.add_argument(
        "--include", metavar="REGEX",
        help="只处理染色体名称匹配该正则表达式的记录",
    )
    parser.add_argument(
        "--exclude", metavar="REGEX",
        help="排除染色体名称匹配该正则表达式的记录",
    )
    parser.add_argument(
        "--exact", action="store_true",
        help="要求 pyBigWig 精确计算；更准确但大型文件可能更慢",
    )
    parser.add_argument(
        "--plot-type", choices=("histogram", "line"), default="histogram",
        help="生成配置使用的 Circos 轨道类型",
    )
    parser.add_argument("--r0", default="0.50r", help="轨道内半径，写入配置文件")
    parser.add_argument("--r1", default="0.65r", help="轨道外半径，写入配置文件")
    parser.add_argument(
        "--color", default="circle_color_01",
        help="轨道颜色名称，需已在 Circos 颜色配置中定义",
    )
    parser.add_argument(
        "--display-min", type=float,
        help="配置中的显示下限；默认根据数据自动计算",
    )
    parser.add_argument(
        "--display-max", type=float,
        help="配置中的显示上限；默认使用指定百分位数",
    )
    parser.add_argument(
        "--max-percentile", type=float, default=99.0,
        help="未指定 --display-max 时使用的自动上限百分位数",
    )
    parser.add_argument(
        "--config-data-path",
        help="配置文件中 file 参数使用的路径；默认 data/<输出数据文件名>",
    )
    parser.add_argument(
        "--overwrite", action="store_true",
        help="允许覆盖已有的同名数据和配置文件",
    )
    return parser


def make_track_config(
    data_reference: str,
    statistic: str,
    plot_type: str,
    r0: str,
    r1: str,
    color: str,
    display_min: float,
    display_max: float,
    window_size: int,
    step_size: int,
) -> str:
    """生成可以手动 include 或复制到主配置中的 Circos 配置文本。"""
    fill_parameters = ""
    if plot_type == "histogram":
        fill_parameters = f"""# 柱形轨道填充色；名称必须在 colors.conf 中存在。
fill_color = {color}
# 关闭柱形边框，避免窗口很多时出现密集黑线。
thickness = 0p
# 不将不相邻的窗口延伸并连接起来。
extend_bin = no
"""
    else:
        fill_parameters = """# 折线宽度，p 表示像素。
thickness = 2p
"""

    legend_symbol = "line" if plot_type == "line" else "box"
    return f"""# bigWig 覆盖度轨道配置
# 原始 bigWig 已按 {window_size} bp 窗口、{step_size} bp 步长汇总。
# 窗口统计方法：{statistic}
# 本文件只包含一个 <plot> 块，应放进 config/tracks.conf 的 <plots> 内。
# legend_label = Sequencing coverage
# legend_symbol = {legend_symbol}
<plot>
# 可选 histogram（柱形）或 line（折线）。
type = {plot_type}
# 路径相对于运行 circos 命令时的当前项目目录，而非本配置文件的位置。
file = {data_reference}
# r0 和 r1 分别是轨道内、外半径；r 表示相对于 Circos 圆半径。
r0 = {r0}
r1 = {r1}
# min/max 是绘图显示范围，超出 max 的极高值可能被截断。
min = {display_min:.8g}
max = {display_max:.8g}
# 主轨道颜色；需要在 colors.conf 中定义。
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


def main() -> int:
    args = make_parser().parse_args()
    bigwig_path = args.bigwig.resolve()
    output_dir = args.output_dir.resolve()
    step_size = args.step_size or args.window_size

    if not bigwig_path.is_file():
        print(f"错误：找不到 bigWig 文件：{bigwig_path}", file=sys.stderr)
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

    try:
        import pyBigWig  # type: ignore[import-not-found]
    except ImportError:
        print("错误：未安装 pyBigWig。", file=sys.stderr)
        print("可运行：conda install -c bioconda pybigwig", file=sys.stderr)
        return 1

    prefix = args.prefix or safe_prefix(bigwig_path)
    if not re.fullmatch(r"[A-Za-z0-9_.-]+", prefix):
        print("错误：--prefix 只能包含字母、数字、点、下划线和连字符。", file=sys.stderr)
        return 2

    data_path = output_dir / f"{prefix}_coverage.txt"
    config_path = output_dir / f"{prefix}_coverage_plot.conf"
    existing = [path for path in (data_path, config_path) if path.exists()]
    if existing and not args.overwrite:
        names = "、".join(path.name for path in existing)
        print(f"错误：{names} 已存在；如需覆盖请添加 --overwrite。", file=sys.stderr)
        return 2

    output_dir.mkdir(parents=True, exist_ok=True)
    values_written: list[float] = []
    chromosomes_written = 0
    windows_written = 0
    bigwig: Any = None

    try:
        bigwig = pyBigWig.open(str(bigwig_path))
        if bigwig is None or not bigwig.isBigWig():
            raise ValueError("输入文件不是有效的 bigWig")
        chromosomes = bigwig.chroms()

        with data_path.open("w", encoding="utf-8", newline="\n") as output:
            for chromosome, chromosome_length in chromosomes.items():
                if include_pattern and not include_pattern.search(chromosome):
                    continue
                if exclude_pattern and exclude_pattern.search(chromosome):
                    continue
                chromosomes_written += 1

                for start in range(0, chromosome_length, step_size):
                    end = min(start + args.window_size, chromosome_length)
                    result = bigwig.stats(
                        chromosome,
                        start,
                        end,
                        type=args.statistic,
                        nBins=1,
                        exact=args.exact,
                    )[0]
                    if result is None or not math.isfinite(result):
                        if args.missing == "skip":
                            continue
                        result = 0.0
                    output.write(f"{chromosome}\t{start}\t{end}\t{result:.8g}\n")
                    values_written.append(float(result))
                    windows_written += 1
    except (OSError, RuntimeError, ValueError) as error:
        data_path.unlink(missing_ok=True)
        config_path.unlink(missing_ok=True)
        print(f"错误：处理 bigWig 失败：{error}", file=sys.stderr)
        return 1
    finally:
        if bigwig is not None:
            bigwig.close()

    if not values_written:
        data_path.unlink(missing_ok=True)
        print("错误：筛选后没有可输出的有效窗口。", file=sys.stderr)
        return 1

    observed_min = min(values_written)
    display_min = (
        args.display_min if args.display_min is not None else min(0.0, observed_min)
    )
    display_max = (
        args.display_max
        if args.display_max is not None
        else percentile(values_written, args.max_percentile)
    )
    if display_max <= display_min:
        display_max = display_min + (1.0 if display_min == 0 else abs(display_min) * 0.05)

    data_reference = (
        args.config_data_path
        if args.config_data_path
        else f"data/{data_path.name}"
    )
    config_text = make_track_config(
        data_reference=data_reference,
        statistic=args.statistic,
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

    print(f"完成：处理 {chromosomes_written} 条染色体、输出 {windows_written} 个窗口")
    print(f"数据文件：{data_path}")
    print(f"配置文件：{config_path}")
    print(f"配置显示范围：{display_min:.8g} 至 {display_max:.8g}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
