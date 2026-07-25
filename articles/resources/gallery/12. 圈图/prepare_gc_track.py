#!/usr/bin/env python3
"""从 FASTA 生成 Circos GC 数据和单个 <plot> 配置片段。

本脚本输出：
1. <prefix>_gc.txt：四列 GC 数据；
2. <prefix>_gc_plot.conf：不含 <plots> 外壳的一个 <plot> 块。
"""

from __future__ import annotations

import argparse
import gzip
import math
import re
import sys
from pathlib import Path
from typing import Iterator, TextIO


def positive_integer(text: str) -> int:
    value = int(text)
    if value <= 0:
        raise argparse.ArgumentTypeError("必须是大于 0 的整数")
    return value


def nonnegative_integer(text: str) -> int:
    value = int(text)
    if value < 0:
        raise argparse.ArgumentTypeError("必须是大于或等于 0 的整数")
    return value


def open_fasta(path: Path) -> TextIO:
    if path.suffix.lower() == ".gz":
        return gzip.open(path, "rt", encoding="utf-8")
    return path.open("rt", encoding="utf-8")


def read_fasta(path: Path) -> Iterator[tuple[str, str]]:
    """按输入顺序返回标题首字段和大写序列。"""
    name: str | None = None
    parts: list[str] = []
    names_seen: set[str] = set()
    with open_fasta(path) as handle:
        for line_number, raw_line in enumerate(handle, start=1):
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(parts).upper()
                header = line[1:].strip()
                if not header:
                    raise ValueError(f"第 {line_number} 行的 FASTA 标题为空")
                name = header.split()[0]
                if name in names_seen:
                    raise ValueError(f"FASTA 中存在重复序列名称：{name}")
                names_seen.add(name)
                parts = []
            else:
                if name is None:
                    raise ValueError("序列出现在第一个 FASTA 标题之前")
                parts.append("".join(line.split()))
    if name is not None:
        yield name, "".join(parts).upper()


def gc_value(sequence: str, denominator: str) -> float | None:
    """计算 GC；acgt 模式不把 N 等非标准字符计入分母。"""
    gc = sequence.count("G") + sequence.count("C")
    total = (
        gc + sequence.count("A") + sequence.count("T")
        if denominator == "acgt"
        else len(sequence)
    )
    return gc / total if total else None


def percentile(values: list[float], percentage: float) -> float:
    ordered = sorted(values)
    if len(ordered) == 1:
        return ordered[0]
    position = (len(ordered) - 1) * percentage / 100
    lower, upper = math.floor(position), math.ceil(position)
    if lower == upper:
        return ordered[lower]
    fraction = position - lower
    return ordered[lower] * (1 - fraction) + ordered[upper] * fraction


def safe_prefix(path: Path) -> str:
    name = path.name
    if name.lower().endswith(".gz"):
        name = name[:-3]
    for suffix in (".fasta", ".fna", ".fa"):
        if name.lower().endswith(suffix):
            name = name[: -len(suffix)]
            break
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", name).strip("._") or "genome"


def make_plot(
    data_reference: str,
    plot_type: str,
    r0: str,
    r1: str,
    color: str,
    display_min: float,
    display_max: float,
    window: int,
    step: int,
    scale: str,
) -> str:
    """生成不带 <plots> 容器的单个 GC <plot>。"""
    style = (
        f"fill_color = {color}\nthickness = 0p\nextend_bin = no"
        if plot_type == "histogram"
        else "thickness = 2p"
    )
    return f"""# GC 含量轨道：窗口 {window} bp，步长 {step} bp，尺度 {scale}
# 将本文件内容放入 config/tracks.conf 的 <plots>...</plots> 内。
# legend_label = GC content
# legend_symbol = {plot_type if plot_type == "line" else "box"}
<plot>
type = {plot_type}
file = {data_reference}
r0 = {r0}
r1 = {r1}
min = {display_min:.8g}
max = {display_max:.8g}
color = {color}
{style}
<axes>
<axis>
spacing = 0.2r
color = lgrey
thickness = 1p
</axis>
</axes>
</plot>
"""


def make_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="从 FASTA 生成 GC 数据和单个 Circos plot 配置。",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("fasta", type=Path, help="输入 FASTA 或 .gz FASTA")
    parser.add_argument(
        "-o", "--output-dir", type=Path, default=Path("circos_gc"),
        help="GC 数据和 plot 配置的输出目录",
    )
    parser.add_argument("--prefix", help="输出前缀；默认来自 FASTA 文件名")
    parser.add_argument("-w", "--window-size", type=positive_integer, default=100_000)
    parser.add_argument("-s", "--step-size", type=positive_integer)
    parser.add_argument(
        "--min-sequence-length", type=nonnegative_integer, default=0
    )
    parser.add_argument("--include", metavar="REGEX")
    parser.add_argument("--exclude", metavar="REGEX")
    parser.add_argument(
        "--gc-denominator", choices=("acgt", "all"), default="acgt"
    )
    parser.add_argument(
        "--gc-scale", choices=("fraction", "percent"), default="fraction"
    )
    parser.add_argument(
        "--plot-type", choices=("line", "histogram"), default="line"
    )
    parser.add_argument("--r0", default="0.68r")
    parser.add_argument("--r1", default="0.84r")
    parser.add_argument("--color", default="circle_color_03")
    parser.add_argument("--display-min", type=float)
    parser.add_argument("--display-max", type=float)
    parser.add_argument("--min-percentile", type=float, default=1.0)
    parser.add_argument("--max-percentile", type=float, default=99.0)
    parser.add_argument(
        "--config-data-path",
        help="plot 中的 file 路径；默认 data/<数据文件名>",
    )
    parser.add_argument("--overwrite", action="store_true")
    return parser


def main() -> int:
    args = make_parser().parse_args()
    fasta = args.fasta.resolve()
    output_dir = args.output_dir.resolve()
    step = args.step_size or args.window_size
    if not fasta.is_file():
        print(f"错误：找不到 FASTA：{fasta}", file=sys.stderr)
        return 2
    if not 0 <= args.min_percentile < args.max_percentile <= 100:
        print("错误：百分位数需满足 0 <= min < max <= 100。", file=sys.stderr)
        return 2
    try:
        include = re.compile(args.include) if args.include else None
        exclude = re.compile(args.exclude) if args.exclude else None
    except re.error as error:
        print(f"错误：正则表达式无效：{error}", file=sys.stderr)
        return 2

    prefix = args.prefix or safe_prefix(fasta)
    if not re.fullmatch(r"[A-Za-z0-9_.-]+", prefix):
        print("错误：--prefix 含有不安全字符。", file=sys.stderr)
        return 2
    data_path = output_dir / f"{prefix}_gc.txt"
    plot_path = output_dir / f"{prefix}_gc_plot.conf"
    if not args.overwrite and (data_path.exists() or plot_path.exists()):
        print("错误：输出文件已存在；请添加 --overwrite。", file=sys.stderr)
        return 2

    output_dir.mkdir(parents=True, exist_ok=True)
    values: list[float] = []
    sequence_count = 0
    try:
        with data_path.open("w", encoding="utf-8", newline="\n") as output:
            for name, sequence in read_fasta(fasta):
                if len(sequence) < args.min_sequence_length:
                    continue
                if include and not include.search(name):
                    continue
                if exclude and exclude.search(name):
                    continue
                sequence_count += 1
                for start in range(0, len(sequence), step):
                    end = min(start + args.window_size, len(sequence))
                    value = gc_value(sequence[start:end], args.gc_denominator)
                    if value is None:
                        continue
                    if args.gc_scale == "percent":
                        value *= 100
                    output.write(f"{name}\t{start}\t{end}\t{value:.8g}\n")
                    values.append(value)
    except (OSError, UnicodeError, ValueError) as error:
        data_path.unlink(missing_ok=True)
        print(f"错误：生成 GC 数据失败：{error}", file=sys.stderr)
        return 1
    if not values:
        data_path.unlink(missing_ok=True)
        print("错误：没有可输出的有效 GC 窗口。", file=sys.stderr)
        return 1

    display_min = (
        args.display_min
        if args.display_min is not None
        else percentile(values, args.min_percentile)
    )
    display_max = (
        args.display_max
        if args.display_max is not None
        else percentile(values, args.max_percentile)
    )
    if display_max <= display_min:
        limit = 100.0 if args.gc_scale == "percent" else 1.0
        padding = 1.0 if args.gc_scale == "percent" else 0.01
        display_min = max(0.0, display_min - padding)
        display_max = min(limit, display_max + padding)
        if display_max <= display_min:
            display_max = display_min + padding

    data_reference = args.config_data_path or f"data/{data_path.name}"
    try:
        plot_path.write_text(
            make_plot(
                data_reference,
                args.plot_type,
                args.r0,
                args.r1,
                args.color,
                display_min,
                display_max,
                args.window_size,
                step,
                args.gc_scale,
            ),
            encoding="utf-8",
        )
    except OSError as error:
        data_path.unlink(missing_ok=True)
        print(f"错误：写入 plot 配置失败：{error}", file=sys.stderr)
        return 1
    print(f"完成：处理 {sequence_count} 条序列，输出 {len(values)} 个窗口")
    print(f"GC 数据：{data_path}")
    print(f"plot 配置：{plot_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
