#!/usr/bin/env python3
"""从 FASTA 生成 Circos 核型数据 karyotype.txt。

本脚本只输出一个核型数据文件，不计算 GC，也不生成任何配置文件。
"""

from __future__ import annotations

import argparse
import gzip
import re
import sys
from pathlib import Path
from typing import Iterator, TextIO


# 与 create_circos_project.py 生成的 config/colors.conf 名称一致。
COLORS = tuple(f"circle_color_{index:02d}" for index in range(1, 13))


def nonnegative_integer(text: str) -> int:
    value = int(text)
    if value < 0:
        raise argparse.ArgumentTypeError("必须是大于或等于 0 的整数")
    return value


def open_fasta(path: Path) -> TextIO:
    """支持普通 FASTA 和 gzip 压缩 FASTA。"""
    if path.suffix.lower() == ".gz":
        return gzip.open(path, "rt", encoding="utf-8")
    return path.open("rt", encoding="utf-8")


def fasta_lengths(path: Path) -> Iterator[tuple[str, int]]:
    """按 FASTA 输入顺序返回标题首字段和序列长度。"""
    name: str | None = None
    length = 0
    names_seen: set[str] = set()
    with open_fasta(path) as handle:
        for line_number, raw_line in enumerate(handle, start=1):
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    yield name, length
                header = line[1:].strip()
                if not header:
                    raise ValueError(f"第 {line_number} 行的 FASTA 标题为空")
                name = header.split()[0]
                if name in names_seen:
                    raise ValueError(f"FASTA 中存在重复序列名称：{name}")
                names_seen.add(name)
                length = 0
            else:
                if name is None:
                    raise ValueError("序列出现在第一个 FASTA 标题之前")
                length += len("".join(line.split()))
    if name is not None:
        yield name, length


def make_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="从 FASTA 生成 Circos karyotype.txt。",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("fasta", type=Path, help="输入 .fa/.fasta 或 .gz FASTA")
    parser.add_argument(
        "-o", "--output", type=Path, default=Path("karyotype.txt"),
        help="输出核型数据文件",
    )
    parser.add_argument(
        "--min-sequence-length", type=nonnegative_integer, default=0,
        help="忽略短于该长度的序列，单位 bp",
    )
    parser.add_argument("--include", metavar="REGEX", help="只保留名称匹配的序列")
    parser.add_argument("--exclude", metavar="REGEX", help="排除名称匹配的序列")
    parser.add_argument(
        "--label-prefix", default="",
        help="仅为显示标签添加前缀，不改变用于匹配数据轨道的序列 ID",
    )
    parser.add_argument("--overwrite", action="store_true", help="覆盖已有输出文件")
    return parser


def main() -> int:
    args = make_parser().parse_args()
    fasta = args.fasta.resolve()
    output = args.output.resolve()
    if not fasta.is_file():
        print(f"错误：找不到 FASTA：{fasta}", file=sys.stderr)
        return 2
    if output.exists() and not args.overwrite:
        print(f"错误：输出已存在：{output}；请添加 --overwrite。", file=sys.stderr)
        return 2
    try:
        include = re.compile(args.include) if args.include else None
        exclude = re.compile(args.exclude) if args.exclude else None
    except re.error as error:
        print(f"错误：正则表达式无效：{error}", file=sys.stderr)
        return 2

    output.parent.mkdir(parents=True, exist_ok=True)
    count = 0
    try:
        with output.open("w", encoding="utf-8", newline="\n") as handle:
            for sequence_id, length in fasta_lengths(fasta):
                if length < args.min_sequence_length:
                    continue
                if include and not include.search(sequence_id):
                    continue
                if exclude and exclude.search(sequence_id):
                    continue
                color = COLORS[count % len(COLORS)]
                label = f"{args.label_prefix}{sequence_id}"
                # 七列含义：类型、分隔符、ID、显示标签、起点、终点、颜色名。
                handle.write(
                    f"chr - {sequence_id} {label} 0 {length} {color}\n"
                )
                count += 1
    except (OSError, UnicodeError, ValueError) as error:
        output.unlink(missing_ok=True)
        print(f"错误：生成核型数据失败：{error}", file=sys.stderr)
        return 1
    if count == 0:
        output.unlink(missing_ok=True)
        print("错误：筛选后没有可输出的序列。", file=sys.stderr)
        return 1
    print(f"完成：输出 {count} 条核型记录")
    print(f"核型数据：{output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
