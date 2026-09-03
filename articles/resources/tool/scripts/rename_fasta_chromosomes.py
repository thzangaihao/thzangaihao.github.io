#!/usr/bin/env python3
"""Rename every FASTA record using a common prefix and sequential numbers."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path


def positive_integer(value: str) -> int:
    """Return a positive integer for argparse."""
    try:
        number = int(value)
    except ValueError as exc:
        raise argparse.ArgumentTypeError("序号长度必须是整数") from exc
    if number < 1:
        raise argparse.ArgumentTypeError("序号长度必须大于或等于 1")
    return number


def rename_fasta(input_path: Path, output_path: Path, prefix: str, width: int) -> int:
    """Rename FASTA headers and return the number of records written."""
    if not input_path.is_file():
        raise ValueError(f"输入文件不存在或不是普通文件：{input_path}")
    if input_path.resolve() == output_path.resolve():
        raise ValueError("输入文件和输出文件不能是同一个文件")

    output_path.parent.mkdir(parents=True, exist_ok=True)
    record_count = 0
    saw_sequence_before_header = False

    try:
        with input_path.open("r", encoding="utf-8-sig", newline="") as source, \
                output_path.open("w", encoding="utf-8", newline="") as destination:
            for line_number, line in enumerate(source, start=1):
                if line.startswith(">"):
                    record_count += 1
                    destination.write(f">{prefix}{record_count:0{width}d}\n")
                else:
                    if line.strip() and record_count == 0:
                        saw_sequence_before_header = True
                        raise ValueError(
                            f"第 {line_number} 行出现在首个 FASTA 序列头之前，文件格式无效"
                        )
                    destination.write(line)
    except Exception:
        # Do not leave a misleading partial result when parsing or writing fails.
        if output_path.exists():
            output_path.unlink()
        raise

    if record_count == 0:
        output_path.unlink(missing_ok=True)
        if saw_sequence_before_header:
            raise ValueError("输入文件不是有效的 FASTA 文件")
        raise ValueError("输入文件中没有找到 FASTA 序列头（以 > 开头的行）")

    return record_count


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="按出现顺序统一重命名 FASTA 中的所有序列头。"
    )
    parser.add_argument("input", type=Path, help="输入 FASTA 文件")
    parser.add_argument("prefix", help="新序列头的前缀，例如 chr 或 c")
    parser.add_argument(
        "width",
        type=positive_integer,
        help="序号最小长度，例如 2 会生成 01、02、03……",
    )
    parser.add_argument(
        "-o",
        "--output",
        type=Path,
        help="输出 FASTA 文件（默认：输入文件名.renamed.原扩展名）",
    )
    return parser


def default_output_path(input_path: Path) -> Path:
    if input_path.suffix:
        return input_path.with_name(
            f"{input_path.stem}.renamed{input_path.suffix}"
        )
    return input_path.with_name(f"{input_path.name}.renamed.fasta")


def main() -> int:
    parser = build_parser()
    args = parser.parse_args()
    output_path = args.output or default_output_path(args.input)

    try:
        count = rename_fasta(args.input, output_path, args.prefix, args.width)
    except (OSError, UnicodeError, ValueError) as exc:
        print(f"错误：{exc}", file=sys.stderr)
        return 1

    print(f"完成：已重命名 {count} 条序列")
    print(f"输出文件：{output_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
