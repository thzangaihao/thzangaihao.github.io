#!/usr/bin/env python3
"""Generate overlapping k-mers from a FASTA file."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path
from typing import Iterator, TextIO


def positive_integer(value: str) -> int:
    try:
        number = int(value)
    except ValueError as exc:
        raise argparse.ArgumentTypeError("k 必须是正整数") from exc
    if number < 1:
        raise argparse.ArgumentTypeError("k 必须是正整数")
    return number


def read_fasta(handle: TextIO) -> Iterator[tuple[str, str]]:
    header: str | None = None
    parts: list[str] = []
    for line_number, raw_line in enumerate(handle, 1):
        line = raw_line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if header is not None:
                yield header, "".join(parts)
            header = line[1:].strip()
            if not header:
                raise ValueError(f"第 {line_number} 行的 FASTA 序列头为空")
            parts = []
        else:
            if header is None:
                raise ValueError(f"第 {line_number} 行的序列出现在 FASTA 序列头之前")
            parts.append("".join(line.split()))
    if header is not None:
        yield header, "".join(parts)


def generate_kmers(input_path: Path, output_path: Path, k: int) -> int:
    if not input_path.is_file():
        raise ValueError(f"输入文件不存在：{input_path}")
    if input_path.resolve() == output_path.resolve():
        raise ValueError("输入文件和输出文件不能相同")

    written = 0
    with input_path.open("r", encoding="utf-8-sig") as source, output_path.open(
        "w", encoding="utf-8", newline="\n"
    ) as destination:
        for header, sequence in read_fasta(source):
            count = max(0, len(sequence) - k + 1)
            width = max(5, len(str(count)))
            for index in range(count):
                destination.write(f">{header}|{index + 1:0{width}d}\n")
                destination.write(sequence[index : index + k] + "\n")
            written += count
    return written


def main() -> int:
    parser = argparse.ArgumentParser(
        description="从 FASTA 生成逐碱基滑动、互相重叠的 k-mer FASTA 文件。"
    )
    parser.add_argument("input", type=Path, help="输入 FASTA 文件")
    parser.add_argument("output", type=Path, help="输出 FASTA 文件")
    parser.add_argument("-k", "--kmer-size", type=positive_integer, required=True, help="k-mer 长度")
    args = parser.parse_args()
    try:
        count = generate_kmers(args.input, args.output, args.kmer_size)
    except (OSError, UnicodeError, ValueError) as exc:
        parser.exit(1, f"错误：{exc}\n")
    print(f"已生成 {count} 条 k-mer：{args.output}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
