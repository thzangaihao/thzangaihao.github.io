#!/usr/bin/env python3
"""Rename FASTA sequences and, optionally, matching GFF/GTF/GFF3 seqids."""

from __future__ import annotations

import argparse
import os
import sys
import tempfile
from pathlib import Path
from typing import Callable, TextIO


def positive_integer(value: str) -> int:
    try:
        number = int(value)
    except ValueError as exc:
        raise argparse.ArgumentTypeError("序号长度必须是整数") from exc
    if number < 1:
        raise argparse.ArgumentTypeError("序号长度必须大于或等于 1")
    return number


def _check_input_output(input_path: Path, output_path: Path) -> None:
    if not input_path.is_file():
        raise ValueError(f"输入文件不存在或不是普通文件：{input_path}")
    if input_path.resolve() == output_path.resolve():
        raise ValueError("输入文件和输出文件不能是同一个文件")


def _atomic_write(output_path: Path, writer: Callable[[TextIO], None]) -> None:
    """Write to a temporary file before replacing the destination."""
    output_path.parent.mkdir(parents=True, exist_ok=True)
    temporary_name: str | None = None
    try:
        with tempfile.NamedTemporaryFile(
            "w", encoding="utf-8", newline="", delete=False,
            dir=output_path.parent, prefix=f".{output_path.name}.", suffix=".tmp"
        ) as destination:
            temporary_name = destination.name
            writer(destination)
        os.replace(temporary_name, output_path)
    except Exception:
        if temporary_name:
            Path(temporary_name).unlink(missing_ok=True)
        raise


def build_fasta_mapping(input_path: Path, prefix: str, width: int) -> dict[str, str]:
    """Build the shared old-id to new-id mapping in FASTA record order."""
    if not input_path.is_file():
        raise ValueError(f"输入文件不存在或不是普通文件：{input_path}")
    mapping: dict[str, str] = {}
    with input_path.open("r", encoding="utf-8-sig", newline="") as source:
        for line_number, line in enumerate(source, start=1):
            if line.startswith(">"):
                value = line[1:].strip()
                if not value:
                    raise ValueError(f"FASTA 第 {line_number} 行的序列名称为空")
                old_id = value.split(maxsplit=1)[0]
                if old_id in mapping:
                    raise ValueError(f"FASTA 序列名称重复：{old_id}")
                mapping[old_id] = f"{prefix}{len(mapping) + 1:0{width}d}"
            elif line.strip() and not mapping:
                raise ValueError(f"第 {line_number} 行出现在首个 FASTA 序列头之前")
    if not mapping:
        raise ValueError("输入文件中没有找到 FASTA 序列头（以 > 开头的行）")
    return mapping


def rename_fasta(input_path: Path, output_path: Path, prefix: str, width: int,
                 mapping: dict[str, str] | None = None) -> int:
    """Rename FASTA headers and return the number of records written."""
    _check_input_output(input_path, output_path)
    mapping = mapping or build_fasta_mapping(input_path, prefix, width)

    def write(destination: TextIO) -> None:
        with input_path.open("r", encoding="utf-8-sig", newline="") as source:
            for line in source:
                if line.startswith(">"):
                    old_id = line[1:].strip().split(maxsplit=1)[0]
                    destination.write(f">{mapping[old_id]}\n")
                else:
                    destination.write(line)

    _atomic_write(output_path, write)
    return len(mapping)


def _annotation_seqids(input_path: Path) -> set[str]:
    seqids: set[str] = set()
    in_fasta = False
    with input_path.open("r", encoding="utf-8-sig", newline="") as source:
        for line_number, line in enumerate(source, start=1):
            stripped = line.rstrip("\r\n")
            if stripped == "##FASTA":
                in_fasta = True
            elif in_fasta and line.startswith(">"):
                value = line[1:].strip()
                if not value:
                    raise ValueError(f"注释文件第 {line_number} 行的 FASTA 名称为空")
                seqids.add(value.split(maxsplit=1)[0])
            elif not in_fasta and stripped.startswith("##sequence-region"):
                fields = stripped.split()
                if len(fields) >= 2:
                    seqids.add(fields[1])
            elif not in_fasta and stripped and not line.startswith("#"):
                fields = stripped.split("\t")
                if len(fields) < 9:
                    raise ValueError(f"注释文件第 {line_number} 行不足 9 列")
                seqids.add(fields[0])
    return seqids


def rename_annotation(input_path: Path, output_path: Path,
                      mapping: dict[str, str]) -> tuple[int, int]:
    """Rename annotation seqids; return (feature records, distinct seqids)."""
    _check_input_output(input_path, output_path)
    unknown = sorted(_annotation_seqids(input_path) - mapping.keys())
    if unknown:
        preview = ", ".join(unknown[:10]) + ("……" if len(unknown) > 10 else "")
        raise ValueError(f"注释文件含有 FASTA 中不存在的序列名称：{preview}")
    record_count = 0
    changed_seqids: set[str] = set()

    def write(destination: TextIO) -> None:
        nonlocal record_count
        in_fasta = False
        with input_path.open("r", encoding="utf-8-sig", newline="") as source:
            for line in source:
                ending = "\r\n" if line.endswith("\r\n") else "\n" if line.endswith("\n") else ""
                stripped = line.rstrip("\r\n")
                if stripped == "##FASTA":
                    in_fasta = True
                    destination.write(line)
                elif in_fasta and line.startswith(">"):
                    value = line[1:].strip()
                    old_id, *description = value.split(maxsplit=1)
                    changed_seqids.add(old_id)
                    tail = f" {description[0]}" if description else ""
                    destination.write(f">{mapping[old_id]}{tail}{ending}")
                elif not in_fasta and stripped.startswith("##sequence-region"):
                    fields = stripped.split()
                    old_id = fields[1]
                    changed_seqids.add(old_id)
                    fields[1] = mapping[old_id]
                    destination.write(" ".join(fields) + ending)
                elif not in_fasta and stripped and not line.startswith("#"):
                    fields = stripped.split("\t")
                    old_id = fields[0]
                    changed_seqids.add(old_id)
                    fields[0] = mapping[old_id]
                    destination.write("\t".join(fields) + ending)
                    record_count += 1
                else:
                    destination.write(line)

    _atomic_write(output_path, write)
    return record_count, len(changed_seqids)


def default_output_path(input_path: Path) -> Path:
    if input_path.suffix:
        return input_path.with_name(f"{input_path.stem}.renamed{input_path.suffix}")
    return input_path.with_name(f"{input_path.name}.renamed.fasta")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="按 FASTA 顺序重命名序列，并可同步修改 GFF/GTF/GFF3。"
    )
    parser.add_argument("input", type=Path, help="输入 FASTA 文件")
    parser.add_argument("prefix", help="新序列名称前缀，例如 chr 或 c")
    parser.add_argument("width", type=positive_integer, help="序号最小位数，例如 2 生成 01、02")
    parser.add_argument("-o", "--output", type=Path, help="输出 FASTA 文件（默认添加 .renamed）")
    parser.add_argument("-a", "--annotation", type=Path, help="配套 GFF、GTF 或 GFF3 文件")
    parser.add_argument("--annotation-output", type=Path,
                        help="注释输出文件（默认添加 .renamed）")
    return parser


def main() -> int:
    parser = build_parser()
    args = parser.parse_args()
    if args.annotation_output and not args.annotation:
        parser.error("--annotation-output 必须与 --annotation 一起使用")
    fasta_output = args.output or default_output_path(args.input)
    annotation_output = args.annotation_output or (
        default_output_path(args.annotation) if args.annotation else None
    )
    try:
        if args.annotation and annotation_output:
            input_paths = {args.input.resolve(), args.annotation.resolve()}
            output_paths = {fasta_output.resolve(), annotation_output.resolve()}
            if len(output_paths) != 2:
                raise ValueError("FASTA 与注释文件不能使用同一个输出路径")
            if input_paths & output_paths:
                raise ValueError("输出路径不能覆盖任一输入文件")
        mapping = build_fasta_mapping(args.input, args.prefix, args.width)
        # Validate both inputs before writing either member of the output pair.
        if args.annotation:
            _check_input_output(args.annotation, annotation_output)
            unknown = sorted(_annotation_seqids(args.annotation) - mapping.keys())
            if unknown:
                raise ValueError("注释文件含有 FASTA 中不存在的序列名称：" + ", ".join(unknown[:10]))
        count = rename_fasta(args.input, fasta_output, args.prefix, args.width, mapping)
        print(f"完成：已重命名 {count} 条 FASTA 序列")
        print(f"FASTA 输出：{fasta_output}")
        if args.annotation and annotation_output:
            records, seqids = rename_annotation(args.annotation, annotation_output, mapping)
            print(f"完成：已更新 {records} 条注释记录，涉及 {seqids} 个序列名称")
            print(f"注释输出：{annotation_output}")
    except (OSError, UnicodeError, ValueError) as exc:
        print(f"错误：{exc}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
