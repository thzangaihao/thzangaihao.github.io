#!/usr/bin/env python3
"""Remap selected 1-based transcript-relative links and bedGraph tracks."""

from __future__ import annotations

import argparse
import re
import sys
from dataclasses import dataclass
from pathlib import Path


@dataclass(frozen=True)
class MappingBlock:
    relative_start: int  # 1-based, inclusive
    relative_end: int  # 1-based, inclusive
    genomic_start: int  # 1-based, inclusive
    genomic_end: int  # 1-based, inclusive


def read_mapping(path: Path) -> list[MappingBlock]:
    """Read ordered, closed genomic intervals and build a spliced coordinate map."""
    blocks: list[MappingBlock] = []
    next_relative_start = 1
    with path.open("r", encoding="utf-8-sig") as handle:
        for line_number, line in enumerate(handle, start=1):
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            fields = [field for field in re.split(r"[,\t ]+", line) if field]
            if len(fields) != 2:
                raise ValueError(f"mapping.txt 第 {line_number} 行必须只有起点和终点")
            try:
                genomic_start, genomic_end = map(int, fields)
            except ValueError as exc:
                raise ValueError(f"mapping.txt 第 {line_number} 行包含非整数坐标") from exc
            if genomic_start < 1 or genomic_end < genomic_start:
                raise ValueError(
                    f"mapping.txt 第 {line_number} 行不是有效的左闭右闭区间"
                )
            length = genomic_end - genomic_start + 1
            relative_end = next_relative_start + length - 1
            blocks.append(
                MappingBlock(
                    next_relative_start,
                    relative_end,
                    genomic_start,
                    genomic_end,
                )
            )
            next_relative_start = relative_end + 1
    if not blocks:
        raise ValueError("mapping.txt 中没有有效区间")
    return blocks


def map_interval(
    start: int, end: int, blocks: list[MappingBlock]
) -> list[tuple[int, int]]:
    """Map a relative half-open interval, splitting it at mapping boundaries."""
    if start < 1 or end <= start:
        return []
    mapped: list[tuple[int, int]] = []
    for block in blocks:
        # A closed relative block [a,b] is [a,b+1) in track coordinates.
        overlap_start = max(start, block.relative_start)
        overlap_end = min(end, block.relative_end + 1)
        if overlap_start >= overlap_end:
            continue
        # mapping.txt is 1-based closed, while links/bedGraph output is BED-like
        # 0-based half-open. Thus relative base 1 -> [genomic_start-1, genomic_start).
        genomic_start = block.genomic_start - 1 + overlap_start - block.relative_start
        genomic_end = block.genomic_start - 1 + overlap_end - block.relative_start
        mapped.append((genomic_start, genomic_end))
    return mapped


def remap_links(path: Path, output: Path, blocks: list[MappingBlock]) -> tuple[int, int]:
    written = skipped = 0
    output_lines: list[str] = []
    with path.open("r", encoding="utf-8-sig") as handle:
        for line_number, line in enumerate(handle, start=1):
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\r\n").split("\t")
            if len(fields) < 6:
                print(f"警告：{path.name} 第 {line_number} 行少于 6 列，已跳过", file=sys.stderr)
                skipped += 1
                continue
            try:
                start1, end1 = int(fields[1]), int(fields[2])
                start2, end2 = int(fields[4]), int(fields[5])
            except ValueError:
                print(f"警告：{path.name} 第 {line_number} 行坐标无效，已跳过", file=sys.stderr)
                skipped += 1
                continue
            mapped1 = map_interval(start1, end1, blocks)
            mapped2 = map_interval(start2, end2, blocks)
            # One links anchor cannot represent a discontinuous genomic interval.
            if len(mapped1) != 1 or len(mapped2) != 1:
                skipped += 1
                continue
            fields[1], fields[2] = map(str, mapped1[0])
            fields[4], fields[5] = map(str, mapped2[0])
            output_lines.append("\t".join(fields))
            written += 1
    output.write_text(
        "\n".join(output_lines) + ("\n" if output_lines else ""),
        encoding="utf-8",
        newline="\n",
    )
    return written, skipped


def remap_bedgraph(path: Path, output: Path, blocks: list[MappingBlock]) -> tuple[int, int]:
    written = skipped = 0
    output_lines: list[str] = []
    with path.open("r", encoding="utf-8-sig") as handle:
        for line_number, line in enumerate(handle, start=1):
            if not line.strip() or line.startswith(("#", "track", "browser")):
                continue
            fields = line.rstrip("\r\n").split("\t")
            if len(fields) < 4:
                print(f"警告：{path.name} 第 {line_number} 行少于 4 列，已跳过", file=sys.stderr)
                skipped += 1
                continue
            try:
                start, end = int(fields[1]), int(fields[2])
            except ValueError:
                print(f"警告：{path.name} 第 {line_number} 行坐标无效，已跳过", file=sys.stderr)
                skipped += 1
                continue
            mapped_parts = map_interval(start, end, blocks)
            if not mapped_parts:
                skipped += 1
                continue
            for genomic_start, genomic_end in mapped_parts:
                output_lines.append(
                    "\t".join(
                        (fields[0], str(genomic_start), str(genomic_end), *fields[3:])
                    )
                )
                written += 1
    output.write_text(
        "\n".join(output_lines) + ("\n" if output_lines else ""),
        encoding="utf-8",
        newline="\n",
    )
    return written, skipped


def track_type(path: Path) -> str:
    suffix = path.suffix.lower()
    if suffix == ".links":
        return "links"
    if suffix == ".bedgraph":
        return "bedgraph"
    raise ValueError(
        f"不支持的文件类型：{path.name}；仅支持 .links、.bedGraph 和 .bedgraph"
    )


def output_path_for(path: Path) -> Path:
    return path.with_name(f"{path.stem}_mapped{path.suffix}")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "根据有序的左闭右闭 mapping 区间，映射指定的 .links、.bedGraph 或 .bedgraph 文件。"
        )
    )
    parser.add_argument("mapping", type=Path, help="每行格式为 genomic_start,genomic_end")
    parser.add_argument(
        "tracks",
        type=Path,
        nargs="+",
        help="一个或多个需要映射的 .links/.bedGraph/.bedgraph 文件",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    if not args.mapping.is_file():
        print(f"错误：找不到映射文件：{args.mapping}", file=sys.stderr)
        return 1
    try:
        blocks = read_mapping(args.mapping)
        total_length = blocks[-1].relative_end
        print(f"映射区间：{len(blocks)} 段；相对坐标范围：1-{total_length}")
        for path in args.tracks:
            if not path.is_file():
                raise ValueError(f"找不到待映射文件：{path}")
            file_type = track_type(path)
            output = output_path_for(path)
            if file_type == "links":
                written, skipped = remap_links(path, output, blocks)
            else:
                written, skipped = remap_bedgraph(path, output, blocks)
            print(f"{path.name} -> {output.name}：输出 {written} 行，跳过 {skipped} 行")
    except (OSError, ValueError) as exc:
        print(f"错误：{exc}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
