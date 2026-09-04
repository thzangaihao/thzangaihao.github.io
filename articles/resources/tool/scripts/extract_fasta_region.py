#!/usr/bin/env python3
"""Interactively extract one chromosome region from a FASTA file."""

from __future__ import annotations

import argparse
import gzip
import re
import sys
from pathlib import Path
from typing import Iterator, TextIO


def open_text(path: Path) -> TextIO:
    """Open a plain-text or gzip-compressed file for reading."""
    if path.suffix.lower() == ".gz":
        return gzip.open(path, "rt", encoding="utf-8")
    return path.open("r", encoding="utf-8")


def read_fasta(path: Path) -> Iterator[tuple[str, str, str]]:
    """Yield (record ID, full header, sequence) from a FASTA file."""
    record_id: str | None = None
    header = ""
    parts: list[str] = []

    with open_text(path) as handle:
        for line_number, raw_line in enumerate(handle, start=1):
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if record_id is not None:
                    yield record_id, header, "".join(parts).upper()
                header = line[1:].strip()
                if not header:
                    raise ValueError(f"FASTA 第 {line_number} 行缺少序列名称")
                record_id = header.split()[0]
                parts = []
            else:
                if record_id is None:
                    raise ValueError(
                        f"FASTA 第 {line_number} 行的序列出现在标题行之前"
                    )
                parts.append("".join(line.split()))

    if record_id is not None:
        yield record_id, header, "".join(parts).upper()


def choose_record(records: list[tuple[str, str, int]]) -> int:
    """Ask the user to choose a record by number or ID."""
    id_to_index = {record_id: index for index, (record_id, _, _) in enumerate(records)}
    while True:
        answer = input("\n请输入染色体编号或名称：").strip()
        if answer in id_to_index:
            return id_to_index[answer]
        if answer.isdigit() and 1 <= int(answer) <= len(records):
            return int(answer) - 1
        print("输入无效，请输入表格中的编号或完整染色体名称。", file=sys.stderr)


def choose_region(sequence_length: int) -> tuple[int, int]:
    """Ask for a 1-based inclusive region and return zero-based slice bounds."""
    pattern = re.compile(r"^\s*([0-9][0-9,]*)\s*[-:]\s*([0-9][0-9,]*)\s*$")
    while True:
        answer = input(
            f"请输入提取范围（1-{sequence_length}，格式如 1001-5000；输入 all 提取整条）："
        ).strip()
        if answer.lower() in {"all", "全部"}:
            return 0, sequence_length
        match = pattern.match(answer)
        if not match:
            print("范围格式错误，请使用 起点-终点，例如 1001-5000。", file=sys.stderr)
            continue
        start = int(match.group(1).replace(",", ""))
        end = int(match.group(2).replace(",", ""))
        if start < 1 or end < start or end > sequence_length:
            print(
                f"范围必须满足 1 <= 起点 <= 终点 <= {sequence_length}。",
                file=sys.stderr,
            )
            continue
        return start - 1, end


def safe_filename(record_id: str, start: int, end: int) -> str:
    safe_id = re.sub(r"[^A-Za-z0-9_.-]+", "_", record_id).strip("._") or "sequence"
    return f"{safe_id}_{start}-{end}.fasta"


def write_fasta(path: Path, header: str, sequence: str, width: int = 80) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="\n") as handle:
        handle.write(f">{header}\n")
        for offset in range(0, len(sequence), width):
            handle.write(sequence[offset : offset + width] + "\n")


def scan_fasta(path: Path) -> list[tuple[str, str, int]]:
    """Return record metadata without retaining genome sequences in memory."""
    metadata: list[tuple[str, str, int]] = []
    for record_id, header, sequence in read_fasta(path):
        metadata.append((record_id, header, len(sequence)))
    return metadata


def load_record(path: Path, wanted_id: str) -> tuple[str, str]:
    for record_id, header, sequence in read_fasta(path):
        if record_id == wanted_id:
            return header, sequence
    raise ValueError(f"没有找到染色体：{wanted_id}")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="列出 FASTA 中的染色体，并交互式提取指定坐标范围。"
    )
    parser.add_argument("fasta", type=Path, help="输入 FASTA 或 FASTA.GZ 文件")
    parser.add_argument(
        "-o", "--output", type=Path,
        help="输出 FASTA 路径；不指定时在交互过程中询问",
    )
    parser.add_argument(
        "--line-width", type=int, default=80,
        help="输出序列每行碱基数（默认：80）",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    input_path = args.fasta.expanduser().resolve()
    if not input_path.is_file():
        raise FileNotFoundError(f"输入文件不存在：{input_path}")
    if args.line_width < 1:
        raise ValueError("--line-width 必须大于 0")

    records = scan_fasta(input_path)
    if not records:
        raise ValueError("输入文件中没有 FASTA 记录")
    ids = [record_id for record_id, _, _ in records]
    seen: set[str] = set()
    duplicates = sorted({record_id for record_id in ids if record_id in seen or seen.add(record_id)})
    if duplicates:
        raise ValueError(f"FASTA 中存在重复序列名称：{', '.join(duplicates)}")

    print(f"\n输入文件：{input_path}")
    print(f"共发现 {len(records)} 条序列：\n")
    index_width = len(str(len(records)))
    name_width = max(6, max(len(record_id) for record_id, _, _ in records))
    print(f"{'编号':>{index_width + 2}}  {'染色体':<{name_width}}  {'长度(bp)':>15}")
    print(f"{'-' * (index_width + 2)}  {'-' * name_width}  {'-' * 15}")
    for index, (record_id, _, sequence_length) in enumerate(records, start=1):
        print(f"{index:>{index_width + 2}}  {record_id:<{name_width}}  {sequence_length:>15,}")

    selected = choose_record(records)
    record_id, _, sequence_length = records[selected]
    slice_start, slice_end = choose_region(sequence_length)
    _, sequence = load_record(input_path, record_id)
    start = slice_start + 1
    end = slice_end
    extracted = sequence[slice_start:slice_end]

    default_output = Path(safe_filename(record_id, start, end))
    if args.output is None:
        answer = input(f"输出 FASTA 路径（直接回车使用 {default_output}）：").strip()
        output_path = Path(answer) if answer else default_output
    else:
        output_path = args.output
    output_path = output_path.expanduser().resolve()

    output_header = f"{record_id}:{start}-{end} source={input_path.name}"
    write_fasta(output_path, output_header, extracted, args.line_width)
    print(f"\n提取完成：{record_id}:{start}-{end}（{len(extracted):,} bp）")
    print(f"输出文件：{output_path}")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (EOFError, KeyboardInterrupt):
        print("\n操作已取消。", file=sys.stderr)
        raise SystemExit(130)
    except (OSError, ValueError) as error:
        print(f"错误：{error}", file=sys.stderr)
        raise SystemExit(2)
