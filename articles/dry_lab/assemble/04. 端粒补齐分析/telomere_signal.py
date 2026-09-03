#!/usr/bin/env python3
"""Calculate telomere signals and write an RIdeogram heatmap table."""

from __future__ import annotations

import argparse
import csv
import math
import sys
from pathlib import Path
from typing import Iterable, Iterator, Sequence, TextIO


DNA_COMPLEMENT = str.maketrans("ACGT", "TGCA")
VALID_DNA = frozenset("ACGT")


def reverse_complement(sequence: str) -> str:
    return sequence.translate(DNA_COMPLEMENT)[::-1]


def read_fasta(handle: TextIO) -> Iterator[tuple[str, str]]:
    """Yield (record ID, sequence) pairs from a FASTA file."""
    record_id: str | None = None
    sequence_parts: list[str] = []

    for line_number, raw_line in enumerate(handle, start=1):
        line = raw_line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if record_id is not None:
                yield record_id, "".join(sequence_parts).upper()
            header = line[1:].strip()
            if not header:
                raise ValueError(f"FASTA 第 {line_number} 行缺少序列名称")
            record_id = header.split()[0]
            sequence_parts = []
        else:
            if record_id is None:
                raise ValueError(f"FASTA 第 {line_number} 行的序列出现在标题行之前")
            sequence_parts.append("".join(line.split()))

    if record_id is not None:
        yield record_id, "".join(sequence_parts).upper()


def motif_intervals(sequence: str, motifs: Iterable[str]) -> list[tuple[int, int]]:
    """Return merged, zero-based half-open intervals covered by motif matches."""
    intervals: list[tuple[int, int]] = []
    for motif in motifs:
        start = sequence.find(motif)
        while start != -1:
            intervals.append((start, start + len(motif)))
            start = sequence.find(motif, start + 1)  # retain overlapping matches

    if not intervals:
        return []

    intervals.sort()
    merged: list[tuple[int, int]] = []
    current_start, current_end = intervals[0]
    for start, end in intervals[1:]:
        if start <= current_end:
            current_end = max(current_end, end)
        else:
            merged.append((current_start, current_end))
            current_start, current_end = start, end
    merged.append((current_start, current_end))
    return merged


def block_signals(
    sequence: str, motifs: Iterable[str], block_size: int
) -> Iterator[tuple[int, int, int, float]]:
    """Yield block ID, 1-based inclusive start/end, and signal."""
    intervals = motif_intervals(sequence, motifs)
    interval_index = 0

    for block_id, block_start in enumerate(range(0, len(sequence), block_size), start=1):
        block_end = min(block_start + block_size, len(sequence))
        while interval_index < len(intervals) and intervals[interval_index][1] <= block_start:
            interval_index += 1

        covered = 0
        scan_index = interval_index
        while scan_index < len(intervals) and intervals[scan_index][0] < block_end:
            start, end = intervals[scan_index]
            covered += max(0, min(end, block_end) - max(start, block_start))
            scan_index += 1

        signal = covered / (block_end - block_start)
        yield block_id, block_start + 1, block_end, signal


def plot_signals(
    chromosome_signals: dict[str, Sequence[tuple[int, int, float]]], output: Path
) -> bool:
    """Plot all chromosomes in one SVG; return False if matplotlib is absent."""
    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        from matplotlib.ticker import MaxNLocator
    except ImportError:
        print(
            "提示：未安装 matplotlib，已跳过折线图绘制；热图 TSV 不受影响。",
            file=sys.stderr,
        )
        return False

    chromosome_count = len(chromosome_signals)
    column_count = 1 if chromosome_count <= 6 else 2
    row_count = math.ceil(chromosome_count / column_count)
    figure, axes = plt.subplots(
        row_count,
        column_count,
        figsize=(8.5 * column_count, max(2.2 * row_count, 3.5)),
        squeeze=False,
        sharey=True,
        constrained_layout=True,
    )

    for axis, (chromosome_id, points) in zip(
        axes.flat, chromosome_signals.items()
    ):
        positions_mb = [((start + end) / 2) / 1_000_000 for start, end, _ in points]
        values = [value for _, _, value in points]
        axis.plot(positions_mb, values, color="#D73027", linewidth=0.8)
        axis.set_title(chromosome_id, fontsize=9, loc="left")
        axis.set_ylim(-0.02, 1.02)
        axis.set_xlabel("Position (Mb)", fontsize=8)
        axis.set_ylabel("Telomere signal", fontsize=8)
        axis.grid(axis="y", color="#D9D9D9", linewidth=0.5, alpha=0.8)
        axis.tick_params(labelsize=7)
        axis.xaxis.set_major_locator(MaxNLocator(nbins=6))

    for axis in list(axes.flat)[chromosome_count:]:
        axis.remove()

    figure.suptitle("Telomere signal by chromosome", fontsize=12)
    figure.savefig(output, format="svg", bbox_inches="tight")
    plt.close(figure)
    return True


def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="将 FASTA 中每条染色体划分为不重叠小块并计算端粒信号强度。"
    )
    parser.add_argument("fasta", type=Path, help="输入 FASTA 文件")
    parser.add_argument("telomere_sequence", help="端粒重复序列，例如 TTAGGG")
    parser.add_argument(
        "-o", "--output", type=Path, default=Path("telomere_signal_RIdeogram.tsv"),
        help="输出 RIdeogram 热图 TSV 路径（默认：telomere_signal_RIdeogram.tsv）",
    )
    parser.add_argument(
        "-b", "--block-size", type=int, default=10_000,
        help="小块长度，单位 bp（默认：10000）",
    )
    parser.add_argument(
        "--forward-only", action="store_true",
        help="只检测给定序列，不检测其反向互补序列",
    )
    parser.add_argument(
        "--precision", type=int, default=6,
        help="信号强度的小数位数（默认：6）",
    )
    parser.add_argument(
        "--plot-output", type=Path,
        help="折线图 SVG 路径（默认：在 TSV 文件名后添加 _lineplot.svg）",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_arguments()
    motif = args.telomere_sequence.upper().replace(" ", "")

    if args.block_size <= 0:
        raise ValueError("--block-size 必须是正整数")
    if args.precision < 0:
        raise ValueError("--precision 不能小于 0")
    if not motif or any(base not in VALID_DNA for base in motif):
        raise ValueError("端粒序列不能为空，且只能包含 A、C、G、T")

    motifs = {motif}
    if not args.forward_only:
        motifs.add(reverse_complement(motif))

    record_count = 0
    seen_ids: set[str] = set()
    chromosome_signals: dict[str, list[tuple[int, int, float]]] = {}
    with args.fasta.open("r", encoding="utf-8") as fasta_handle, args.output.open(
        "w", encoding="utf-8", newline=""
    ) as output_handle:
        writer = csv.writer(output_handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["Chr", "Start", "End", "Value"])
        for chromosome_id, sequence in read_fasta(fasta_handle):
            if chromosome_id in seen_ids:
                raise ValueError(f"FASTA 中存在重复的序列编号：{chromosome_id}")
            seen_ids.add(chromosome_id)
            record_count += 1
            chromosome_signals[chromosome_id] = []
            for _, start, end, signal in block_signals(sequence, motifs, args.block_size):
                writer.writerow([chromosome_id, start, end, f"{signal:.{args.precision}f}"])
                chromosome_signals[chromosome_id].append((start, end, signal))

    if record_count == 0:
        raise ValueError("输入文件中没有 FASTA 记录")
    print(f"完成：已分析 {record_count} 条序列，结果写入 {args.output}")

    plot_output = args.plot_output or args.output.with_name(
        f"{args.output.stem}_lineplot.svg"
    )
    if plot_signals(chromosome_signals, plot_output):
        print(f"折线图已写入 {plot_output}")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (OSError, ValueError) as error:
        print(f"错误：{error}", file=sys.stderr)
        raise SystemExit(1)
