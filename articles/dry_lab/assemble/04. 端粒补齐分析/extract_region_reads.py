#!/usr/bin/env python3
"""Interactively select genome regions and extract overlapping reads."""

from __future__ import annotations

import argparse
import gzip
import shutil
import subprocess
import sys
from pathlib import Path
from typing import TextIO


def open_text(path: Path) -> TextIO:
    if path.suffix.lower() == ".gz":
        return gzip.open(path, "rt", encoding="utf-8")
    return path.open("r", encoding="utf-8")


def fasta_lengths(path: Path) -> list[tuple[str, int]]:
    records: list[tuple[str, int]] = []
    name: str | None = None
    length = 0
    seen: set[str] = set()
    with open_text(path) as handle:
        for line_number, raw_line in enumerate(handle, start=1):
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    records.append((name, length))
                header = line[1:].strip()
                if not header:
                    raise ValueError(f"FASTA 第 {line_number} 行缺少序列名称")
                name = header.split()[0]
                if name in seen:
                    raise ValueError(f"FASTA 中存在重复序列名称：{name}")
                seen.add(name)
                length = 0
            else:
                if name is None:
                    raise ValueError(f"FASTA 第 {line_number} 行的序列出现在标题行之前")
                length += len("".join(line.split()))
    if name is not None:
        records.append((name, length))
    if not records:
        raise ValueError("输入文件中没有 FASTA 记录")
    return records


def show_records(records: list[tuple[str, int]]) -> None:
    number_width = len(str(len(records)))
    name_width = max(6, max(len(name) for name, _ in records))
    print(f"\n共发现 {len(records)} 条染色体/序列：\n")
    print(f"{'编号':>{number_width + 2}}  {'染色体':<{name_width}}  {'长度(bp)':>15}")
    print(f"{'-' * (number_width + 2)}  {'-' * name_width}  {'-' * 15}")
    for number, (name, length) in enumerate(records, start=1):
        print(f"{number:>{number_width + 2}}  {name:<{name_width}}  {length:>15,}")


def ask_chromosome(records: list[tuple[str, int]]) -> tuple[str, int]:
    by_name = {name: (name, length) for name, length in records}
    while True:
        answer = input("\n请输入染色体编号或名称：").strip()
        if answer in by_name:
            return by_name[answer]
        if answer.isdigit() and 1 <= int(answer) <= len(records):
            return records[int(answer) - 1]
        print("输入无效，请输入表格中的编号或完整名称。", file=sys.stderr)


def ask_positive_integer(prompt: str, default: int | None = None) -> int:
    while True:
        suffix = f"（默认 {default:,}）" if default is not None else ""
        answer = input(f"{prompt}{suffix}：").strip().replace(",", "")
        if not answer and default is not None:
            return default
        if answer.isdigit() and int(answer) > 0:
            return int(answer)
        print("请输入大于 0 的整数。", file=sys.stderr)


def ask_region(length: int, default_end_size: int) -> tuple[int, int, str]:
    while True:
        mode = input("选择区域：左端(L)、右端(R)或自定义(C) [L/R/C]：").strip().lower()
        if mode in {"l", "left", "左", "左端"}:
            size = min(ask_positive_integer("左端长度 bp", default_end_size), length)
            return 0, size, "left"
        if mode in {"r", "right", "右", "右端"}:
            size = min(ask_positive_integer("右端长度 bp", default_end_size), length)
            return length - size, length, "right"
        if mode in {"c", "custom", "自定义"}:
            start = ask_positive_integer(f"起点（1-{length:,}，1-based）")
            end = ask_positive_integer(f"终点（{start:,}-{length:,}，两端包含）")
            if 1 <= start <= end <= length:
                return start - 1, end, "custom"
            print(f"坐标必须满足 1 <= 起点 <= 终点 <= {length:,}。", file=sys.stderr)
            continue
        print(f"无法识别区域类型：{mode or '<空>'}", file=sys.stderr)


def choose_regions(records: list[tuple[str, int]], end_size: int) -> list[tuple[str, int, int, str]]:
    selected: list[tuple[str, int, int, str]] = []
    while True:
        chromosome, length = ask_chromosome(records)
        start, end, label = ask_region(length, end_size)
        region = (chromosome, start, end, label)
        if region not in selected:
            selected.append(region)
            print(f"已添加：{chromosome}:{start + 1}-{end} ({label})")
        else:
            print("该区域已经添加。")
        if input("继续添加区域吗？[y/N]：").strip().lower() not in {"y", "yes", "是"}:
            return selected


def resolve_program(value: str, label: str) -> str:
    found = shutil.which(value)
    if found:
        return found
    candidate = Path(value).expanduser()
    if candidate.is_file():
        return str(candidate.resolve())
    raise FileNotFoundError(f"找不到 {label}：{value}")


def write_bed(path: Path, regions: list[tuple[str, int, int, str]]) -> None:
    with path.open("w", encoding="utf-8", newline="\n") as handle:
        for chromosome, start, end, label in regions:
            handle.write(f"{chromosome}\t{start}\t{end}\t{chromosome}_{label}\n")


def extract_names(samtools: str, bam: Path, bed: Path, output: Path, threads: int) -> int:
    command = [samtools, "view", "-@", str(threads), "-L", str(bed), str(bam)]
    names: set[str] = set()
    process = subprocess.Popen(command, stdout=subprocess.PIPE, text=True, encoding="utf-8")
    assert process.stdout is not None
    for line in process.stdout:
        name = line.split("\t", 1)[0]
        if name:
            names.add(name)
    return_code = process.wait()
    if return_code != 0:
        raise RuntimeError(f"samtools view 运行失败，退出码：{return_code}")
    with output.open("w", encoding="utf-8", newline="\n") as handle:
        for name in sorted(names):
            handle.write(name + "\n")
    return len(names)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="交互式选择基因组区域，通过 samtools 和 seqkit 提取完整 reads。"
    )
    parser.add_argument("genome", type=Path, help="用于比对的参考基因组 FASTA(.gz)")
    parser.add_argument("bam", type=Path, help="坐标排序后的 BAM 文件")
    parser.add_argument("reads", type=Path, help="生成 BAM 时使用的原始 FASTQ(.gz)")
    parser.add_argument("-o", "--output-dir", type=Path, default=Path("terminal_reads"))
    parser.add_argument("--prefix", default="terminal", help="输出文件前缀")
    parser.add_argument("--end-size", type=int, default=50_000, help="默认末端长度（bp）")
    parser.add_argument("-t", "--threads", type=int, default=8)
    parser.add_argument("--samtools", default="samtools", help="samtools 命令或路径")
    parser.add_argument("--seqkit", default="seqkit", help="seqkit 命令或路径")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    genome = args.genome.expanduser().resolve()
    bam = args.bam.expanduser().resolve()
    reads = args.reads.expanduser().resolve()
    for path, label in ((genome, "参考基因组"), (bam, "BAM"), (reads, "FASTQ")):
        if not path.is_file():
            raise FileNotFoundError(f"{label}文件不存在：{path}")
    if args.threads < 1 or args.end_size < 1:
        raise ValueError("线程数和 --end-size 必须大于 0")

    samtools = resolve_program(args.samtools, "samtools")
    seqkit = resolve_program(args.seqkit, "seqkit")
    records = fasta_lengths(genome)
    show_records(records)
    regions = choose_regions(records, args.end_size)

    output_dir = args.output_dir.expanduser().resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    bed = output_dir / f"{args.prefix}.regions.bed"
    names = output_dir / f"{args.prefix}.read_names.txt"
    output_fastq = output_dir / f"{args.prefix}.reads.fastq.gz"
    write_bed(bed, regions)
    print(f"\nBED 已写入：{bed}")

    count = extract_names(samtools, bam, bed, names, args.threads)
    print(f"找到 {count:,} 条不重复的候选 reads；名称已写入：{names}")
    if count == 0:
        print("没有提取 FASTQ：所选区域在 BAM 中没有比对记录。", file=sys.stderr)
        return 3

    command = [seqkit, "grep", "--threads", str(args.threads), "-f", str(names), str(reads), "-o", str(output_fastq)]
    print("正在从原始 FASTQ 提取完整 reads……")
    result = subprocess.run(command, check=False)
    if result.returncode != 0:
        raise RuntimeError(f"seqkit grep 运行失败，退出码：{result.returncode}")
    print(f"提取完成：{output_fastq}")
    print("提示：这些是与区域重叠的候选 reads，尚未按向外悬垂长度和方向筛选。")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (EOFError, KeyboardInterrupt):
        print("\n操作已取消。", file=sys.stderr)
        raise SystemExit(130)
    except (OSError, ValueError, RuntimeError) as error:
        print(f"错误：{error}", file=sys.stderr)
        raise SystemExit(2)
