#!/usr/bin/env python3
"""A small, dependency-free wrapper for running TGS-GapCloser.

The wrapper accepts long reads in FASTQ/FASTQ.GZ, converts them to FASTA (the
format required by TGS-GapCloser), and then starts TGS-GapCloser.  It is meant
for assemblies that contain explicit N gaps.  TGS-GapCloser cannot perform
true one-sided extension beyond a chromosome end; see --help for details.
"""

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


def fastq_to_fasta(fastq: Path, fasta: Path) -> int:
    """Convert FASTQ/FASTQ.GZ (including wrapped records) to FASTA."""
    count = 0
    with open_text(fastq) as source, fasta.open("w", encoding="utf-8", newline="\n") as target:
        while True:
            header = source.readline()
            if not header:
                break
            if not header.startswith("@"):
                raise ValueError(f"FASTQ 第 {count + 1} 条记录缺少 @ 标题行")
            sequence_parts: list[str] = []
            while True:
                line = source.readline()
                if not line:
                    raise ValueError(f"FASTQ 第 {count + 1} 条记录不完整")
                if line.startswith("+"):
                    break
                sequence_parts.append(line.strip())
            sequence = "".join(sequence_parts)
            quality_length = 0
            while quality_length < len(sequence):
                quality = source.readline()
                if not quality:
                    raise ValueError(f"FASTQ 第 {count + 1} 条记录缺少质量值")
                quality_length += len(quality.rstrip("\r\n"))
            if not sequence or quality_length != len(sequence):
                raise ValueError(f"FASTQ 第 {count + 1} 条记录的序列与质量长度不一致")
            name = header[1:].strip().split()[0]
            if not name:
                raise ValueError(f"FASTQ 第 {count + 1} 条记录缺少 read 名称")
            target.write(f">{name}\n{sequence}\n")
            count += 1
    if count == 0:
        raise ValueError("FASTQ 中没有读段")
    return count


def decompress_gzip(source: Path, target: Path) -> None:
    with gzip.open(source, "rb") as compressed, target.open("wb") as plain:
        shutil.copyfileobj(compressed, plain)


def validate_fasta(path: Path, label: str) -> None:
    with open_text(path) as handle:
        first = next((line.strip() for line in handle if line.strip()), "")
    if not first.startswith(">"):
        raise ValueError(f"{label} 不是 FASTA 格式：{path}")


def detect_reads_format(path: Path) -> str:
    with open_text(path) as handle:
        first = next((line.strip() for line in handle if line.strip()), "")
    if first.startswith("@"):
        return "fastq"
    if first.startswith(">"):
        return "fasta"
    raise ValueError(f"无法识别测序数据格式：{path}")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="将长读段 FASTQ 自动转成 FASTA，然后运行 TGS-GapCloser。",
        epilog=(
            "注意：TGS-GapCloser 只能填补 scaffold 内部由 N 表示的 gap"
        ),
    )
    parser.add_argument("reads", type=Path, help="ONT/PacBio 长读段 FASTQ(.gz) 或 FASTA(.gz)")
    parser.add_argument("chromosomes", type=Path, help="主染色体 FASTA")
    parser.add_argument("-o", "--output-dir", type=Path, default=Path("tgs_gapcloser_out"))
    parser.add_argument("--prefix", default="gapclosed", help="输出前缀（默认：gapclosed）")
    parser.add_argument("-t", "--threads", type=int, default=16)
    parser.add_argument("--read-type", choices=("ont", "pb", "hifi"), default="ont")
    parser.add_argument(
        "--executable", default=None,
        help="tgsgapcloser/tgsgapcloser2 路径；默认从 PATH 自动寻找",
    )
    parser.add_argument(
        "--racon", default=None,
        help="用 Racon 校正原始长读段；不指定时使用 --ne（建议输入已校正 reads）",
    )
    parser.add_argument("--min-identity", type=float, default=None)
    parser.add_argument("--min-match", type=int, default=None)
    parser.add_argument("--min-reads", type=int, default=1, help="支持 gap 的最少 reads 数")
    parser.add_argument("--dry-run", action="store_true", help="只准备输入并显示命令")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    reads = args.reads.expanduser().resolve()
    chromosomes = args.chromosomes.expanduser().resolve()
    if not reads.is_file() or not chromosomes.is_file():
        raise FileNotFoundError("测序数据或染色体 FASTA 不存在")
    if args.threads < 1 or args.min_reads < 1:
        raise ValueError("线程数和 --min-reads 必须大于 0")
    if args.min_identity is not None and not 0 < args.min_identity <= 1:
        raise ValueError("--min-identity 必须在 (0, 1] 范围内")
    validate_fasta(chromosomes, "主染色体文件")

    out_dir = args.output_dir.expanduser().resolve()
    out_dir.mkdir(parents=True, exist_ok=True)
    if chromosomes.suffix.lower() == ".gz":
        plain_chromosomes = out_dir / "chromosomes.for_tgsgapcloser.fasta"
        decompress_gzip(chromosomes, plain_chromosomes)
        chromosomes = plain_chromosomes
        print(f"已解压主染色体 FASTA：{chromosomes}")
    reads_format = detect_reads_format(reads)
    if reads_format == "fastq":
        tgs_reads = out_dir / "long_reads.for_tgsgapcloser.fasta"
        count = fastq_to_fasta(reads, tgs_reads)
        print(f"已将 {count:,} 条长读段转换为 {tgs_reads}")
    else:
        validate_fasta(reads, "长读段文件")
        if reads.suffix.lower() == ".gz":
            tgs_reads = out_dir / "long_reads.for_tgsgapcloser.fasta"
            decompress_gzip(reads, tgs_reads)
            print(f"已解压长读段 FASTA：{tgs_reads}")
        else:
            tgs_reads = reads

    executable = args.executable
    if executable is None:
        executable = shutil.which("tgsgapcloser2") or shutil.which("tgsgapcloser")
    elif not Path(executable).is_file() and shutil.which(executable) is None:
        raise FileNotFoundError(f"找不到 TGS-GapCloser：{executable}")
    if executable is None:
        raise FileNotFoundError(
            "PATH 中找不到 tgsgapcloser2 或 tgsgapcloser；请先安装，或用 --executable 指定路径"
        )
    if args.read_type == "hifi" and "tgsgapcloser2" not in Path(executable).name.lower():
        raise ValueError("--read-type hifi 需要 TGS-GapCloser2；旧版仅支持 ont/pb")

    output_prefix = out_dir / args.prefix
    command = [
        str(executable), "--scaff", str(chromosomes), "--reads", str(tgs_reads),
        "--output", str(output_prefix), "--tgstype", args.read_type,
        "--thread", str(args.threads), "--min_nread", str(args.min_reads),
    ]
    if args.racon:
        racon = shutil.which(args.racon) or (str(Path(args.racon).resolve()) if Path(args.racon).is_file() else None)
        if racon is None:
            raise FileNotFoundError(f"找不到 Racon：{args.racon}")
        command.extend(["--racon", racon])
    else:
        command.append("--ne")
        print("警告：未指定 --racon，将跳过纠错；原始 ONT/PacBio reads 建议提供 Racon。", file=sys.stderr)
    if args.min_identity is not None:
        command.extend(["--min_idy", str(args.min_identity)])
    if args.min_match is not None:
        command.extend(["--min_match", str(args.min_match)])

    print("运行命令：", subprocess.list2cmdline(command))
    if args.dry_run:
        return 0
    with (out_dir / "pipeline.stdout.log").open("w", encoding="utf-8") as stdout, \
         (out_dir / "pipeline.stderr.log").open("w", encoding="utf-8") as stderr:
        result = subprocess.run(command, cwd=out_dir, stdout=stdout, stderr=stderr, check=False)
    if result.returncode != 0:
        print(f"TGS-GapCloser 失败（退出码 {result.returncode}），请查看 {out_dir / 'pipeline.stderr.log'}", file=sys.stderr)
        return result.returncode
    print(f"运行完成；结果和日志位于 {out_dir}")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (OSError, ValueError) as error:
        print(f"错误：{error}", file=sys.stderr)
        raise SystemExit(2)
