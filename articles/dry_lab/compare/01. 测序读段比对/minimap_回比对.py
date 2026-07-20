#!/usr/bin/env python3
"""使用 minimap2 将多种格式的测序读段回比对至参考序列。

支持 FASTQ/FASTA（可 gzip 压缩）、SAM、BAM 和 CRAM。SAM/BAM/CRAM 会由
samtools fastq 流式转换后交给 minimap2，不会产生临时 FASTQ 文件。
"""

from __future__ import annotations

import argparse
import os
import shutil
import subprocess
import sys
from datetime import datetime
from pathlib import Path


SEQUENCING_PRESETS = {
    "hifi": "map-hifi",
    "ont": "map-ont",
    "clr": "map-pb",
    "sr": "sr",
}
SEQUENCE_SUFFIXES = {
    ".fa", ".fasta", ".fna", ".fas", ".fq", ".fastq",
}
ALIGNMENT_SUFFIXES = {".sam", ".bam", ".cram"}


def log_info(message: str) -> None:
    print(f"[{datetime.now().strftime('%H:%M:%S')}] {message}", flush=True)


def fail(message: str) -> "None":
    raise SystemExit(f"错误：{message}")


def file_suffix(path: Path) -> str:
    """识别普通或 gzip 压缩序列文件的实际扩展名。"""
    suffixes = [suffix.lower() for suffix in path.suffixes]
    if suffixes and suffixes[-1] in {".gz", ".bgz"} and len(suffixes) > 1:
        return suffixes[-2]
    return suffixes[-1] if suffixes else ""


def sample_stem(path: Path) -> str:
    name = path.name
    for suffix in (".fasta.gz", ".fastq.gz", ".fna.gz", ".fq.gz", ".fa.gz"):
        if name.lower().endswith(suffix):
            return name[: -len(suffix)]
    return path.stem


def split_paths(value: str) -> list[str]:
    """交互模式下允许用英文逗号输入双端 reads。"""
    return [item.strip().strip('"').strip("'") for item in value.split(",") if item.strip()]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="将 FASTQ/FASTA/SAM/BAM/CRAM 测序读段回比对并生成排序 BAM。"
    )
    parser.add_argument("-r", "--reference", help="参考 FASTA 路径")
    parser.add_argument(
        "-i", "--reads", nargs="+",
        help="一个 reads 文件，或两个双端 FASTQ/FASTA 文件",
    )
    parser.add_argument(
        "-x", "--technology", choices=SEQUENCING_PRESETS,
        help="测序类型：hifi、ont、clr 或 sr",
    )
    parser.add_argument("-t", "--threads", type=int, help="线程数")
    parser.add_argument("-o", "--output", help="输出 BAM 路径")
    parser.add_argument(
        "--no-index", action="store_true", help="不为输出 BAM 创建 BAI/CSI 索引"
    )
    return parser.parse_args()


def complete_interactively(args: argparse.Namespace) -> argparse.Namespace:
    print("=" * 60)
    print("       测序读段多格式回比对与覆盖度准备工具")
    print("=" * 60)
    if not args.reference:
        args.reference = input("参考 FASTA 路径 >>> ").strip().strip('"')
    if not args.reads:
        value = input("Reads 路径（双端文件用英文逗号分隔）>>> ").strip()
        args.reads = split_paths(value)
    if not args.technology:
        value = input("测序类型 [hifi/ont/clr/sr，默认 hifi] >>> ").strip().lower()
        args.technology = value or "hifi"
    if args.threads is None:
        default_threads = min(16, os.cpu_count() or 1)
        value = input(f"线程数 [默认 {default_threads}] >>> ").strip()
        try:
            args.threads = int(value) if value else default_threads
        except ValueError:
            fail("线程数必须是整数")
    return args


def validate(args: argparse.Namespace) -> tuple[Path, list[Path], Path]:
    if args.technology not in SEQUENCING_PRESETS:
        fail("测序类型必须是 hifi、ont、clr 或 sr")
    if not args.threads or args.threads < 1:
        fail("线程数必须大于 0")

    reference = Path(args.reference).expanduser().resolve()
    reads = [Path(item).expanduser().resolve() for item in (args.reads or [])]
    if not reference.is_file():
        fail(f"参考序列不存在：{reference}")
    if file_suffix(reference) not in {".fa", ".fasta", ".fna", ".fas"}:
        fail("参考序列必须是 FASTA 格式（.fa/.fasta/.fna/.fas，可使用 .gz）")
    if not reads:
        fail("至少需要一个 reads 文件")
    if len(reads) > 2:
        fail("最多支持两个双端 reads 文件")
    for path in reads:
        if not path.is_file():
            fail(f"reads 文件不存在：{path}")
        if file_suffix(path) not in SEQUENCE_SUFFIXES | ALIGNMENT_SUFFIXES:
            fail(f"不支持的 reads 格式：{path.name}")

    alignment_inputs = [p for p in reads if file_suffix(p) in ALIGNMENT_SUFFIXES]
    if alignment_inputs and len(reads) != 1:
        fail("SAM/BAM/CRAM 每次只能输入一个；双端信息将由 samtools 自动提取")
    if len(reads) == 2 and args.technology != "sr":
        fail("两个 reads 文件仅适用于双端短读段，请使用 --technology sr")

    output = (
        Path(args.output).expanduser()
        if args.output
        else Path.cwd() / f"{sample_stem(reference)}.sorted.bam"
    )
    if output.suffix.lower() != ".bam":
        output = output.with_suffix(".bam")
    output = output.resolve()
    if output in reads or output == reference:
        fail("输出路径不能覆盖输入文件")
    output.parent.mkdir(parents=True, exist_ok=True)
    return reference, reads, output


def require_tools() -> None:
    missing = [name for name in ("minimap2", "samtools") if shutil.which(name) is None]
    if missing:
        fail("找不到外部程序：" + ", ".join(missing) + "。请先安装并加入 PATH。")


def minimap_command(reference: Path, reads: list[Path], args: argparse.Namespace) -> list[str]:
    return [
        "minimap2", "-ax", SEQUENCING_PRESETS[args.technology],
        "-t", str(args.threads), str(reference), *(str(path) for path in reads),
    ]


def run_pipeline(reference: Path, reads: list[Path], output: Path, args: argparse.Namespace) -> None:
    """运行转换/比对/排序管道，并正确传播任一子进程的错误。"""
    input_format = file_suffix(reads[0])
    converter: subprocess.Popen[bytes] | None = None
    mapper_stdin = None

    if input_format in ALIGNMENT_SUFFIXES:
        conversion_cmd = ["samtools", "fastq", "-@", str(args.threads)]
        if input_format == ".cram":
            conversion_cmd.extend(["-T", str(reference)])
        conversion_cmd.append(str(reads[0]))
        log_info(f"正在流式读取 {input_format[1:].upper()} 并转换为 FASTQ...")
        converter = subprocess.Popen(conversion_cmd, stdout=subprocess.PIPE)
        mapper_stdin = converter.stdout
        mapper_reads = [Path("-")]
    else:
        mapper_reads = reads

    map_cmd = minimap_command(reference, mapper_reads, args)
    view_cmd = ["samtools", "view", "-@", str(args.threads), "-u", "-"]
    sort_cmd = [
        "samtools", "sort", "-@", str(args.threads), "-o", str(output), "-",
    ]
    log_info(f"开始 minimap2 比对（preset: {SEQUENCING_PRESETS[args.technology]}）...")

    mapper = subprocess.Popen(map_cmd, stdin=mapper_stdin, stdout=subprocess.PIPE)
    if converter and converter.stdout:
        converter.stdout.close()
    viewer = subprocess.Popen(view_cmd, stdin=mapper.stdout, stdout=subprocess.PIPE)
    if mapper.stdout:
        mapper.stdout.close()
    sorter = subprocess.Popen(sort_cmd, stdin=viewer.stdout)
    if viewer.stdout:
        viewer.stdout.close()

    sort_code = sorter.wait()
    view_code = viewer.wait()
    map_code = mapper.wait()
    convert_code = converter.wait() if converter else 0
    codes = {
        "samtools fastq": convert_code,
        "minimap2": map_code,
        "samtools view": view_code,
        "samtools sort": sort_code,
    }
    failed = [f"{name} (退出码 {code})" for name, code in codes.items() if code]
    if failed:
        fail("处理管道运行失败：" + "；".join(failed))


def main() -> None:
    args = complete_interactively(parse_args())
    reference, reads, output = validate(args)
    require_tools()

    log_info("输入格式：" + file_suffix(reads[0]).lstrip(".").upper())
    run_pipeline(reference, reads, output, args)
    log_info(f"排序后的 BAM 已生成：{output}")

    if not args.no_index:
        log_info("正在创建 BAM 索引...")
        result = subprocess.run(
            ["samtools", "index", "-@", str(args.threads), str(output)],
            check=False,
        )
        if result.returncode:
            fail(f"BAM 已生成，但索引创建失败（退出码 {result.returncode}）")
        log_info(f"索引已生成：{output}.bai")

    print("-" * 60)
    log_info("完成。可将参考 FASTA 和输出 BAM 导入 IGV。")
    print(f"参考序列：{reference}")
    print(f"比对结果：{output}")


if __name__ == "__main__":
    main()
