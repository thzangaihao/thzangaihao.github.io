#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""用 MUMmer4 对两个基因组 FASTA 进行一键式共线性分析。"""

from __future__ import annotations

import argparse
import csv
import glob
import json
import os
import re
import shutil
import subprocess
import sys
import time
from datetime import datetime
from pathlib import Path
from typing import Callable, TypeVar


FASTA_PATTERNS = ("*.fa", "*.fasta", "*.fna", "*.fas")
T = TypeVar("T")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="输入两个基因组 FASTA，运行 MUMmer4 核心分析并输出可移植的绘图数据。"
    )
    parser.add_argument("reference", nargs="?", help="参考基因组 FASTA")
    parser.add_argument("query", nargs="?", help="查询基因组 FASTA")
    parser.add_argument("-o", "--outdir", help="结果目录（默认自动创建带时间戳的目录）")
    parser.add_argument(
        "-t", "--threads", "--cores", dest="threads", type=int,
        default=max(1, os.cpu_count() or 1), help="运算核心数（默认使用全部可用核心）",
    )
    parser.add_argument("--min-match", type=int, default=100, help="nucmer 最小精确匹配长度，默认 100")
    parser.add_argument("--min-cluster", type=int, default=500, help="nucmer 最小聚类长度，默认 500")
    parser.add_argument("--min-identity", type=float, default=90.0, help="过滤后的最小一致性百分比，默认 90")
    parser.add_argument("--min-length", type=int, default=1000, help="过滤后的最小比对长度，默认 1000")
    parser.add_argument(
        "--filter-mode", choices=("one-to-one", "many-to-many"), default="one-to-one",
        help="一对一共线性（-1）或保留重复区域（-m），默认 one-to-one",
    )
    return parser.parse_args()


def find_fastas(root: Path) -> list[Path]:
    files: set[Path] = set()
    for pattern in FASTA_PATTERNS:
        files.update(Path(p).resolve() for p in glob.glob(str(root / "**" / pattern), recursive=True))
    return sorted(files, key=lambda p: str(p).lower())


def choose_fasta(files: list[Path], label: str, excluded: Path | None = None) -> Path:
    choices = [p for p in files if p != excluded]
    if not choices:
        raise SystemExit(f"错误：没有找到可用的{label} FASTA（支持 .fa/.fasta/.fna/.fas）。")
    print(f"\n扫描到 {len(choices)} 个可用文件，请选择{label}：")
    for index, path in enumerate(choices, 1):
        print(f"  [{index}] {path}")
    while True:
        answer = input("请输入编号（q 退出）：").strip().lower()
        if answer == "q":
            raise SystemExit(0)
        try:
            return choices[int(answer) - 1]
        except (ValueError, IndexError):
            print("输入无效，请重新输入。")


def checked_fasta(path_text: str, label: str) -> Path:
    path = Path(path_text).expanduser().resolve()
    if not path.is_file():
        raise SystemExit(f"错误：{label}不存在或不是文件：{path}")
    with path.open("r", encoding="utf-8", errors="replace") as handle:
        if not any(line.startswith(">") for line in handle):
            raise SystemExit(f"错误：{label}不像 FASTA 文件（未找到 > 开头的序列标题）：{path}")
    return path


def safe_name(path: Path) -> str:
    name = re.sub(r"\.(fa|fasta|fna|fas)(\.gz)?$", "", path.name, flags=re.IGNORECASE)
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", name).strip("._-") or "genome"


def require_programs() -> None:
    required = ["nucmer", "delta-filter", "show-coords", "show-snps"]
    missing = [name for name in required if shutil.which(name) is None]
    if missing:
        raise SystemExit(
            "错误：缺少程序：" + ", ".join(missing)
            + "\n请先安装 MUMmer4。R 与绘图组件不属于核心分析依赖。"
        )


def run(command: list[str], stdout: Path | None = None, cwd: Path | None = None) -> None:
    print("  $ " + " ".join(command))
    if stdout is None:
        subprocess.run(command, check=True, cwd=cwd)
    else:
        with stdout.open("w", encoding="utf-8", newline="") as handle:
            subprocess.run(command, check=True, stdout=handle, cwd=cwd)


def write_plot_files(
    coords: Path,
    plot_table: Path,
    metadata: Path,
    reference: Path,
    query: Path,
    args: argparse.Namespace,
) -> None:
    """把 MUMmer 坐标转换为跨电脑可直接绘图的自描述文件。"""
    columns = [
        "ref_start", "ref_end", "query_start", "query_end", "ref_align_length",
        "query_align_length", "identity", "ref_id", "query_id",
    ]
    row_count = 0
    with coords.open("r", encoding="utf-8", errors="replace", newline="") as source, \
            plot_table.open("w", encoding="utf-8", newline="") as target:
        writer = csv.writer(target, delimiter="\t", lineterminator="\n")
        writer.writerow(columns)
        for line_number, line in enumerate(source, 1):
            fields = line.rstrip("\r\n").split("\t")
            if len(fields) < 9:
                print(f"警告：跳过坐标文件第 {line_number} 行（列数不足）。", file=sys.stderr)
                continue
            writer.writerow(fields[:9])
            row_count += 1

    info = {
        "format": "mummer-synteny-plot-v1",
        "reference_name": safe_name(reference),
        "query_name": safe_name(query),
        "reference_file": reference.name,
        "query_file": query.name,
        "alignment_count": row_count,
        "min_identity": args.min_identity,
        "min_length": args.min_length,
        "filter_mode": args.filter_mode,
        "plot_table": plot_table.name,
    }
    metadata.write_text(json.dumps(info, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")


def prompt_value(
    label: str,
    default: T,
    converter: Callable[[str], T],
    validator: Callable[[T], bool],
) -> T:
    """读取一个带默认值的参数，直到用户输入合法值。"""
    while True:
        text = input(f"{label} [{default}]：").strip()
        if not text:
            return default
        try:
            value = converter(text)
            if validator(value):
                return value
        except ValueError:
            pass
        print("输入无效，请重新输入。")


def configure_interactively(args: argparse.Namespace, base_dir: Path) -> None:
    """将所有运行参数集中在一个交互阶段完成设置。"""
    available_cores = max(1, os.cpu_count() or 1)
    print("\n" + "-" * 64)
    print("参数配置（直接回车采用方括号中的默认值）")
    print("-" * 64)
    args.threads = prompt_value(
        f"运算核心数（本机检测到 {available_cores} 核）", min(args.threads, available_cores),
        int, lambda value: 1 <= value <= available_cores,
    )
    args.min_match = prompt_value("最小精确匹配长度", args.min_match, int, lambda value: value > 0)
    args.min_cluster = prompt_value("最小聚类长度", args.min_cluster, int, lambda value: value > 0)
    args.min_identity = prompt_value(
        "最小一致性（%）", args.min_identity, float, lambda value: 0 <= value <= 100
    )
    args.min_length = prompt_value("最小比对长度", args.min_length, int, lambda value: value > 0)

    mode_default = "1" if args.filter_mode == "one-to-one" else "2"
    mode = prompt_value(
        "过滤模式（1=一对一共线性，2=保留重复区域）", mode_default,
        str, lambda value: value in {"1", "2"},
    )
    args.filter_mode = "one-to-one" if mode == "1" else "many-to-many"

    default_outdir = "自动创建"
    outdir = input(f"结果目录 [{default_outdir}]：").strip()
    if outdir:
        args.outdir = str((base_dir / Path(outdir).expanduser()).resolve())


def print_configuration(args: argparse.Namespace, reference: Path, query: Path, outdir: Path) -> None:
    mode = "一对一共线性" if args.filter_mode == "one-to-one" else "保留重复区域"
    print("\n" + "-" * 64)
    print("运行配置摘要")
    print("-" * 64)
    print(f"参考基因组：{reference}")
    print(f"查询基因组：{query}")
    print(f"结果目录：  {outdir}")
    print(f"运算核心数：{args.threads}")
    print(f"匹配/聚类长度：{args.min_match} / {args.min_cluster}")
    print(f"一致性/比对长度：{args.min_identity}% / {args.min_length}")
    print(f"过滤模式：{mode}")
    print("-" * 64)


def main() -> None:
    args = parse_args()

    print("=" * 64)
    print(" MUMmer4 双基因组共线性一键分析")
    print("=" * 64)
    base_dir = Path.cwd()
    interactive = not (args.reference and args.query)
    files = find_fastas(base_dir) if interactive else []
    reference = checked_fasta(args.reference, "参考基因组") if args.reference else choose_fasta(files, "参考基因组")
    query = checked_fasta(args.query, "查询基因组") if args.query else choose_fasta(files, "查询基因组", reference)
    if reference == query:
        raise SystemExit("错误：参考基因组和查询基因组不能是同一个文件。")

    if interactive:
        configure_interactively(args, base_dir)
    if args.threads < 1 or args.min_match < 1 or args.min_cluster < 1 or args.min_length < 1:
        raise SystemExit("错误：核心数和长度阈值必须为正整数。")
    if not 0 <= args.min_identity <= 100:
        raise SystemExit("错误：--min-identity 必须在 0 到 100 之间。")

    require_programs()
    ref_name, qry_name = safe_name(reference), safe_name(query)
    pair_name = f"{ref_name}_vs_{qry_name}"
    outdir = Path(args.outdir).expanduser().resolve() if args.outdir else (
        base_dir / f"MUMmer_Synteny_{pair_name}_{datetime.now():%Y%m%d_%H%M%S}"
    )
    print_configuration(args, reference, query, outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    prefix = outdir / pair_name
    raw_delta = Path(f"{prefix}.delta")
    filtered_delta = Path(f"{prefix}.filtered.delta")
    coords = Path(f"{prefix}.coords.tsv")
    plot_table = Path(f"{prefix}.plot.tsv")
    plot_metadata = Path(f"{prefix}.plot.json")
    snps = Path(f"{prefix}.snps.tsv")
    start = time.time()

    print(f"\n参考：{reference}\n查询：{query}\n结果：{outdir}\n")
    print("[1/4] 全基因组比对")
    run([
        "nucmer", "--threads", str(args.threads), "--maxmatch",
        "-l", str(args.min_match), "-c", str(args.min_cluster),
        "-p", str(prefix), str(reference), str(query),
    ])

    print("\n[2/4] 筛选共线性比对")
    mode = "-1" if args.filter_mode == "one-to-one" else "-m"
    run([
        "delta-filter", mode, "-i", str(args.min_identity), "-l", str(args.min_length),
        str(raw_delta),
    ], filtered_delta)

    print("\n[3/4] 导出坐标与变异表")
    run(["show-coords", "-THr", str(filtered_delta)], coords)
    run(["show-snps", "-THClr", str(filtered_delta)], snps)

    print("\n[4/4] 生成可移植绘图核心文件")
    write_plot_files(coords, plot_table, plot_metadata, reference, query, args)

    print("\n" + "=" * 64)
    print(f"分析完成，用时 {(time.time() - start) / 60:.2f} 分钟。")
    print(f"过滤比对：{filtered_delta}")
    print(f"坐标表：  {coords}")
    print(f"变异表：  {snps}")
    print(f"绘图数据：{plot_table}")
    print(f"绘图信息：{plot_metadata}")
    print("=" * 64)


if __name__ == "__main__":
    try:
        main()
    except KeyboardInterrupt:
        print("\n用户终止运行。", file=sys.stderr)
        raise SystemExit(130)
    except subprocess.CalledProcessError as exc:
        print(f"\n错误：命令执行失败（退出码 {exc.returncode}）。", file=sys.stderr)
        raise SystemExit(exc.returncode)
