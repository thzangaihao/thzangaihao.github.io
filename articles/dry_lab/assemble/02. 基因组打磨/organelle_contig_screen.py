#!/usr/bin/env python3
"""批量评估组装 contig 与多个线粒体/叶绿体参考基因组的相似性。"""

import csv
import glob
import gzip
import math
import os
import re
import shutil
import subprocess
import sys
import time
from collections import defaultdict
from datetime import datetime
from pathlib import Path


FASTA_EXTENSIONS = (".fasta.gz", ".fa.gz", ".fna.gz", ".fas.gz", ".fasta", ".fa", ".fna", ".fas")
DEFAULT_THREADS = min(os.cpu_count() or 8, 32)
OUTFMT_FIELDS = (
    "qseqid", "qlen", "sseqid", "slen", "pident", "length", "mismatch",
    "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore",
)


def log(message):
    print(f"[{time.strftime('%Y-%m-%d %H:%M:%S')}] {message}", flush=True)


def ask_text(prompt, default=None):
    suffix = f"（回车使用 {default}）" if default is not None else ""
    value = input(f"{prompt}{suffix} >>> ").strip()
    if value:
        return value
    if default is not None:
        return str(default)
    raise ValueError("输入不能为空。")


def ask_yes_no(prompt, default=True):
    hint = "Y/n" if default else "y/N"
    while True:
        value = input(f"{prompt} [{hint}] >>> ").strip().lower()
        if not value:
            return default
        if value in {"y", "yes", "是"}:
            return True
        if value in {"n", "no", "否"}:
            return False
        print("请输入 y 或 n。")


def safe_name(path):
    name = Path(path).name
    for extension in FASTA_EXTENSIONS:
        if name.lower().endswith(extension):
            name = name[: -len(extension)]
            break
    return re.sub(r"[^0-9A-Za-z._-]+", "_", name).strip("._") or "sample"


def unique_labels(paths):
    counts = defaultdict(int)
    labels = []
    for path in paths:
        base = safe_name(path)
        counts[base] += 1
        labels.append(base if counts[base] == 1 else f"{base}_{counts[base]}")
    return labels


def scan_fastas():
    found = []
    for extension in FASTA_EXTENSIONS:
        found.extend(glob.glob(f"**/*{extension}", recursive=True))
    return sorted({path for path in found if not Path(path).name.startswith(".")})


def select_files(files, prompt, excluded=()):
    excluded = {os.path.normcase(os.path.abspath(path)) for path in excluded}
    candidates = [path for path in files if os.path.normcase(os.path.abspath(path)) not in excluded]
    if not candidates:
        raise ValueError("没有可选择的 FASTA 文件。")
    print(f"\n{prompt}")
    for index, path in enumerate(candidates, 1):
        print(f"  [{index}] {path}")
    print("输入编号（如 1,3 或 1-4），或输入 all 选择全部。")
    choice = input("你的选择 >>> ").strip().lower()
    if choice == "all":
        return candidates
    indices = []
    try:
        for part in choice.split(","):
            part = part.strip()
            if "-" in part:
                start, end = (int(item) for item in part.split("-", 1))
                indices.extend(range(start - 1, end))
            else:
                indices.append(int(part) - 1)
    except ValueError as error:
        raise ValueError("选择格式错误。") from error
    if not indices or any(index < 0 or index >= len(candidates) for index in indices):
        raise ValueError("选择中包含无效编号。")
    return [candidates[index] for index in dict.fromkeys(indices)]


def classify_reference(path):
    lower = Path(path).name.lower()
    inferred = "mitochondrion" if any(x in lower for x in ("mito", "mt")) else (
        "chloroplast" if any(x in lower for x in ("chloro", "plast", "cp")) else "other"
    )
    choices = {"m": "mitochondrion", "c": "chloroplast", "o": "other"}
    default = {value: key for key, value in choices.items()}[inferred]
    while True:
        value = input(
            f"参考 {path} 的类型：m=线粒体，c=叶绿体，o=其它（回车使用 {default}）>>> "
        ).strip().lower() or default
        if value in choices:
            return choices[value]
        print("请输入 m、c 或 o。")


def fasta_lengths(path):
    opener = gzip.open if str(path).lower().endswith(".gz") else open
    lengths = {}
    current = None
    with opener(path, "rt", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if line.startswith(">"):
                current = line[1:].strip().split()[0]
                if not current:
                    raise ValueError(f"空 FASTA 序列名：{path}")
                if current in lengths:
                    raise ValueError(f"FASTA 序列名重复：{current}（{path}）")
                lengths[current] = 0
            elif current is not None:
                lengths[current] += len("".join(line.split()))
    if not lengths:
        raise ValueError(f"没有在文件中读到 FASTA 序列：{path}")
    return lengths


def prepare_fasta(path, staging_dir, label):
    """BLAST+ 不直接依赖 gzip 支持；压缩输入统一解压到本次输出目录。"""
    if not str(path).lower().endswith(".gz"):
        return os.path.abspath(path)
    staging_dir.mkdir(parents=True, exist_ok=True)
    destination = staging_dir / f"{label}.fasta"
    log(f"解压暂存：{path} -> {destination}")
    with gzip.open(path, "rb") as source, open(destination, "wb") as target:
        shutil.copyfileobj(source, target, length=16 * 1024 * 1024)
    return str(destination)


def merge_length(intervals):
    if not intervals:
        return 0
    total = 0
    start, end = sorted(intervals)[0]
    for next_start, next_end in sorted(intervals)[1:]:
        if next_start <= end + 1:
            end = max(end, next_end)
        else:
            total += end - start + 1
            start, end = next_start, next_end
    return total + end - start + 1


def organelle_score(query_coverage_pct, identity_pct):
    """启发式证据分数，并非经过统计校准的概率。"""
    evidence = (query_coverage_pct / 100.0) * (identity_pct / 100.0) ** 2
    return 1.0 - math.exp(-6.0 * evidence)


def build_databases(references, labels, database_dir, staging_dir):
    databases = []
    for path, label in zip(references, labels):
        target_dir = database_dir / label
        target_dir.mkdir(parents=True, exist_ok=True)
        prefix = target_dir / "reference"
        log(f"建库：{label}")
        prepared_path = prepare_fasta(path, staging_dir, label)
        subprocess.run(
            ["makeblastdb", "-in", prepared_path, "-dbtype", "nucl", "-out", str(prefix)],
            check=True,
        )
        databases.append(prefix)
    return databases


def run_blast(assembly, database, output_file, threads):
    command = [
        "blastn", "-query", os.path.abspath(assembly), "-db", str(database),
        "-task", "blastn", "-evalue", "1e-10", "-max_target_seqs", "1000",
        "-num_threads", str(threads), "-outfmt", "6 " + " ".join(OUTFMT_FIELDS),
        "-out", str(output_file),
    ]
    subprocess.run(command, check=True)


def summarize_blast(path, contig_lengths):
    hits = defaultdict(list)
    with open(path, encoding="utf-8") as handle:
        reader = csv.DictReader(handle, fieldnames=OUTFMT_FIELDS, delimiter="\t")
        for row in reader:
            hits[row["qseqid"]].append(row)
    summaries = {}
    for contig, rows in hits.items():
        intervals = [(min(int(r["qstart"]), int(r["qend"])), max(int(r["qstart"]), int(r["qend"]))) for r in rows]
        aligned = merge_length(intervals)
        hsp_bases = sum(int(r["length"]) for r in rows)
        identity = sum(float(r["pident"]) * int(r["length"]) for r in rows) / hsp_bases
        coverage = 100.0 * aligned / contig_lengths[contig]
        summaries[contig] = {
            "aligned_bp": aligned,
            "query_coverage_pct": coverage,
            "identity_pct": identity,
            "bitscore": sum(float(r["bitscore"]) for r in rows),
            "evalue": min(float(r["evalue"]) for r in rows),
            "hsp_count": len(rows),
            "organelle_score": organelle_score(coverage, identity),
        }
    return summaries


def write_matrix(path, assembly, contigs, labels, reference_types, results):
    metrics = ("aligned_bp", "query_coverage_pct", "identity_pct", "bitscore", "evalue", "hsp_count", "organelle_score")
    fields = ["assembly", "contig", "contig_length", "best_reference", "best_reference_type", "organelle_score"]
    fields.extend(f"{label}__{metric}" for label in labels for metric in metrics)
    with open(path, "w", newline="", encoding="utf-8-sig") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        for contig, length in contigs.items():
            candidates = [(results[label].get(contig, {}).get("organelle_score", 0.0), label) for label in labels]
            best_score, best_label = max(candidates)
            row = {
                "assembly": os.path.abspath(assembly), "contig": contig, "contig_length": length,
                "best_reference": best_label if best_score else "",
                "best_reference_type": reference_types[best_label] if best_score else "",
                "organelle_score": f"{best_score:.6f}",
            }
            for label in labels:
                summary = results[label].get(contig, {})
                for metric in metrics:
                    value = summary.get(metric, 0)
                    if metric in {"query_coverage_pct", "identity_pct", "bitscore"}:
                        value = f"{value:.4f}"
                    elif metric == "organelle_score":
                        value = f"{value:.6f}"
                    row[f"{label}__{metric}"] = value
            writer.writerow(row)


def main():
    print("=" * 72)
    print("细胞器 contig 批量筛查（BLAST+）")
    print("=" * 72)
    missing = [name for name in ("makeblastdb", "blastn") if shutil.which(name) is None]
    if missing:
        raise SystemExit(f"找不到命令：{', '.join(missing)}。请先安装 NCBI BLAST+ 并加入 PATH。")
    files = scan_fastas()
    if len(files) < 2:
        raise SystemExit("当前目录及子目录中至少需要两个 FASTA 文件。")
    try:
        references = select_files(files, "请选择一个或多个细胞器参考基因组：")
        assemblies = select_files(files, "请选择一个或多个组装结果：", excluded=references)
        reference_labels = unique_labels(references)
        reference_types = {label: classify_reference(path) for path, label in zip(references, reference_labels)}
        threads = int(ask_text("BLAST 使用的线程数", DEFAULT_THREADS))
        if threads < 1:
            raise ValueError("线程数必须大于 0。")
        default_output = f"organelle_screen_{datetime.now().strftime('%Y%m%d_%H%M%S')}"
        output_dir = Path(ask_text("输出文件夹", default_output)).resolve()
    except (ValueError, EOFError) as error:
        raise SystemExit(f"参数设置失败：{error}") from error
    if output_dir.exists() and any(output_dir.iterdir()):
        raise SystemExit(f"输出文件夹已存在且非空：{output_dir}")
    print("\n运行配置：")
    print(f"  参考基因组：{len(references)} 个")
    print(f"  组装结果：{len(assemblies)} 个")
    print(f"  输出文件夹：{output_dir}")
    if not ask_yes_no("确认先建库、再开始批量比对", True):
        return
    database_dir = output_dir / "00_blast_databases"
    database_dir.mkdir(parents=True, exist_ok=True)
    databases = build_databases(references, reference_labels, database_dir, output_dir / "00_staged_inputs" / "references")
    assembly_labels = unique_labels(assemblies)
    for index, (assembly, assembly_label) in enumerate(zip(assemblies, assembly_labels), 1):
        log(f"处理组装 [{index}/{len(assemblies)}]：{assembly_label}")
        sample_dir = output_dir / assembly_label
        raw_dir = sample_dir / "01_blast_raw"
        result_dir = sample_dir / "02_results"
        raw_dir.mkdir(parents=True, exist_ok=True)
        result_dir.mkdir(parents=True, exist_ok=True)
        contigs = fasta_lengths(assembly)
        prepared_assembly = prepare_fasta(assembly, sample_dir / "00_staged_input", assembly_label)
        results = {}
        for reference_label, database in zip(reference_labels, databases):
            raw_file = raw_dir / f"{reference_label}.blast.tsv"
            log(f"  比对到：{reference_label}")
            run_blast(prepared_assembly, database, raw_file, threads)
            results[reference_label] = summarize_blast(raw_file, contigs)
        matrix = result_dir / f"{assembly_label}.organelle_matrix.tsv"
        write_matrix(matrix, assembly, contigs, reference_labels, reference_types, results)
        log(f"  矩阵完成：{matrix}")
    log(f"全部任务完成：{output_dir}")
    print("提示：organelle_score 是覆盖度和一致性的启发式评分，不是统计学概率；请结合 contig 长度、重复序列和生物学背景人工判断。")


if __name__ == "__main__":
    try:
        main()
    except subprocess.CalledProcessError as error:
        raise SystemExit(f"外部命令运行失败（退出码 {error.returncode}）：{' '.join(map(str, error.cmd))}") from error
    except KeyboardInterrupt:
        raise SystemExit("用户中断运行。")
