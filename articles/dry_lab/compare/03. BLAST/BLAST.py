#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""交互式 NCBI BLAST+ 全流程：递归选择 FASTA、建库、比对并输出 TSV。"""

from __future__ import annotations

import hashlib
import multiprocessing
import os
import re
import shutil
import subprocess
import sys
import time
from pathlib import Path


FASTA_SUFFIXES = {".fa", ".fasta", ".fna", ".ffn", ".faa", ".fas", ".pep"}
PROGRAMS = {
    "1": ("blastn", "nucl", "核酸查询 vs 核酸库"),
    "2": ("blastp", "prot", "蛋白查询 vs 蛋白库"),
    "3": ("blastx", "prot", "核酸查询（翻译）vs 蛋白库"),
    "4": ("tblastn", "nucl", "蛋白查询 vs 核酸库（翻译）"),
    "5": ("tblastx", "nucl", "核酸查询（翻译）vs 核酸库（翻译）"),
}
FIELDS = [
    "qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
    "qstart", "qend", "sstart", "send", "evalue", "bitscore",
    "qlen", "slen", "qcovs", "stitle",
]


def base_dir() -> Path:
    """始终以脚本所在目录为搜索起点，而非命令执行目录。"""
    return Path(sys.executable if getattr(sys, "frozen", False) else __file__).resolve().parent


def find_fastas(root: Path) -> list[Path]:
    return sorted(
        (p for p in root.rglob("*") if p.is_file() and p.suffix.lower() in FASTA_SUFFIXES),
        key=lambda p: str(p.relative_to(root)).lower(),
    )


def parse_selection(text: str, count: int, multiple: bool) -> list[int]:
    text = text.strip().lower()
    if text in {"q", "quit", "exit"}:
        raise KeyboardInterrupt
    if multiple and text in {"all", "a"}:
        return list(range(count))
    selected: set[int] = set()
    for item in text.split(","):
        item = item.strip()
        if not item:
            raise ValueError("存在空的选择项")
        if "-" in item:
            if not multiple or item.count("-") != 1:
                raise ValueError("单选时不能使用范围")
            start, end = (int(x.strip()) for x in item.split("-"))
            if start > end:
                raise ValueError("范围起始编号不能大于结束编号")
            selected.update(range(start - 1, end))
        else:
            selected.add(int(item) - 1)
    if not selected or min(selected) < 0 or max(selected) >= count:
        raise ValueError(f"编号必须在 1-{count} 之间")
    if not multiple and len(selected) != 1:
        raise ValueError("这里只能选择一个文件")
    return sorted(selected)


def choose_files(files: list[Path], root: Path, description: str, multiple: bool) -> list[Path]:
    if not files:
        raise RuntimeError(f"脚本同级目录及子目录下没有找到{description}")
    print(f"\n找到 {len(files)} 个可用 FASTA（{description}）：")
    for number, path in enumerate(files, 1):
        print(f"  [{number:>3}] {path.relative_to(root)}")
    hint = "1,3,5-8 或 all" if multiple else "一个编号"
    while True:
        try:
            chosen = parse_selection(input(f"请选择{description}（{hint}，q 退出）："), len(files), multiple)
            return [files[i] for i in chosen]
        except ValueError as error:
            print(f"输入无效：{error}")


def validate_fasta(path: Path) -> int:
    """做基础格式检查并返回记录数；不依赖 Biopython。"""
    records = 0
    has_sequence = False
    try:
        with path.open("r", encoding="utf-8", errors="replace") as handle:
            for line_no, raw in enumerate(handle, 1):
                line = raw.strip()
                if not line:
                    continue
                if line.startswith(">"):
                    if len(line) == 1:
                        raise ValueError(f"第 {line_no} 行的 FASTA 标题为空")
                    if records and not has_sequence:
                        raise ValueError(f"第 {line_no} 行之前的记录没有序列")
                    records += 1
                    has_sequence = False
                else:
                    if records == 0:
                        raise ValueError(f"第 {line_no} 行在第一个 '>' 标题之前出现序列")
                    compact = re.sub(r"\s+", "", line)
                    if not compact or not re.fullmatch(r"[A-Za-z*.-]+", compact):
                        raise ValueError(f"第 {line_no} 行包含非法字符")
                    has_sequence = True
    except OSError as error:
        raise ValueError(f"无法读取：{error}") from error
    if records == 0:
        raise ValueError("未发现以 '>' 开头的 FASTA 记录")
    if not has_sequence:
        raise ValueError("最后一条记录没有序列")
    return records


def ask_program() -> tuple[str, str, str]:
    print("\n请选择分析类型：")
    for key, (program, _, description) in PROGRAMS.items():
        print(f"  [{key}] {program:<8} {description}")
    while True:
        choice = input("请输入编号（默认 1/blastn）：").strip() or "1"
        if choice in PROGRAMS:
            return PROGRAMS[choice]
        print("输入无效，请输入 1-5。")


def ask_number(prompt: str, default: str, kind: type, minimum: float = 0):
    while True:
        raw = input(f"{prompt}（默认 {default}）：").strip() or default
        try:
            value = kind(raw)
            if value <= minimum:
                raise ValueError
            return value
        except ValueError:
            print(f"请输入大于 {minimum:g} 的有效数值。")


def safe_stem(path: Path) -> str:
    stem = re.sub(r"[^A-Za-z0-9._-]+", "_", path.stem).strip("._-") or "database"
    digest = hashlib.sha1(str(path.resolve()).encode("utf-8")).hexdigest()[:8]
    return f"{stem}_{digest}"


def run(command: list[str], description: str, log_handle) -> None:
    print(f"\n[运行] {description}")
    print("       " + " ".join(command))
    log_handle.write("\n$ " + " ".join(command) + "\n")
    log_handle.flush()
    subprocess.run(command, check=True, stdout=log_handle, stderr=subprocess.STDOUT)


def add_header(raw_path: Path, final_path: Path) -> None:
    with final_path.open("w", encoding="utf-8", newline="\n") as output:
        output.write("\t".join(FIELDS) + "\n")
        with raw_path.open("r", encoding="utf-8", errors="replace") as source:
            shutil.copyfileobj(source, output)
    raw_path.unlink()


def main() -> int:
    print("=" * 68)
    print(" NCBI BLAST+ 交互式建库与批量比对流程")
    print("=" * 68)
    root = base_dir()
    program, dbtype, description = ask_program()
    for executable in ("makeblastdb", program):
        if shutil.which(executable) is None:
            print(f"错误：PATH 中未找到 {executable}。可先安装：conda install -c bioconda blast")
            return 1

    files = find_fastas(root)
    try:
        query = choose_files(files, root, "查询序列文件（单选）", multiple=False)[0]
        databases = choose_files(files, root, "数据库源文件（可多选）", multiple=True)
    except (RuntimeError, KeyboardInterrupt) as error:
        print(f"操作结束{(': ' + str(error)) if str(error) else '。'}")
        return 1

    print("\n正在检查 FASTA 格式……")
    for path in [query, *databases]:
        try:
            count = validate_fasta(path)
            print(f"  通过：{path.relative_to(root)}（{count} 条序列）")
        except ValueError as error:
            print(f"错误：{path.relative_to(root)}：{error}")
            return 1

    default_threads = max(1, int(multiprocessing.cpu_count() * 0.8))
    threads = ask_number("线程数", str(default_threads), int)
    evalue = ask_number("E-value 阈值", "1e-5", float)
    max_targets = ask_number("每条查询最多保留的目标序列数", "10", int)
    out_name = input("结果目录名（将在脚本目录下创建，默认 blast_results）：").strip() or "blast_results"
    if Path(out_name).name != out_name or out_name in {".", ".."}:
        print("错误：结果目录名只能是单个目录名称，不能包含路径。")
        return 1
    output_dir = root / out_name
    database_dir = output_dir / "databases"

    print("\n任务确认：")
    print(f"  模式：{program}（{description}）")
    print(f"  查询：{query.relative_to(root)}")
    print(f"  数据库：{len(databases)} 个；线程：{threads}；E-value：{evalue:g}；max_target_seqs：{max_targets}")
    if input("确认开始？[Y/n]：").strip().lower() in {"n", "no"}:
        print("操作已取消。")
        return 0

    database_dir.mkdir(parents=True, exist_ok=True)
    log_path = output_dir / "blast_pipeline.log"
    successes: list[Path] = []
    failures: list[tuple[Path, str]] = []
    started = time.time()
    with log_path.open("a", encoding="utf-8") as log:
        log.write(f"\n=== Run {time.strftime('%Y-%m-%d %H:%M:%S')} | {program} ===\n")
        for index, fasta in enumerate(databases, 1):
            tag = safe_stem(fasta)
            db_prefix = database_dir / tag
            raw_result = output_dir / f"{query.stem}_vs_{tag}.raw.tsv"
            final_result = output_dir / f"{query.stem}_vs_{tag}.blast.tsv"
            print(f"\n{'-' * 68}\n[{index}/{len(databases)}] 数据库：{fasta.relative_to(root)}")
            try:
                run([
                    "makeblastdb", "-in", str(fasta), "-dbtype", dbtype,
                    "-parse_seqids", "-out", str(db_prefix),
                ], "建立 BLAST 数据库", log)
                run([
                    program, "-query", str(query), "-db", str(db_prefix),
                    "-out", str(raw_result), "-outfmt", "6 " + " ".join(FIELDS),
                    "-evalue", str(evalue), "-max_target_seqs", str(max_targets),
                    "-num_threads", str(threads),
                ], "运行序列比对", log)
                add_header(raw_result, final_result)
                successes.append(final_result)
                print(f"[完成] {final_result.relative_to(root)}")
            except (subprocess.CalledProcessError, OSError) as error:
                raw_result.unlink(missing_ok=True)
                failures.append((fasta, str(error)))
                print(f"[失败] 详见日志：{log_path.relative_to(root)}")

    print("\n" + "=" * 68)
    print(f"流程结束，耗时 {time.time() - started:.1f} 秒；成功 {len(successes)}，失败 {len(failures)}。")
    for path in successes:
        print(f"  结果：{path.relative_to(root)}")
    for path, error in failures:
        print(f"  失败：{path.relative_to(root)}（{error}）")
    print(f"  日志：{log_path.relative_to(root)}")
    return 1 if failures else 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except KeyboardInterrupt:
        print("\n用户中止。")
        raise SystemExit(130)
