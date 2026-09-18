#!/usr/bin/env python3
"""Build BLAST databases from FASTA files, individually or combined by type."""

from __future__ import annotations

import argparse
import shutil
import subprocess
import sys
import tempfile
import time
from pathlib import Path


NUCLEOTIDE_EXTENSIONS = {".fna", ".ffn", ".frn"}
PROTEIN_EXTENSIONS = {".faa", ".pep"}
FASTA_EXTENSIONS = NUCLEOTIDE_EXTENSIONS | PROTEIN_EXTENSIONS | {
    ".fa", ".fasta", ".fas", ".fsa"
}
NUCLEOTIDE_LETTERS = set("ACGTUNRYKMSWBDHV.-*")


def detect_dbtype(path: Path) -> str:
    if path.suffix.lower() in NUCLEOTIDE_EXTENSIONS:
        return "nucl"
    if path.suffix.lower() in PROTEIN_EXTENSIONS:
        return "prot"
    with path.open("r", encoding="utf-8-sig") as handle:
        for line in handle:
            if line.startswith(">") or not line.strip():
                continue
            if any(char.upper() not in NUCLEOTIDE_LETTERS for char in line if not char.isspace()):
                return "prot"
    return "nucl"


def combine_fasta(paths: list[Path], destination: Path, dbtype: str) -> int:
    """Stream all FASTA records into one input for makeblastdb."""
    records = 0
    with destination.open("w", encoding="utf-8", newline="\n") as output:
        for file_index, path in enumerate(paths, 1):
            print(f"[{dbtype}] 合并文件 {file_index}/{len(paths)}：{path.name}", flush=True)
            has_header = False
            bytes_read = 0
            last_update = time.monotonic()
            with path.open("r", encoding="utf-8-sig") as source:
                for line_number, raw_line in enumerate(source, 1):
                    bytes_read += len(raw_line)
                    line = raw_line.strip()
                    if not line:
                        continue
                    if line.startswith(">"):
                        if not line[1:].strip():
                            raise ValueError(f"{path.name} 第 {line_number} 行的序列头为空")
                        has_header = True
                        records += 1
                    elif not has_header:
                        raise ValueError(f"{path.name} 第 {line_number} 行的序列出现在序列头之前")
                    output.write(line + "\n")
                    if line_number % 100000 == 0 and time.monotonic() - last_update >= 5:
                        print(f"[{dbtype}] {path.name} 已读取约 {bytes_read / 1024**2:.1f} MiB，"
                              f"累计 {records:,} 条序列", flush=True)
                        last_update = time.monotonic()
            if not has_header:
                raise ValueError(f"{path.name} 中没有 FASTA 序列")
            print(f"[{dbtype}] 已合并 {path.name}；累计 {records:,} 条序列", flush=True)
    return records


def build_database(makeblastdb: str, paths: list[Path], dbtype: str,
                   prefix: Path) -> tuple[Path, int]:
    temporary_path: Path | None = None
    try:
        with tempfile.NamedTemporaryFile(
            mode="w", suffix=".fasta", prefix=f".{prefix.name}.",
            dir=prefix.parent, delete=False,
        ) as temporary:
            temporary_path = Path(temporary.name)
        records = combine_fasta(paths, temporary_path, dbtype)
        print(f"[{dbtype}] 合并完成，共 {records:,} 条序列；开始运行 makeblastdb。",
              flush=True)
        result = subprocess.run(
            [makeblastdb, "-in", str(temporary_path), "-dbtype", dbtype,
             "-out", str(prefix), "-title", prefix.name], check=False,
        )
        if result.returncode:
            raise RuntimeError(f"makeblastdb 退出码 {result.returncode}")
        return prefix, records
    finally:
        if temporary_path is not None:
            temporary_path.unlink(missing_ok=True)


def build_one_database(makeblastdb: str, path: Path, dbtype: str, prefix: Path) -> None:
    print(f"[{dbtype}] 开始建库：{path.name} -> {prefix}", flush=True)
    result = subprocess.run(
        [makeblastdb, "-in", str(path), "-dbtype", dbtype,
         "-out", str(prefix), "-title", path.name], check=False,
    )
    if result.returncode:
        raise RuntimeError(f"{path.name}：makeblastdb 退出码 {result.returncode}")


def main() -> int:
    parser = argparse.ArgumentParser(description="将 FASTA 逐个建库，或按核酸/蛋白类型合并建库。")
    parser.add_argument(
        "--input-dir", type=Path, default=Path(__file__).resolve().parent,
        help="输入目录；默认是本脚本所在目录",
    )
    parser.add_argument(
        "--output-dir", type=Path, help="输出目录；默认为输入目录下的 blastdb",
    )
    parser.add_argument(
        "--prefix", help="数据库前缀；一对一模式下用作各文件库名的前缀",
    )
    parser.add_argument(
        "--mode", choices=("one-to-one", "many-to-one"), default="many-to-one",
        help="建库模式；默认 many-to-one（同类型文件合并为一个库）",
    )
    parser.add_argument(
        "--dbtype", choices=("auto", "nucl", "prot"), default="auto",
        help="指定所有输入的序列类型；auto 按扩展名和内容判断",
    )
    parser.add_argument("-y", "--yes", action="store_true", help="跳过建库前的确认")
    args = parser.parse_args()
    if args.prefix is not None and (
        not args.prefix.strip() or Path(args.prefix).name != args.prefix
        or args.prefix in {".", ".."}
    ):
        parser.error("--prefix 必须是非空文件名，不能包含目录路径")

    input_dir = args.input_dir.resolve()
    if not input_dir.is_dir():
        parser.exit(1, f"错误：输入目录不存在：{input_dir}\n")
    makeblastdb = shutil.which("makeblastdb")
    if makeblastdb is None:
        parser.exit(1, "错误：找不到 makeblastdb，请先安装 NCBI BLAST+ 并加入 PATH。\n")

    files = sorted(
        (path for path in input_dir.iterdir()
         if path.is_file() and path.suffix.lower() in FASTA_EXTENSIONS),
        key=lambda path: path.name.lower(),
    )
    if not files:
        parser.exit(1, f"错误：{input_dir} 下没有找到 FASTA 文件。\n")

    groups: dict[str, list[Path]] = {"nucl": [], "prot": []}
    try:
        for path in files:
            dbtype = args.dbtype if args.dbtype != "auto" else detect_dbtype(path)
            groups[dbtype].append(path)
        sizes = {path: path.stat().st_size for path in files}
    except (OSError, UnicodeError) as exc:
        parser.exit(1, f"错误：读取 FASTA 文件失败：{exc}\n")

    output_dir = (args.output_dir or input_dir / "blastdb").resolve()
    print("找到以下 FASTA 文件：")
    for index, path in enumerate(files, 1):
        dbtype = "nucl" if path in groups["nucl"] else "prot"
        print(f"  {index}. {path.name}  [{dbtype}]  {sizes[path]:,} 字节")
    print(f"统计：共 {len(files)} 个文件；核酸 {len(groups['nucl'])} 个，"
          f"蛋白 {len(groups['prot'])} 个；输入总大小 {sum(sizes.values()):,} 字节。")
    active_groups = [(dbtype, paths) for dbtype, paths in groups.items() if paths]
    if args.mode == "one-to-one":
        jobs = [
            (path, dbtype, output_dir / (f"{args.prefix}_{path.name}" if args.prefix else path.name))
            for dbtype, paths in active_groups for path in paths
        ]
        print(f"建库模式：一对一；计划生成 {len(jobs)} 个数据库。")
        for path, _, prefix in jobs:
            print(f"计划生成：{path.name} -> {prefix}")
    else:
        prefixes = {
            dbtype: output_dir / (
                f"{args.prefix}_{dbtype}" if args.prefix and len(active_groups) > 1
                else args.prefix or f"combined_{dbtype}"
            )
            for dbtype, _ in active_groups
        }
        print(f"建库模式：多对一；计划生成 {len(active_groups)} 个数据库。")
        for dbtype, paths in active_groups:
            print(f"计划生成：{prefixes[dbtype]}（{len(paths)} 个输入文件）")

    if not args.yes:
        try:
            answer = input("确认开始建库？输入 y 或 yes 继续：").strip().lower()
        except EOFError:
            answer = ""
        if answer not in {"y", "yes"}:
            print("已取消建库。")
            return 0

    try:
        print("已确认，开始建库。", flush=True)
        output_dir.mkdir(parents=True, exist_ok=True)
        if args.mode == "one-to-one":
            for index, (path, dbtype, prefix) in enumerate(jobs, 1):
                print(f"进度：{index}/{len(jobs)}", flush=True)
                build_one_database(makeblastdb, path, dbtype, prefix)
                print(f"完成：{path.name}；BLAST 数据库前缀：{prefix}", flush=True)
        else:
            for dbtype, paths in active_groups:
                prefix, records = build_database(makeblastdb, paths, dbtype, prefixes[dbtype])
                print(f"完成：{dbtype} 总库，{records} 条序列；"
                      f"BLAST 数据库前缀：{prefix}", flush=True)
    except (OSError, UnicodeError, ValueError, RuntimeError) as exc:
        parser.exit(1, f"建库失败：{exc}\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())
