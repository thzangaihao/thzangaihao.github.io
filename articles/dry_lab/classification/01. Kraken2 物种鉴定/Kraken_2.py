#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Kraken2 交互式自建库 + 物种鉴定流水线

功能：
1. 交互式输入参考基因组 FASTA/FNA 路径，用于自建 Kraken2 数据库
2. 交互式输入 NCBI taxdump.tar.gz 路径，并自动解压到数据库 taxonomy 目录
3. 交互式选择 TaxID 解析方式：
   - 手动输入每个参考 FASTA 对应的 NCBI TaxID
   - 使用 nucl_gb.accession2taxid / nucl_wgs.accession2taxid 等映射文件
   - 直接使用已经包含 kraken:taxid|TaxID 的 FASTA header
4. 自动调用 kraken2-build 完成 add-to-library 和 build
5. 交互式输入待鉴定 fastq/fasta 路径，自动识别双端样本并输出 Kraken2 结果

运行：
    python Kraken_2.py

必要外部数据：
- taxdump.tar.gz：NCBI taxonomy 数据，必须包含 names.dmp 和 nodes.dmp
- 参考 FASTA/FNA：用于建库的基因组序列
- 待鉴定 fastq/fasta：需要做物种鉴定的测序数据或序列

说明：
如果选择手动 TaxID 模式，脚本会按你输入的 TaxID 改写 FASTA header：
    >contig_1 kraken:taxid|12345
如果选择 accession2taxid 模式，脚本会把映射文件复制或解压到数据库 taxonomy 目录。
"""

from __future__ import annotations

import gzip
import os
import re
import shutil
import subprocess
import sys
import tarfile
import time
from pathlib import Path


FASTA_SUFFIXES = {
    ".fa",
    ".fasta",
    ".fna",
    ".ffn",
    ".fas",
    ".fa.gz",
    ".fasta.gz",
    ".fna.gz",
    ".ffn.gz",
}
QUERY_SUFFIXES = FASTA_SUFFIXES | {
    ".fq",
    ".fastq",
    ".fq.gz",
    ".fastq.gz",
}
DB_REQUIRED_FILES = ("hash.k2d", "opts.k2d", "taxo.k2d")


def base_dir() -> Path:
    if getattr(sys, "frozen", False):
        return Path(sys.executable).resolve().parent
    return Path(__file__).resolve().parent


def print_header(title: str) -> None:
    print("\n" + "=" * 70)
    print(title)
    print("=" * 70)


def ask(prompt: str, default: str | None = None) -> str:
    suffix = f" [{default}]" if default else ""
    value = input(f"{prompt}{suffix}: ").strip()
    return value or (default or "")


def ask_yes_no(prompt: str, default: bool = True) -> bool:
    default_text = "Y/n" if default else "y/N"
    while True:
        value = input(f"{prompt} ({default_text}): ").strip().lower()
        if not value:
            return default
        if value in {"y", "yes"}:
            return True
        if value in {"n", "no"}:
            return False
        print("请输入 y 或 n。")


def ask_int(prompt: str, default: int) -> int:
    while True:
        value = ask(prompt, str(default))
        try:
            number = int(value)
        except ValueError:
            print("请输入整数。")
            continue
        if number > 0:
            return number
        print("数值必须大于 0。")


def ask_taxid(prompt: str, default: str | None = None) -> str:
    while True:
        value = ask(prompt, default).strip()
        if re.fullmatch(r"\d+", value):
            return value
        print("TaxID 必须是纯数字，例如 562。")


def path_suffix(path: Path) -> str:
    name = path.name.lower()
    for suffix in sorted(QUERY_SUFFIXES | {".tar.gz", ".gz"}, key=len, reverse=True):
        if name.endswith(suffix):
            return suffix
    return path.suffix.lower()


def is_fasta(path: Path) -> bool:
    return path_suffix(path) in FASTA_SUFFIXES


def is_query_file(path: Path) -> bool:
    return path_suffix(path) in QUERY_SUFFIXES


def glob_like(pattern: str) -> list[str]:
    import glob

    return glob.glob(os.path.expanduser(os.path.expandvars(pattern)), recursive=True)


def expand_paths(text: str, allowed_checker, recursive: bool = True) -> list[Path]:
    """
    支持输入单个文件、逗号分隔的多个路径、目录、glob 通配符。
    目录会自动递归搜索。
    """
    results: list[Path] = []
    for raw_item in re.split(r"[,;]", text):
        item = raw_item.strip().strip('"').strip("'")
        if not item:
            continue

        expanded = Path(os.path.expanduser(os.path.expandvars(item)))
        if any(ch in item for ch in "*?[]"):
            matched = [Path(p) for p in glob_like(item)]
        elif expanded.is_dir():
            pattern = "**/*" if recursive else "*"
            matched = [p for p in expanded.glob(pattern) if p.is_file()]
        else:
            matched = [expanded]

        for path in matched:
            resolved = path.resolve()
            if resolved.is_file() and allowed_checker(resolved):
                results.append(resolved)

    return sorted(set(results))


def require_tools() -> None:
    missing = [tool for tool in ("kraken2", "kraken2-build") if shutil.which(tool) is None]
    if missing:
        print("错误：未找到以下命令：", ", ".join(missing))
        print("请先安装 Kraken2，并激活包含 kraken2/kraken2-build 的 conda 环境。")
        sys.exit(1)


def open_text(path: Path):
    if path.name.lower().endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8", errors="ignore")
    return path.open("r", encoding="utf-8", errors="ignore")


def ensure_taxdump(taxdump_path: Path, db_dir: Path) -> None:
    if not taxdump_path.is_file():
        raise FileNotFoundError(f"taxdump.tar.gz 不存在：{taxdump_path}")
    if not taxdump_path.name.endswith(".tar.gz"):
        raise ValueError("taxonomy 文件应为 NCBI taxdump.tar.gz。")

    taxonomy_dir = db_dir / "taxonomy"
    taxonomy_dir.mkdir(parents=True, exist_ok=True)

    print(f"\n正在解压 taxonomy 到：{taxonomy_dir}")
    with tarfile.open(taxdump_path, "r:gz") as tar:
        wanted = {"names.dmp", "nodes.dmp", "merged.dmp", "delnodes.dmp"}
        members = [m for m in tar.getmembers() if Path(m.name).name in wanted]
        tar.extractall(taxonomy_dir, members=members)

    missing = [name for name in ("names.dmp", "nodes.dmp") if not (taxonomy_dir / name).exists()]
    if missing:
        raise RuntimeError(f"taxdump 解压后缺少必要文件：{', '.join(missing)}")


def validate_taxids_exist(taxids: set[str], db_dir: Path) -> None:
    nodes_file = db_dir / "taxonomy" / "nodes.dmp"
    found: set[str] = set()
    with nodes_file.open("r", encoding="utf-8", errors="ignore") as handle:
        for line in handle:
            taxid = line.split("|", 1)[0].strip()
            if taxid in taxids:
                found.add(taxid)
                if found == taxids:
                    break

    missing = sorted(taxids - found, key=int)
    if missing:
        raise ValueError(
            "以下 TaxID 不存在于 taxdump 的 nodes.dmp 中："
            + ", ".join(missing)
            + "。请检查 TaxID 或更新 taxdump.tar.gz。"
        )


def fasta_has_kraken_taxid(fasta: Path) -> bool:
    try:
        with open_text(fasta) as handle:
            for line in handle:
                if line.startswith(">"):
                    return "kraken:taxid|" in line
    except OSError:
        return False
    return False


def sanitize_filename(name: str) -> str:
    name = re.sub(r"\.(fa|fasta|fna|ffn|fas)(\.gz)?$", "", name, flags=re.IGNORECASE)
    name = re.sub(r"[^A-Za-z0-9._-]+", "_", name)
    return name.strip("._-") or "reference"


def rewrite_fasta_with_taxid(source: Path, target: Path, taxid: str) -> None:
    target.parent.mkdir(parents=True, exist_ok=True)
    with open_text(source) as inp, target.open("w", encoding="utf-8", newline="\n") as out:
        for line in inp:
            if line.startswith(">"):
                header = line.rstrip("\r\n")
                if "kraken:taxid|" in header:
                    header = re.sub(r"kraken:taxid\|\d+", f"kraken:taxid|{taxid}", header)
                else:
                    header = f"{header} kraken:taxid|{taxid}"
                out.write(header + "\n")
            else:
                out.write(line)


def prepare_taxid_fastas(fasta_taxids: dict[Path, str], db_dir: Path) -> list[Path]:
    prepared_dir = db_dir / "library_with_taxid"
    prepared_files: list[Path] = []
    used_names: set[str] = set()

    print_header("生成带 kraken:taxid 的参考 FASTA")
    for index, (source, taxid) in enumerate(fasta_taxids.items(), start=1):
        stem = sanitize_filename(source.name)
        target_name = f"{stem}.taxid_{taxid}.fasta"
        counter = 2
        while target_name in used_names:
            target_name = f"{stem}.{counter}.taxid_{taxid}.fasta"
            counter += 1
        used_names.add(target_name)

        target = prepared_dir / target_name
        print(f"[{index}/{len(fasta_taxids)}] {source.name} -> {target.name}")
        rewrite_fasta_with_taxid(source, target, taxid)
        prepared_files.append(target)

    return prepared_files


def looks_like_accession_map(path: Path) -> bool:
    name = path.name.lower()
    return "accession2taxid" in name


def copy_accession_maps(map_paths: list[Path], db_dir: Path) -> None:
    taxonomy_dir = db_dir / "taxonomy"
    taxonomy_dir.mkdir(parents=True, exist_ok=True)

    print_header("准备 accession2taxid 映射文件")
    for source in map_paths:
        if source.name.lower().endswith(".gz"):
            target = taxonomy_dir / source.name[:-3]
            print(f"解压：{source.name} -> {target.name}")
            with gzip.open(source, "rb") as inp, target.open("wb") as out:
                shutil.copyfileobj(inp, out)
        else:
            target = taxonomy_dir / source.name
            if source.resolve() == target.resolve():
                print(f"已在 taxonomy 目录中：{source.name}")
                continue
            print(f"复制：{source.name}")
            shutil.copy2(source, target)


def run_command(cmd: list[str], log_file: Path | None = None) -> None:
    print("\n$ " + " ".join(quote_arg(part) for part in cmd))
    if log_file:
        log_file.parent.mkdir(parents=True, exist_ok=True)
        with log_file.open("a", encoding="utf-8") as log:
            log.write("\n$ " + " ".join(cmd) + "\n")
            subprocess.run(cmd, check=True, stdout=log, stderr=subprocess.STDOUT)
    else:
        subprocess.run(cmd, check=True)


def quote_arg(value: str) -> str:
    return f'"{value}"' if any(ch.isspace() for ch in value) else value


def build_database(
    db_dir: Path,
    fasta_files: list[Path],
    taxdump_path: Path,
    threads: int,
    build_mode: str,
    fasta_taxids: dict[Path, str] | None = None,
    accession_maps: list[Path] | None = None,
) -> None:
    db_dir.mkdir(parents=True, exist_ok=True)
    log_file = db_dir / "build.log"

    ensure_taxdump(taxdump_path, db_dir)

    if build_mode == "manual_taxid":
        if not fasta_taxids:
            raise ValueError("手动 TaxID 模式需要 FASTA 与 TaxID 对应关系。")
        validate_taxids_exist(set(fasta_taxids.values()), db_dir)
        prepared_fastas = prepare_taxid_fastas(fasta_taxids, db_dir)
    elif build_mode == "accession2taxid":
        if not accession_maps:
            raise ValueError("accession2taxid 模式需要至少一个 accession2taxid 映射文件。")
        copy_accession_maps(accession_maps, db_dir)
        prepared_fastas = fasta_files
    elif build_mode == "header_taxid":
        missing_taxid = [path for path in fasta_files if not fasta_has_kraken_taxid(path)]
        if missing_taxid:
            names = ", ".join(path.name for path in missing_taxid[:10])
            raise ValueError(f"以下 FASTA 未检测到 kraken:taxid| 标记：{names}")
        prepared_fastas = fasta_files
    else:
        raise ValueError(f"未知建库模式：{build_mode}")

    print_header("添加参考 FASTA 到 Kraken2 library")
    for index, fasta in enumerate(prepared_fastas, start=1):
        print(f"[{index}/{len(prepared_fastas)}] {fasta.name}")
        run_command(
            ["kraken2-build", "--add-to-library", str(fasta), "--db", str(db_dir)],
            log_file=log_file,
        )

    print_header("构建 Kraken2 数据库")
    run_command(
        ["kraken2-build", "--build", "--threads", str(threads), "--db", str(db_dir)],
        log_file=log_file,
    )

    missing = [name for name in DB_REQUIRED_FILES if not (db_dir / name).exists()]
    if missing:
        raise RuntimeError(f"数据库构建完成后仍缺少：{', '.join(missing)}")


def collect_reference_fastas() -> list[Path]:
    while True:
        text = ask("\n请输入参考 FASTA/FNA 文件路径、目录或 glob，多个路径用逗号分隔")
        files = expand_paths(text, is_fasta)
        if files:
            print(f"找到 {len(files)} 个参考 FASTA/FNA：")
            for path in files[:20]:
                mark = "已有 kraken:taxid|" if fasta_has_kraken_taxid(path) else "需要输入 TaxID"
                print(f"  - {path} ({mark})")
            if len(files) > 20:
                print(f"  ... 还有 {len(files) - 20} 个文件")
            if ask_yes_no("确认使用这些参考序列建库吗", True):
                return files
        else:
            print("没有找到 FASTA/FNA 文件，请重新输入。")


def collect_taxdump() -> Path:
    while True:
        value = ask("\n请输入 NCBI taxdump.tar.gz 路径")
        path = Path(os.path.expanduser(os.path.expandvars(value))).resolve()
        if path.is_file() and path.name.endswith(".tar.gz"):
            return path
        print("路径无效，或文件名不是 taxdump.tar.gz。")


def collect_build_mode() -> str:
    print_header("选择参考序列 TaxID 解析方式")
    print("[1] 手动输入 TaxID，并由脚本生成 kraken:taxid| 标记 FASTA")
    print("[2] 使用 nucl_gb.accession2taxid / nucl_wgs.accession2taxid 等映射文件")
    print("[3] 参考 FASTA header 已经包含 kraken:taxid|TaxID")

    choices = {
        "1": "manual_taxid",
        "2": "accession2taxid",
        "3": "header_taxid",
    }
    while True:
        value = ask("请选择建库方式", "1")
        if value in choices:
            return choices[value]
        print("请输入 1、2 或 3。")


def collect_accession_maps() -> list[Path]:
    print_header("输入 accession2taxid 映射文件")
    print("可输入文件、目录或 glob；多个路径用逗号分隔。")
    print("常见文件包括 nucl_gb.accession2taxid、nucl_wgs.accession2taxid，也支持 .gz。")

    while True:
        text = ask("请输入 accession2taxid 路径")
        files = expand_paths(text, looks_like_accession_map)
        if files:
            print(f"找到 {len(files)} 个 accession2taxid 文件：")
            for path in files[:20]:
                print(f"  - {path}")
            if len(files) > 20:
                print(f"  ... 还有 {len(files) - 20} 个文件")
            if ask_yes_no("确认使用这些映射文件吗", True):
                return files
        else:
            print("没有找到 accession2taxid 文件，请重新输入。")


def collect_reference_taxids(fasta_files: list[Path]) -> dict[Path, str]:
    print_header("设置参考序列 TaxID")
    print("TaxID 必须来自 NCBI taxonomy，并存在于你提供的 taxdump.tar.gz 中。")
    print("例如 Escherichia coli 的 TaxID 是 562。")

    fasta_taxids: dict[Path, str] = {}
    if len(fasta_files) > 1 and ask_yes_no("这些参考 FASTA 是否都属于同一个 TaxID", False):
        taxid = ask_taxid("请输入统一使用的 TaxID")
        return {path: taxid for path in fasta_files}

    last_taxid: str | None = None
    for index, fasta in enumerate(fasta_files, start=1):
        print(f"\n[{index}/{len(fasta_files)}] {fasta.name}")
        if last_taxid:
            print(f"上一个 TaxID：{last_taxid}")
        taxid = ask_taxid("请输入该 FASTA 对应的 TaxID", last_taxid)
        fasta_taxids[fasta] = taxid
        last_taxid = taxid

    return fasta_taxids


def collect_db_dir() -> Path:
    default = str(base_dir() / f"Custom_Kraken2_DB_{time.strftime('%Y%m%d_%H%M%S')}")
    value = ask("\n请输入数据库输出目录", default)
    return Path(os.path.expanduser(os.path.expandvars(value))).resolve()


def is_complete_db(db_dir: Path) -> bool:
    return db_dir.is_dir() and all((db_dir / name).exists() for name in DB_REQUIRED_FILES)


def sample_name_from_read(path: Path) -> tuple[str, str | None]:
    name = path.name
    suffix = path_suffix(path)
    if suffix:
        name = name[: -len(suffix)]
    match = re.match(r"(.+?)(?:[._-]R?([12])(?:[._-]\d+)?)$", name, re.IGNORECASE)
    if match:
        return match.group(1), match.group(2)
    return name, None


def collect_query_samples() -> dict[str, dict[str, Path | None]]:
    while True:
        value = ask("\n请输入待鉴定 fastq/fasta 文件路径、目录或 glob，多个路径用逗号分隔")
        files = expand_paths(value, is_query_file)
        if not files:
            print("没有找到待鉴定序列文件，请重新输入。")
            continue

        samples: dict[str, dict[str, Path | None]] = {}
        for file_path in files:
            sample, read = sample_name_from_read(file_path)
            item = samples.setdefault(sample, {"r1": None, "r2": None, "single": None})
            if read == "1":
                item["r1"] = file_path
            elif read == "2":
                item["r2"] = file_path
            else:
                item["single"] = file_path

        valid = {
            name: paths
            for name, paths in samples.items()
            if paths["single"] or (paths["r1"] and paths["r2"])
        }
        if valid:
            print(f"找到 {len(valid)} 个可分析样本：")
            for index, (name, paths) in enumerate(sorted(valid.items()), start=1):
                mode = "paired" if paths["r1"] and paths["r2"] else "single"
                print(f"  [{index}] {name} ({mode})")
            if ask_yes_no("确认分析这些样本吗", True):
                return dict(sorted(valid.items()))

        print("没有找到完整样本。双端数据需要同时有 R1/R2；单端 fasta/fastq 可直接分析。")


def run_kraken2(db_dir: Path, samples: dict[str, dict[str, Path | None]], threads: int) -> Path:
    timestamp = time.strftime("%Y%m%d_%H%M%S")
    out_dir = base_dir() / f"03_Kraken2_鉴定结果_{timestamp}"
    out_dir.mkdir(parents=True, exist_ok=True)

    print_header("运行 Kraken2 鉴定")
    print(f"结果目录：{out_dir}")

    for index, (sample_name, paths) in enumerate(samples.items(), start=1):
        print(f"\n[{index}/{len(samples)}] {sample_name}")
        report_file = out_dir / f"{sample_name}_report.txt"
        output_file = out_dir / f"{sample_name}_kraken.out"
        log_file = out_dir / f"{sample_name}.log"

        cmd = [
            "kraken2",
            "--db",
            str(db_dir),
            "--threads",
            str(threads),
            "--use-names",
            "--report",
            str(report_file),
        ]
        if paths["r1"] and paths["r2"]:
            cmd.extend(["--paired", str(paths["r1"]), str(paths["r2"])])
        else:
            cmd.append(str(paths["single"]))

        print("$ " + " ".join(quote_arg(part) for part in cmd) + f" > {quote_arg(str(output_file))}")
        with output_file.open("w", encoding="utf-8") as out, log_file.open("w", encoding="utf-8") as log:
            subprocess.run(cmd, check=True, stdout=out, stderr=log)
        print(f"完成：{report_file.name}")

    return out_dir


def main() -> None:
    print_header("Kraken2 交互式自建库 + 物种鉴定流水线")
    require_tools()

    if ask_yes_no("是否使用已有 Kraken2 数据库", False):
        db_dir = Path(ask("请输入已有 Kraken2 数据库目录")).resolve()
        if not is_complete_db(db_dir):
            print("错误：该目录不是完整 Kraken2 数据库，至少需要 hash.k2d、opts.k2d、taxo.k2d。")
            sys.exit(1)
    else:
        db_dir = collect_db_dir()
        fasta_files = collect_reference_fastas()
        taxdump_path = collect_taxdump()
        build_mode = collect_build_mode()
        fasta_taxids = None
        accession_maps = None
        if build_mode == "manual_taxid":
            fasta_taxids = collect_reference_taxids(fasta_files)
        elif build_mode == "accession2taxid":
            accession_maps = collect_accession_maps()
        threads = ask_int("\n请输入建库线程数", max(1, os.cpu_count() or 1))
        build_database(
            db_dir=db_dir,
            fasta_files=fasta_files,
            taxdump_path=taxdump_path,
            threads=threads,
            build_mode=build_mode,
            fasta_taxids=fasta_taxids,
            accession_maps=accession_maps,
        )

    samples = collect_query_samples()
    threads = ask_int("\n请输入鉴定线程数", max(1, os.cpu_count() or 1))

    print_header("最终确认")
    print(f"Kraken2 数据库：{db_dir}")
    print(f"待分析样本数：{len(samples)}")
    for name in samples:
        print(f"  - {name}")
    if not ask_yes_no("开始运行鉴定吗", True):
        print("已取消。")
        return

    out_dir = run_kraken2(db_dir, samples, threads)
    print_header("全部完成")
    print(f"Kraken2 结果目录：{out_dir}")
    print("主要输出：*_report.txt 为分类报告，*_kraken.out 为每条序列的分类结果。")


if __name__ == "__main__":
    try:
        main()
    except KeyboardInterrupt:
        print("\n用户中断，已退出。")
        sys.exit(130)
    except subprocess.CalledProcessError as exc:
        print(f"\n命令运行失败，退出码：{exc.returncode}")
        print("请查看数据库目录 build.log 或结果目录中的 *.log 获取 Kraken2 原始报错。")
        sys.exit(exc.returncode)
    except Exception as exc:
        print(f"\n错误：{exc}")
        sys.exit(1)
