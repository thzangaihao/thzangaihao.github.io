#!/usr/bin/env python3
"""交互式选择 BAM/FASTQ，并用 minimap2 回比对后合并为一个排序 BAM。"""

from __future__ import annotations

import argparse
import os
import re
import shutil
import subprocess
import tempfile
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path


SEQUENCING_PRESETS = {
    "hifi": "map-hifi",
    "ont": "map-ont",
    "clr": "map-pb",
    "sr": "sr",
}
FASTA_SUFFIXES = {".fa", ".fasta", ".fna", ".fas"}
FASTQ_SUFFIXES = (".fastq.gz", ".fq.gz", ".fastq", ".fq")


@dataclass(frozen=True)
class InputUnit:
    """一个可独立比对的输入：BAM、单端 FASTQ 或双端 FASTQ。"""

    kind: str
    files: tuple[Path, ...]
    label: str


def log_info(message: str) -> None:
    print(f"[{datetime.now().strftime('%H:%M:%S')}] {message}", flush=True)


def fail(message: str) -> "None":
    raise SystemExit(f"错误：{message}")


def file_suffix(path: Path) -> str:
    suffixes = [suffix.lower() for suffix in path.suffixes]
    if suffixes and suffixes[-1] in {".gz", ".bgz"} and len(suffixes) > 1:
        return suffixes[-2]
    return suffixes[-1] if suffixes else ""


def strip_fastq_suffix(name: str) -> str:
    lower_name = name.lower()
    for suffix in FASTQ_SUFFIXES:
        if lower_name.endswith(suffix):
            return name[: -len(suffix)]
    return Path(name).stem


def sample_stem(path: Path) -> str:
    name = path.name
    for suffix in (".fasta.gz", ".fna.gz", ".fa.gz", ".fasta", ".fna", ".fa", ".fas"):
        if name.lower().endswith(suffix):
            return name[: -len(suffix)]
    return path.stem


def pair_signature(path: Path) -> tuple[str, int] | None:
    """返回去掉 R1/R2 标志后的配对键及 read 方向。"""
    stem = strip_fastq_suffix(path.name)
    patterns = (
        r"^(.*?)(?:[._-]R)([12])(?:[._-]?00[12])?$",
        r"^(.*?)(?:[._-])([12])$",
    )
    for pattern in patterns:
        match = re.match(pattern, stem, flags=re.IGNORECASE)
        if match and match.group(1):
            return match.group(1).lower(), int(match.group(2))
    return None


def scan_inputs(roots: list[Path]) -> tuple[list[Path], list[InputUnit], list[InputUnit]]:
    """递归扫描并将 FASTQ 自动分为双端和单端。"""
    seen: set[Path] = set()
    bam_files: list[Path] = []
    fastq_files: list[Path] = []

    for root in roots:
        if not root.is_dir():
            continue
        for path in root.rglob("*"):
            if not path.is_file():
                continue
            resolved = path.resolve()
            if resolved in seen:
                continue
            seen.add(resolved)
            lower_name = path.name.lower()
            if lower_name.endswith(".bam"):
                bam_files.append(resolved)
            elif lower_name.endswith(FASTQ_SUFFIXES):
                fastq_files.append(resolved)

    pair_groups: dict[tuple[Path, str], dict[int, list[Path]]] = {}
    unmatched: set[Path] = set(fastq_files)
    for path in fastq_files:
        signature = pair_signature(path)
        if signature:
            key, mate = signature
            pair_groups.setdefault((path.parent, key), {1: [], 2: []})[mate].append(path)

    paired_units: list[InputUnit] = []
    for (_, key), mates in sorted(pair_groups.items(), key=lambda item: str(item[0])):
        r1s, r2s = sorted(mates[1]), sorted(mates[2])
        for r1, r2 in zip(r1s, r2s):
            unmatched.discard(r1)
            unmatched.discard(r2)
            paired_units.append(InputUnit("paired", (r1, r2), key))

    single_units = [
        InputUnit("fastq", (path,), strip_fastq_suffix(path.name))
        for path in sorted(unmatched)
    ]
    return sorted(bam_files), paired_units, single_units


def relative_display(path: Path, roots: list[Path]) -> str:
    for root in roots:
        try:
            return str(path.relative_to(root))
        except ValueError:
            pass
    return str(path)


def parse_number_selection(text: str, maximum: int) -> list[int]:
    """解析 1,3-5 或 all/a 格式的多选输入。"""
    text = text.strip().lower()
    if text in {"a", "all", "*"}:
        return list(range(maximum))
    selected: set[int] = set()
    try:
        for part in text.replace("，", ",").split(","):
            part = part.strip()
            if not part:
                continue
            if "-" in part:
                start, end = (int(value) for value in part.split("-", 1))
                if start > end:
                    start, end = end, start
                selected.update(range(start - 1, end))
            else:
                selected.add(int(part) - 1)
    except ValueError:
        fail("选择格式无效，请输入如 1,3-5 或 all")
    if not selected or min(selected) < 0 or max(selected) >= maximum:
        fail(f"选择编号必须在 1-{maximum} 之间")
    return sorted(selected)


def interactive_input_units() -> list[InputUnit]:
    script_dir = Path(__file__).resolve().parent
    cwd = Path.cwd().resolve()
    roots = list(dict.fromkeys((script_dir, cwd)))
    print("\n正在递归搜索以下目录：")
    for root in roots:
        print(f"  - {root}")

    bam_files, paired_units, single_units = scan_inputs(roots)
    categories: list[tuple[str, list[InputUnit]]] = [
        ("BAM", [InputUnit("bam", (path,), path.stem) for path in bam_files]),
        ("单端/未配对 FASTQ", single_units),
        ("双端 FASTQ", paired_units),
    ]
    print("\n发现的输入类型：")
    for index, (name, units) in enumerate(categories, 1):
        print(f"  [{index}] {name}: {len(units)} 组")
    available = [index for index, (_, units) in enumerate(categories, 1) if units]
    if not available:
        fail("脚本同级目录及当前工作目录下未找到 BAM/FASTQ 文件")

    category_text = input(
        "选择输入类型（可多选，如 1,2,3；直接回车选择所有已发现类型）>>> "
    ).strip()
    category_indexes = (
        parse_number_selection(category_text, len(categories))
        if category_text
        else [index - 1 for index in available]
    )

    candidates: list[InputUnit] = []
    for category_index in category_indexes:
        candidates.extend(categories[category_index][1])
    if not candidates:
        fail("所选类型下没有可用文件")

    print("\n可选择的比对单元：")
    for index, unit in enumerate(candidates, 1):
        if unit.kind == "paired":
            r1 = relative_display(unit.files[0], roots)
            r2 = relative_display(unit.files[1], roots)
            print(f"  [{index}] [双端] {r1}  <->  {r2}")
        else:
            print(
                f"  [{index}] [{unit.kind.upper()}] "
                f"{relative_display(unit.files[0], roots)}"
            )
    selection = input(
        "选择一个或多个比对单元（如 1,3-5；all 表示全部）>>> "
    )
    return [candidates[index] for index in parse_number_selection(selection, len(candidates))]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="交互选择 BAM/FASTQ，回比对、合并并生成坐标排序 BAM。"
    )
    parser.add_argument("-r", "--reference", help="参考 FASTA 路径")
    parser.add_argument(
        "-i", "--reads", nargs="+",
        help="非交互输入文件；可混合多个 BAM/FASTQ，FASTQ 将自动识别双端",
    )
    parser.add_argument(
        "-x", "--technology", choices=SEQUENCING_PRESETS,
        help="测序类型：hifi、ont、clr 或 sr",
    )
    parser.add_argument("-t", "--threads", type=int, help="线程数")
    parser.add_argument("-o", "--output", help="输出 BAM 路径")
    parser.add_argument("--no-index", action="store_true", help="不创建 BAM 索引")
    return parser.parse_args()


def units_from_paths(values: list[str]) -> list[InputUnit]:
    paths = [Path(value).expanduser().resolve() for value in values]
    for path in paths:
        if not path.is_file():
            fail(f"输入文件不存在：{path}")
    unsupported = [
        path for path in paths
        if not path.name.lower().endswith((".bam",) + FASTQ_SUFFIXES)
    ]
    if unsupported:
        fail("不支持的输入格式：" + ", ".join(path.name for path in unsupported))

    bam_units = [InputUnit("bam", (path,), path.stem) for path in paths if path.suffix.lower() == ".bam"]
    fastqs = [path for path in paths if path.suffix.lower() != ".bam"]
    groups: dict[tuple[Path, str], dict[int, list[Path]]] = {}
    unmatched = set(fastqs)
    for path in fastqs:
        signature = pair_signature(path)
        if signature:
            key, mate = signature
            groups.setdefault((path.parent, key), {1: [], 2: []})[mate].append(path)
    paired: list[InputUnit] = []
    for (_, key), mates in sorted(groups.items(), key=lambda item: str(item[0])):
        for r1, r2 in zip(sorted(mates[1]), sorted(mates[2])):
            unmatched.discard(r1)
            unmatched.discard(r2)
            paired.append(InputUnit("paired", (r1, r2), key))
    singles = [InputUnit("fastq", (path,), strip_fastq_suffix(path.name)) for path in sorted(unmatched)]
    return bam_units + paired + singles


def complete_interactively(args: argparse.Namespace) -> tuple[argparse.Namespace, list[InputUnit]]:
    print("=" * 66)
    print("       BAM / FASTQ 多文件回比对、合并与自动排序工具")
    print("=" * 66)
    if not args.reference:
        args.reference = input("参考 FASTA 路径 >>> ").strip().strip('"').strip("'")
    units = units_from_paths(args.reads) if args.reads else interactive_input_units()
    if not args.technology:
        default = "sr" if any(unit.kind == "paired" for unit in units) else "hifi"
        value = input(
            f"测序类型 [hifi/ont/clr/sr，默认 {default}] >>> "
        ).strip().lower()
        args.technology = value or default
    if args.threads is None:
        default_threads = min(16, os.cpu_count() or 1)
        value = input(f"线程数 [默认 {default_threads}] >>> ").strip()
        try:
            args.threads = int(value) if value else default_threads
        except ValueError:
            fail("线程数必须是整数")
    return args, units


def validate(
    args: argparse.Namespace, units: list[InputUnit]
) -> tuple[Path, list[InputUnit], Path]:
    if args.technology not in SEQUENCING_PRESETS:
        fail("测序类型必须是 hifi、ont、clr 或 sr")
    if not args.threads or args.threads < 1:
        fail("线程数必须大于 0")
    if not units:
        fail("至少需要选择一个比对单元")
    if any(unit.kind == "paired" for unit in units) and args.technology != "sr":
        fail("双端 FASTQ 必须使用短读段模式（--technology sr）")

    reference = Path(args.reference).expanduser().resolve()
    if not reference.is_file():
        fail(f"参考序列不存在：{reference}")
    if file_suffix(reference) not in FASTA_SUFFIXES:
        fail("参考序列必须是 .fa/.fasta/.fna/.fas 格式（可 gzip 压缩）")

    output = (
        Path(args.output).expanduser()
        if args.output
        else Path.cwd() / f"{sample_stem(reference)}.sorted.bam"
    )
    if output.suffix.lower() != ".bam":
        output = output.with_suffix(".bam")
    output = output.resolve()
    input_paths = {path for unit in units for path in unit.files}
    if output in input_paths or output == reference:
        fail("输出路径不能覆盖输入文件")
    output.parent.mkdir(parents=True, exist_ok=True)
    return reference, units, output


def require_tools() -> None:
    missing = [name for name in ("minimap2", "samtools") if shutil.which(name) is None]
    if missing:
        fail("找不到外部程序：" + ", ".join(missing) + "。请先安装并加入 PATH。")


def run_mapping_unit(
    reference: Path,
    unit: InputUnit,
    temporary_bam: Path,
    args: argparse.Namespace,
) -> None:
    preset = SEQUENCING_PRESETS[args.technology]
    collator: subprocess.Popen[bytes] | None = None
    converter: subprocess.Popen[bytes] | None = None
    mapper_stdin = None
    if unit.kind == "bam":
        collate_cmd = [
            "samtools", "collate", "-u", "-O", "-@", str(args.threads),
            str(unit.files[0]),
        ]
        conversion_cmd = [
            "samtools", "fastq", "-@", str(args.threads), "-n", "-"
        ]
        log_info(f"按 read 名称整理并流式读取 BAM：{unit.files[0].name}")
        collator = subprocess.Popen(collate_cmd, stdout=subprocess.PIPE)
        converter = subprocess.Popen(
            conversion_cmd, stdin=collator.stdout, stdout=subprocess.PIPE
        )
        if collator.stdout:
            collator.stdout.close()
        mapper_stdin = converter.stdout
        mapper_reads = ["-"]
    else:
        mapper_reads = [str(path) for path in unit.files]

    map_cmd = [
        "minimap2", "-ax", preset, "-t", str(args.threads),
        str(reference), *mapper_reads,
    ]
    if unit.kind == "bam" and args.technology == "sr":
        map_cmd[1:1] = ["--frag=yes"]
    view_cmd = [
        "samtools", "view", "-@", str(args.threads), "-b",
        "-o", str(temporary_bam), "-",
    ]
    log_info(f"开始比对 [{unit.kind}] {unit.label}")
    mapper = subprocess.Popen(map_cmd, stdin=mapper_stdin, stdout=subprocess.PIPE)
    if converter and converter.stdout:
        converter.stdout.close()
    viewer = subprocess.Popen(view_cmd, stdin=mapper.stdout)
    if mapper.stdout:
        mapper.stdout.close()

    view_code = viewer.wait()
    map_code = mapper.wait()
    convert_code = converter.wait() if converter else 0
    collate_code = collator.wait() if collator else 0
    codes = {
        "samtools collate": collate_code,
        "samtools fastq": convert_code,
        "minimap2": map_code,
        "samtools view": view_code,
    }
    failed = [f"{name}（退出码 {code}）" for name, code in codes.items() if code]
    if failed:
        fail(f"比对单元 {unit.label} 处理失败：" + "；".join(failed))


def merge_and_sort(temporary_bams: list[Path], output: Path, threads: int) -> None:
    sort_cmd = ["samtools", "sort", "-@", str(threads), "-o", str(output), "-"]
    if len(temporary_bams) == 1:
        command = [
            "samtools", "sort", "-@", str(threads), "-o", str(output),
            str(temporary_bams[0]),
        ]
        result = subprocess.run(command, check=False)
        if result.returncode:
            fail(f"BAM 排序失败（退出码 {result.returncode}）")
        return

    merge_cmd = [
        "samtools", "merge", "-u", "-@", str(threads), "-",
        *(str(path) for path in temporary_bams),
    ]
    merger = subprocess.Popen(merge_cmd, stdout=subprocess.PIPE)
    sorter = subprocess.Popen(sort_cmd, stdin=merger.stdout)
    if merger.stdout:
        merger.stdout.close()
    sort_code = sorter.wait()
    merge_code = merger.wait()
    if merge_code or sort_code:
        fail(
            f"BAM 合并/排序失败（samtools merge={merge_code}, "
            f"samtools sort={sort_code}）"
        )


def main() -> None:
    args, units = complete_interactively(parse_args())
    reference, units, output = validate(args, units)
    require_tools()

    print("\n已选择的比对单元：")
    for index, unit in enumerate(units, 1):
        print(f"  {index}. [{unit.kind}] " + " <-> ".join(str(path) for path in unit.files))
    print(f"最终输出：{output}\n")

    with tempfile.TemporaryDirectory(prefix="minimap_align_", dir=output.parent) as temp_dir:
        temporary_bams: list[Path] = []
        for index, unit in enumerate(units, 1):
            temporary_bam = Path(temp_dir) / f"unit_{index:04d}.bam"
            run_mapping_unit(reference, unit, temporary_bam, args)
            temporary_bams.append(temporary_bam)
        log_info(f"正在合并 {len(temporary_bams)} 个比对单元并按坐标排序...")
        merge_and_sort(temporary_bams, output, args.threads)

    log_info(f"排序 BAM 已生成：{output}")
    if not args.no_index:
        log_info("正在创建 BAM 索引...")
        result = subprocess.run(
            ["samtools", "index", "-@", str(args.threads), str(output)],
            check=False,
        )
        if result.returncode:
            fail(f"BAM 已生成，但索引创建失败（退出码 {result.returncode}）")
        log_info(f"索引已生成：{output}.bai")

    print("-" * 66)
    log_info("完成。中间 FASTQ/SAM/BAM 已自动清理。")
    print(f"参考序列：{reference}")
    print(f"比对结果：{output}")


if __name__ == "__main__":
    main()
