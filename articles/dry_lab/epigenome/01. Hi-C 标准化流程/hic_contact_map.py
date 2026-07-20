#!/usr/bin/env python3
"""Interactive scaffold-level Hi-C contact-map pipeline.

Workflow
--------
paired FASTQ -> BWA-MEM -> pairtools parse/sort -> merge -> dedup/select
             -> cooler (.cool) -> balanced multi-resolution cooler (.mcool)

The script is standalone, does not modify input files, and resumes completed
steps when the same output directory is reused.
"""

from __future__ import annotations

import argparse
import os
import re
import shlex
import shutil
import subprocess
import sys
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path


FASTQ_SUFFIXES = (".fastq.gz", ".fq.gz", ".fastq", ".fq")
PAIR_PATTERNS = (
    re.compile(r"^(.*?)([_\.-])R([12])([_\.-].*)?$", re.I),
    re.compile(r"^(.*?)([_\.-])READ([12])([_\.-].*)?$", re.I),
    re.compile(r"^(.*?)([_\.-])([12])([_\.-].*)?$", re.I),
)


@dataclass(frozen=True)
class ReadPair:
    label: str
    r1: Path
    r2: Path


def strip_fastq_suffix(name: str) -> str:
    lower = name.lower()
    for suffix in FASTQ_SUFFIXES:
        if lower.endswith(suffix):
            return name[: -len(suffix)]
    return name


def pair_key(path: Path) -> tuple[str, int] | None:
    stem = strip_fastq_suffix(path.name)
    for pattern in PAIR_PATTERNS:
        match = pattern.match(stem)
        if match:
            tail = match.group(4) or ""
            key = f"{match.group(1)}{match.group(2)}{{READ}}{tail}"
            return key.lower(), int(match.group(3))
    return None


def discover_pairs(root: Path) -> tuple[list[ReadPair], list[Path], list[str]]:
    files = sorted(
        path.resolve()
        for path in root.rglob("*")
        if path.is_file() and path.name.lower().endswith(FASTQ_SUFFIXES)
    )
    grouped: dict[tuple[Path, str], dict[int, list[Path]]] = {}
    unmatched: list[Path] = []
    for path in files:
        parsed = pair_key(path)
        if parsed is None:
            unmatched.append(path)
            continue
        key, read_number = parsed
        grouped.setdefault((path.parent, key), {1: [], 2: []})[read_number].append(path)

    pairs: list[ReadPair] = []
    problems: list[str] = []
    for (parent, key), reads in sorted(grouped.items(), key=lambda item: str(item[0])):
        if len(reads[1]) == 1 and len(reads[2]) == 1:
            label = strip_fastq_suffix(reads[1][0].name)
            label = re.sub(
                r"(?i)([_\.-])(R|READ)?1(?=([_\.-]|$))", "", label, count=1
            )
            pairs.append(ReadPair(label, reads[1][0], reads[2][0]))
        else:
            problems.append(f"{parent / key}: R1={len(reads[1])}, R2={len(reads[2])}")
    return pairs, unmatched, problems


def choose_pairs(pairs: list[ReadPair]) -> list[ReadPair]:
    if not pairs:
        raise SystemExit("没有发现可配对 FASTQ；文件名需包含 R1/R2 或 _1/_2。")
    print("\n发现以下 FASTQ 配对：")
    for number, pair in enumerate(pairs, 1):
        print(f"  [{number}] {pair.label}")
        print(f"      R1: {pair.r1}")
        print(f"      R2: {pair.r2}")
    while True:
        answer = input("\n选择配对 [all，或 1,3-4；默认 all]: ").strip()
        if not answer or answer.lower() in {"all", "a", "全部"}:
            return pairs
        selected: set[int] = set()
        try:
            for token in answer.replace("，", ",").split(","):
                token = token.strip()
                if "-" in token:
                    start, end = map(int, token.split("-", 1))
                    selected.update(range(start, end + 1))
                else:
                    selected.add(int(token))
            if selected and min(selected) >= 1 and max(selected) <= len(pairs):
                return [pair for i, pair in enumerate(pairs, 1) if i in selected]
        except ValueError:
            pass
        print("输入无效，请重新输入。", file=sys.stderr)


def ask_file(prompt: str) -> Path:
    while True:
        value = input(f"{prompt}: ").strip().strip("'\"")
        path = Path(value).expanduser()
        if path.is_file():
            return path.resolve()
        print("文件不存在，请重新输入。", file=sys.stderr)


def safe_label(text: str, number: int) -> str:
    cleaned = re.sub(r"[^A-Za-z0-9_.-]+", "_", text).strip("_.-")
    return cleaned or f"lane{number}"


def require_tools(names: list[str]) -> None:
    missing = [name for name in names if shutil.which(name) is None]
    if missing:
        raise SystemExit(
            "缺少必要软件: " + ", ".join(missing)
            + "\n请先激活 hic-contact；仍缺失时安装："
            + "mamba install -n hic-contact -c conda-forge -c bioconda "
            + "pairtools cooler"
        )


def ready(path: Path) -> bool:
    return path.is_file() and path.stat().st_size > 0


def done(path: Path) -> bool:
    return ready(path)


def write_done(path: Path, dry_run: bool) -> None:
    if not dry_run:
        path.write_text("complete\n", encoding="utf-8")


class Runner:
    def __init__(self, log_path: Path, dry_run: bool):
        self.log_path = log_path
        self.dry_run = dry_run

    def log(self, message: str) -> None:
        print(message, flush=True)
        with self.log_path.open("a", encoding="utf-8") as handle:
            handle.write(message + "\n")

    def run(self, command: list[str], stdout_path: Path | None = None) -> None:
        shown = shlex.join([str(part) for part in command])
        if stdout_path:
            shown += f" > {shlex.quote(str(stdout_path))}"
        self.log(f"[{datetime.now().isoformat(timespec='seconds')}] $ {shown}")
        if self.dry_run:
            return
        with self.log_path.open("a", encoding="utf-8") as log:
            if stdout_path:
                partial = Path(str(stdout_path) + ".part")
                partial.unlink(missing_ok=True)
                try:
                    with partial.open("wb") as output:
                        subprocess.run(command, stdout=output, stderr=log, check=True)
                    partial.replace(stdout_path)
                except Exception:
                    partial.unlink(missing_ok=True)
                    raise
            else:
                subprocess.run(
                    command, stdout=log, stderr=subprocess.STDOUT, check=True
                )

    def pipe(self, commands: list[list[str]], expected_output: Path) -> None:
        shown = " | ".join(shlex.join([str(part) for part in cmd]) for cmd in commands)
        self.log(f"[{datetime.now().isoformat(timespec='seconds')}] $ {shown}")
        if self.dry_run:
            return
        processes: list[subprocess.Popen] = []
        with self.log_path.open("a", encoding="utf-8") as log:
            previous = None
            try:
                for index, command in enumerate(commands):
                    process = subprocess.Popen(
                        command,
                        stdin=previous.stdout if previous else None,
                        stdout=(
                            subprocess.PIPE
                            if index < len(commands) - 1
                            else subprocess.DEVNULL
                        ),
                        stderr=log,
                    )
                    if previous and previous.stdout:
                        previous.stdout.close()
                    processes.append(process)
                    previous = process
                codes = [process.wait() for process in processes]
                if any(code != 0 for code in codes):
                    expected_output.unlink(missing_ok=True)
                    raise subprocess.CalledProcessError(
                        next(code for code in codes if code != 0), shown
                    )
            except Exception:
                for process in processes:
                    if process.poll() is None:
                        process.terminate()
                raise


def prepare_reference(source: Path, destination: Path) -> None:
    if destination.exists():
        return
    try:
        destination.symlink_to(source)
    except OSError:
        print("无法创建 FASTA 软链接，改为复制参考序列。")
        shutil.copy2(source, destination)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="交互式 scaffold 层级 Hi-C 互作矩阵分析"
    )
    parser.add_argument("--search-dir", type=Path, default=Path.cwd())
    parser.add_argument("--fasta", type=Path, help="YaHS 最终 scaffold FASTA")
    parser.add_argument("--outdir", type=Path)
    parser.add_argument(
        "--threads", type=int, default=max(1, min(16, os.cpu_count() or 1))
    )
    parser.add_argument("--mapq", type=int, default=30, help="唯一比对最低 MAPQ")
    parser.add_argument("--bin-size", type=int, default=10000, help="基础分辨率/bp")
    parser.add_argument(
        "--zoom-resolutions",
        default="10000N",
        help="cooler zoomify 分辨率（默认 10000N）",
    )
    parser.add_argument("--assembly-name", default="hic_scaffold")
    parser.add_argument("--skip-fastqc", action="store_true")
    parser.add_argument("--no-balance", action="store_true")
    parser.add_argument("--dry-run", action="store_true")
    args = parser.parse_args()

    search_dir = args.search_dir.expanduser().resolve()
    if not search_dir.is_dir():
        raise SystemExit(f"FASTQ 搜索目录不存在: {search_dir}")
    pairs, unmatched, problems = discover_pairs(search_dir)
    if unmatched:
        print("\n警告：以下 FASTQ 无法识别配对，将跳过：", file=sys.stderr)
        for path in unmatched:
            print(f"  {path}", file=sys.stderr)
    if problems:
        print("\n警告：以下分组不是一对一配对，将跳过：", file=sys.stderr)
        for problem in problems:
            print(f"  {problem}", file=sys.stderr)
    selected = choose_pairs(pairs)

    fasta = args.fasta.expanduser().resolve() if args.fasta else ask_file(
        "参考基因组 FASTA 路径"
    )
    if not fasta.is_file():
        raise SystemExit(f"FASTA 不存在: {fasta}")
    outdir = args.outdir
    if outdir is None:
        outdir = Path(input("输出目录 [默认 ./hic_contact_result]: ").strip() or "hic_contact_result")
    outdir = outdir.expanduser().resolve()

    if sys.stdin.isatty() and not any(
        arg == "--threads" or arg.startswith("--threads=") for arg in sys.argv[1:]
    ):
        value = input(f"线程数 [默认 {args.threads}]: ").strip()
        if value:
            args.threads = int(value)
    if sys.stdin.isatty() and not any(
        arg == "--mapq" or arg.startswith("--mapq=") for arg in sys.argv[1:]
    ):
        value = input(f"最低 MAPQ [默认 {args.mapq}]: ").strip()
        if value:
            args.mapq = int(value)
    if sys.stdin.isatty() and not any(
        arg == "--bin-size" or arg.startswith("--bin-size=") for arg in sys.argv[1:]
    ):
        value = input(f"基础分辨率/bp [默认 {args.bin_size}]: ").strip()
        if value:
            args.bin_size = int(value)

    if args.threads < 1 or args.bin_size < 1 or args.mapq < 0:
        raise SystemExit("线程数和 bin size 必须大于 0，MAPQ 不能小于 0。")
    require_tools(
        ["bwa", "samtools", "pairtools", "cooler"]
        + ([] if args.skip_fastqc else ["fastqc", "multiqc"])
    )

    directories = {
        "reference": outdir / "reference",
        "qc": outdir / "qc",
        "pairs": outdir / "pairs",
        "cool": outdir / "cool",
        "tmp": outdir / "tmp",
    }
    for directory in (outdir, *directories.values()):
        directory.mkdir(parents=True, exist_ok=True)
    runner = Runner(outdir / "pipeline.log", args.dry_run)

    reference = directories["reference"] / "scaffolds.fa"
    prepare_reference(fasta, reference)
    fai = Path(str(reference) + ".fai")
    chromsizes = directories["reference"] / "scaffolds.chrom.sizes"

    print("\n分析配置：")
    print(f"  FASTQ 配对数: {len(selected)}")
    print(f"  参考基因组: {fasta}")
    print(f"  输出目录:       {outdir}")
    print(f"  线程数:         {args.threads}")
    print(f"  最低 MAPQ:      {args.mapq}")
    print(f"  基础分辨率:     {args.bin_size:,} bp")
    print(f"  多分辨率规则:   {args.zoom_resolutions}")
    if sys.stdin.isatty() and input("开始运行？[Y/n]: ").strip().lower() in {"n", "no"}:
        raise SystemExit("已取消。")

    if not ready(fai):
        runner.run(["samtools", "faidx", str(reference)])
    if not ready(chromsizes):
        runner.run(
            ["cut", "-f1,2", str(fai)], stdout_path=chromsizes
        )
    if not ready(Path(str(reference) + ".bwt")):
        runner.run(["bwa", "index", str(reference)])

    if not args.skip_fastqc:
        qc_done = directories["qc"] / ".multiqc_complete"
        if not done(qc_done):
            runner.run(
                [
                    "fastqc", "-t", str(args.threads),
                    "-o", str(directories["qc"]),
                    *[str(path) for pair in selected for path in (pair.r1, pair.r2)],
                ]
            )
            runner.run(
                [
                    "multiqc", "-f", "-o", str(directories["qc"]),
                    str(directories["qc"]),
                ]
            )
            write_done(qc_done, args.dry_run)

    lane_pairs: list[Path] = []
    for number, pair in enumerate(selected, 1):
        label = f"{number:02d}_{safe_label(pair.label, number)}"
        output = directories["pairs"] / f"{label}.sorted.pairs.gz"
        lane_pairs.append(output)
        lane_done = directories["pairs"] / f".{label}.complete"
        if done(lane_done) and ready(output):
            print(f"跳过已完成 lane: {output}")
            continue
        output.unlink(missing_ok=True)
        read_group = (
            f"@RG\\tID:{label}\\tSM:HIC\\tLB:HIC\\tPL:ILLUMINA\\tPU:{label}"
        )
        runner.pipe(
            [
                [
                    "bwa", "mem", "-5SP", "-t", str(args.threads),
                    "-R", read_group, str(reference), str(pair.r1), str(pair.r2),
                ],
                [
                    "pairtools", "parse",
                    "--chroms-path", str(chromsizes),
                    "--assembly", args.assembly_name,
                    "--min-mapq", str(args.mapq),
                    "--walks-policy", "5unique",
                    "--drop-sam",
                ],
                [
                    "pairtools", "sort",
                    "--nproc", str(args.threads),
                    "--tmpdir", str(directories["tmp"]),
                    "--output", str(output),
                ],
            ],
            output,
        )
        if not args.dry_run and not ready(output):
            raise RuntimeError(f"lane pairs 未生成: {output}")
        write_done(lane_done, args.dry_run)

    merged = directories["pairs"] / "all.sorted.pairs.gz"
    merge_done = directories["pairs"] / ".merge_complete"
    if not (done(merge_done) and ready(merged)):
        merged.unlink(missing_ok=True)
        if len(lane_pairs) == 1:
            if not args.dry_run:
                shutil.copy2(lane_pairs[0], merged)
        else:
            runner.run(
                [
                    "pairtools", "merge", "--nproc", str(args.threads),
                    "--output", str(merged), *[str(path) for path in lane_pairs],
                ]
            )
        write_done(merge_done, args.dry_run)

    dedup = directories["pairs"] / "all.dedup.pairs.gz"
    dups = directories["pairs"] / "all.dups.pairs.gz"
    unmapped = directories["pairs"] / "all.unmapped.pairs.gz"
    dedup_stats = directories["pairs"] / "dedup.stats.txt"
    dedup_done = directories["pairs"] / ".dedup_complete"
    if not (done(dedup_done) and ready(dedup)):
        for path in (dedup, dups, unmapped, dedup_stats):
            path.unlink(missing_ok=True)
        runner.run(
            [
                "pairtools", "dedup", "--max-mismatch", "3",
                "--output", str(dedup),
                "--output-dups", str(dups),
                "--output-unmapped", str(unmapped),
                "--output-stats", str(dedup_stats),
                str(merged),
            ]
        )
        write_done(dedup_done, args.dry_run)

    valid = directories["pairs"] / "all.valid.UU.pairs.gz"
    valid_stats = directories["pairs"] / "valid.stats.txt"
    select_done = directories["pairs"] / ".select_complete"
    if not (done(select_done) and ready(valid)):
        valid.unlink(missing_ok=True)
        valid_stats.unlink(missing_ok=True)
        runner.run(
            [
                "pairtools", "select", '(pair_type == "UU")',
                "--output", str(valid), str(dedup),
            ]
        )
        runner.run(["pairtools", "stats", "--output", str(valid_stats), str(valid)])
        write_done(select_done, args.dry_run)

    base_cool = directories["cool"] / f"contacts.{args.bin_size}.cool"
    cool_done = directories["cool"] / ".cload_complete"
    if not (done(cool_done) and ready(base_cool)):
        base_cool.unlink(missing_ok=True)
        runner.run(
            [
                "cooler", "cload", "pairs",
                "--assembly", args.assembly_name,
                "-c1", "2", "-p1", "3", "-c2", "4", "-p2", "5",
                f"{chromsizes}:{args.bin_size}", str(valid), str(base_cool),
            ]
        )
        write_done(cool_done, args.dry_run)

    if not args.no_balance:
        balance_done = directories["cool"] / ".balance_complete"
        if not done(balance_done):
            runner.run(
                ["cooler", "balance", "--nproc", str(args.threads), str(base_cool)]
            )
            write_done(balance_done, args.dry_run)

    mcool = directories["cool"] / "contacts.mcool"
    zoom_done = directories["cool"] / ".zoomify_complete"
    if not (done(zoom_done) and ready(mcool)):
        mcool.unlink(missing_ok=True)
        command = [
            "cooler", "zoomify", "--nproc", str(args.threads),
            "--resolutions", args.zoom_resolutions,
            "--out", str(mcool),
        ]
        if not args.no_balance:
            command.append("--balance")
        command.append(str(base_cool))
        runner.run(command)
        write_done(zoom_done, args.dry_run)

    runner.run(["cooler", "info", str(base_cool)], stdout_path=directories["cool"] / "contacts.info.txt")
    runner.run(["cooler", "ls", str(mcool)], stdout_path=directories["cool"] / "contacts.resolutions.txt")

    print("\n流程完成。")
    print(f"有效互作 pairs: {valid}")
    print(f"基础矩阵:       {base_cool}")
    print(f"多分辨率矩阵:   {mcool}")
    print(f"互作统计:       {valid_stats}")
    print(f"运行日志:       {outdir / 'pipeline.log'}")


if __name__ == "__main__":
    try:
        main()
    except subprocess.CalledProcessError as error:
        print(
            f"\n命令运行失败（退出码 {error.returncode}），请查看 pipeline.log。",
            file=sys.stderr,
        )
        raise SystemExit(error.returncode)
    except KeyboardInterrupt:
        print("\n用户中止。", file=sys.stderr)
        raise SystemExit(130)
