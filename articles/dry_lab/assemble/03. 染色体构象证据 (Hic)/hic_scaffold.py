#!/usr/bin/env python3
"""Interactive Hi-C scaffolding pipeline: FASTQ -> BWA -> BAM -> YaHS.

Designed for paired-end Hi-C lanes named like *_R1/*_R2 or *_1/*_2.
The script never modifies input FASTQ files. It creates a reference link/copy,
logs every command, and supports resuming completed steps.
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
    low = name.lower()
    for suffix in FASTQ_SUFFIXES:
        if low.endswith(suffix):
            return name[: -len(suffix)]
    return name


def pair_key(path: Path) -> tuple[str, int] | None:
    stem = strip_fastq_suffix(path.name)
    for pattern in PAIR_PATTERNS:
        match = pattern.match(stem)
        if match:
            # Keep text after R1/R2 in the key, e.g. R1_001 pairs with R2_001.
            tail = match.group(4) or ""
            key = f"{match.group(1)}{match.group(2)}{{READ}}{tail}"
            return key.lower(), int(match.group(3))
    return None


def discover_pairs(root: Path) -> tuple[list[ReadPair], list[Path], list[str]]:
    files = sorted(
        p.resolve()
        for p in root.rglob("*")
        if p.is_file() and p.name.lower().endswith(FASTQ_SUFFIXES)
    )
    grouped: dict[tuple[Path, str], dict[int, list[Path]]] = {}
    unmatched: list[Path] = []
    for path in files:
        parsed = pair_key(path)
        if parsed is None:
            unmatched.append(path)
            continue
        key, read_no = parsed
        grouped.setdefault((path.parent, key), {1: [], 2: []})[read_no].append(path)

    pairs: list[ReadPair] = []
    problems: list[str] = []
    for (parent, key), reads in sorted(grouped.items(), key=lambda x: str(x[0])):
        if len(reads[1]) == 1 and len(reads[2]) == 1:
            label = strip_fastq_suffix(reads[1][0].name)
            label = re.sub(r"(?i)([_\.-])(R|READ)?1(?=([_\.-]|$))", "", label, count=1)
            pairs.append(ReadPair(label=label, r1=reads[1][0], r2=reads[2][0]))
        else:
            problems.append(
                f"{parent / key}: R1={len(reads[1])}, R2={len(reads[2])}"
            )
    return pairs, unmatched, problems


def choose_pairs(pairs: list[ReadPair]) -> list[ReadPair]:
    if not pairs:
        raise SystemExit("没有发现可配对的 FASTQ。请检查命名是否包含 R1/R2 或 _1/_2。")
    print("\n发现以下 FASTQ 配对：")
    for i, pair in enumerate(pairs, 1):
        print(f"  [{i}] {pair.label}")
        print(f"      R1: {pair.r1}")
        print(f"      R2: {pair.r2}")
    while True:
        answer = input("\n选择配对 [all，或编号如 1,3-4；默认 all]: ").strip()
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


def ask_path(prompt: str, default: Path | None = None) -> Path:
    suffix = f" [默认 {default}]" if default else ""
    while True:
        value = input(f"{prompt}{suffix}: ").strip().strip("'\"")
        path = Path(value).expanduser() if value else default
        if path and path.is_file():
            return path.resolve()
        print("文件不存在，请重新输入。", file=sys.stderr)


def require_tools(names: list[str]) -> None:
    missing = [name for name in names if shutil.which(name) is None]
    if missing:
        raise SystemExit("缺少必要软件: " + ", ".join(missing))


class Runner:
    def __init__(self, log_path: Path, dry_run: bool = False):
        self.log_path = log_path
        self.dry_run = dry_run

    def _write(self, text: str) -> None:
        print(text, flush=True)
        with self.log_path.open("a", encoding="utf-8") as handle:
            handle.write(text + "\n")

    def run(
        self,
        command: list[str],
        stdout_path: Path | None = None,
        env: dict[str, str] | None = None,
    ) -> None:
        shown = shlex.join(map(str, command))
        if env:
            shown = " ".join(f"{key}={shlex.quote(value)}" for key, value in env.items()) + " " + shown
        if stdout_path:
            shown += f" > {shlex.quote(str(stdout_path))}"
        self._write(f"[{datetime.now().isoformat(timespec='seconds')}] $ {shown}")
        if self.dry_run:
            return
        with self.log_path.open("a", encoding="utf-8") as log:
            child_env = os.environ.copy()
            if env:
                child_env.update(env)
            if stdout_path:
                with stdout_path.open("wb") as output:
                    subprocess.run(command, stdout=output, stderr=log, check=True, env=child_env)
            else:
                subprocess.run(
                    command, stdout=log, stderr=subprocess.STDOUT,
                    check=True, env=child_env,
                )

    def pipe(self, commands: list[list[str]], output_path: Path) -> None:
        shown = " | ".join(shlex.join(map(str, cmd)) for cmd in commands)
        self._write(
            f"[{datetime.now().isoformat(timespec='seconds')}] $ {shown}"
        )
        if self.dry_run:
            return
        processes: list[subprocess.Popen] = []
        with self.log_path.open("a", encoding="utf-8") as log:
            previous = None
            try:
                for i, command in enumerate(commands):
                    proc = subprocess.Popen(
                        command,
                        stdin=previous.stdout if previous else None,
                        stdout=subprocess.PIPE if i < len(commands) - 1 else subprocess.DEVNULL,
                        stderr=log,
                    )
                    if previous and previous.stdout:
                        previous.stdout.close()
                    processes.append(proc)
                    previous = proc
                return_codes = [proc.wait() for proc in processes]
                if any(code != 0 for code in return_codes):
                    output_path.unlink(missing_ok=True)
                    raise subprocess.CalledProcessError(
                        next(code for code in return_codes if code != 0), shown
                    )
            except Exception:
                for proc in processes:
                    if proc.poll() is None:
                        proc.terminate()
                raise


def ready(path: Path) -> bool:
    return path.is_file() and path.stat().st_size > 0


def safe_label(text: str, number: int) -> str:
    cleaned = re.sub(r"[^A-Za-z0-9_.-]+", "_", text).strip("_.-")
    return cleaned or f"lane{number}"


def prepare_reference(source: Path, destination: Path) -> None:
    if destination.exists():
        return
    try:
        destination.symlink_to(source)
    except OSError:
        print("无法创建参考序列软链接，改为复制 FASTA。")
        shutil.copy2(source, destination)


def main() -> None:
    parser = argparse.ArgumentParser(description="交互式 Hi-C + YaHS 挂载流程")
    parser.add_argument("--search-dir", type=Path, default=Path.cwd(), help="FASTQ 搜索目录")
    parser.add_argument("--fasta", type=Path, help="待挂载的组装 FASTA")
    parser.add_argument("--outdir", type=Path, help="输出目录")
    parser.add_argument("--threads", type=int, default=max(1, min(16, os.cpu_count() or 1)))
    parser.add_argument(
        "--java-memory", default="64g",
        help="Picard 最大 Java 堆内存（默认 64g，例如 16g、32g、128g）",
    )
    parser.add_argument("--motif", help="限制酶基序，如 GATC；未知或 Omni-C 时不填")
    parser.add_argument("--skip-fastqc", action="store_true")
    parser.add_argument("--dry-run", action="store_true", help="只显示命令，不运行")
    args = parser.parse_args()

    search_dir = args.search_dir.expanduser().resolve()
    if not search_dir.is_dir():
        raise SystemExit(f"FASTQ 搜索目录不存在: {search_dir}")
    pairs, unmatched, problems = discover_pairs(search_dir)
    if unmatched:
        print("\n警告：以下 FASTQ 未识别配对，将跳过：", file=sys.stderr)
        for path in unmatched:
            print(f"  {path}", file=sys.stderr)
    if problems:
        print("\n警告：以下分组不是一对一配对，将跳过：", file=sys.stderr)
        for problem in problems:
            print(f"  {problem}", file=sys.stderr)
    selected = choose_pairs(pairs)

    fasta = args.fasta.expanduser().resolve() if args.fasta else ask_path("输入组装 FASTA 路径")
    if not fasta.is_file():
        raise SystemExit(f"FASTA 不存在: {fasta}")
    outdir = (args.outdir or Path(input("输出目录 [默认 ./hic_yahs_result]: ").strip() or "hic_yahs_result"))
    outdir = outdir.expanduser().resolve()
    if args.threads < 1:
        raise SystemExit("线程数必须大于 0")
    if sys.stdin.isatty() and "--threads" not in sys.argv:
        value = input(f"线程数 [默认 {args.threads}]: ").strip()
        if value:
            args.threads = int(value)
    java_memory_given = any(
        arg == "--java-memory" or arg.startswith("--java-memory=")
        for arg in sys.argv[1:]
    )
    if sys.stdin.isatty() and not java_memory_given:
        value = input(
            f"Picard Java 最大堆内存 [默认 {args.java_memory}，例如 32g/64g]: "
        ).strip()
        if value:
            args.java_memory = value
    if not re.fullmatch(r"(?i)[1-9][0-9]*[mg]", args.java_memory):
        raise SystemExit("Java 内存格式错误，请使用 8g、16g、32000m 等格式")
    motif = args.motif
    if motif is None and sys.stdin.isatty():
        motif = input("限制酶基序 [如 GATC；未知/Omni-C 留空]: ").strip() or None

    require_tools(["bwa", "samtools", "yahs", "picard"] + ([] if args.skip_fastqc else ["fastqc", "multiqc"]))
    for directory in (
        outdir, outdir / "reference", outdir / "qc", outdir / "bam",
        outdir / "yahs", outdir / "tmp",
    ):
        directory.mkdir(parents=True, exist_ok=True)
    log_path = outdir / "pipeline.log"
    runner = Runner(log_path, args.dry_run)

    reference = outdir / "reference" / "assembly.fa"
    prepare_reference(fasta, reference)
    print("\n分析配置：")
    print(f"  FASTQ 配对数: {len(selected)}")
    print(f"  FASTA:        {fasta}")
    print(f"  输出目录:     {outdir}")
    print(f"  线程数:       {args.threads}")
    print(f"  Picard 堆内存:{args.java_memory}")
    print(f"  酶切基序:     {motif or '未指定'}")
    if sys.stdin.isatty() and input("开始运行？[Y/n]: ").strip().lower() in {"n", "no"}:
        raise SystemExit("已取消。")

    if not ready(Path(str(reference) + ".fai")):
        runner.run(["samtools", "faidx", str(reference)])
    if not ready(Path(str(reference) + ".bwt")):
        runner.run(["bwa", "index", str(reference)])

    if not args.skip_fastqc:
        qc_flag = outdir / "qc" / ".multiqc_complete"
        if not ready(qc_flag):
            runner.run(["fastqc", "-t", str(args.threads), "-o", str(outdir / "qc")] + [str(p) for pair in selected for p in (pair.r1, pair.r2)])
            runner.run(["multiqc", "-f", "-o", str(outdir / "qc"), str(outdir / "qc")])
            if not args.dry_run:
                qc_flag.write_text("complete\n", encoding="utf-8")

    lane_bams: list[Path] = []
    for i, pair in enumerate(selected, 1):
        label = f"{i:02d}_{safe_label(pair.label, i)}"
        bam = outdir / "bam" / f"{label}.sorted.bam"
        lane_bams.append(bam)
        if ready(bam):
            print(f"跳过已完成比对: {bam}")
            continue
        read_group = f"@RG\\tID:{label}\\tSM:HIC\\tLB:HIC\\tPL:ILLUMINA\\tPU:{label}"
        runner.pipe(
            [
                ["bwa", "mem", "-5SP", "-t", str(args.threads), "-R", read_group, str(reference), str(pair.r1), str(pair.r2)],
                ["samtools", "view", "-@", str(args.threads), "-bh", "-F", "0x904", "-"],
                ["samtools", "sort", "-@", str(args.threads), "-o", str(bam), "-"],
            ],
            bam,
        )

    merged = outdir / "bam" / "hic.merged.sorted.bam"
    if not ready(merged):
        if len(lane_bams) == 1:
            if not args.dry_run:
                shutil.copy2(lane_bams[0], merged)
        else:
            runner.run(["samtools", "merge", "-@", str(args.threads), "-f", str(merged)] + [str(p) for p in lane_bams])

    dedup = outdir / "bam" / "hic.dedup.bam"
    metrics = outdir / "bam" / "hic.duplicate_metrics.txt"
    dedup_done = outdir / "bam" / ".hic_dedup_complete"
    if not ready(dedup_done):
        # A previous out-of-memory failure can leave a non-empty but truncated BAM.
        for partial in (dedup, dedup.with_suffix(".bai"), Path(str(dedup) + ".bai"), metrics):
            partial.unlink(missing_ok=True)
        runner.run([
            "picard", "MarkDuplicates",
            f"I={merged}", f"O={dedup}", f"M={metrics}",
            "REMOVE_DUPLICATES=false", "ASSUME_SORTED=true",
            "VALIDATION_STRINGENCY=LENIENT", "CREATE_INDEX=true",
            f"TMP_DIR={outdir / 'tmp'}", "MAX_RECORDS_IN_RAM=250000",
        # Bioconda's Picard launcher may append its own small -Xmx value after
        # JAVA_TOOL_OPTIONS. HotSpot processes _JAVA_OPTIONS last, so this
        # reliably overrides the launcher's default heap limit.
        ], env={"_JAVA_OPTIONS": f"-Xmx{args.java_memory}"})
        if not args.dry_run:
            dedup_done.write_text("complete\n", encoding="utf-8")
    if not ready(Path(str(dedup) + ".bai")) and not ready(dedup.with_suffix(".bai")):
        runner.run(["samtools", "index", "-@", str(args.threads), str(dedup)])
    runner.run(["samtools", "flagstat", "-@", str(args.threads), str(dedup)], outdir / "bam" / "hic.dedup.flagstat.txt")

    prefix = outdir / "yahs" / "hic"
    final_fasta = Path(str(prefix) + "_scaffolds_final.fa")
    if not ready(final_fasta):
        command = ["yahs", str(reference), str(dedup), "-o", str(prefix)]
        if motif:
            command.extend(["-e", motif])
        runner.run(command)

    print("\n流程完成。")
    print(f"最终 FASTA: {final_fasta}")
    print(f"最终 AGP:   {prefix}_scaffolds_final.agp")
    print(f"运行日志:   {log_path}")
    print(f"质控报告:   {outdir / 'qc' / 'multiqc_report.html'}")


if __name__ == "__main__":
    try:
        main()
    except subprocess.CalledProcessError as exc:
        print(f"\n命令运行失败（退出码 {exc.returncode}），请查看 pipeline.log。", file=sys.stderr)
        raise SystemExit(exc.returncode)
    except KeyboardInterrupt:
        print("\n用户中止。", file=sys.stderr)
        raise SystemExit(130)
