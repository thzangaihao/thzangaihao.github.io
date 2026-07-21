#!/usr/bin/env python3
"""Run AUGUSTUS with optional RNA-seq evidence from BAM or FASTQ files."""

import glob
import os
import re
import shutil
import subprocess
import sys
import time
from datetime import datetime


FASTQ_PATTERNS = (
    "**/*.fastq", "**/*.fq", "**/*.fastq.gz", "**/*.fq.gz",
    "**/*.fastq.bz2", "**/*.fq.bz2",
)


def log_info(message):
    print(f"[{time.strftime('%Y-%m-%d %H:%M:%S')}] {message}")


def require_commands(commands):
    missing = [command for command in commands if shutil.which(command) is None]
    if missing:
        raise SystemExit("错误：未找到命令：" + ", ".join(missing))


def find_files(patterns):
    files = []
    for pattern in patterns:
        files.extend(glob.glob(pattern, recursive=True))
    return sorted({f for f in files if not os.path.basename(f).startswith(".")})


def select_file(prompt, patterns, allow_skip=False):
    files = find_files(patterns)
    if not files:
        if allow_skip:
            return None
        raise SystemExit(f"未找到相关文件：{', '.join(patterns)}")

    print(f"\n{prompt}")
    if allow_skip:
        print("  [0] 跳过")
    for number, path in enumerate(files, 1):
        size = os.path.getsize(path) / 1024 ** 2
        print(f"  [{number}] {path} ({size:.1f} MB)")

    while True:
        choice = input("请输入编号 >>> ").strip()
        if allow_skip and choice == "0":
            return None
        try:
            index = int(choice) - 1
        except ValueError:
            print("请输入有效编号。")
            continue
        if 0 <= index < len(files):
            return files[index]
        print("请输入有效编号。")


def select_multiple(prompt, items, formatter=str, default_all=True):
    """Select one or more items using comma-separated numbers or ``all``."""
    if not items:
        return []
    print(f"\n{prompt}")
    for number, item in enumerate(items, 1):
        print(f"  [{number}] {formatter(item)}")
    default = "all" if default_all else ""
    while True:
        choice = input(f"输入编号（逗号分隔）或 all [默认 {default or '无'}] >>> ").strip().lower()
        if not choice and default_all:
            return list(items)
        if choice == "all":
            return list(items)
        try:
            indexes = [int(value.strip()) - 1 for value in choice.split(",")]
        except ValueError:
            print("请输入有效编号，例如 1,3,5 或 all。")
            continue
        if indexes and len(set(indexes)) == len(indexes) and all(0 <= i < len(items) for i in indexes):
            return [items[i] for i in indexes]
        print("编号重复或超出范围，请重新输入。")


def ask_choice(prompt, choices, default):
    labels = "/".join(choices)
    while True:
        answer = input(f"{prompt} [{labels}，默认 {default}] >>> ").strip().lower()
        if not answer:
            return default
        if answer in choices:
            return answer
        print("请输入：" + "、".join(choices))


def strip_fastq_suffix(path):
    return re.sub(r"(?i)(?:\.fastq|\.fq)(?:\.gz|\.bz2)?$", "", path)


def infer_fastq_pair(selected_read):
    """Return an ordered (R1, R2) pair using only the final number before the suffix.

    Examples: leaf_1_1.fq.gz <-> leaf_1_2.fq.gz and sample_R1.fastq
    <-> sample_R2.fastq. Numbers earlier in the sample name are left untouched.
    """
    directory, filename = os.path.split(selected_read)
    match = re.fullmatch(
        r"(?i)(?P<prefix>.+)(?P<separator>[_\.-])(?P<read>R?[12])"
        r"(?P<suffix>\.(?:fastq|fq)(?:\.(?:gz|bz2))?)",
        filename,
    )
    if not match:
        return None

    read_token = match.group("read")
    uses_r = read_token.lower().startswith("r")
    read_number = read_token[-1]
    mate_number = "2" if read_number == "1" else "1"
    mate_token = (read_token[0] if uses_r else "") + mate_number
    mate_name = "".join((
        match.group("prefix"), match.group("separator"), mate_token,
        match.group("suffix"),
    ))
    mate = os.path.join(directory, mate_name)
    if not os.path.isfile(mate):
        return None
    return (selected_read, mate) if read_number == "1" else (mate, selected_read)


def discover_fastq_datasets():
    """Discover paired and unpaired FASTQ datasets and print the pairing result."""
    files = find_files(FASTQ_PATTERNS)
    used = set()
    datasets = []
    for path in files:
        key = os.path.normcase(os.path.abspath(path))
        if key in used:
            continue
        pair = infer_fastq_pair(path)
        if pair:
            read1, read2 = pair
            used.update(os.path.normcase(os.path.abspath(p)) for p in pair)
            datasets.append((read1, read2))
        else:
            used.add(key)
            datasets.append((path, None))

    print("\nFASTQ 自动配对结果：")
    for number, (read1, read2) in enumerate(datasets, 1):
        if read2:
            print(f"  [{number}] 双端\n       R1: {read1}\n       R2: {read2}")
        else:
            print(f"  [{number}] 单端/未配对\n       SE: {read1}")
    paired = sum(read2 is not None for _, read2 in datasets)
    print(f"共发现 {len(datasets)} 个数据集：{paired} 个双端，{len(datasets) - paired} 个单端/未配对。")
    return datasets


def choose_fastq_datasets():
    datasets = discover_fastq_datasets()
    if not datasets:
        raise SystemExit("错误：未发现 FASTQ 文件。")
    return select_multiple(
        "请选择用于整合预测的转录组数据集：",
        datasets,
        formatter=lambda reads: (
            f"双端: {reads[0]} + {reads[1]}" if reads[1] else f"单端: {reads[0]}"
        ),
    )


def run_checked(command, **kwargs):
    log_info("执行：" + " ".join(map(str, command)))
    subprocess.run(command, check=True, **kwargs)


def build_and_align(reference, read1, read2, aligner, out_dir, threads):
    """Build an index and stream alignments directly into coordinate-sorted BAM."""
    require_commands([aligner, f"{aligner}-build", "samtools"])
    index_dir = os.path.join(out_dir, f"{aligner}_index")
    os.makedirs(index_dir, exist_ok=True)
    prefix = os.path.join(index_dir, "genome")
    bam = os.path.join(out_dir, f"{os.path.basename(strip_fastq_suffix(read1))}.{aligner}.sorted.bam")

    index_suffixes = (".1.ht2", ".1.ht2l") if aligner == "hisat2" else (".1.bt2", ".1.bt2l")
    if not any(os.path.isfile(prefix + suffix) for suffix in index_suffixes):
        build = [f"{aligner}-build"]
        if aligner == "hisat2":
            build += ["-p", str(threads)]
        build += [reference, prefix]
        run_checked(build)
    else:
        log_info(f"复用已建立的 {aligner} 索引：{prefix}")

    align = [aligner, "-p", str(threads), "-x", prefix]
    if aligner == "hisat2":
        align.append("--dta")
    if read2:
        align += ["-1", read1, "-2", read2]
    else:
        align += ["-U", read1]

    sort = ["samtools", "sort", "-@", str(threads), "-o", bam, "-"]
    log_info("开始比对并直接生成坐标排序 BAM。")
    align_process = subprocess.Popen(align, stdout=subprocess.PIPE)
    try:
        sort_process = subprocess.run(sort, stdin=align_process.stdout)
    finally:
        if align_process.stdout:
            align_process.stdout.close()
    align_code = align_process.wait()
    if align_code or sort_process.returncode:
        if os.path.exists(bam):
            os.remove(bam)
        raise subprocess.CalledProcessError(align_code or sort_process.returncode, align if align_code else sort)
    run_checked(["samtools", "quickcheck", "-v", bam])
    run_checked(["samtools", "index", "-@", str(threads), bam])
    log_info(f"已生成 RNA-seq BAM：{bam}")
    return bam


def choose_rnaseq_evidence(reference, out_dir, threads):
    bam_files = find_files(("**/*.bam",))
    fastq_files = find_files(FASTQ_PATTERNS)
    available = ["none"]
    if bam_files:
        available.insert(0, "bam")
    if fastq_files:
        available.insert(0, "fastq")
    mode = ask_choice("RNA-seq 证据类型", available, available[0])
    if mode == "none":
        return []
    if mode == "bam":
        return select_multiple(
            "请选择用于整合预测的 BAM：",
            bam_files,
            formatter=lambda path: f"{path} ({os.path.getsize(path) / 1024 ** 2:.1f} MB)",
        )

    datasets = choose_fastq_datasets()
    aligner = ask_choice("序列比对工具（转录组推荐 hisat2）", ["hisat2", "bowtie2"], "hisat2")
    if aligner == "bowtie2":
        log_info("警告：bowtie2 不支持跨内含子剪接比对，真核 RNA-seq 通常应选择 hisat2。")
    bam_files = []
    for number, (read1, read2) in enumerate(datasets, 1):
        log_info(f"处理转录组数据集 {number}/{len(datasets)}：{read1}")
        bam_files.append(build_and_align(reference, read1, read2, aligner, out_dir, threads))
    return bam_files


def build_combined_hints(bam_files, out_dir, sample, threads):
    """Merge all RNA-seq BAMs, then extract one integrated hints file."""
    evidence_bam = bam_files[0]
    if len(bam_files) > 1:
        require_commands(["samtools"])
        evidence_bam = os.path.join(out_dir, f"{sample}.rnaseq_merged.sorted.bam")
        log_info(f"正在合并 {len(bam_files)} 个 RNA-seq BAM 用于整合预测。")
        run_checked([
            "samtools", "merge", "-@", str(threads), "-f", evidence_bam,
            *bam_files,
        ])
        run_checked(["samtools", "quickcheck", "-v", evidence_bam])
        run_checked(["samtools", "index", "-@", str(threads), evidence_bam])
        log_info(f"整合 BAM：{evidence_bam}")

    combined = os.path.join(out_dir, f"{sample}.combined.hints.gff")
    run_checked(["bam2hints", f"--in={evidence_bam}", f"--out={combined}"])
    log_info(f"已从 {len(bam_files)} 组 RNA-seq 数据生成整合 hints：{combined}")
    return combined, evidence_bam


def main():
    print("=" * 64)
    print("AUGUSTUS 基因预测（支持 BAM 或 FASTQ RNA-seq 证据）")
    print("=" * 64)
    require_commands(["augustus", "bam2hints"])

    reference = select_file(
        "请选择待预测的基因组 FASTA：",
        ("**/*.fasta", "**/*.fa", "**/*.fna", "**/*.p_ctg.fasta"),
    )
    sample = os.path.splitext(os.path.basename(reference))[0]
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    out_dir = os.path.abspath(f"augustus_{sample}_{timestamp}")
    os.makedirs(out_dir, exist_ok=True)

    try:
        threads = int(input(f"线程数 [默认 {os.cpu_count() or 1}] >>> ").strip() or (os.cpu_count() or 1))
        if threads < 1:
            raise ValueError
    except ValueError:
        raise SystemExit("错误：线程数必须是正整数。")

    bam_files = choose_rnaseq_evidence(reference, out_dir, threads)
    species = input("AUGUSTUS 物种模型（如 aspergillus_fumigatus）[必填] >>> ").strip()
    if not species:
        raise SystemExit("错误：物种模型不能为空。")

    output = os.path.join(out_dir, f"{sample}.augustus.gff")
    hints = None
    evidence_bam = None
    if bam_files:
        print("\n参与整合预测的 RNA-seq BAM：")
        for number, bam in enumerate(bam_files, 1):
            print(f"  [{number}] {bam}")
        hints, evidence_bam = build_combined_hints(bam_files, out_dir, sample, threads)

    command = ["augustus", f"--species={species}"]
    if hints:
        command += [
            f"--hintsfile={hints}",
            "--allow_hinted_splicesites=atac",
            "--extrinsicCfgFile=extrinsic.M.RM.E.W.cfg",
        ]
    command.append(reference)
    log_info("开始 AUGUSTUS 预测。")
    with open(output, "w", encoding="utf-8") as handle:
        run_checked(command, stdout=handle, text=True)
    log_info(f"预测完成：{output}")
    if bam_files:
        log_info(f"共整合 {len(bam_files)} 组 RNA-seq 证据。")
        log_info(f"用于提取 Hints 的 BAM：{evidence_bam}")
        log_info(f"整合 Hints 文件：{hints}")


if __name__ == "__main__":
    try:
        main()
    except subprocess.CalledProcessError as error:
        raise SystemExit(f"命令执行失败（退出码 {error.returncode}）：{error.cmd}")
