import glob
import gzip
import os
import shutil
import subprocess
import sys
import time
from datetime import datetime


THREADS = 128
READ_EXTENSIONS = (".fastq.gz", ".fq.gz", ".fastq", ".fq", ".bam")


def log_info(message):
    print(f"[{time.strftime('%Y-%m-%d %H:%M:%S')}] {message}", flush=True)


def select_threads(default=THREADS):
    """在终端选择本次运行使用的核心数。"""
    detected = os.cpu_count()
    detected_text = str(detected) if detected is not None else "未知"
    print(f"\n检测到当前系统有 {detected_text} 个逻辑核心。")

    while True:
        try:
            choice = input(f"请输入本次运行使用的核心数（直接回车使用 {default}）>>> ").strip()
        except EOFError:
            log_info(f"未检测到交互式输入，使用默认核心数：{default}")
            return default

        if not choice:
            return default

        try:
            threads = int(choice)
        except ValueError:
            print("输入无效：核心数必须是正整数，请重新输入。")
            continue

        if threads < 1:
            print("输入无效：核心数必须大于 0，请重新输入。")
            continue

        if detected is not None and threads > detected:
            print(f"提示：选择的核心数 ({threads}) 超过检测到的逻辑核心数 ({detected})。")
        return threads


def check_dependencies():
    missing = [command for command in ("samtools", "hifiasm") if shutil.which(command) is None]
    if missing:
        log_info(f"严重错误：找不到命令：{', '.join(missing)}。请先安装并加入 PATH。")
        sys.exit(1)


def scan_read_files():
    files = []
    for extension in READ_EXTENSIONS:
        files.extend(glob.glob(f"**/*{extension}", recursive=True))
    return sorted({path for path in files if not os.path.basename(path).startswith(".")})


def select_files(files, prompt, allow_none=False, excluded=None):
    excluded = set(excluded or [])
    candidates = [path for path in files if path not in excluded]
    if not candidates:
        if allow_none:
            return []
        raise RuntimeError("没有可选择的测序文件。")

    print("\n" + prompt)
    for number, path in enumerate(candidates, 1):
        print(f"  [{number}] {path}")
    suffix = "；输入 none 跳过" if allow_none else ""
    print(f"请输入编号（可用逗号分隔）或 all 选择全部{suffix}。")
    choice = input("你的选择 >>> ").strip().lower()
    if allow_none and choice in {"", "none", "no", "n", "0"}:
        return []
    if choice == "all":
        return candidates
    try:
        indices = [int(value.strip()) - 1 for value in choice.split(",")]
    except ValueError as error:
        raise ValueError("输入格式错误，应为编号、逗号分隔编号或 all。") from error
    if not indices or any(index < 0 or index >= len(candidates) for index in indices):
        raise ValueError("选择中包含无效编号。")
    return [candidates[index] for index in dict.fromkeys(indices)]


def sample_stem(path):
    name = os.path.basename(path)
    for extension in READ_EXTENSIONS:
        if name.lower().endswith(extension):
            return name[: -len(extension)]
    return os.path.splitext(name)[0]


def bam_to_fastq(bam_file, fastq_file, threads):
    log_info(f"  -> BAM 转 FASTQ：{bam_file}")
    with open(fastq_file, "wb") as output:
        subprocess.run(
            ["samtools", "fastq", "-@", str(threads), bam_file],
            stdout=output,
            check=True,
        )
    return fastq_file


def prepare_read(path, output_dir, threads, label=""):
    """BAM 转 FASTQ；FASTQ/FQ（可 gzip）直接使用绝对路径。"""
    if path.lower().endswith(".bam"):
        filename = f"{sample_stem(path)}{label}.fastq"
        return bam_to_fastq(path, os.path.join(output_dir, filename), threads)
    return os.path.abspath(path)


def combine_ultralong_reads(paths, output_file):
    """将多个超长读长 FASTQ/FQ（含 gzip）流式合并成 gzip FASTQ。"""
    if len(paths) == 1:
        return paths[0]
    log_info(f"  -> 合并 {len(paths)} 个超长读长文件：{output_file}")
    with gzip.open(output_file, "wb") as destination:
        for path in paths:
            opener = gzip.open if path.lower().endswith(".gz") else open
            with opener(path, "rb") as source:
                shutil.copyfileobj(source, destination, length=16 * 1024 * 1024)
    return output_file


def run_hifiasm(hifi_file, prefix, threads, ploidy_mode, ultralong_file=None):
    cmd = ["hifiasm", "-o", prefix, "-t", str(threads)]
    if ploidy_mode in {"1", "2"}:
        cmd.append("--primary")
    if ultralong_file:
        cmd.extend(["--ul", ultralong_file])
    cmd.append(hifi_file)

    log_file = f"{prefix}.log"
    log_info(f"  -> 执行命令：{' '.join(cmd)}")
    log_info(f"  -> hifiasm 日志：{log_file}")
    with open(log_file, "w", encoding="utf-8") as log:
        subprocess.run(cmd, stdout=log, stderr=subprocess.STDOUT, check=True)


def extract_fasta(gfa_file, fasta_file):
    if not os.path.exists(gfa_file):
        log_info(f"  -> 警告：找不到图文件 {gfa_file}，跳过提取。")
        return
    count = 0
    with open(gfa_file, encoding="utf-8") as gfa, open(fasta_file, "w", encoding="utf-8") as fasta:
        for line in gfa:
            if line.startswith("S\t"):
                parts = line.rstrip().split("\t")
                if len(parts) >= 3:
                    fasta.write(f">{parts[1]}\n{parts[2]}\n")
                    count += 1
    log_info(f"  -> 从 {os.path.basename(gfa_file)} 提取 {count} 条 contig。")


def main():
    print("=" * 64)
    print("HiFi 基因组组装（支持 hifiasm 超长 ONT 读长证据）")
    print("=" * 64)
    check_dependencies()
    threads = select_threads()
    log_info(f"本次运行使用核心数：{threads}")

    log_info("正在扫描当前目录及子目录中的 BAM/FASTQ/FQ 文件……")
    all_files = scan_read_files()
    if not all_files:
        log_info("未找到测序文件，脚本退出。")
        return

    try:
        hifi_files = select_files(all_files, "请选择作为组装主输入的 PacBio HiFi 文件：")
        ultralong_files = select_files(
            all_files,
            "请选择用于证据支持的超长 ONT 读长文件（可不选）：",
            allow_none=True,
            excluded=hifi_files,
        )
    except (ValueError, RuntimeError) as error:
        log_info(f"选择失败：{error}")
        sys.exit(1)

    print("\n请选择组装策略：")
    print("  [1] 单倍体/纯合模式（输出 primary contig）")
    print("  [2] 主-副组装模式（输出 primary/alternate contig）")
    print("  [3] 二倍体分型模式（输出 hap1/hap2）")
    ploidy_mode = input("策略选择 (1/2/3) >>> ").strip()
    if ploidy_mode not in {"1", "2", "3"}:
        log_info("输入无效，使用 [3] 二倍体分型模式。")
        ploidy_mode = "3"

    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    work_dir = os.path.abspath(f"assemble_{timestamp}")
    os.makedirs(work_dir, exist_ok=True)
    log_info(f"总输出目录：{work_dir}")

    # 超长读长只准备一次，可作为所有所选 HiFi 样本的共同证据。
    prepared_ul = None
    if ultralong_files:
        ul_dir = os.path.join(work_dir, "00_ultralong_reads")
        os.makedirs(ul_dir, exist_ok=True)
        prepared = [prepare_read(path, ul_dir, threads, f".ul{index}")
                    for index, path in enumerate(ultralong_files, 1)]
        prepared_ul = combine_ultralong_reads(
            prepared, os.path.join(ul_dir, "ultralong.combined.fastq.gz")
        )
        log_info(f"超长读长证据已启用：{prepared_ul}")
    else:
        log_info("本次不使用超长读长证据。")

    for index, hifi_file in enumerate(hifi_files, 1):
        sample_name = sample_stem(hifi_file)
        print("-" * 64)
        log_info(f"任务 [{index}/{len(hifi_files)}]：{sample_name}")
        sample_dir = os.path.join(work_dir, sample_name)
        fastq_dir = os.path.join(sample_dir, "01_fastq")
        hifiasm_dir = os.path.join(sample_dir, "02_hifiasm")
        fasta_dir = os.path.join(sample_dir, "03_fasta")
        for directory in (fastq_dir, hifiasm_dir, fasta_dir):
            os.makedirs(directory, exist_ok=True)
        prefix = os.path.join(hifiasm_dir, sample_name)

        try:
            hifi_input = prepare_read(hifi_file, fastq_dir, threads)
            run_hifiasm(hifi_input, prefix, threads, ploidy_mode, prepared_ul)
            targets = {
                # --primary 产生 prefix.p_ctg.gfa / prefix.a_ctg.gfa（无 bp. 前缀）。
                "1": [("p_ctg.gfa", "p_ctg.fasta")],
                "2": [("p_ctg.gfa", "p_ctg.fasta"), ("a_ctg.gfa", "a_ctg.fasta")],
                "3": [("bp.hap1.p_ctg.gfa", "hap1.fasta"), ("bp.hap2.p_ctg.gfa", "hap2.fasta")],
            }[ploidy_mode]
            for gfa_suffix, fasta_suffix in targets:
                extract_fasta(
                    f"{prefix}.{gfa_suffix}",
                    os.path.join(fasta_dir, f"{sample_name}.{fasta_suffix}"),
                )
            log_info(f"样本 {sample_name} 完成：{fasta_dir}")
        except Exception as error:
            log_info(f"样本 {sample_name} 失败，继续下一样本。错误：{error}")

    print("-" * 64)
    log_info(f"全部任务结束：{work_dir}")


if __name__ == "__main__":
    main()
