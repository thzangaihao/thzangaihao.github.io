#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import datetime
import os
import re
import shutil
import subprocess
import sys

# ==========================================================
#                      关键参数配置区
# ==========================================================
THREADS = 32               # 每个任务使用的核心数
HISAT2_BIN = "hisat2"     # HISAT2 可执行文件
SAMTOOLS_BIN = "samtools" # Samtools 可执行文件
DTA_MODE = True            # 是否开启 --dta（用于后续 StringTie）
# ==========================================================

FASTA_SUFFIXES = (".fa", ".fasta", ".fna", ".fas")
FASTQ_SUFFIXES = (".fq", ".fastq", ".fq.gz", ".fastq.gz")


def get_base_dir():
    """扫描脚本所在目录；将脚本复制到项目目录即可自动扫描整个项目。"""
    return os.path.dirname(os.path.abspath(__file__))


def find_files(base_dir, suffixes):
    """递归查找指定后缀文件，并跳过本脚本生成的输出目录。"""
    found = []
    for root, dirs, files in os.walk(base_dir):
        dirs[:] = [d for d in dirs if not d.startswith("hisat2_mapping_")]
        for name in files:
            if name.lower().endswith(suffixes):
                found.append(os.path.join(root, name))
    return sorted(found)


def get_r2_path(r1_path):
    """根据常见命名规则推导 R2；无法识别时返回 None。"""
    directory, name = os.path.split(r1_path)
    r2_name, count = re.subn(
        r"(?i)([._-])R1(?=([._-]|$))", r"\1R2", name, count=1
    )
    if count:
        return os.path.join(directory, r2_name)
    r2_name, count = re.subn(
        r"(?i)([._-])1(?=\.(fastq|fq)(\.gz)?$)", r"\g<1>2", name, count=1
    )
    if count:
        return os.path.join(directory, r2_name)
    return None


def sample_name_from_r1(r1_path):
    """从 R1 文件名提取用于输出目录的样本名。"""
    name = os.path.basename(r1_path)
    name = re.sub(r"(?i)\.(fastq|fq)(\.gz)?$", "", name)
    name = re.split(r"(?i)[._-]R1(?=([._-]|$))", name, maxsplit=1)[0]
    name = re.sub(r"(?i)[._-]1$", "", name)
    return name.rstrip("._-") or "sample"


def discover_pairs(base_dir):
    """发现并验证所有双端 FASTQ 配对。"""
    pairs = []
    seen = set()
    for path in find_files(base_dir, FASTQ_SUFFIXES):
        r2_path = get_r2_path(path)
        if not r2_path:
            continue
        key = (os.path.normcase(path), os.path.normcase(r2_path))
        if key in seen:
            continue
        if os.path.isfile(r2_path):
            pairs.append((sample_name_from_r1(path), path, r2_path))
            seen.add(key)
        else:
            print(f"⚠️ 跳过未配对的 R1: {os.path.relpath(path, base_dir)}")
            print(f"   未找到对应 R2: {os.path.relpath(r2_path, base_dir)}")
    return pairs


def parse_selection(text, item_count):
    """解析 1、1,3、1-3 或 all。"""
    if text.strip().lower() == "all":
        return list(range(item_count))
    selected = set()
    for part in text.replace(" ", "").split(","):
        if not part:
            continue
        if "-" in part:
            start, end = map(int, part.split("-", 1))
            if start > end:
                start, end = end, start
            selected.update(range(start - 1, end))
        else:
            selected.add(int(part) - 1)
    if not selected or any(i < 0 or i >= item_count for i in selected):
        raise ValueError
    return sorted(selected)


def choose_reference(base_dir):
    references = find_files(base_dir, FASTA_SUFFIXES)
    if not references:
        print("❌ 未在脚本所在目录及子目录中找到参考基因组 FASTA 文件。")
        return None
    print(f"\n🧬 找到 {len(references)} 个候选参考基因组：")
    for i, path in enumerate(references, 1):
        print(f"  [{i}] {os.path.relpath(path, base_dir)}")
    while True:
        choice = input("\n👉 请选择一个参考基因组编号（退出:q）: ").strip().lower()
        if choice == "q":
            return None
        try:
            index = int(choice) - 1
            if 0 <= index < len(references):
                print(f"✅ 参考基因组: {references[index]}")
                return references[index]
        except ValueError:
            pass
        print("⚠️ 请输入列表中的有效编号。")


def choose_pairs(base_dir):
    pairs = discover_pairs(base_dir)
    if not pairs:
        print("❌ 未找到完整的双端 FASTQ 配对。")
        return []
    print(f"\n📂 找到 {len(pairs)} 个双端测序任务：")
    for i, (sample, r1, r2) in enumerate(pairs, 1):
        print(f"\n  [{i}] 样本: {sample}")
        print(f"      R1: {os.path.relpath(r1, base_dir)}")
        print(f"      R2: {os.path.relpath(r2, base_dir)}")
    while True:
        choice = input(
            "\n👉 请选择任务（单选:1，多选:1,3，范围:1-3，全部:all，退出:q）: "
        ).strip().lower()
        if choice == "q":
            return []
        try:
            selected = parse_selection(choice, len(pairs))
            result = [pairs[i] for i in selected]
            print(f"✅ 已选择 {len(result)} 个双端测序任务。")
            return result
        except (ValueError, TypeError):
            print("⚠️ 输入格式或编号无效，请重新输入。")


def check_tools():
    tools = [f"{HISAT2_BIN}-build", HISAT2_BIN, SAMTOOLS_BIN]
    missing = [tool for tool in tools if shutil.which(tool) is None]
    if missing:
        print(f"❌ 未找到软件: {', '.join(missing)}")
        print("请确认已安装 HISAT2 和 Samtools，并激活正确的 Conda 环境。")
        return False
    return True


def run_command(cmd, log_file=None):
    """安全执行命令；可将 stderr 写入日志。"""
    print("\n[执行命令] " + " ".join(cmd))
    if log_file:
        with open(log_file, "w", encoding="utf-8") as log:
            subprocess.run(cmd, check=True, stderr=log)
    else:
        subprocess.run(cmd, check=True)


def run_mapping(sample, r1, r2, index_prefix, batch_dir, task_no, total):
    """运行一个双端样本；失败后返回 False，让队列继续。"""
    safe_sample = re.sub(r"[^A-Za-z0-9._-]+", "_", sample)
    sample_dir = os.path.join(batch_dir, safe_sample)
    if os.path.exists(sample_dir):
        sample_dir = os.path.join(batch_dir, f"{safe_sample}_{task_no}")
    os.makedirs(sample_dir, exist_ok=True)
    bam_out = os.path.join(sample_dir, f"{safe_sample}.sorted.bam")
    hisat2_log = os.path.join(sample_dir, f"{safe_sample}.hisat2.log")

    print("\n" + "=" * 70)
    print(f"🚀 [任务 {task_no}/{total}] 样本: {sample}")
    print(f"   R1: {os.path.abspath(r1)}")
    print(f"   R2: {os.path.abspath(r2)}")
    print(f"   输出: {sample_dir}")
    print("=" * 70)

    hisat2_cmd = [HISAT2_BIN, "-p", str(THREADS)]
    if DTA_MODE:
        hisat2_cmd.append("--dta")
    hisat2_cmd.extend(["-x", index_prefix, "-1", r1, "-2", r2])
    sort_cmd = [SAMTOOLS_BIN, "sort", "-@", str(THREADS), "-o", bam_out, "-"]

    print("\n--- HISAT2 比对并转换为排序 BAM ---")
    print("[执行管道] " + " ".join(hisat2_cmd) + " | " + " ".join(sort_cmd))
    try:
        with open(hisat2_log, "w", encoding="utf-8") as log:
            hisat2_process = subprocess.Popen(
                hisat2_cmd, stdout=subprocess.PIPE, stderr=log
            )
            sort_process = subprocess.Popen(sort_cmd, stdin=hisat2_process.stdout)
            hisat2_process.stdout.close()
            sort_code = sort_process.wait()
            hisat2_code = hisat2_process.wait()
        if hisat2_code != 0 or sort_code != 0:
            raise subprocess.CalledProcessError(
                hisat2_code or sort_code, "HISAT2 | samtools sort"
            )
        print("\n--- 建立 BAM 索引 ---")
        run_command([SAMTOOLS_BIN, "index", "-@", str(THREADS), bam_out])
        print(f"✅ {sample} 完成: {bam_out}")
        return True
    except (subprocess.CalledProcessError, OSError) as error:
        print(f"❌ {sample} 运行失败: {error}")
        print(f"   请检查日志: {hisat2_log}")
        return False


def main():
    print("=== HISAT2 双端 RNA-seq 批量比对脚本 ===")
    if not check_tools():
        return 1
    base_dir = get_base_dir()
    print(f"📁 扫描目录: {base_dir}")
    reference = choose_reference(base_dir)
    if not reference:
        return 0
    pairs = choose_pairs(base_dir)
    if not pairs:
        return 0

    timestamp = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
    batch_dir = os.path.join(base_dir, f"hisat2_mapping_{timestamp}")
    index_prefix = os.path.join(batch_dir, "genome_index", "genome")
    os.makedirs(os.path.dirname(index_prefix), exist_ok=True)
    print(f"\n📁 批次输出目录: {batch_dir}")
    print("\n--- 建立一次 HISAT2 索引，供全部样本共用 ---")
    try:
        run_command([
            f"{HISAT2_BIN}-build", "-p", str(THREADS), reference, index_prefix
        ])
    except (subprocess.CalledProcessError, OSError) as error:
        print(f"❌ 参考基因组索引建立失败: {error}")
        return 1

    success = 0
    for task_no, (sample, r1, r2) in enumerate(pairs, 1):
        if run_mapping(
            sample, r1, r2, index_prefix, batch_dir, task_no, len(pairs)
        ):
            success += 1

    print("\n" + "=" * 70)
    print(f"🎉 队列结束：成功 {success}/{len(pairs)}，失败 {len(pairs) - success}")
    print(f"📁 结果目录: {batch_dir}")
    print("=" * 70)
    return 0 if success == len(pairs) else 1


if __name__ == "__main__":
    try:
        sys.exit(main())
    except KeyboardInterrupt:
        print("\n🛑 用户中断运行。")
        sys.exit(130)
