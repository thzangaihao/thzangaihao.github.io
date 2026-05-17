#!/usr/bin/env python3
import os
import subprocess
import datetime
import sys

# ==========================================================
#                      关键参数配置区
# ==========================================================
THREADS = 32              # 使用的核心数
HISAT2_BIN = "hisat2"    # HISAT2 可执行文件路径
SAMTOOLS_BIN = "samtools" # Samtools 可执行文件路径
DTA_MODE = True          # 是否开启 --dta (用于后续 StringTie 组装)
# ==========================================================

def run_command(cmd):
    """运行 shell 命令并检查错误"""
    print(f"\n[执行命令]: {cmd}")
    try:
        subprocess.run(cmd, shell=True, check=True)
    except subprocess.CalledProcessError as e:
        print(f"错误: 命令执行失败。\n{e}")
        sys.exit(1)

def main():
    print("=== HISAT2 RNA-seq Mapping 自动化脚本 ===")

    # 1. 交互式输入路径
    ref_fasta = input("请输入参考基因组 FASTA 文件路径: ").strip()
    read1 = input("请输入双端测序 R1 FASTQ 文件路径: ").strip()
    read2 = input("请输入双端测序 R2 FASTQ 文件路径: ").strip()

    # 路径检查
    for f in [ref_fasta, read1, read2]:
        if not os.path.exists(f):
            print(f"错误: 找不到文件 {f}")
            return

    # 2. 创建带时间戳的输出文件夹
    timestamp = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
    out_dir = os.path.abspath(f"mapping_output_{timestamp}")
    os.makedirs(out_dir, exist_ok=True)
    
    index_prefix = os.path.join(out_dir, "genome_index")
    sam_out = os.path.join(out_dir, "aligned.sam")
    bam_out = os.path.join(out_dir, "aligned_sorted.bam")

    print(f"\n结果将存放在: {out_dir}")

    # 3. 建立索引 (hisat2-build)
    print("\n--- 正在建立 HISAT2 索引 ---")
    build_cmd = f"{HISAT2_BIN}-build -p {THREADS} {ref_fasta} {index_prefix}"
    run_command(build_cmd)

    # 4. 执行 Mapping
    print("\n--- 正在进行 Mapping ---")
    dta_flag = "--dta" if DTA_MODE else ""
    mapping_cmd = (
        f"{HISAT2_BIN} -p {THREADS} {dta_flag} "
        f"-x {index_prefix} "
        f"-1 {read1} -2 {read2} "
        f"-S {sam_out}"
    )
    run_command(mapping_cmd)

    # 5. SAM 转 BAM 并排序 (samtools)
    print("\n--- 正在转换并排序 BAM 文件 ---")
    sort_cmd = f"{SAMTOOLS_BIN} sort -@ {THREADS} -o {bam_out} {sam_out}"
    run_command(sort_cmd)

    # 6. 建立 BAM 索引
    print("\n--- 正在为 BAM 文件建立索引 ---")
    index_bam_cmd = f"{SAMTOOLS_BIN} index {bam_out}"
    run_command(index_bam_cmd)

    # 可选：删除巨大的 SAM 中间文件
    if os.path.exists(bam_out):
        os.remove(sam_out)
        print(f"\n清理完成：已移除临时 SAM 文件。")

    print(f"\n任务完成！最终文件: {bam_out}")

if __name__ == "__main__":
    main()