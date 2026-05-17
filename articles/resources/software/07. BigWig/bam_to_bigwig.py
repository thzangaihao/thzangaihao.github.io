#!/usr/bin/env python3
import os
import subprocess
import datetime
import sys

# ==========================================================
#                      关键参数配置区
# ==========================================================
THREADS = 8                # 使用的核心数
BAM_COVERAGE_BIN = "bamCoverage" # deepTools 的命令名称
BIN_SIZE = 1               # 分辨率：1 表示单碱基分辨率，数值越大文件越小但越模糊
NORMALIZATION = "RPKM"     # 标准化方法: RPKM, CPM, BPM, RPGC, None
SAMTOOLS_BIN = "samtools"  # 用于检查和建立索引
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
    print("=== BAM 转换为 BigWig 自动化脚本 ===")

    # 1. 交互式输入路径
    bam_file = input("请输入 BAM 文件路径: ").strip()

    # 路径检查
    if not os.path.exists(bam_file):
        print(f"错误: 找不到文件 {bam_file}")
        return

    # 2. 创建带时间戳的输出文件夹
    timestamp = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
    out_dir = os.path.abspath(f"bigwig_output_{timestamp}")
    os.makedirs(out_dir, exist_ok=True)
    
    # 定义输出文件名
    base_name = os.path.basename(bam_file).replace(".bam", "")
    bw_out = os.path.join(out_dir, f"{base_name}.bw")

    print(f"\n结果将存放在: {out_dir}")

    # 3. 检查并建立 BAM 索引 (.bai)
    # bamCoverage 要求输入的 BAM 文件必须经过排序且有索引
    bai_file = bam_file + ".bai"
    if not os.path.exists(bai_file):
        print("\n--- 未检测到索引文件，正在建立索引 ---")
        run_command(f"{SAMTOOLS_BIN} index {bam_file}")

    # 4. 执行 bamCoverage 转换
    print(f"\n--- 正在转换 BAM 为 BigWig (标准化方法: {NORMALIZATION}) ---")
    
    # 构建命令
    norm_flag = f"--normalizeUsing {NORMALIZATION}" if NORMALIZATION else ""
    bw_cmd = (
        f"{BAM_COVERAGE_BIN} -p {THREADS} "
        f"-b {bam_file} "
        f"-o {bw_out} "
        f"--binSize {BIN_SIZE} "
        f"{norm_flag}"
    )
    
    run_command(bw_cmd)

    print(f"\n任务完成！")
    print(f"生成的 BigWig 文件: {bw_out}")
    print("您现在可以将该文件拖入 IGV 进行可视化。")

if __name__ == "__main__":
    main()