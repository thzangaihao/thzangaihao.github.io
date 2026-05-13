#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import sys
import glob
import subprocess
from datetime import datetime
import concurrent.futures

'''
诱变突变检测 (全样本版)
功能：自动化执行 SNP/Indel (bcftools) 与 SV (Delly) 检测。
特性：取消了对 CK 样本的屏蔽，可用于提取背景噪音。
'''

def get_base_dir():
    return os.path.dirname(os.path.abspath(__file__))

def check_env():
    tools = ['bcftools', 'delly', 'samtools']
    for t in tools:
        if subprocess.call(['which', t], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL) != 0:
            print(f"❌ 错误: 未找到 {t}！请执行: conda install -c bioconda bcftools delly samtools")
            sys.exit(1)

def interactive_select(files, desc):
    """通用的交互式选择逻辑"""
    if not files:
        print(f"⚠️ 未找到任何 {desc}！")
        return []
    
    print(f"\n📂 扫描到以下 {len(files)} 个 {desc}:")
    for i, f in enumerate(files, 1):
        print(f"  [{i}] {os.path.relpath(f, get_base_dir())}")
        
    while True:
        choice = input(f"\n👉 请选择 {desc} (单选:1, 多选:1,3, 范围:1-3, 全部:all, 退出:q): ").strip().lower()
        if choice == 'q': sys.exit(0)
        if choice == 'all': return files
        
        try:
            selected = []
            parts = choice.replace(' ', '').split(',')
            for part in parts:
                if '-' in part:
                    start, end = map(int, part.split('-'))
                    selected.extend(files[start-1:end])
                else:
                    selected.append(files[int(part)-1])
            return list(set(selected))
        except:
            print("⚠️ 输入格式错误，请重新选择。")

def get_run_configs():
    """交互式获取过滤阈值与多线程资源配置"""
    print("\n🔬 --- 变异过滤参数设置 ---")
    try:
        qual = input("👉 最小质量分 QUAL (默认 30): ") or "30"
        dp = input("👉 最小测序深度 DP (默认 10): ") or "10"
    except:
        qual, dp = "30", "10"
        
    print("\n💻 --- 运算资源与多线程配置 ---")
    print("⚠️ 提示: 总CPU占用 = 同时处理的样本数 × 每个样本分配的核心数")
    try:
        max_workers = int(input("👉 同时处理几个样本？(默认 2): ") or "2")
        threads_per_task = int(input("👉 每个样本分配几个CPU核心？(默认 2): ") or "2")
    except:
        max_workers, threads_per_task = 2, 2
        
    return qual, dp, max_workers, threads_per_task

def process_single_bam(bam, ref_fasta, out_root, qual_th, dp_th, threads):
    """处理单个 BAM 文件的核心逻辑 (供线程池调用)"""
    sample_name = os.path.basename(bam).replace('.bam', '')
    sample_dir = os.path.join(out_root, sample_name)
    os.makedirs(sample_dir, exist_ok=True)
    
    print(f"\n🚀 开始处理样本: {sample_name}")

    # --- SNP/Indel 检测 (注入 --threads 参数加速) ---
    raw_vcf = os.path.join(sample_dir, f"{sample_name}_raw.vcf")
    filtered_vcf = os.path.join(sample_dir, f"{sample_name}_SNP_filtered.vcf")
    
    # mpileup 和 call 都开启多线程
    snp_cmd = (f"bcftools mpileup --threads {threads} -Ou -f {ref_fasta} {bam} | "
               f"bcftools call --threads {threads} -mv -Ov -o {raw_vcf}")
    subprocess.run(snp_cmd, shell=True, check=True, stderr=subprocess.DEVNULL)

    # 过滤操作也极度消耗 IO，开启多线程
    filter_cmd = (f"bcftools filter --threads {threads} -e 'QUAL<{qual_th} || DP<{dp_th}' -s LOWQUAL -m + {raw_vcf} | "
                  f"bcftools view --threads {threads} -f PASS > {filtered_vcf}")
    subprocess.run(filter_cmd, shell=True, stderr=subprocess.DEVNULL)

    # --- SV 结构变异检测 (利用 OMP_NUM_THREADS 唤醒 Delly 的 OpenMP 多线程) ---
    sv_vcf = os.path.join(sample_dir, f"{sample_name}_SV_delly.vcf")
    # 注意命令开头的 OMP_NUM_THREADS
    sv_cmd = f"OMP_NUM_THREADS={threads} delly call -g {ref_fasta} -o {sv_vcf} {bam}"
    subprocess.run(sv_cmd, shell=True, stderr=subprocess.DEVNULL)

    print(f"✅ 样本 {sample_name} 分析完成！")
    return sample_name

def run_pipeline():
    base_dir = get_base_dir()
    
    # 1. 选择参考基因组
    ref_list = glob.glob(os.path.join(base_dir, "05_Cleaned_Ref_*", "*.fasta"))
    selected_refs = interactive_select(ref_list, "参考基因组 (.fasta)")
    if not selected_refs: return
    ref_fasta = selected_refs[0]

    # 2. 选择待检测的 BAM 文件 (已解除对 CK 的屏蔽)
    bam_list = glob.glob(os.path.join(base_dir, "**", "*.bam"), recursive=True)
    selected_bams = interactive_select(bam_list, "待检测 BAM 样本")
    if not selected_bams: return

    # 3. 设置阈值与资源配置
    qual_th, dp_th, max_workers, threads_per_task = get_run_configs()

    # 4. 创建输出目录
    timestamp = datetime.now().strftime("%Y%m%d_%H%M")
    out_root = os.path.join(base_dir, f"06_Variant_Results_{timestamp}")
    os.makedirs(out_root, exist_ok=True)

    # 5. 准备索引
    if not os.path.exists(ref_fasta + ".fai"):
        print("\n⚙️ 正在建立参考基因组索引...")
        subprocess.run(f"samtools faidx {ref_fasta}", shell=True)

    # 6. 启动多线程并发池
    print(f"\n🔥 启动多线程引擎: 同时处理 {max_workers} 个样本，单样本分配 {threads_per_task} 核...")
    
    with concurrent.futures.ThreadPoolExecutor(max_workers=max_workers) as executor:
        futures = [executor.submit(process_single_bam, bam, ref_fasta, out_root, qual_th, dp_th, threads_per_task) for bam in selected_bams]
        
        for future in concurrent.futures.as_completed(futures):
            try:
                future.result() 
            except Exception as exc:
                print(f"❌ 某个样本处理时发生错误: {exc}")

    print(f"\n🎉 所有变异检测已完成！结果存放在: {out_root}")

if __name__ == "__main__":
    try:
        check_env()
        run_pipeline()
    except KeyboardInterrupt:
        print("\n🛑 用户强制退出。")