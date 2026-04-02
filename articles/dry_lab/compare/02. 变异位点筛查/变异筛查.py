#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import sys
import glob
import subprocess
from datetime import datetime
import concurrent.futures

'''
诱变突变检测 (全样本版) - 倍性适配与深度过滤增强版
功能：自动化执行 SNP/Indel (bcftools) 与 SV (Delly) 检测。
特性：
1. 支持交互式选择倍性（单倍体/二倍体），并自动调整 AF 默认阈值。
2. 强制提取输出 GQ 标签，彻底修复单倍体模式下过滤表达式找不到变量的报错。
3. 释放了子进程的报错日志，方便后续查错。
'''

def get_base_dir():
    """获取脚本自身所在的文件夹作为检索基准"""
    return os.path.dirname(os.path.abspath(__file__))

def check_env():
    tools = ['bcftools', 'delly', 'samtools']
    for t in tools:
        if subprocess.call(['which', t], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL) != 0:
            print(f"❌ 错误: 未找到 {t}！请执行: conda install -c bioconda bcftools delly samtools")
            sys.exit(1)

def interactive_select(files, desc, display_root):
    """通用的交互式选择逻辑 (路径展示基于 display_root 进行相对化)"""
    if not files:
        print(f"⚠️ 未找到任何 {desc}！")
        return []
    
    print(f"\n📂 在当前目录及其子目录下，扫描到以下 {len(files)} 个 {desc}:")
    for i, f in enumerate(files, 1):
        rel_path = os.path.relpath(f, display_root)
        print(f"  [{i}] {rel_path}")
        
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
    """交互式获取倍性、多维度过滤阈值与多线程资源配置"""
    print("\n🧬 --- 样本倍性设置 ---")
    ploidy_input = input("👉 请输入物种倍性 (单倍体输入1, 二倍体输入2) [默认 1]: ").strip()
    ploidy = "2" if ploidy_input == "2" else "1"
    
    # 根据倍性动态设定默认的 AF 阈值
    default_af = "0.25" if ploidy == "2" else "0.8"

    print("\n🔬 --- 变异硬过滤 (Hard Filtering) 参数设置 ---")
    try:
        qual = input("👉 最小质量分 QUAL (默认 30): ") or "30"
        dp   = input("👉 最小测序深度 DP (默认 30): ") or "30"
        af   = input(f"👉 最小等位基因频率 AF (默认 {default_af}): ") or default_af
        gq   = input("👉 最小基因型质量 GQ (默认 20): ") or "20"
        mq   = input("👉 最小比对质量 MQ (建议 20, 默认 20): ") or "20"
        sp   = input("👉 最大链偏好性 SP (相当于 SOR/FS) (默认 60): ") or "60"
    except:
        qual, dp, af, gq, mq, sp = "30", "10", default_af, "20", "20", "60"
        
    print("\n💻 --- 运算资源与多线程配置 ---")
    try:
        max_workers = int(input("👉 同时处理几个样本？(默认 2): ") or "2")
        threads_per_task = int(input("👉 每个样本分配几个CPU核心？(默认 2): ") or "2")
    except:
        max_workers, threads_per_task = 2, 2
        
    return ploidy, (qual, dp, af, gq, mq, sp), max_workers, threads_per_task

def process_single_bam(bam, ref_fasta, out_root, ploidy, thresholds, threads):
    """处理单个 BAM 文件的核心逻辑"""
    qual_th, dp_th, af_th, gq_th, mq_th, sp_th = thresholds
    sample_name = os.path.basename(bam).replace('.bam', '')
    sample_dir = os.path.join(out_root, sample_name)
    os.makedirs(sample_dir, exist_ok=True)
    
    print(f"\n🚀 开始处理样本: {sample_name}")

    # --- 1. SNP/Indel 检测 ---
    raw_vcf = os.path.join(sample_dir, f"{sample_name}_raw.vcf")
    filtered_vcf = os.path.join(sample_dir, f"{sample_name}_SNP_filtered.vcf")
    
    # 核心修复 1: mpileup 保留必要格式, call 阶段加入 --ploidy 设定，并通过 -f GQ 强制输出基因型质量
    snp_cmd = (f"bcftools mpileup --threads {threads} -a FORMAT/AD,FORMAT/DP,FORMAT/SP -Ou -f {ref_fasta} {bam} | "
               f"bcftools call --threads {threads} --ploidy {ploidy} -f GQ -mv -Ov -o {raw_vcf}")
    # 移除 stderr=subprocess.DEVNULL，让真实报错可以直接在终端显示
    subprocess.run(snp_cmd, shell=True, check=True)

    # 核心修复 2: 更稳健的过滤表达式，防止因为格式偏差导致直接报错
    filter_expr = (
        f"QUAL < {qual_th} || "
        f"INFO/DP < {dp_th} || "
        f"FORMAT/GQ < {gq_th} || "
        f"MQ < {mq_th} || "
        f"FORMAT/SP > {sp_th}"
    )
    
    # 仅当 AF 大于 0 时才加入 AF 过滤逻辑，提高灵活性
    if float(af_th) > 0:
        filter_expr += f" || (FORMAT/AD[0:1]/FORMAT/DP) < {af_th}"

    filter_cmd = (f"bcftools filter --threads {threads} -e '{filter_expr}' -s LOWQUAL -m + {raw_vcf} | "
                  f"bcftools view --threads {threads} -f PASS > {filtered_vcf}")
    subprocess.run(filter_cmd, shell=True, check=True)

    # --- 2. SV 结构变异检测 ---
    sv_vcf = os.path.join(sample_dir, f"{sample_name}_SV_delly.vcf")
    sv_cmd = f"OMP_NUM_THREADS={threads} delly call -g {ref_fasta} -o {sv_vcf} {bam}"
    subprocess.run(sv_cmd, shell=True, stderr=subprocess.DEVNULL)

    print(f"✅ 样本 {sample_name} 分析完成！")
    return sample_name

def run_pipeline():
    base_dir = get_base_dir()
    
    # 1. 扫描 FASTA
    print("\n🔍 正在当前目录及其子目录下扫描参考基因组与样本文件，请稍候...")
    ref_list = glob.glob(os.path.join(base_dir, "**", "*.fasta"), recursive=True)
    selected_refs = interactive_select(ref_list, "参考基因组 (.fasta)", base_dir)
    if not selected_refs: return
    ref_fasta = selected_refs[0]

    # 2. 扫描 BAM
    bam_list = glob.glob(os.path.join(base_dir, "**", "*.bam"), recursive=True)
    selected_bams = interactive_select(bam_list, "待检测 BAM 样本", base_dir)
    if not selected_bams: return

    # 3. 获取交互配置
    ploidy, thresholds, max_workers, threads_per_task = get_run_configs()

    # 4. 创建输出目录
    timestamp = datetime.now().strftime("%Y%m%d_%H%M")
    out_root = os.path.join(base_dir, f"06_Variant_Results_{timestamp}")
    os.makedirs(out_root, exist_ok=True)

    if not os.path.exists(ref_fasta + ".fai"):
        print("\n⚙️ 正在建立参考基因组索引...")
        subprocess.run(f"samtools faidx {ref_fasta}", shell=True)

    print(f"\n🔥 启动多线程引擎: {max_workers} 个样本并发，单样本 {threads_per_task} 核，设定倍性: {ploidy}...")
    
    with concurrent.futures.ThreadPoolExecutor(max_workers=max_workers) as executor:
        futures = [executor.submit(process_single_bam, bam, ref_fasta, out_root, ploidy, thresholds, threads_per_task) for bam in selected_bams]
        
        for future in concurrent.futures.as_completed(futures):
            try:
                future.result() 
            except subprocess.CalledProcessError as e:
                print(f"❌ 进程执行失败，返回码 {e.returncode}")
            except Exception as exc:
                print(f"❌ 某个样本处理时发生错误: {exc}")

    print(f"\n🎉 所有变异检测已完成！结果存放在: {out_root}")

if __name__ == "__main__":
    try:
        check_env()
        run_pipeline()
    except KeyboardInterrupt:
        print("\n🛑 用户强制退出。")