#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import sys
import glob
import subprocess
import time
from datetime import datetime

'''
通用 De novo 组装流水线 (SPAdes 批量排队版)
功能：自动化完成 Fastp 质控 -> SPAdes 组装 -> QUAST 评估
特性：全动态扫描双端文件、交互式多选排队、结果统一存放在带时间戳的主文件夹及各样本独立的子文件夹下。
'''

def get_base_dir():
    return os.path.dirname(os.path.abspath(__file__))

def check_env():
    tools = ['fastp', 'spades.py', 'quast.py']
    for tool in tools:
        if subprocess.call(['which', tool], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL) != 0:
            print(f"❌ 错误: 未找到 {tool}，请确认已在 Conda 环境中安装并激活！")
            sys.exit(1)

def choose_samples():
    """交互式选择测序样本 (支持任意双端 fastq/fq)"""
    base_dir = get_base_dir()
    
    # 扩大扫描范围：支持 .fq.gz, .fastq.gz 以及未压缩的 .fq, .fastq
    patterns = [
        '*_R1.fq.gz', '*_1.fq.gz', '*_R1.fastq.gz', '*_1.fastq.gz',
        '*_R1.fq', '*_1.fq', '*_R1.fastq', '*_1.fastq'
    ]
    
    raw_files = []
    for pat in patterns:
        raw_files.extend(glob.glob(os.path.join(base_dir, '**', pat), recursive=True))
    
    valid_samples = []
    for f in raw_files:
        basename = os.path.basename(f)
        # 防污染机制 (排除已经清洗过的 clean 数据)
        if 'clean' not in basename.lower():
            valid_samples.append(f)
            
    valid_samples = sorted(list(set(valid_samples)))
    
    if not valid_samples:
        print("⚠️ 提示：未在当前目录及任何子目录中找到原始双端测序样本 (R1文件)！")
        return []

    print(f"\n📂 自动扫描到以下 {len(valid_samples)} 个测序样本:")
    for i, f in enumerate(valid_samples, 1):
        rel_path = os.path.relpath(f, base_dir)
        print(f"  [{i}] {rel_path}")

    while True:
        try:
            choice = input("\n👉 请输入要组装的样本编号 (单选:1, 多选:1,3, 范围:1-3, 全部:all, 退出:q): ").strip().lower()
            if choice == 'q':
                sys.exit(0)
            
            selected_paths = []
            
            if choice == 'all':
                selected_paths = valid_samples
            else:
                parts = choice.replace(' ', '').split(',')
                indices = set()
                for part in parts:
                    if '-' in part:
                        start, end = map(int, part.split('-'))
                        indices.update(range(start, end + 1))
                    else:
                        indices.add(int(part))
                
                for idx in sorted(list(indices)):
                    if 1 <= idx <= len(valid_samples):
                        selected_paths.append(valid_samples[idx-1])
                    else:
                        print(f"⚠️ 警告: 编号 {idx} 超出范围，已忽略。")
            
            if not selected_paths:
                print("⚠️ 未选中任何有效样本，请重新输入。")
                continue

            # 提取信息并自动推导对应的 R2 文件
            final_samples = []
            for r1_path in selected_paths:
                # 智能替换，寻找 R2
                r2_path = r1_path.replace('_R1.fq', '_R2.fq')\
                                 .replace('_1.fq', '_2.fq')\
                                 .replace('_R1.fastq', '_R2.fastq')\
                                 .replace('_1.fastq', '_2.fastq')
                
                if not os.path.exists(r2_path):
                    print(f"❌ 跳过: 找到了 R1，但找不到对应的 R2 文件: {os.path.basename(r2_path)}")
                    continue
                
                # 智能提取纯粹的样本名
                sample_name = os.path.basename(r1_path)
                for suffix in ['_R1.fq.gz', '_1.fq.gz', '_R1.fastq.gz', '_1.fastq.gz', 
                               '_R1.fq', '_1.fq', '_R1.fastq', '_1.fastq']:
                    sample_name = sample_name.replace(suffix, '')
                
                final_samples.append((sample_name, r1_path, r2_path))
                
            print(f"\n✅ 已成功选中 {len(final_samples)} 个样本进入排队序列。")
            return final_samples

        except ValueError:
            print("⚠️ 输入格式错误！请确保输入的是数字编号、范围或 'all'。")

def get_system_resources():
    print("\n💻 --- 运算资源配置 (将应用于每一个排队的任务) ---")
    while True:
        try:
            threads = input("👉 请输入使用的 CPU 核心数 (建议填 4-8): ").strip()
            threads = int(threads)
            if threads > 0: break
            else: print("⚠️ 线程数必须大于0！")
        except ValueError: print("⚠️ 请输入有效的整数。")
            
    while True:
        try:
            memory = input("👉 请输入最大内存 (单位 GB, 建议填 16-32): ").strip()
            memory = int(memory)
            if memory > 0: break
            else: print("⚠️ 内存必须大于0！")
        except ValueError: print("⚠️ 请输入有效的整数。")
            
    return threads, memory

def run_pipeline(sample_name, r1, r2, threads, memory, current_idx, total_samples, batch_out_dir):
    # 为当前样本在时间戳主目录下创建独立的子文件夹
    out_base = os.path.join(batch_out_dir, sample_name)
    os.makedirs(out_base, exist_ok=True)
    
    print("\n" + "="*60)
    print(f"🚀 [任务队列 {current_idx}/{total_samples}] 开始处理: {sample_name}")
    print(f"🕒 当前任务开始时间: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print("="*60)
    
    try:
        # 1. Fastp 数据清洗
        print(f"\n--- [1/3] 执行 Fastp 数据清洗 ---")
        clean_r1 = os.path.join(out_base, f"{sample_name}_clean_R1.fq.gz")
        clean_r2 = os.path.join(out_base, f"{sample_name}_clean_R2.fq.gz")
        fastp_html = os.path.join(out_base, f"{sample_name}_fastp.html")
        fastp_cmd = f"fastp -i {r1} -I {r2} -o {clean_r1} -O {clean_r2} -w {min(threads, 16)} -h {fastp_html}"
        subprocess.run(fastp_cmd, shell=True, check=True)

        # 2. SPAdes 从头组装
        print(f"\n--- [2/3] 启动 SPAdes 从头组装 ---")
        spades_out = os.path.join(out_base, "spades_run")
        spades_cmd = f"spades.py -1 {clean_r1} -2 {clean_r2} -o {spades_out} -t {threads} -m {memory}"
        subprocess.run(spades_cmd, shell=True, check=True)

        # 3. 提取结果与 QUAST 评估
        print(f"\n--- [3/3] 提取组装结果并生成评估报告 ---")
        final_fasta = os.path.join(out_base, f"{sample_name}_ref_draft.fasta")
        
        # 将 spades 核心结果拷贝出来重命名，方便查看
        scaffolds_path = os.path.join(spades_out, "scaffolds.fasta")
        if os.path.exists(scaffolds_path):
            subprocess.run(f"cp {scaffolds_path} {final_fasta}", shell=True, check=True)
            
            # 运行 QUAST
            quast_out = os.path.join(out_base, "quast_report")
            quast_cmd = f"quast.py {final_fasta} -o {quast_out} -t {threads}"
            subprocess.run(quast_cmd, shell=True, check=True)
            print(f"\n✅ {sample_name} 组装及评估完成！")
        else:
            print(f"⚠️ 警告: SPAdes 运行结束，但未找到 scaffolds.fasta 文件，跳过 QUAST。")

        print(f"📁 样本结果保存在: {out_base}")

    except subprocess.CalledProcessError as e:
        print(f"\n❌ {sample_name} 任务在运行过程中出错: {e}")
        print("⚠️ 脚本将跳过此样本，继续处理队列中的下一个样本。")

if __name__ == "__main__":
    try:
        check_env()
        samples_to_run = choose_samples()
        
        if samples_to_run:
            threads, memory = get_system_resources()
            
            # 【核心修改】创建带有时间戳的主输出文件夹
            timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
            batch_out_dir = os.path.join(get_base_dir(), f"02_DeNovo_Assembly_{timestamp}")
            os.makedirs(batch_out_dir, exist_ok=True)
            
            print("\n" + "#"*60)
            print(f"📁 本次批量排队任务的主输出目录为:\n   {batch_out_dir}")
            print("#"*60)
            
            start_time = time.time()
            total = len(samples_to_run)
            
            # 这里的 for 循环就是你的“队列”，它会按顺序逐个完成流水线任务
            for i, (sample_name, r1, r2) in enumerate(samples_to_run, 1):
                run_pipeline(sample_name, r1, r2, threads, memory, i, total, batch_out_dir)
            
            end_time = time.time()
            hours, rem = divmod(end_time - start_time, 3600)
            minutes, seconds = divmod(rem, 60)
            
            print("\n" + "🎉 "*15)
            print(f"🌟 全部 {total} 个排队样本处理完毕！")
            print(f"⏱️ 总计耗时: {int(hours)}小时 {int(minutes)}分钟 {int(seconds)}秒")
            print("🎉 "*15 + "\n")
            
    except KeyboardInterrupt:
        print("\n\n🛑 程序被手动中断 (Ctrl+C)。退出运行。")
        sys.exit(0)