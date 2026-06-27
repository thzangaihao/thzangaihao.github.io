#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import sys
import glob
import subprocess
import time
import multiprocessing
from collections import OrderedDict

'''
🧬 RSEM 转录组定量批量队列自动化流程 (Batch Queue Version)
修复：完美适配 .gz 文件，利用 bash 进程替换动态传入数据流，防崩溃。
升级：采用严格的尾部后缀匹配算法，彻底解决样本名自带 _1 或 .R2 被误认的 Bug。
'''

def get_base_dir():
    return os.path.dirname(os.path.abspath(__file__))

def check_dependency():
    deps = ['rsem-prepare-reference', 'rsem-calculate-expression', 'bowtie2', 'gunzip']
    missing = [dep for dep in deps if subprocess.call(['which', dep], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL) != 0]
    if missing:
        print(f"❌ 错误：未找到以下依赖软件：{', '.join(missing)}")
        sys.exit(1)

def interactive_select(all_files, desc, valid_exts):
    valid_files = [f for f in all_files if f.lower().endswith(valid_exts)]
    if not valid_files:
        path = input(f"👉 未扫描到{desc}，请手动输入绝对路径 (q 退出): ").strip()
        if path.lower() == 'q': sys.exit(0)
        if os.path.isfile(path): return path
        print("❌ 文件不存在。")
        sys.exit(1)

    print(f"\n📂 扫描到以下 {len(valid_files)} 个候选 {desc}:")
    for i, f in enumerate(valid_files, 1):
        print(f"  [{i}] {os.path.relpath(f, get_base_dir())}")
        
    while True:
        choice = input(f"👉 请选择 {desc} (输入编号或绝对路径，q 退出): ").strip()
        if choice.lower() == 'q': sys.exit(0)
        if choice.isdigit() and 1 <= int(choice) <= len(valid_files):
            return valid_files[int(choice)-1]
        if os.path.isfile(choice): return choice
        print("⚠️ 输入无效，请重新输入。")

def build_sample_queue(fastq_dir, is_paired):
    """核心队列逻辑：扫描文件夹并严格配对样本"""
    all_fastqs = glob.glob(os.path.join(fastq_dir, "**", "*.fq"), recursive=True) + \
                 glob.glob(os.path.join(fastq_dir, "**", "*.fastq"), recursive=True) + \
                 glob.glob(os.path.join(fastq_dir, "**", "*.fq.gz"), recursive=True) + \
                 glob.glob(os.path.join(fastq_dir, "**", "*.fastq.gz"), recursive=True)
                 
    sample_queue = OrderedDict()
    
    # 扩展的后缀匹配列表，优先匹配最长的后缀
    valid_extensions = ['.raw.fastq.gz', '.raw.fq.gz', '.fastq.gz', '.fq.gz', '.raw.fastq', '.raw.fq', '.fastq', '.fq']
    
    if not is_paired:
        for fq in sorted(all_fastqs):
            basename = os.path.basename(fq)
            for ext in valid_extensions:
                if basename.endswith(ext):
                    sample_name = basename[:-len(ext)]
                    sample_queue[sample_name] = [fq]
                    break
        return sample_queue

    # 双端测序严格配对逻辑 (修复名称中含 _1 的 bug)
    for fq in sorted(all_fastqs):
        basename = os.path.basename(fq)
        
        # 1. 匹配并去除文件扩展名
        matched_ext = ""
        for ext in valid_extensions:
            if basename.endswith(ext):
                matched_ext = ext
                break
        
        if not matched_ext:
            continue
            
        name_without_ext = basename[:-len(matched_ext)]
        
        # 2. 严格检查 Read 1 标识符（必须在去除后缀后的最末尾）
        r1_suffix, r2_suffix = "", ""
        if name_without_ext.endswith('.R1'):
            r1_suffix, r2_suffix = '.R1', '.R2'
        elif name_without_ext.endswith('_R1'):
            r1_suffix, r2_suffix = '_R1', '_R2'
        elif name_without_ext.endswith('_1'):
            r1_suffix, r2_suffix = '_1', '_2'
        else:
            # 如果不是 R1 文件（例如是 R2 文件，或者未识别格式），直接跳过
            # 我们只需要找到 R1，然后去反推 R2 的位置即可，这样可以避免 R2 文件报错
            continue
            
        # 3. 提取干净的样本名并寻找 R2
        sample_name = name_without_ext[:-len(r1_suffix)]
        expected_r2_basename = sample_name + r2_suffix + matched_ext
        expected_r2_path = os.path.join(os.path.dirname(fq), expected_r2_basename)
        
        if os.path.exists(expected_r2_path):
            sample_queue[sample_name] = [fq, expected_r2_path]
        else:
            print(f"⚠️ 警告: 找不到 {basename} 对应的 Read 2 文件 ({expected_r2_basename})，已跳过该样本。")
            
    return sample_queue

def run_cmd(cmd, desc):
    """执行系统命令，强制使用 bash 以支持高级语法"""
    print(f"\n▶️ [{time.strftime('%H:%M:%S')}] 正在执行: {desc}")
    print(f"   💻 命令: {cmd}") 
    try:
        subprocess.run(cmd, shell=True, check=True, executable='/bin/bash')
        print(f"✅ {desc} 完成！")
    except subprocess.CalledProcessError:
        print(f"\n❌ 错误：{desc} 失败。请检查上方的报错信息。")
        sys.exit(1)

def main():
    print("=" * 60)
    print(" 🧬 RSEM 批量队列自动化流程 (严格配对版)")
    print("=" * 60)

    check_dependency()
    base_dir = get_base_dir()
    
    print("🔍 正在扫描当前目录文件...")
    all_files = glob.glob(os.path.join(base_dir, "**", "*.*"), recursive=True)

    # 1. 配置参考基因组
    print("\n" + "-"*40)
    print(" 🛠️ 步骤 1：配置参考基因组与注释")
    genome_fasta = interactive_select(all_files, "参考基因组文件", ('.fa', '.fasta', '.fna'))
    annotation_gtf = interactive_select(all_files, "基因组注释文件", ('.gtf', '.gff', '.gff3'))

    ref_basename = os.path.basename(genome_fasta).rsplit('.', 1)[0]
    ref_dir = os.path.join(base_dir, f"{ref_basename}_rsem_idx")
    ref_prefix = os.path.join(ref_dir, "ref")

    # 2. 配置测序数据目录与队列
    print("\n" + "-"*40)
    print(" 🛠️ 步骤 2：配置测序数据与队列")
    fastq_dir = input(f"👉 请输入存放测序数据的文件夹路径 (直接回车扫描当前目录): ").strip()
    if not fastq_dir: fastq_dir = base_dir
    
    is_paired = input("👉 数据是双端测序 (Paired-end) 吗？(y/n, 默认 y): ").strip().lower() != 'n'
    
    # 构建任务队列
    sample_queue = build_sample_queue(fastq_dir, is_paired)
    
    if not sample_queue:
        print("❌ 在指定目录中未找到任何匹配的测序文件。")
        sys.exit(1)

    print(f"\n📋 成功构建任务队列，共扫描到 \033[92m{len(sample_queue)}\033[0m 个样本：")
    for i, (sample, files) in enumerate(sample_queue.items(), 1):
        print(f"  [{i}] {sample} -> {len(files)} 个文件")

    # 3. 性能设置与输出目录
    print("\n" + "-"*40)
    threads = max(1, int(multiprocessing.cpu_count() * 0.8))
    user_threads = input(f"👉 请输入 RSEM 调用的线程数 (直接回车使用 {threads}): ").strip()
    if user_threads.isdigit(): threads = int(user_threads)

    out_dir = os.path.join(base_dir, "RSEM_Results")
    if not os.path.exists(out_dir): os.makedirs(out_dir)

    print("\n" + "="*60)
    print("🚀 队列执行即将开始")
    print("="*60)

    # 4. 构建参考库
    if not os.path.exists(ref_dir):
        os.makedirs(ref_dir)
        gtf_flag = "--gff3" if annotation_gtf.lower().endswith(('.gff', '.gff3')) else "--gtf"
        build_cmd = f"rsem-prepare-reference --num-threads {threads} {gtf_flag} '{annotation_gtf}' --bowtie2 '{genome_fasta}' '{ref_prefix}'"
        run_cmd(build_cmd, "构建 RSEM 参考基因组库")
    else:
        print(f"\n✅ 检测到参考基因组索引已存在，跳过建库步骤。")

    # 5. 执行任务队列
    for idx, (sample_name, files) in enumerate(sample_queue.items(), 1):
        output_prefix = os.path.join(out_dir, sample_name)
        result_file = f"{output_prefix}.genes.results"
        
        print("\n" + "-"*40)
        print(f"⏳ 正在处理队列任务 [{idx}/{len(sample_queue)}]: \033[93m{sample_name}\033[0m")
        
        if os.path.exists(result_file):
            print(f"⏭️ 检测到 {sample_name} 的结果文件已存在，跳过该样本。")
            continue

        rsem_params = f"--num-threads {threads} --bowtie2"
        if is_paired: rsem_params += " --paired-end"
        
        if files[0].endswith('.gz'):
            fq_inputs = " ".join([f"<(gunzip -c '{f}')" for f in files])
        else:
            fq_inputs = " ".join([f"'{f}'" for f in files])
            
        quant_cmd = f"rsem-calculate-expression {rsem_params} {fq_inputs} '{ref_prefix}' '{output_prefix}'"
        
        run_cmd(quant_cmd, f"定量计算: {sample_name}")

    print("\n" + "="*60)
    print("🎉 所有队列任务处理完毕！")
    print(f"📂 结果文件保存在: {out_dir}/")
    print("="*60 + "\n")

if __name__ == "__main__":
    main()