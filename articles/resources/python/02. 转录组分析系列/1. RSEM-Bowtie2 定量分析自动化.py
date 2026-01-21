#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
简化版本：仅保留
1. 文件查找与样本选择功能
2. RSEM + Bowtie2 构建转录本索引
3. RSEM 表达量预测（输出 .genes.results / .isoforms.results）
4. 不包含 HISAT2 基因组比对，不生成 BAM
"""
import os
import sys
import glob
import subprocess
import shutil
import re
import time

# 获取脚本路径

def get_base_dir():
    if getattr(sys, 'frozen', False):
        return os.path.dirname(sys.executable)
    return os.path.dirname(os.path.abspath(__file__))

# 查找文件

def find_files(ext, path=None):
    if path is None:
        path = get_base_dir()
    return glob.glob(os.path.join(path, '**', f'*{ext}'), recursive=True)

# 交互式选择文件

def choose_file(files, desc="文件"):
    if not files:
        return None
    if len(files) == 1:
        print(f"自动选择 {desc}: {files[0]}")
        return files[0]
    print(f"找到 {len(files)} 个 {desc}:")
    for i, f in enumerate(files, 1):
        print(f" [{i}] {f}")
    while True:
        c = input(f"请选择 {desc} 编号: ").strip()
        try:
            return files[int(c)-1]
        except:
            print("无效输入")

# 扫描 fastq

def scan_fastq(ext, path=None):
    if path is None:
        path = get_base_dir()
    fastqs = find_files(ext, path)
    pairs = {}
    for f in fastqs:
        b = os.path.basename(f)
        if '_R1' in b or '_1' in b:
            key = re.split(r'_R1|_1', b)[0]
            pairs.setdefault(key, {})['r1'] = f
        elif '_R2' in b or '_2' in b:
            key = re.split(r'_R2|_2', b)[0]
            pairs.setdefault(key, {})['r2'] = f
    samples = []
    for k, d in pairs.items():
        if 'r1' in d and 'r2' in d:
            samples.append({'name': k, 'r1': d['r1'], 'r2': d['r2']})
    return samples

# 选择样本

def select_samples(samples):
    print("检测到以下样本：")
    for i, s in enumerate(samples, 1):
        print(f" {i}. {s['name']}")
    sel = input("选择样本 (all 或 1,3 或 1-4): ").strip()
    if sel == 'all':
        return samples
    if '-' in sel:
        a, b = map(int, sel.split('-'))
        return samples[a-1:b]
    if ',' in sel:
        idxs = [int(x) for x in sel.split(',')]
        return [samples[i-1] for i in idxs]
    return [samples[int(sel)-1]]

# 构建 RSEM 转录本索引

def build_rsem_reference(annotation, genome, outdir, threads=8):
    # 不需要这一步了，RSEM 会自己生成 transcripts.fa
    # tx_fa = os.path.join(outdir, 'transcripts.fa')
    # subprocess.check_call(['gffread', annotation, '-g', genome, '-w', tx_fa])
    
    rsem_ref = os.path.join(outdir, 'rsem_ref')
    print("构建 RSEM 转录本索引 (Bowtie2)...")
    
    # 关键修改：使用 --gtf (或 --gff3) 选项，并直接传入基因组文件
    cmd = [
        'rsem-prepare-reference', 
        '--bowtie2', 
        '--gtf', annotation,  # 指明注释文件，RSEM 会从中读取基因-转录本对应关系
        genome,               # 这里传入基因组 FASTA，而不是转录本 FASTA
        rsem_ref
    ]
    
    # 注意：如果是 GFF3 格式，需要把 --gtf 换成 --gff3
    if annotation.endswith('.gff3') or annotation.endswith('.gff'):
         cmd[2] = '--gff3'

    subprocess.check_call(cmd)
    return rsem_ref

# RSEM 定量

def rsem_quant(sample, ref, outdir, threads=8):
    name = sample['name']
    r1 = sample['r1']
    r2 = sample['r2']
    sdir = os.path.join(outdir, name)
    os.makedirs(sdir, exist_ok=True)
    print(f"运行 RSEM: {name}")
    cmd = [
        'rsem-calculate-expression',
        '--paired-end', '--bowtie2',
        '-p', str(threads),
        r1, r2,
        ref,
        os.path.join(sdir, name)
    ]
    subprocess.check_call(cmd)

# 主函数

def main():
    base = get_base_dir()
    outdir = os.path.join(base, 'rsem_out_' + time.strftime('%Y%m%d_%H%M%S'))
    os.makedirs(outdir, exist_ok=True)

    # 选择 genome fasta
    fasta = choose_file(find_files('.fa') + find_files('.fasta'), "参考基因组FASTA")

    # 选择注释
    annotation = choose_file(find_files('.gtf') + find_files('.gff') + find_files('.gff3'), "注释文件GTF/GFF")

    # 选择 fastq 格式
    print("选择 fastq 格式：1=.fastq  2=.fastq.gz  3=.fq  4=.fq.gz")
    ext = { '1':'.fastq', '2':'.fastq.gz', '3':'.fq', '4':'.fq.gz' }[ input("输入 1-4: ") ]

    samples = scan_fastq(ext)
    selected = select_samples(samples)

    threads = int(input("线程数（默认8）: ") or 8)

    # 构建 RSEM 索引
    ref = build_rsem_reference(annotation, fasta, outdir, threads)

    # RSEM 定量
    for s in selected:
        rsem_quant(s, ref, outdir, threads)

    print("全部完成！结果位于：", outdir)

if __name__ == '__main__':
    main()
