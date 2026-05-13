#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import sys
import glob
from datetime import datetime
from collections import defaultdict

'''
诱变全样本变异基因汇总工具 (富集分析预处理版)
功能：批量读取多个样本的 HIGH_MODERATE.vcf，提取所有受损基因。
输出：生成包含基因名、突变样本数、样本来源的 TSV 制表符分隔大表，无缝对接 Excel 和富集软件。
'''

def get_base_dir():
    return os.path.dirname(os.path.abspath(__file__))

def parse_snpeff_vcf(vcf_file):
    """解析单份 VCF，提取发生 HIGH/MODERATE 突变的基因"""
    mutated_genes = defaultdict(list)
    
    with open(vcf_file, 'r') as f:
        for line in f:
            if line.startswith('#'): continue
            
            parts = line.strip().split('\t')
            if len(parts) < 8: continue
            
            info = parts[7]
            if 'ANN=' in info:
                ann_str = [x for x in info.split(';') if x.startswith('ANN=')][0]
                annotations = ann_str.replace('ANN=', '').split(',')
                
                for ann in annotations:
                    ann_parts = ann.split('|')
                    if len(ann_parts) >= 4:
                        impact = ann_parts[2]
                        gene_name = ann_parts[3]
                        
                        # 只提取 HIGH 和 MODERATE 级别的破坏性突变
                        if impact in ['HIGH', 'MODERATE'] and gene_name:
                            mutated_genes[gene_name].append(impact)
    
    return mutated_genes

def interactive_select(files, desc):
    """通用的交互式多选逻辑"""
    if not files:
        print(f"⚠️ 未找到任何 {desc}！请确保目录下有 *_HIGH_MODERATE.vcf 文件。")
        return []
    
    print(f"\n📂 扫描到以下 {len(files)} 个 {desc}:")
    for i, f in enumerate(files, 1):
        print(f"  [{i}] {os.path.relpath(f, get_base_dir())}")
        
    while True:
        choice = input(f"\n👉 请选择要汇总的样本 (多选如1-3,5 或all，输入q退出): ").strip().lower()
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
            print("⚠️ 输入无效，请重新输入。")

def run_extraction():
    base_dir = get_base_dir()
    print("\n" + "="*50)
    print(" 📊 全样本变异基因汇总大表生成系统 (GO/KEGG 专用)")
    print("="*50)
    
    # 获取所有的提纯后 VCF 文件
    vcf_list = glob.glob(os.path.join(base_dir, "**", "*_HIGH_MODERATE.vcf"), recursive=True)
    selected_vcfs = interactive_select(vcf_list, "待汇总的变异文件 (_HIGH_MODERATE.vcf)")
    
    if not selected_vcfs:
        return

    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    out_dir = os.path.join(base_dir, f"10_Enrichment_Table_{timestamp}")
    os.makedirs(out_dir, exist_ok=True)
    
    # 使用 .tsv 后缀，Excel 可以直接完美解析分列
    report_file = os.path.join(out_dir, "Mutated_Genes_Summary_for_Enrichment.tsv")

    print(f"\n🚀 正在从 {len(selected_vcfs)} 个样本中提取并汇总基因数据...")
    
    # 结构: { gene_name: set(sample_names) }
    gene_to_samples = defaultdict(set)
    
    for vcf in selected_vcfs:
        sample_name = os.path.basename(vcf).replace('_HIGH_MODERATE.vcf', '')
        mutated_genes = parse_snpeff_vcf(vcf)
        
        for gene in mutated_genes.keys():
            gene_to_samples[gene].add(sample_name)

    # 按照基因命中的样本数量降序排列
    sorted_genes = sorted(gene_to_samples.items(), key=lambda x: len(x[1]), reverse=True)
    total_samples = len(selected_vcfs)
    
    # 写入大表 (Tab 分隔，完美适配 Excel)
    with open(report_file, 'w', encoding='utf-8') as f:
        # 写入表头
        f.write("Gene_ID\tHit_Count\tTotal_Samples\tMutated_In_Samples\n")
        
        for gene, samples in sorted_genes:
            hit_count = len(samples)
            sample_list_str = ", ".join(sorted(list(samples)))
            
            # 写入每一行数据
            f.write(f"{gene}\t{hit_count}\t{total_samples}\t{sample_list_str}\n")

    print("\n" + "="*50)
    print(f"🎉 汇总大表生成成功！共提取到 {len(sorted_genes)} 个独立的突变基因。")
    print(f"📂 报告保存在: {report_file}")
    print("\n💡 下一步操作指南：")
    print("1. 请用 Excel 打开这个 .tsv 文件。")
    print("2. ⚠️ 强烈建议：剔除那些 Hit_Count 等于总样本数（例如 8/8）的基因，它们大概率是参考基因组的背景噪音！")
    print("3. 复制剩余干净的 'Gene_ID' 这一整列。")
    print("4. 将其粘贴到 KOBAS 或 DAVID 等富集分析网站中，开启你的通路挖掘之旅！")

if __name__ == "__main__":
    try:
        run_extraction()
    except KeyboardInterrupt:
        print("\n🛑 用户强制退出。")