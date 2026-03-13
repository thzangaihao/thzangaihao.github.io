#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import sys
import glob
from datetime import datetime
from collections import defaultdict

'''
诱变核心靶标基因联合分析工具 (自动保存报告版)
功能：读取多个样本的 HIGH_MODERATE.vcf，寻找共有突变基因。
更新：新增自动创建时间戳文件夹，并将完整结果导出为 TXT 报告。
'''

def get_base_dir():
    return os.path.dirname(os.path.abspath(__file__))

def parse_snpeff_vcf(vcf_file):
    """解析单份 VCF，提取发生 HIGH/MODERATE 突变的基因及其突变详情"""
    mutated_genes = defaultdict(list)
    
    with open(vcf_file, 'r') as f:
        for line in f:
            if line.startswith('#'): continue
            
            parts = line.strip().split('\t')
            if len(parts) < 8: continue
            
            chrom = parts[0]
            pos = parts[1]
            ref = parts[3]
            alt = parts[4]
            info = parts[7]
            
            if 'ANN=' in info:
                ann_str = [x for x in info.split(';') if x.startswith('ANN=')][0]
                annotations = ann_str.replace('ANN=', '').split(',')
                
                for ann in annotations:
                    ann_parts = ann.split('|')
                    if len(ann_parts) >= 4:
                        impact = ann_parts[2] # HIGH 或 MODERATE
                        gene_name = ann_parts[3] # 基因名
                        mutation_type = ann_parts[1] # 突变类型
                        
                        if impact in ['HIGH', 'MODERATE'] and gene_name:
                            detail = f"{chrom}:{pos} ({ref}->{alt}, {mutation_type}, {impact})"
                            mutated_genes[gene_name].append(detail)
    
    return mutated_genes

def interactive_select(files, desc):
    """交互式多选"""
    if not files:
        print(f"⚠️ 未找到任何 {desc}！")
        return []
    
    print(f"\n📂 扫描到以下 {len(files)} 个 {desc}:")
    for i, f in enumerate(files, 1):
        print(f"  [{i}] {os.path.relpath(f, get_base_dir())}")
        
    while True:
        choice = input(f"\n👉 请选择你要联合分析的样本 (多选如1,2,3 或all，输入q退出): ").strip().lower()
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

def run_intersection_analysis():
    base_dir = get_base_dir()
    print("\n" + "="*50)
    print(" 🎯 核心靶标基因联合分析系统 (Venn Intersection)")
    print("="*50)
    
    vcf_list = glob.glob(os.path.join(base_dir, "**", "*_HIGH_MODERATE.vcf"), recursive=True)
    selected_vcfs = interactive_select(vcf_list, "待联合分析的提纯变异文件 (.vcf)")
    
    if len(selected_vcfs) < 2:
        print("\n⚠️ 联合分析至少需要选择 2 个以上的独立样本！")
        return

    # === 新增：创建带时间戳的输出目录 ===
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    out_dir = os.path.join(base_dir, f"08_Core_Genes_{timestamp}")
    os.makedirs(out_dir, exist_ok=True)
    report_file = os.path.join(out_dir, "Candidate_Genes_Report.txt")

    print(f"\n🚀 正在从 {len(selected_vcfs)} 个样本中提取破坏性突变基因...")
    
    gene_to_samples = defaultdict(lambda: defaultdict(list))
    
    for vcf in selected_vcfs:
        sample_name = os.path.basename(vcf).replace('_HIGH_MODERATE.vcf', '')
        mutated_genes = parse_snpeff_vcf(vcf)
        
        for gene, details in mutated_genes.items():
            gene_to_samples[gene][sample_name].extend(details)

    sorted_genes = sorted(gene_to_samples.items(), key=lambda x: len(x[1]), reverse=True)
    total_samples = len(selected_vcfs)
    
    # 准备报告内容
    report_lines = []
    report_lines.append("="*50)
    report_lines.append(" 候选耐盐关键基因")
    report_lines.append(f" 分析时间: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    report_lines.append(f" 参与分析的样本数: {total_samples}")
    report_lines.append("="*50)
    
    found_core_genes = False
    
    for gene, sample_data in sorted_genes:
        hit_count = len(sample_data)
        
        # 只记录在至少 2 个样本中共同突变的基因
        if hit_count >= 2:
            found_core_genes = True
            hit_ratio = f"{hit_count}/{total_samples}"
            
            stars = "*" * int((hit_count / total_samples) * 5)
            if hit_count == total_samples: stars = "*****"
            
            report_lines.append(f"\n[{stars}] 基因名称: {gene} (命中样本数: {hit_ratio})")
            
            for sample, details in sample_data.items():
                report_lines.append(f"  ├── 样本 [{sample}]")
                for d in set(details):
                    report_lines.append(f"  │    └── {d}")
                    
    if not found_core_genes:
        report_lines.append("\n⚠️ 在选中的样本中，没有发现任何一个基因发生过共有的破坏性突变。")
    else:
        report_lines.append("\n" + "="*50)
        report_lines.append("💡 下一步湿实验建议：")
        report_lines.append("请排查上述文件中坐标完全一致的 '假阳性' 突变。")
        report_lines.append("寻找 '突变在同一个基因上，但突变坐标不同' 的结果，去设计 sgRNA 或同源重组片段！")

    # 将内容同时输出到屏幕和文件
    with open(report_file, 'w', encoding='utf-8') as f:
        for line in report_lines:
            print(line)
            f.write(line + "\n")

    print(f"\n分析完成！完整的候选基因报告已永久保存至: ")
    print(f"📂 {report_file}")

if __name__ == "__main__":
    try:
        run_intersection_analysis()
    except KeyboardInterrupt:
        print("\n用户强制退出。")