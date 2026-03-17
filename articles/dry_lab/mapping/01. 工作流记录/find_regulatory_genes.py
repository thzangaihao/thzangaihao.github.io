#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import sys
import glob
from datetime import datetime
from collections import defaultdict

'''
诱变调控区 (Promoter/UTR) 靶标联合分析工具 (自动保存报告版)
功能：读取全量 _annotated.vcf，专门提取上游调控区发生突变的共有基因。
更新：新增自动创建时间戳文件夹，并将完整结果导出为 TXT 报告。
'''

def get_base_dir():
    return os.path.dirname(os.path.abspath(__file__))

def parse_snpeff_full_vcf(vcf_file):
    """解析全量 VCF，专门提取调控区突变"""
    mutated_genes = defaultdict(list)
    
    # 核心目标：只关注可能影响基因表达量的调控区突变
    regulatory_types = ['upstream_gene_variant', '5_prime_UTR_variant', 'promoter']
    
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
                        mutation_type = ann_parts[1] 
                        gene_name = ann_parts[3] 
                        
                        if any(reg_type in mutation_type for reg_type in regulatory_types) and gene_name:
                            detail = f"{chrom}:{pos} ({ref}->{alt}, {mutation_type})"
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
        choice = input(f"\n👉 请选择你要分析的样本 (多选如1,2,3 或all，输入q退出): ").strip().lower()
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

def run_regulatory_analysis():
    base_dir = get_base_dir()
    print("\n" + "="*50)
    print(" 🎛️ 调控区 (Promoter/UTR) 靶标基因挖掘系统")
    print("="*50)
    
    # ⚠️ 强制扫描全量注释文件
    vcf_list = glob.glob(os.path.join(base_dir, "**", "*_annotated.vcf"), recursive=True)
    selected_vcfs = interactive_select(vcf_list, "待分析的全量注释文件 (_annotated.vcf)")
    
    if len(selected_vcfs) < 2:
        print("\n⚠️ 联合分析至少需要选择 2 个以上的独立样本！")
        return

    # === 创建带时间戳的输出目录 ===
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    out_dir = os.path.join(base_dir, f"09_Regulatory_Genes_{timestamp}")
    os.makedirs(out_dir, exist_ok=True)
    report_file = os.path.join(out_dir, "Regulatory_Genes_Report.txt")

    print(f"\n🚀 正在从 {len(selected_vcfs)} 个样本中提取调控区突变...")
    
    gene_to_samples = defaultdict(lambda: defaultdict(list))
    
    for vcf in selected_vcfs:
        sample_name = os.path.basename(vcf).replace('_annotated.vcf', '')
        mutated_genes = parse_snpeff_full_vcf(vcf)
        
        for gene, details in mutated_genes.items():
            gene_to_samples[gene][sample_name].extend(details)

    sorted_genes = sorted(gene_to_samples.items(), key=lambda x: len(x[1]), reverse=True)
    total_samples = len(selected_vcfs)
    
    # 准备写入报告的内容
    report_lines = []
    report_lines.append("="*50)
    report_lines.append(" 🎛️ 候选耐盐调控基因 (上游启动子区) 排行榜")
    report_lines.append(f" 分析时间: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    report_lines.append(f" 参与分析的样本数: {total_samples}")
    report_lines.append("="*50)
    
    found_core_genes = False
    
    for gene, sample_data in sorted_genes:
        hit_count = len(sample_data)
        
        # 只保留至少在 2 个样本中同时出现调控区突变的基因
        if hit_count >= 2:
            found_core_genes = True
            hit_ratio = f"{hit_count}/{total_samples}"
            stars = "⭐" * int((hit_count / total_samples) * 5)
            if hit_count == total_samples: stars = "🌟🌟🌟🌟🌟 (完美交集)"
            
            report_lines.append(f"\n[{stars}] 基因名称: {gene} (调控区命中样本数: {hit_ratio})")
            
            for sample, details in sample_data.items():
                report_lines.append(f"  ├── 样本 [{sample}]")
                for d in set(details):
                    report_lines.append(f"  │    └── {d}")
                    
    if not found_core_genes:
        report_lines.append("\n⚠️ 遗憾：没有发现共有的调控区突变。")
    else:
        report_lines.append("\n" + "="*50)
        report_lines.append("💡 湿实验建议：")
        report_lines.append("排查假阳性后，若发现真实的共同调控靶标，该基因的蛋白质序列可能完好。")
        report_lines.append("建议后续提取 RNA，通过 RT-qPCR 验证其在野生型和突变株间的表达量差异！")

    # 同时输出到屏幕和文件
    with open(report_file, 'w', encoding='utf-8') as f:
        for line in report_lines:
            print(line)
            f.write(line + "\n")

    # 结束提示语（参照 Core 脚本）
    print(f"\n🎉 分析完成！完整的候选调控基因报告已永久保存至: ")
    print(f"📂 {report_file}")

if __name__ == "__main__":
    try:
        run_regulatory_analysis()
    except KeyboardInterrupt:
        print("\n🛑 用户强制退出。")