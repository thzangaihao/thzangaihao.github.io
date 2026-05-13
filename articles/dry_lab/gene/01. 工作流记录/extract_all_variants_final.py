#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import sys
import glob
from datetime import datetime
from collections import defaultdict

'''
诱变全变异类型汇总工具 (SNP + SV 终极版)
功能：一键汇总 3 个菌株所有样本的 SNP 和 SV 变异基因。
输出：生成包含基因名、突变类型、影响等级和命中样本数的 TSV 大表。
'''

def get_base_dir():
    return os.path.dirname(os.path.abspath(__file__))

def parse_any_vcf(vcf_file):
    """通用解析器：从 SnpEff 注释后的 VCF 中提取基因及其突变属性"""
    mutated_info = defaultdict(lambda: {'impacts': set(), 'types': set()})
    
    with open(vcf_file, 'r') as f:
        for line in f:
            if line.startswith('#'): continue
            parts = line.strip().split('\t')
            if len(parts) < 8: continue
            
            info = parts[7]
            if 'ANN=' in info:
                # 提取 SnpEff 的 ANN 字段
                ann_str = [x for x in info.split(';') if x.startswith('ANN=')][0]
                annotations = ann_str.replace('ANN=', '').split(',')
                for ann in annotations:
                    ann_parts = ann.split('|')
                    if len(ann_parts) >= 4:
                        m_type = ann_parts[1]   # 突变类型 (如 missense_variant)
                        impact = ann_parts[2]   # 影响级别 (如 HIGH, MODERATE, MODIFIER)
                        gene_name = ann_parts[3] # 基因 ID (如 TGAM01_v210528)
                        
                        if gene_name:
                            mutated_info[gene_name]['impacts'].add(impact)
                            mutated_info[gene_name]['types'].add(m_type)
    return mutated_info

def interactive_select(files, desc):
    """交互式多选逻辑"""
    if not files:
        print(f"⚠️ 未找到任何 {desc}！请确保目录下有 *_annotated.vcf 文件。")
        return []
    
    print(f"\n📂 扫描到以下 {len(files)} 个 {desc}:")
    for i, f in enumerate(files, 1):
        print(f"  [{i}] {os.path.relpath(f, get_base_dir())}")
        
    while True:
        choice = input(f"\n👉 请选择要汇总的样本 (多选如1,3-5 或all，输入q退出): ").strip().lower()
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

def run_final_extraction():
    base_dir = get_base_dir()
    print("\n" + "="*50)
    print(" 🚀 全样本变异基因终极汇总系统 (SNP + SV)")
    print("="*50)
    
    # 扫描所有的注释后 VCF (包含 SNP 和 SV)
    vcf_list = glob.glob(os.path.join(base_dir, "**", "*_annotated.vcf"), recursive=True)
    selected_vcfs = interactive_select(vcf_list, "全量注释文件 (_annotated.vcf)")
    
    if not selected_vcfs:
        return

    # 创建带时间戳的输出目录
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    out_dir = os.path.join(base_dir, f"10_Final_Enrichment_List_{timestamp}")
    os.makedirs(out_dir, exist_ok=True)
    report_file = os.path.join(out_dir, "Final_Combined_Gene_List.tsv")

    print(f"\n🚀 正在从 {len(selected_vcfs)} 个 VCF 文件中提取汇总数据...")
    
    # 核心结构: { gene_id: { 'samples': set(), 'impacts': set(), 'types': set() } }
    grand_summary = defaultdict(lambda: {'samples': set(), 'impacts': set(), 'types': set()})

    for vcf in selected_vcfs:
        # 自动提取样本名 (去除后缀)
        file_name = os.path.basename(vcf)
        s_name = file_name.split('_')[0] 
        
        # 标记是来自 SNP 还是 SV
        v_type_tag = "SV" if "SV" in file_name.upper() else "SNP"
        sample_tag = f"{s_name}({v_type_tag})"
        
        v_data = parse_any_vcf(vcf)
        for gene, info in v_data.items():
            grand_summary[gene]['samples'].add(sample_tag)
            grand_summary[gene]['impacts'].update(info['impacts'])
            grand_summary[gene]['types'].update(info['types'])

    # 按照命中样本数排序
    sorted_genes = sorted(grand_summary.items(), key=lambda x: len(x[1]['samples']), reverse=True)

    # 写入大表
    with open(report_file, 'w', encoding='utf-8') as f:
        f.write("Gene_ID\tHit_Samples_Count\tMax_Impacts\tMutation_Types\tDetailed_Samples\n")
        for gene, data in sorted_genes:
            hit_count = len(data['samples'])
            impact_str = ",".join(sorted(list(data['impacts'])))
            type_str = ",".join(sorted(list(data['types'])))
            sample_str = ",".join(sorted(list(data['samples'])))
            
            f.write(f"{gene}\t{hit_count}\t{impact_str}\t{type_str}\t{sample_str}\n")

    print("\n" + "="*50)
    print(f"🎉 汇总完成！共提取到 {len(grand_summary)} 个独立的变异基因。")
    print(f"📂 最终表格已保存至: {report_file}")
    print("\n💡 建议：打开 Excel 后，先剔除背景假阳性，再进行富集分析。")

if __name__ == "__main__":
    try:
        run_final_extraction()
    except KeyboardInterrupt:
        print("\n🛑 用户强制退出。")