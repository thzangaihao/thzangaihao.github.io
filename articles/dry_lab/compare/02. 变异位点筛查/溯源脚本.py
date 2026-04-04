#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import sys
import glob
import csv
import re
from collections import defaultdict
from datetime import datetime

'''
CDS变异终极汇总大表生成器 (精准 ID 匹配版)
功能：整合 VCF、GFF3/GTF 注释与 InterProScan TSV 结果。
修复：彻底解决了 GFF3 (含 gene- 前缀) 与 InterProScan TSV (含复杂 rna-gnl| 管道符前缀) ID 无法匹配的问题。
'''

def get_base_dir():
    return os.path.dirname(os.path.abspath(__file__))

def interactive_select(files, desc, display_root):
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

def clean_feature_id(raw_id):
    """
    【核心修复】：ID 卸妆液。
    剥离不同软件强加的前缀，提取真正的核心基因 ID。
    """
    if not raw_id: return "Unknown"
    
    # 1. 切除复杂的管道符前缀 (如 rna-gnl|WGS:JPDN|mrna.TGAM01_v203339 -> mrna.TGAM01_v203339)
    cleaned = raw_id.split('|')[-1]
    
    # 2. 切除常见的 GFF/NCBI 前缀 (如 gene-TGAM01, mrna.TGAM01, rna-TGAM01, cds-)
    cleaned = re.sub(r'^(gene-|rna-|mrna\.|transcript-|cds-)', '', cleaned, flags=re.IGNORECASE)
    
    return cleaned

def parse_gff(gff_file):
    print(f"\n⚙️ 正在解析基因注释文件: {os.path.basename(gff_file)}...")
    gff_map = defaultdict(list)
    
    with open(gff_file, 'r') as f:
        for line in f:
            if line.startswith('#') or not line.strip(): continue
            parts = line.strip().split('\t')
            if len(parts) < 9: continue
            
            chrom, feature, start, end, attr = parts[0], parts[2], int(parts[3]), int(parts[4]), parts[8]
            
            attr_dict = {}
            for item in attr.strip().split(';'):
                if '=' in item:
                    k, v = item.split('=', 1)
                    attr_dict[k.strip()] = v.strip(' "')
                elif ' ' in item:
                    k, v = item.split(' ', 1)
                    attr_dict[k.strip()] = v.strip(' "')

            gene_id = "Unknown"
            
            if feature in ['gene', 'mRNA', 'transcript']:
                gene_id = attr_dict.get('ID') or attr_dict.get('Name') or attr_dict.get('locus_tag') or "Unknown"
            elif feature == 'CDS':
                parent = attr_dict.get('Parent')
                if parent:
                    gene_id = parent.split(',')[0]
                else:
                    gene_id = attr_dict.get('ID') or "Unknown"

            if gene_id != "Unknown":
                # 存入字典前，彻底清洗 ID
                clean_gene = clean_feature_id(gene_id)
                gff_map[chrom].append((start, end, clean_gene))
    
    print(f"✅ GFF 解析完毕，共建立 {sum(len(v) for v in gff_map.values())} 个基因组坐标映射。")
    return gff_map

def parse_interpro(tsv_file):
    print(f"⚙️ 正在解析 InterProScan 注释文件: {os.path.basename(tsv_file)}...")
    ipr_map = defaultdict(set)
    
    with open(tsv_file, 'r') as f:
        for line in f:
            if not line.strip(): continue
            parts = line.split('\t')
            
            # 存入字典前，用同样的规则清洗 IPR TSV 里的 ID
            protein_id = clean_feature_id(parts[0].strip())
            desc = ""
            
            if len(parts) > 12 and parts[12].strip() != '-':
                desc = parts[12].strip()
            elif len(parts) > 5 and parts[5].strip() != '-':
                desc = parts[5].strip()
                
            if desc:
                ipr_map[protein_id].add(desc)
                
    final_map = {k: " | ".join(v) for k, v in ipr_map.items()}
    print(f"✅ InterProScan 解析完毕，共提取 {len(final_map)} 个蛋白的注释信息。")
    return final_map

def extract_dp(info_str, format_str, sample_str):
    if format_str and sample_str and "DP" in format_str:
        try:
            dp_idx = format_str.split(':').index('DP')
            return sample_str.split(':')[dp_idx]
        except:
            pass
            
    for item in info_str.split(';'):
        if item.startswith('DP='):
            return item.split('=')[1]
            
    return "N/A"

def get_gene_and_anno(chrom, pos, gff_map, ipr_map):
    pos = int(pos)
    matched_gene = "Intergenic"
    
    if chrom in gff_map:
        for start, end, gene_id in gff_map[chrom]:
            if start <= pos <= end:
                matched_gene = gene_id
                break
                
    anno = ipr_map.get(matched_gene, "无功能注释")
    
    # 后缀模糊匹配容错
    if anno == "无功能注释" and matched_gene != "Intergenic":
        possible_bases = [
            matched_gene.split('.')[0],
            matched_gene.rsplit('_', 1)[0]
        ]
        for base in possible_bases:
            if base in ipr_map:
                anno = ipr_map[base]
                break
                
    return matched_gene, anno

def run_pipeline():
    base_dir = get_base_dir()
    
    vcf_list = glob.glob(os.path.join(base_dir, "**", "*_CDS.vcf"), recursive=True)
    selected_vcfs = interactive_select(vcf_list, "CDS 变异文件 (*_CDS.vcf)", base_dir)
    if not selected_vcfs: return
    
    anno_list = glob.glob(os.path.join(base_dir, "**", "*.gff*"), recursive=True) + glob.glob(os.path.join(base_dir, "**", "*.gtf"), recursive=True)
    selected_anno = interactive_select(anno_list, "基因注释文件 (.gff3 / .gtf)", base_dir)
    if not selected_anno: return
    
    tsv_list = glob.glob(os.path.join(base_dir, "**", "*.tsv"), recursive=True)
    if not tsv_list:
        print("\n⚠️ 未扫描到 InterProScan TSV 文件，将跳过功能注释列。")
        selected_tsv = []
    else:
        selected_tsv = interactive_select(tsv_list, "InterProScan 注释文件 (.tsv) [可按 q 退出或直接选]", base_dir)
    
    gff_map = parse_gff(selected_anno[0])
    ipr_map = parse_interpro(selected_tsv[0]) if selected_tsv else {}
    
    timestamp = datetime.now().strftime("%Y%m%d_%H%M")
    out_csv = os.path.join(base_dir, f"08_Mutations_Summary_{timestamp}.csv")
    
    headers = [
        "Sample_Name", "Chromosome", "Position", "Mutation_Type", 
        "REF_Allele", "ALT_Allele", "Depth_(DP)", "Gene_ID", "InterPro_Annotation"
    ]
    
    print(f"\n🚀 正在整合所有样本的变异数据...")
    total_mutations = 0
    
    with open(out_csv, 'w', newline='', encoding='utf-8-sig') as f_out:
        writer = csv.writer(f_out)
        writer.writerow(headers)
        
        for vcf in selected_vcfs:
            sample_name = os.path.basename(vcf).replace('_CDS.vcf', '').replace('_Gene.vcf', '')
            
            with open(vcf, 'r') as f_in:
                for line in f_in:
                    if line.startswith('#'): continue
                    
                    parts = line.strip().split('\t')
                    if len(parts) < 10: continue
                    
                    chrom = parts[0]
                    pos = parts[1]
                    ref = parts[3]
                    alt = parts[4]
                    info_str = parts[7]
                    format_str = parts[8]
                    sample_str = parts[9]
                    
                    mut_type = "SNP" if len(ref) == 1 and len(alt) == 1 else "Indel"
                    dp = extract_dp(info_str, format_str, sample_str)
                    
                    gene_id, annotation = get_gene_and_anno(chrom, pos, gff_map, ipr_map)
                    
                    writer.writerow([
                        sample_name, chrom, pos, mut_type,
                        ref, alt, dp, gene_id, annotation
                    ])
                    total_mutations += 1

    print(f"\n🎉 汇总大表生成完毕！\n   共写入 {total_mutations} 条变异记录。\n   请查收文件: {out_csv}")

if __name__ == "__main__":
    try:
        run_pipeline()
    except KeyboardInterrupt:
        print("\n🛑 用户强制退出。")