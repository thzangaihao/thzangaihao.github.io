#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import sys
import glob
import subprocess
import re
from datetime import datetime

'''
特定基因组序列精准提取工具
功能：输入多个基因名称，基于 GFF 注释文件在参考基因组中进行精准切片。
特性：
1. 继承交互式文件选择逻辑。
2. 自动处理链的正负性，负链基因自动进行反向互补 (Reverse Complement)。
3. 严格保证输出序列统一为 5' -> 3' 生物学方向。
'''

def get_base_dir():
    return os.path.dirname(os.path.abspath(__file__))

def check_env():
    if subprocess.call(['which', 'samtools'], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL) != 0:
        print("❌ 错误: 未找到 samtools！请执行: conda install -c bioconda samtools")
        sys.exit(1)

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
    """剔除软件强制添加的冗余前缀，保证 ID 干净纯粹"""
    if not raw_id: return "Unknown"
    cleaned = raw_id.split('|')[-1]
    cleaned = re.sub(r'^(gene-|rna-|mrna\.|transcript-|cds-)', '', cleaned, flags=re.IGNORECASE)
    return cleaned

def reverse_complement(seq):
    """DNA 序列反向互补函数"""
    trans = str.maketrans('ATCGatcgNn', 'TAGCtagcNn')
    # 先替换互补碱基，再将序列反转 [::-1]
    return seq.translate(trans)[::-1]

def get_target_genes():
    """获取用户想要提取的基因列表"""
    print("\n🎯 --- 输入目标基因 ---")
    print("请输入你想提取的基因名称，多个基因请用逗号 (,) 分隔。")
    print("例如: TGAM01_v208328, TGAM01_v203339")
    
    genes_input = input("👉 请输入: ").strip()
    if not genes_input:
        print("⚠️ 未输入任何基因，程序退出。")
        sys.exit(0)
        
    # 清洗输入并去重
    targets = {clean_feature_id(g.strip()) for g in genes_input.split(',') if g.strip()}
    return targets

def parse_gff_for_targets(gff_file, target_genes):
    """扫描注释文件，抓取目标基因的坐标与链方向"""
    print(f"\n⚙️ 正在解析注释文件寻找目标... ")
    found_genes = {} # { gene_id: (chrom, start, end, strand) }
    
    with open(gff_file, 'r') as f:
        for line in f:
            if line.startswith('#') or not line.strip(): continue
            parts = line.strip().split('\t')
            if len(parts) < 9: continue
            
            chrom, feature, start, end, strand, attr = parts[0], parts[2], parts[3], parts[4], parts[6], parts[8]
            
            # 我们提取 'gene' 或 'mRNA' 的整体范围
            if feature in ['gene', 'mRNA', 'transcript']:
                attr_dict = {}
                for item in attr.strip().split(';'):
                    if '=' in item:
                        k, v = item.split('=', 1)
                        attr_dict[k.strip()] = v.strip(' "')
                    elif ' ' in item:
                        k, v = item.split(' ', 1)
                        attr_dict[k.strip()] = v.strip(' "')

                gene_id = attr_dict.get('ID') or attr_dict.get('Name') or attr_dict.get('locus_tag') or "Unknown"
                clean_id = clean_feature_id(gene_id)
                
                # 如果这个基因在我们的目标列表里，且尚未被记录
                if clean_id in target_genes and clean_id not in found_genes:
                    found_genes[clean_id] = (chrom, start, end, strand)
                    
    print(f"✅ 共输入 {len(target_genes)} 个目标基因，成功在 GFF 中找到 {len(found_genes)} 个。")
    if len(target_genes) > len(found_genes):
        missing = target_genes - set(found_genes.keys())
        print(f"⚠️ 以下基因未在注释文件中找到: {', '.join(missing)}")
        
    return found_genes

def extract_sequences(fasta_file, found_genes, out_fasta):
    """调用 samtools faidx 提取序列并处理正负链"""
    print(f"\n🧬 正在从参考基因组中提取序列...")
    
    with open(out_fasta, 'w') as f_out:
        for gene_id, (chrom, start, end, strand) in found_genes.items():
            # 组装 samtools 命令
            region = f"{chrom}:{start}-{end}"
            cmd = f"samtools faidx {fasta_file} {region}"
            
            try:
                # 捕获 samtools 输出
                result = subprocess.run(cmd, shell=True, check=True, stdout=subprocess.PIPE, text=True)
                lines = result.stdout.strip().split('\n')
                
                if len(lines) > 1:
                    # 合并跨行的序列
                    raw_seq = "".join(lines[1:])
                    
                    # 核心逻辑：负链基因的反向互补处理
                    if strand == '-':
                        final_seq = reverse_complement(raw_seq)
                        direction_note = "minus_strand_reverse_complemented"
                    else:
                        final_seq = raw_seq
                        direction_note = "plus_strand"
                        
                    # 写入标准 FASTA 格式，Header 中保留完备的位置与链信息
                    header = f">{gene_id} {chrom}:{start}-{end} strand:{strand} info:{direction_note}\n"
                    f_out.write(header)
                    
                    # 格式化序列：每 80 个字符换行，保证 FASTA 文件的可读性
                    for i in range(0, len(final_seq), 80):
                        f_out.write(final_seq[i:i+80] + "\n")
                        
            except subprocess.CalledProcessError:
                print(f"❌ 提取 {gene_id} ({region}) 时发生错误，请检查坐标是否越界。")
                
    print(f"\n🎉 提取完成！所有序列均已校正为 5' -> 3' 方向。\n   结果已保存至: {out_fasta}")

def run_pipeline():
    base_dir = get_base_dir()
    
    # 1. 选择参考基因组
    fasta_list = glob.glob(os.path.join(base_dir, "**", "*.fasta"), recursive=True) + \
                 glob.glob(os.path.join(base_dir, "**", "*.fa"), recursive=True)
    selected_fasta = interactive_select(fasta_list, "参考基因组 (.fasta/.fa)", base_dir)
    if not selected_fasta: return
    ref_fasta = selected_fasta[0]
    
    # 建立/检查索引
    if not os.path.exists(ref_fasta + ".fai"):
        print("\n⚙️ 正在建立参考基因组索引...")
        subprocess.run(f"samtools faidx {ref_fasta}", shell=True)

    # 2. 选择注释文件
    anno_list = glob.glob(os.path.join(base_dir, "**", "*.gff*"), recursive=True) + \
                glob.glob(os.path.join(base_dir, "**", "*.gtf"), recursive=True)
    selected_anno = interactive_select(anno_list, "基因注释文件 (.gff3 / .gtf)", base_dir)
    if not selected_anno: return
    gff_file = selected_anno[0]

    # 3. 输入目标基因
    target_genes = get_target_genes()

    # 4. 查找坐标
    found_genes = parse_gff_for_targets(gff_file, target_genes)
    if not found_genes:
        print("🛑 没有找到任何有效基因记录，程序终止。")
        return

    # 5. 提取并生成文件
    timestamp = datetime.now().strftime("%Y%m%d_%H%M")
    out_fasta = os.path.join(base_dir, f"09_Target_Genes_{timestamp}.fasta")
    
    extract_sequences(ref_fasta, found_genes, out_fasta)

if __name__ == "__main__":
    try:
        check_env()
        run_pipeline()
    except KeyboardInterrupt:
        print("\n🛑 用户强制退出。")