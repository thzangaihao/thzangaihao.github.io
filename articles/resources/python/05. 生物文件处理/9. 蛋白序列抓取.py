#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import sys
import glob
import subprocess
import re
from datetime import datetime

'''
特定基因组序列精准提取与蛋白翻译工具 (单行序列极简版)
功能：输入多个基因名称，基于 GFF 注释精准提取序列。
特性：
1. 交互式选择输出类型：全长核苷酸序列 (DNA) 和/或 CDS翻译的蛋白序列 (Protein)。
2. 自动拼接打断的 CDS 区域（去除内含子）。
3. 自动处理链的正负性，负链自动反向互补。
4. 序列严格按照 "两行格式" (Single-line FASTA) 输出，方便下游正则和脚本处理。
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
    return seq.translate(trans)[::-1]

def translate_dna_to_protein(seq):
    """标准遗传密码翻译引擎"""
    codon_table = {
        'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
        'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
        'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
        'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
        'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
        'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
        'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
        'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
        'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
        'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
        'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
        'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
        'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
        'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
        'TAC':'Y', 'TAT':'Y', 'TAA':'*', 'TAG':'*',
        'TGC':'C', 'TGT':'C', 'TGA':'*', 'TGG':'W',
    }
    seq = seq.upper()
    protein = []
    
    for i in range(0, len(seq) - len(seq) % 3, 3):
        codon = seq[i:i+3]
        protein.append(codon_table.get(codon, 'X'))
        
    return "".join(protein)

def get_target_genes():
    """获取用户想要提取的基因及输出模式"""
    print("\n🎯 --- 输入目标基因 ---")
    print("请输入你想提取的基因名称，多个基因请用逗号 (,) 分隔。")
    genes_input = input("👉 请输入: ").strip()
    if not genes_input:
        print("⚠️ 未输入任何基因，程序退出。")
        sys.exit(0)
        
    targets = {clean_feature_id(g.strip()) for g in genes_input.split(',') if g.strip()}
    
    print("\n📦 --- 选择输出模式 ---")
    print("  [1] 仅输出基因全长核苷酸序列 (DNA, 包含内含子)")
    print("  [2] 仅输出 CDS 翻译的蛋白序列 (Protein, 剔除内含子)")
    print("  [3] 两者都要")
    
    mode = input("👉 请选择 (1/2/3) [默认 3]: ").strip()
    if mode not in ['1', '2', '3']: mode = '3'
    
    return targets, mode

def parse_gff_for_targets(gff_file, target_genes):
    """扫描注释文件，抓取基因全长范围与所有 CDS 片段"""
    print(f"\n⚙️ 正在解析注释文件寻找目标... ")
    
    found_genes = {g: {'chrom': None, 'strand': None, 'gene_coord': None, 'cds_coords': []} for g in target_genes}
    
    with open(gff_file, 'r') as f:
        for line in f:
            if line.startswith('#') or not line.strip(): continue
            parts = line.strip().split('\t')
            if len(parts) < 9: continue
            
            chrom, feature, start, end, strand, attr = parts[0], parts[2], parts[3], parts[4], parts[6], parts[8]
            
            attr_dict = {}
            for item in attr.strip().split(';'):
                if '=' in item:
                    k, v = item.split('=', 1)
                    attr_dict[k.strip()] = v.strip(' "')
                elif ' ' in item:
                    k, v = item.split(' ', 1)
                    attr_dict[k.strip()] = v.strip(' "')

            if feature in ['gene', 'mRNA', 'transcript']:
                gene_id = attr_dict.get('ID') or attr_dict.get('Name') or attr_dict.get('locus_tag') or "Unknown"
                clean_id = clean_feature_id(gene_id)
                
                if clean_id in found_genes:
                    found_genes[clean_id]['chrom'] = chrom
                    found_genes[clean_id]['strand'] = strand
                    found_genes[clean_id]['gene_coord'] = (int(start), int(end))
            
            elif feature == 'CDS':
                parent = attr_dict.get('Parent')
                raw_id = parent.split(',')[0] if parent else attr_dict.get('ID')
                clean_id = clean_feature_id(raw_id)
                
                if clean_id in found_genes:
                    found_genes[clean_id]['cds_coords'].append((int(start), int(end)))
                    if not found_genes[clean_id]['chrom']:
                        found_genes[clean_id]['chrom'] = chrom
                        found_genes[clean_id]['strand'] = strand
                        
    valid_genes = {k: v for k, v in found_genes.items() if v['chrom']}
    
    print(f"✅ 共输入 {len(target_genes)} 个目标基因，成功在 GFF 中定位到 {len(valid_genes)} 个。")
    if len(target_genes) > len(valid_genes):
        missing = target_genes - set(valid_genes.keys())
        print(f"⚠️ 以下基因未在注释文件中找到: {', '.join(missing)}")
        
    return valid_genes

def extract_sequence_chunk(fasta_file, chrom, start, end):
    """底层提取单个坐标块的序列"""
    cmd = f"samtools faidx {fasta_file} {chrom}:{start}-{end}"
    result = subprocess.run(cmd, shell=True, stdout=subprocess.PIPE, text=True, stderr=subprocess.DEVNULL)
    lines = result.stdout.strip().split('\n')
    if len(lines) > 1:
        return "".join(lines[1:])
    return ""

def format_fasta(header, seq):
    """【核心修改】格式化 FASTA 文本，强制单行序列输出"""
    return f"{header}\n{seq}\n"

def process_and_export(fasta_file, found_genes, mode, base_dir):
    """基于模式提取、拼接、反向互补与翻译，并导出文件"""
    timestamp = datetime.now().strftime("%Y%m%d_%H%M")
    
    out_dna = os.path.join(base_dir, f"09_Gene_DNA_{timestamp}.fasta")
    out_prot = os.path.join(base_dir, f"09_CDS_Protein_{timestamp}.fasta")
    
    f_dna = open(out_dna, 'w') if mode in ['1', '3'] else None
    f_prot = open(out_prot, 'w') if mode in ['2', '3'] else None

    print(f"\n🧬 正在从参考基因组中提取并进行生物学转换...")

    for gene_id, info in found_genes.items():
        chrom = info['chrom']
        strand = info['strand']
        
        if f_dna:
            if info['gene_coord']:
                g_start, g_end = info['gene_coord']
            else:
                g_start = min([c[0] for c in info['cds_coords']])
                g_end = max([c[1] for c in info['cds_coords']])
                
            raw_dna = extract_sequence_chunk(fasta_file, chrom, g_start, g_end)
            if raw_dna:
                final_dna = reverse_complement(raw_dna) if strand == '-' else raw_dna
                dir_note = "minus_strand_RC" if strand == '-' else "plus_strand"
                header = f">{gene_id}_DNA {chrom}:{g_start}-{g_end} strand:{strand} info:{dir_note}"
                f_dna.write(format_fasta(header, final_dna))

        if f_prot:
            cds_list = info['cds_coords']
            if not cds_list:
                print(f"⚠️ {gene_id} 没有对应的 CDS 注释，已跳过蛋白提取。")
                continue
                
            cds_list = sorted(cds_list, key=lambda x: x[0])
            
            spliced_dna = ""
            for c_start, c_end in cds_list:
                spliced_dna += extract_sequence_chunk(fasta_file, chrom, c_start, c_end)
                
            if strand == '-':
                spliced_dna = reverse_complement(spliced_dna)
                
            protein_seq = translate_dna_to_protein(spliced_dna)
            
            if protein_seq.endswith('*'):
                protein_seq = protein_seq[:-1]
                
            exon_count = len(cds_list)
            header = f">{gene_id}_Protein {chrom} strand:{strand} exons:{exon_count} length:{len(protein_seq)}aa"
            f_prot.write(format_fasta(header, protein_seq))

    if f_dna: f_dna.close()
    if f_prot: f_prot.close()

    print("\n🎉 提取与翻译全部完成！")
    if mode in ['1', '3']: print(f"   📄 基因核苷酸文件: {out_dna}")
    if mode in ['2', '3']: print(f"   🥩 蛋白序列文件: {out_prot}")

def run_pipeline():
    base_dir = get_base_dir()
    
    fasta_list = glob.glob(os.path.join(base_dir, "**", "*.fasta"), recursive=True) + \
                 glob.glob(os.path.join(base_dir, "**", "*.fa"), recursive=True)
    selected_fasta = interactive_select(fasta_list, "参考基因组 (.fasta/.fa)", base_dir)
    if not selected_fasta: return
    ref_fasta = selected_fasta[0]
    
    if not os.path.exists(ref_fasta + ".fai"):
        print("\n⚙️ 正在建立参考基因组索引...")
        subprocess.run(f"samtools faidx {ref_fasta}", shell=True)

    anno_list = glob.glob(os.path.join(base_dir, "**", "*.gff*"), recursive=True) + \
                glob.glob(os.path.join(base_dir, "**", "*.gtf"), recursive=True)
    selected_anno = interactive_select(anno_list, "基因注释文件 (.gff3 / .gtf)", base_dir)
    if not selected_anno: return
    gff_file = selected_anno[0]

    target_genes, mode = get_target_genes()

    found_genes = parse_gff_for_targets(gff_file, target_genes)
    if not found_genes:
        print("🛑 没有找到任何有效基因记录，程序终止。")
        return

    process_and_export(ref_fasta, found_genes, mode, base_dir)

if __name__ == "__main__":
    try:
        check_env()
        run_pipeline()
    except KeyboardInterrupt:
        print("\n🛑 用户强制退出。")