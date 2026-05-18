#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import sys
import glob
from datetime import datetime

'''
🌌 JCVI ➡️ SynVisio 格式完美转换神器 (V3 防崩溃终极版)
新增功能：
1. [防御升级] 引入 GFF 基因白名单机制，自动剔除 "幽灵基因"，彻底解决解析崩溃报错！
2. [格式对齐] 严格还原 MCScanX 原始空格与 Tab 缩进，迎合 SynVisio 的强迫症解析器。
3. [纯净输出] 移除所有非官方支持的注释行。
'''

def get_base_dir():
    return os.path.dirname(os.path.abspath(__file__))

def interactive_select(files, desc):
    if not files:
        print(f"⚠️ 未扫描到任何 {desc}！请确保文件在当前目录或子目录中。")
        sys.exit(0)
    
    print(f"\n📂 扫描到以下 {len(files)} 个 {desc}:")
    for i, f in enumerate(files, 1):
        print(f"  [{i}] {os.path.relpath(f, get_base_dir())}")
        
    while True:
        choice = input(f"\n👉 请选择 {desc} (输入编号，或 q 退出): ").strip()
        if choice.lower() == 'q': sys.exit(0)
        try:
            return files[int(choice)-1]
        except Exception:
            print("⚠️ 输入无效，请重新输入。")

def parse_bed_and_build_mapping(bed_a, bed_b):
    gene_to_new_chr = {}
    old_to_new_chr = {}
    mapping_records = [] 
    valid_genes = set() # 记录所有合法的基因 ID

    def process_bed(bed_file, prefix_char, species_id):
        chr_set = set()
        genes = []
        with open(bed_file, 'r', encoding='utf-8') as f:
            for line in f:
                if not line.strip() or line.startswith('#'): continue
                parts = line.strip().split()
                if len(parts) >= 4:
                    chrom, gene = parts[0], parts[3]
                    chr_set.add(chrom)
                    genes.append((gene, chrom))
                    valid_genes.add(gene) # 录入白名单
        
        sorted_chrs = sorted(list(chr_set))
        for idx, old_chr in enumerate(sorted_chrs, 1):
            new_chr = f"{prefix_char}{idx}"
            old_to_new_chr[old_chr] = new_chr
            mapping_records.append((species_id, old_chr, new_chr))
            
        for gene, old_chr in genes:
            gene_to_new_chr[gene] = old_to_new_chr[old_chr]

    process_bed(bed_a, "A", "Species_A")
    process_bed(bed_b, "B", "Species_B")
    
    return gene_to_new_chr, old_to_new_chr, mapping_records, valid_genes

def convert_bed_to_gff_mapped(bed_files, out_gff, old_to_new_chr):
    with open(out_gff, 'w', encoding='utf-8') as out:
        for bed in bed_files:
            with open(bed, 'r', encoding='utf-8') as f:
                for line in f:
                    line = line.strip()
                    if not line or line.startswith('#'): continue
                    parts = line.split()
                    if len(parts) >= 4:
                        old_chr = parts[0]
                        new_chr = old_to_new_chr.get(old_chr, old_chr)
                        out.write(f"{new_chr}\t{parts[3]}\t{parts[1]}\t{parts[2]}\n")

def write_alignment_block(out_file, b_idx, block, gene_to_new_chr):
    if not block: return
    g1, g2 = block[0]
    chr1 = gene_to_new_chr.get(g1, "A1")
    chr2 = gene_to_new_chr.get(g2, "B1")
    
    out_file.write(f"## Alignment {b_idx}: score=100.0 e_value=0 N={len(block)} {chr1}&{chr2} plus\n")
    for i, (gene1, gene2) in enumerate(block):
        # 严格模仿 MCScanX 的缩进格式: 两个空格开头，破折号后带空格，Tab 分隔
        out_file.write(f"  {b_idx}-{i:2}:\t{gene1}\t{gene2}\t0\n")

def convert_anchors_to_collinearity_mapped(anchors_in, collin_out, gene_to_new_chr, valid_genes):
    with open(anchors_in, 'r', encoding='utf-8') as f, open(collin_out, 'w', encoding='utf-8') as out:
        # 移除了所有自定义注释行，保持纯净
        block_idx = 0
        current_block = []
        dropped_pairs = 0

        for line in f:
            line = line.strip()
            if not line: continue
            
            if line.startswith('#'):
                if current_block:
                    write_alignment_block(out, block_idx, current_block, gene_to_new_chr)
                    block_idx += 1
                    current_block = []
                continue
            
            parts = line.split()
            if len(parts) >= 2:
                g1, g2 = parts[0], parts[1]
                # 🔥 核心防御：只有当两个基因都在 GFF 中存在时，才加入区块
                if g1 in valid_genes and g2 in valid_genes:
                    current_block.append((g1, g2))
                else:
                    dropped_pairs += 1
        
        if current_block:
            write_alignment_block(out, block_idx, current_block, gene_to_new_chr)
            
        if dropped_pairs > 0:
            print(f"⚠️ 警告: 拦截了 {dropped_pairs} 对'幽灵基因'（在 BED 中不存在）。已自动将其剔除以防 SynVisio 崩溃。")

def export_mapping_table(records, out_file):
    with open(out_file, 'w', encoding='utf-8') as f:
        f.write("Species\tOriginal_Name\tSynVisio_Name\n")
        for rec in records:
            f.write(f"{rec[0]}\t{rec[1]}\t{rec[2]}\n")

def main():
    print("=" * 60)
    print(" 🌌 JCVI ➡️ SynVisio (V3 防崩溃终极版)")
    print("=" * 60)

    base_dir = get_base_dir()
    bed_files = glob.glob(os.path.join(base_dir, "**", "*.bed"), recursive=True)
    anchors_files = glob.glob(os.path.join(base_dir, "**", "*.anchors"), recursive=True)

    if len(bed_files) < 2 or not anchors_files:
        print("❌ 错误：找不到足够的 .bed 或 .anchors 文件！")
        sys.exit(1)

    bed_a = interactive_select(bed_files, "物种 A 的 .bed 文件")
    bed_b = interactive_select([f for f in bed_files if f != bed_a], "物种 B 的 .bed 文件")
    anchors_file = interactive_select(anchors_files, "共线性文件 (.anchors)")

    prefix_a = os.path.basename(bed_a).replace('.bed', '')
    prefix_b = os.path.basename(bed_b).replace('.bed', '')

    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    work_dir = os.path.join(base_dir, f"SynVisio_Ready_{prefix_a}_{prefix_b}_{timestamp}")
    os.makedirs(work_dir, exist_ok=True)

    out_gff = os.path.join(work_dir, "synvisio.gff")
    out_collinearity = os.path.join(work_dir, f"{prefix_a}_{prefix_b}.collinearity")
    out_mapping = os.path.join(work_dir, "chromosome_mapping.tsv")

    try:
        print("\n🛠️ 正在扫描 BED 构建基因白名单...")
        gene_to_new_chr, old_to_new_chr, mapping_records, valid_genes = parse_bed_and_build_mapping(bed_a, bed_b)
        export_mapping_table(mapping_records, out_mapping)

        print("🛠️ 正在整合 GFF 文件...")
        convert_bed_to_gff_mapped([bed_a, bed_b], out_gff, old_to_new_chr)
        
        print("🛠️ 正在清洗并重构 Collinearity 文件 (防崩溃处理中)...")
        convert_anchors_to_collinearity_mapped(anchors_file, out_collinearity, gene_to_new_chr, valid_genes)
        
        print("\n" + "="*60)
        print("🎉 转换完成！幽灵基因已清除，格式已严格对齐！")
        print(f"📂 输出目录: \033[92m{os.path.relpath(work_dir, base_dir)}/\033[0m")
        print("👉 现在请刷新 SynVisio 页面，重新上传这两个文件。")
        print("="*60 + "\n")
    except Exception as e:
        print(f"\n❌ 转换过程中发生错误: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main()