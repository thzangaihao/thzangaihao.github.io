#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import sys
import glob
import subprocess
from datetime import datetime
import concurrent.futures

'''
功能基因组区域变异筛选工具
功能：基于 GFF3/GTF 注释文件，从群体 VCF 文件中筛选出：
      1. 基因区变异 (Gene Region Variants)
      2. 编码区变异 (CDS Region Variants)
依赖：bedtools
'''

def get_base_dir():
    return os.path.dirname(os.path.abspath(__file__))

def check_env():
    if subprocess.call(['which', 'bedtools'], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL) != 0:
        print("❌ 错误: 未找到 bedtools！请执行: conda install -c bioconda bedtools")
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

def parse_annotation_to_bed(anno_file, out_dir):
    """
    解析 GFF3/GTF 文件，生成 gene 和 CDS 的坐标边界文件 (.bed)
    """
    print(f"\n⚙️ 正在解析基因注释文件: {os.path.basename(anno_file)}...")
    genes_bed = os.path.join(out_dir, "temp_genes.bed")
    cds_bed = os.path.join(out_dir, "temp_cds.bed")
    
    gene_count, cds_count = 0, 0

    with open(anno_file, 'r') as f, open(genes_bed, 'w') as g_out, open(cds_bed, 'w') as c_out:
        for line in f:
            if line.startswith('#') or not line.strip(): continue
            parts = line.strip().split('\t')
            if len(parts) < 8: continue
            
            chrom = parts[0]
            feature = parts[2]
            # GFF/GTF 是 1-based, BED 是 0-based，因此 start 需要减 1
            start = int(parts[3]) - 1 
            end = parts[4]
            
            bed_line = f"{chrom}\t{start}\t{end}\n"
            
            # 兼容各种注释文件的命名习惯
            if feature in ['gene', 'mRNA', 'transcript']:
                g_out.write(bed_line)
                gene_count += 1
            elif feature == 'CDS':
                c_out.write(bed_line)
                cds_count += 1

    # 对生成的 BED 文件进行排序，这是 bedtools 稳定运行的前提
    for bed in [genes_bed, cds_bed]:
        subprocess.run(f"bedtools sort -i {bed} > {bed}.sorted && mv {bed}.sorted {bed}", shell=True, check=True)

    print(f"✅ 解析完毕！提取到 {gene_count} 个基因/转录本特征，{cds_count} 个 CDS 特征。")
    return genes_bed, cds_bed

def process_single_vcf(vcf_file, genes_bed, cds_bed, out_root):
    """单样本变异拆分核心逻辑"""
    sample_name = os.path.basename(vcf_file).replace('.vcf', '').replace('_SNP_filtered', '')
    sample_dir = os.path.join(out_root, sample_name)
    os.makedirs(sample_dir, exist_ok=True)
    
    print(f"🚀 正在提取样本: {sample_name}")

    out_gene_vcf = os.path.join(sample_dir, f"{sample_name}_Gene.vcf")
    out_cds_vcf = os.path.join(sample_dir, f"{sample_name}_CDS.vcf")

    # 使用 bedtools intersect 提取交集
    # -a 指定输入的 VCF
    # -b 指定目标的 BED 区间
    # -header 保留 VCF 的表头信息
    # -u (unique) 只要有重叠就输出，且即使匹配到多个目标（如重叠的 CDS）也只输出一次原记录
    
    cmd_gene = f"bedtools intersect -a {vcf_file} -b {genes_bed} -header -u > {out_gene_vcf}"
    cmd_cds  = f"bedtools intersect -a {vcf_file} -b {cds_bed} -header -u > {out_cds_vcf}"

    subprocess.run(cmd_gene, shell=True, check=True)
    subprocess.run(cmd_cds, shell=True, check=True)
    
    return sample_name

def run_pipeline():
    base_dir = get_base_dir()
    
    # 1. 选择注释文件 (GFF3 / GTF)
    anno_list = glob.glob(os.path.join(base_dir, "**", "*.gff3"), recursive=True) + \
                glob.glob(os.path.join(base_dir, "**", "*.gtf"), recursive=True) + \
                glob.glob(os.path.join(base_dir, "**", "*.gff"), recursive=True)
    selected_anno = interactive_select(anno_list, "基因注释文件 (.gff3 / .gtf)", base_dir)
    if not selected_anno: return
    anno_file = selected_anno[0]

    # 2. 选择待处理的已过滤 VCF 文件
    vcf_list = glob.glob(os.path.join(base_dir, "**", "*_SNP_filtered.vcf"), recursive=True)
    selected_vcfs = interactive_select(vcf_list, "待提取的 VCF 文件 (*_SNP_filtered.vcf)", base_dir)
    if not selected_vcfs: return

    # 3. 多线程配置
    print("\n💻 --- 运算资源配置 ---")
    try:
        max_workers = int(input(f"👉 同时处理几个 VCF 文件？(默认 {min(4, len(selected_vcfs))}): ") or str(min(4, len(selected_vcfs))))
    except:
        max_workers = 2

    # 4. 创建输出目录
    timestamp = datetime.now().strftime("%Y%m%d_%H%M")
    out_root = os.path.join(base_dir, f"07_Region_Filtered_{timestamp}")
    os.makedirs(out_root, exist_ok=True)

    # 5. 解析注释文件提取坐标 BED
    genes_bed, cds_bed = parse_annotation_to_bed(anno_file, out_root)

    print(f"\n🔥 启动并发引擎: 开始执行 VCF 区域投影...")
    
    with concurrent.futures.ThreadPoolExecutor(max_workers=max_workers) as executor:
        futures = [executor.submit(process_single_vcf, vcf, genes_bed, cds_bed, out_root) for vcf in selected_vcfs]
        
        for future in concurrent.futures.as_completed(futures):
            try:
                future.result() 
            except Exception as exc:
                print(f"❌ 处理时发生错误: {exc}")

    # 清理临时的 BED 文件保持目录整洁
    if os.path.exists(genes_bed): os.remove(genes_bed)
    if os.path.exists(cds_bed): os.remove(cds_bed)

    print(f"\n🎉 区域提取全部完成！\n   每个样本已生成对应的 _Gene.vcf 和 _CDS.vcf\n   结果存放在: {out_root}")

if __name__ == "__main__":
    try:
        check_env()
        run_pipeline()
    except KeyboardInterrupt:
        print("\n🛑 用户强制退出。")