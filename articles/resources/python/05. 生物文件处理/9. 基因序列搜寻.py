#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import os
import sys
import glob
import pandas as pd
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

# ============= 基础路径获取 =============
def get_base_dir():
    """获取脚本所在目录，作为搜索的起点"""
    if getattr(sys, 'frozen', False):
        return os.path.dirname(sys.executable)
    return os.path.dirname(os.path.abspath(__file__))

def find_files(ext, path=None):
    if not ext.startswith('.'):
        ext = '.' + ext
    if path is None:
        path = get_base_dir()
    search_pattern = os.path.join(path, '**', f'*{ext}')
    return sorted(glob.glob(search_pattern, recursive=True))

# ============= 交互式单一文件选择 =============
def choose_single_file(ext_list, desc):
    """专门为核心输入设计的单选逻辑"""
    files = []
    for ext in ext_list:
        files.extend(find_files(ext))

    if not files:
        print(f"\n提示：在当前目录及其子目录下未找到任何 {desc} ({'/'.join(ext_list)})")
        return None

    print(f"\n>>> 寻找 {desc}:")
    for i, f in enumerate(files, 1):
        print(f"  [{i}] {os.path.basename(f)}")

    while True:
        try:
            user_input = input(f"请输入欲使用的 {desc} 编号: ").strip()
            idx = int(user_input) - 1
            if 0 <= idx < len(files):
                selected_file = files[idx]
                print(f"-> 已成功选择: {os.path.basename(selected_file)}\n")
                return selected_file
            else:
                print("未匹配到有效编号，请重新输入。")
        except ValueError:
            print("输入错误：请确保输入的是数字编号。")
        except KeyboardInterrupt:
            print("\n用户取消操作。")
            sys.exit()

# ============= 核心逻辑：解析 CSV =============
def parse_gene_list(csv_file):
    """解析CSV文件获取目标基因列表"""
    try:
        # 假设第一列就是基因ID，将整列转为字符串并去除可能的前后空格
        df = pd.read_csv(csv_file, header=None)
        gene_list = df[0].astype(str).str.strip().tolist()
        return set(gene_list)
    except Exception as e:
        print(f"读取基因列表 CSV 文件失败: {e}")
        sys.exit(1)

# ============= 核心逻辑：解析 GFF/GTF =============
def parse_annotation(anno_file, target_genes):
    """解析注释文件，返回 {gene_id: (chrom, start, end, strand)}"""
    print(f"正在解析注释文件 (这可能需要几秒钟) ...")
    gene_coords = {}
    
    with open(anno_file, 'r', encoding='utf-8') as f:
        for line in f:
            if line.startswith('#'):
                continue
            parts = line.strip().split('\t')
            if len(parts) < 9:
                continue

            feature_type = parts[2]
            # 锁定特征级别
            if feature_type not in ['gene', 'mRNA', 'transcript']:
                continue

            chrom, start, end, strand, attributes = parts[0], int(parts[3]), int(parts[4]), parts[6], parts[8]

            gene_id = None
            # 兼容 GFF3 (ID=xxx;) 和 GTF (gene_id "xxx";)
            if 'ID=' in attributes:
                for attr in attributes.split(';'):
                    if attr.startswith('ID='):
                        gene_id = attr.split('=')[1]
                        break
            elif 'gene_id' in attributes:
                for attr in attributes.split(';'):
                    attr = attr.strip()
                    if attr.startswith('gene_id'):
                        gene_id = attr.split(' ')[1].replace('"', '')
                        break

            # 如果 ID 存在于我们的 CSV 待搜寻列表中
            if gene_id and gene_id in target_genes:
                # 若同基因有多条转录本，取最长的一条
                if gene_id not in gene_coords or (end - start) > (gene_coords[gene_id][2] - gene_coords[gene_id][1]):
                    gene_coords[gene_id] = (chrom, start, end, strand)

    print(f"在注释文件中成功定位到 {len(gene_coords)} 个目标基因的坐标。")
    return gene_coords

# ============= 核心逻辑：提取序列 =============
def extract_sequences(fasta_file, gene_coords, output_fasta):
    """从 FASTA 提取序列，处理正负链与坐标转换"""
    print(f"正在读取基因组 FASTA 序列并提取... ")
    
    try:
        genome_dict = SeqIO.to_dict(SeqIO.parse(fasta_file, "fasta"))
    except Exception as e:
        print(f"读取 FASTA 文件失败: {e}")
        sys.exit(1)

    records = []
    found_count = 0
    
    for gene_id, (chrom, start, end, strand) in gene_coords.items():
        if chrom in genome_dict:
            # GFF/GTF 为 1-based, Python 索引为 0-based。
            # 起始位置需要 -1，终止位置保持不变 (Python 切片左闭右开)
            raw_seq = genome_dict[chrom].seq[start-1 : end]
            
            # 若基因在负链，提取反向互补序列
            if strand == '-':
                final_seq = raw_seq.reverse_complement()
            else:
                final_seq = raw_seq

            # 构建新的序列记录
            record = SeqRecord(
                final_seq,
                id=gene_id,
                description=f"loc={chrom}:{start}-{end} strand={strand}"
            )
            records.append(record)
            found_count += 1
        else:
            print(f"警告：染色体/Scaffold '{chrom}' 未在 FASTA 中找到 (涉及基因: {gene_id})")

    SeqIO.write(records, output_fasta, "fasta")
    print("-" * 50)
    print(f"提取大功告成！共提取 {found_count} 条序列。")
    print(f"文件已保存至: {output_fasta}")

# ============= 主流程 =============
def main():
    print("=" * 50)
    print(" 基因序列批量提取工具 (Interactive Gene Extractor)")
    print("=" * 50)

    # 1. 选择参考基因组
    fasta_file = choose_single_file(['.fasta', '.fa', '.fna'], "基因组文件 (FASTA)")
    if not fasta_file: return

    # 2. 选择注释文件
    anno_file = choose_single_file(['.gff3', '.gff', '.gtf'], "注释文件 (GFF/GTF)")
    if not anno_file: return

    # 3. 选择基因列表
    csv_file = choose_single_file(['.csv', '.txt'], "待搜寻基因列表 (CSV/TXT)")
    if not csv_file: return

    # 4. 读入 CSV 中的目标基因
    print("\n--- 开始计算 ---")
    target_genes = parse_gene_list(csv_file)
    print(f"从 CSV 文件中共读取到 {len(target_genes)} 个去重后的待搜寻基因 ID。")

    # 5. 解析注释坐标
    gene_coords = parse_annotation(anno_file, target_genes)
    if not gene_coords:
        print("未在注释文件中找到任何匹配的基因坐标，程序结束。请检查基因 ID 格式是否一致。")
        return

    # 6. 确认并提取
    response = input("\n最后确认：是否执行序列提取操作? (y/n): ").strip().lower()
    if response in ['y', 'yes']:
        out_name = input("请输入输出的 FASTA 文件名 (回车默认: extracted_genes.fasta): ").strip()
        if not out_name:
            out_name = "extracted_genes.fasta"
        elif not out_name.endswith('.fasta') and not out_name.endswith('.fa'):
            out_name += ".fasta"
            
        out_path = os.path.join(get_base_dir(), out_name)
        extract_sequences(fasta_file, gene_coords, out_path)
    else:
        print("操作已取消。")

if __name__ == "__main__":
    main()