#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import sys
import glob
import subprocess
from datetime import datetime

'''
🧬 多功能蛋白序列提取系统 (专业版)
功能：
1. 交互式选择：参考基因组、GFF3注释、目标基因列表。
2. 交互式模式：自主选择【仅最长转录本】或【全量转录本】。
3. 智能清洗：自动移除非法字符（* 和 .），完美修复 FASTA 标题行 Bug。
4. 双重过滤：支持全基因组模式（ALL）或精准目标列表过滤。
'''

def get_base_dir():
    return os.path.dirname(os.path.abspath(__file__))

def check_gffread():
    """检查系统是否安装了 gffread"""
    if subprocess.call(['which', 'gffread'], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL) != 0:
        print("❌ 错误：未找到 gffread！请执行: conda install -c bioconda gffread")
        sys.exit(1)

def interactive_select(files, desc, allow_all=False):
    """通用的交互式单选逻辑"""
    if not files and not allow_all:
        print(f"⚠️ 未找到任何 {desc}！")
        return None
    
    print(f"\n📂 扫描到以下 {len(files)} 个 {desc}:")
    for i, f in enumerate(files, 1):
        print(f"  [{i}] {os.path.relpath(f, get_base_dir())}")
        
    if allow_all:
        print(f"  [all] 提取全基因组序列 (不使用目标列表过滤)")
        
    while True:
        prompt = f"\n👉 请选择一个 {desc} (输入编号"
        prompt += "，或输入 all" if allow_all else ""
        prompt += "，q退出): "
        choice = input(prompt).strip().lower()
        
        if choice == 'q': sys.exit(0)
        if allow_all and choice == 'all': return 'ALL'
        
        try:
            return files[int(choice)-1]
        except:
            print("⚠️ 输入无效，请重新输入。")

def parse_target_list(txt_file):
    """读取目标基因列表文件，并进行彻底的格式清洗"""
    target_ids = set()
    with open(txt_file, 'r', encoding='utf-8') as f:
        for line in f:
            line = line.strip().replace('\r', '') 
            if line and not line.startswith('#') and not line.startswith('Gene_ID'):
                gene_id = line.split('\t')[0].strip().replace('gene-', '').replace('rna-', '')
                target_ids.add(gene_id)
    return target_ids

def process_sequences(input_path, output_path, filter_longest=True):
    """读取原始 gffread 输出，进行清洗和(可选的)最长转录本筛选"""
    if filter_longest:
        print(f"⚙️ 正在筛选每个基因的【最长转录本】并清洗非法字符...")
    else:
        print(f"⚙️ 正在处理【全量转录本】并清洗非法字符...")
        
    sequences = []
    current_header = None
    current_seq = []
    
    # 1. 读取临时 FASTA 文件
    with open(input_path, 'r', encoding='utf-8') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith('>'):
                if current_header:
                    sequences.append((current_header, ''.join(current_seq)))
                current_header = line
                current_seq = []
            else:
                current_seq.append(line)
        if current_header:
            sequences.append((current_header, ''.join(current_seq)))
            
    # 2. 核心处理逻辑
    if filter_longest:
        # 最长转录本模式
        processed_data = {}
        for header, seq in sequences:
            header_content = header[1:].strip()
            parts = header_content.split()
            full_id = parts[0]
            
            # 识别基因 ID
            gene_id = None
            for part in parts:
                if part.startswith('gene='):
                    gene_id = part.split('=')[1]
                    break
            if not gene_id:
                gene_id = full_id.replace('gene-', '').replace('rna-', '').replace('cds-', '').split('.')[0]
                
            clean_seq = "".join(filter(str.isalpha, seq))
            seq_len = len(clean_seq)
            
            if gene_id not in processed_data:
                processed_data[gene_id] = (header, clean_seq)
            else:
                if seq_len > len(processed_data[gene_id][1]):
                    processed_data[gene_id] = (header, clean_seq)
                    
        # 写入结果
        with open(output_path, 'w', encoding='utf-8') as out:
            for gene_id, (header, seq) in processed_data.items():
                out.write(f"{header}\n")
                for i in range(0, len(seq), 80):
                    out.write(f"{seq[i:i+80]}\n")
    else:
        # 保留全量转录本模式（仅清洗非法字符）
        with open(output_path, 'w', encoding='utf-8') as out:
            for header, seq in sequences:
                clean_seq = "".join(filter(str.isalpha, seq))
                out.write(f"{header}\n")
                for i in range(0, len(clean_seq), 80):
                    out.write(f"{clean_seq[i:i+80]}\n")

def run_extraction():
    base_dir = get_base_dir()
    check_gffread()
    
    print("\n" + "="*60)
    print(" 🧬 多功能蛋白序列提取系统 (交互组合版)")
    print("="*60)

    # 1. 交互式选择文件
    ref_list = glob.glob(os.path.join(base_dir, "**", "*.fasta"), recursive=True) + \
               glob.glob(os.path.join(base_dir, "**", "*.fna"), recursive=True)
    ref_fasta = interactive_select(ref_list, "参考基因组 (.fasta)")

    gff_list = glob.glob(os.path.join(base_dir, "**", "*.gff*"), recursive=True)
    gff_file = interactive_select(gff_list, "基因注释文件 (.gff3)")

    txt_list = glob.glob(os.path.join(base_dir, "**", "*.txt"), recursive=True) + \
               glob.glob(os.path.join(base_dir, "**", "*.tsv"), recursive=True)
    target_txt = interactive_select(txt_list, "目标基因列表 (.txt / .tsv)", allow_all=True)

    # 2. 新增核心交互：选择转录本提取模式
    print("\n[?] 请选择转录本提取模式:")
    print("  [1] 仅提取【最长转录本】 (推荐用于: 染色体共线性分析、正选择进化分析)")
    print("  [2] 保留【全量转录本】 (推荐用于: 可变剪接分析、全面结构注释)")
    while True:
        mode_choice = input("👉 请选择模式 (输入 1 或 2): ").strip()
        if mode_choice == '1':
            filter_longest = True
            mode_str = "Longest"
            break
        elif mode_choice == '2':
            filter_longest = False
            mode_str = "AllTranscripts"
            break
        else:
            print("⚠️ 输入无效，请重新输入 1 或 2。")

    # 3. 创建输出目录
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    out_dir = os.path.join(base_dir, f"Protein_Output_{timestamp}")
    os.makedirs(out_dir, exist_ok=True)
    
    temp_gffread_faa = os.path.join(out_dir, "temp_gffread_raw.faa")
    temp_processed_faa = os.path.join(out_dir, "temp_processed.faa")

    # 4. 调用 gffread 提取原始序列
    print(f"\n⚙️ 正在调用 gffread 提取基因组基础序列...")
    cmd = f"gffread '{gff_file}' -g '{ref_fasta}' -y '{temp_gffread_faa}'"
    subprocess.run(cmd, shell=True, check=True)

    # 5. 根据用户选择处理序列（清洗 + 可选的过滤最长）
    process_sequences(temp_gffread_faa, temp_processed_faa, filter_longest=filter_longest)
    
    if os.path.exists(temp_gffread_faa):
        os.remove(temp_gffread_faa)

    # 6. 根据目标列表过滤或全量保留
    if target_txt == 'ALL':
        print(f"🔍 检测到 [全基因组模式]，正在跳过列表过滤...")
        final_faa = os.path.join(out_dir, f"Genome_{mode_str}_Proteins.faa")
        os.rename(temp_processed_faa, final_faa)
        
        count = sum(1 for line in open(final_faa, 'r', encoding='utf-8') if line.startswith('>'))
        
        print("\n" + "="*60)
        print(f"🎉 全基因组蛋白序列提取成功！")
        print(f"   🔹 提取模式: {mode_str}")
        print(f"   🔹 总序列数: {count} 条")
        print(f"📂 最终文件路径: {final_faa}")
        
    else:
        final_faa = os.path.join(out_dir, f"Target_{mode_str}_Proteins.faa")
        print(f"🔍 正在匹配目标列表...")
        target_ids = parse_target_list(target_txt)
        
        count = 0
        with open(temp_processed_faa, 'r', encoding='utf-8') as fin, open(final_faa, 'w', encoding='utf-8') as fout:
            write_flag = False
            for line in fin:
                if line.startswith('>'):
                    header_full = line[1:].strip()
                    full_id = header_full.split()[0]
                    clean_id = full_id.replace('gene-', '').replace('rna-', '').replace('cds-', '').split('.')[0]
                    
                    match_found = False
                    for tid in target_ids:
                        if tid == clean_id or tid in full_id:
                            match_found = True
                            break
                    
                    if match_found:
                        write_flag = True
                        fout.write(line) # 完美修复：直接写入包含 '>' 的原始标题行
                        count += 1
                    else:
                        write_flag = False
                else:
                    if write_flag:
                        fout.write(line)

        if os.path.exists(temp_processed_faa):
            os.remove(temp_processed_faa)

        print("\n" + "="*60)
        print(f"🎉 目标列表筛选提取完成！")
        print(f"   🔹 提取模式: {mode_str}")
        print(f"   🔹 成功匹配: {count} 条目标序列")
        print(f"📂 最终文件路径: {final_faa}")

if __name__ == "__main__":
    try:
        run_extraction()
    except KeyboardInterrupt:
        print("\n🛑 用户强制退出。")