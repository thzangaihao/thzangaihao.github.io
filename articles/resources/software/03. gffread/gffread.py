#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import sys
import glob
import subprocess
from datetime import datetime

'''
蛋白序列精准提取工具 (适配 InterProScan 输入)
功能：
1. 交互式选择参考基因组 (FASTA)、注释文件 (GFF) 和目标基因列表 (TXT)。
2. 调用 gffread 提取全量蛋白序列。
3. 根据目标列表过滤并生成最终的 IPRSCAN 输入文件。
'''

def get_base_dir():
    return os.path.dirname(os.path.abspath(__file__))

def check_gffread():
    """检查系统是否安装了 gffread"""
    if subprocess.call(['which', 'gffread'], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL) != 0:
        print("❌ 错误：未找到 gffread！请执行: conda install -c bioconda gffread")
        sys.exit(1)

def interactive_select(files, desc):
    """通用的交互式单选逻辑"""
    if not files:
        print(f"⚠️ 未找到任何 {desc}！")
        return None
    
    print(f"\n📂 扫描到以下 {len(files)} 个 {desc}:")
    for i, f in enumerate(files, 1):
        print(f"  [{i}] {os.path.relpath(f, get_base_dir())}")
        
    while True:
        choice = input(f"\n👉 请选择一个 {desc} (输入编号，q退出): ").strip().lower()
        if choice == 'q': sys.exit(0)
        try:
            return files[int(choice)-1]
        except:
            print("⚠️ 输入无效，请重新输入。")

def parse_target_list(txt_file):
    """读取目标基因列表文件"""
    target_ids = set()
    with open(txt_file, 'r', encoding='utf-8') as f:
        for line in f:
            line = line.strip()
            # 兼容带有表头的 TSV 汇总表，提取第一列
            if line and not line.startswith('#') and not line.startswith('Gene_ID'):
                gene_id = line.split('\t')[0]
                target_ids.add(gene_id)
    return target_ids

def run_extraction():
    base_dir = get_base_dir()
    check_gffread()
    
    print("\n" + "="*50)
    print(" 🧬 目标蛋白序列精准提取系统")
    print("="*50)

    # 1. 交互式选择三个核心文件
    ref_list = glob.glob(os.path.join(base_dir, "**", "*.fasta"), recursive=True) + \
               glob.glob(os.path.join(base_dir, "**", "*.fna"), recursive=True)
    ref_fasta = interactive_select(ref_list, "参考基因组 (.fasta)")

    gff_list = glob.glob(os.path.join(base_dir, "**", "*.gff*"), recursive=True)
    gff_file = interactive_select(gff_list, "基因注释文件 (.gff3)")

    txt_list = glob.glob(os.path.join(base_dir, "**", "*.txt"), recursive=True) + \
               glob.glob(os.path.join(base_dir, "**", "*.tsv"), recursive=True)
    target_txt = interactive_select(txt_list, "目标基因列表 (.txt / .tsv)")

    # 2. 创建输出目录
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    out_dir = os.path.join(base_dir, f"_protein_seqence_{timestamp}")
    os.makedirs(out_dir, exist_ok=True)
    
    temp_all_faa = os.path.join(out_dir, "temp_all_proteins.faa")
    final_faa = os.path.join(out_dir, "Target_Proteins_for_IPRSCAN.faa")

    # 3. 第一步：调用 gffread 提取全量蛋白
    print(f"\n⚙️ 正在调用 gffread 提取全量蛋白序列...")
    # -y 提取蛋白序列
    cmd = f"gffread '{gff_file}' -g '{ref_fasta}' -y '{temp_all_faa}'"
    subprocess.run(cmd, shell=True, check=True)

    # 4. 第二步：根据列表进行二次过滤 (解决 gffread ID 可能带后缀的问题)
    print(f"🔍 正在匹配目标列表...")
    target_ids = parse_target_list(target_txt)
    
    count = 0
    with open(temp_all_faa, 'r') as fin, open(final_faa, 'w') as fout:
        write_flag = False
        for line in fin:
            if line.startswith('>'):
                # 获取 FASTA 标题行的全文，例如 ">gene-TGAM01_v200818 [gene=TGAM...]"
                header_full = line[1:].strip()
                # 提取第一个词 (ID)，去掉常见前缀
                full_id = header_full.split()[0]
                clean_id = full_id.replace('gene-', '').replace('rna-', '').replace('cds-', '').split('.')[0]
                
                # 只要目标 ID 是 clean_id 的一部分，或者完全相等，就提取
                match_found = False
                for tid in target_ids:
                    if tid == clean_id or tid in full_id:
                        match_found = True
                        break
                
                if match_found:
                    write_flag = True
                    fout.write(line)
                    count += 1
                else:
                    write_flag = False
            else:
                if write_flag:
                    fout.write(line)
    # 清理临时文件
    if os.path.exists(temp_all_faa):
        os.remove(temp_all_faa)

    print("\n" + "="*50)
    print(f"🎉 提取完成！共匹配到 {count} 条序列。")
    print(f"📂 最终蛋白序列文件: {final_faa}")
    print("💡 下一步：请将此文件投喂给 InterProScan 进行本地功能注释。")

if __name__ == "__main__":
    try:
        run_extraction()
    except KeyboardInterrupt:
        print("\n🛑 用户强制退出。")