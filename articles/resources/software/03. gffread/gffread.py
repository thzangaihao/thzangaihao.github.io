#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import sys
import glob
import subprocess
from datetime import datetime

'''
蛋白序列提取工具 (支持全基因组模式)
功能：
1. 调用 gffread 提取蛋白序列。
2. 支持精准提取（目标列表过滤）和全量提取（All模式，用于构建背景库）。
'''

def get_base_dir():
    return os.path.dirname(os.path.abspath(__file__))

def check_gffread():
    """检查系统是否安装了 gffread"""
    if subprocess.call(['which', 'gffread'], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL) != 0:
        print("❌ 错误：未找到 gffread！请执行: conda install -c bioconda gffread")
        sys.exit(1)

def interactive_select(files, desc, allow_all=False):
    """通用的交互式单选逻辑，增加 allow_all 参数"""
    if not files and not allow_all:
        print(f"⚠️ 未找到任何 {desc}！")
        return None
    
    print(f"\n📂 扫描到以下 {len(files)} 个 {desc}:")
    for i, f in enumerate(files, 1):
        print(f"  [{i}] {os.path.relpath(f, get_base_dir())}")
        
    if allow_all:
        print(f"  [all] 提取全基因组所有蛋白序列 (用于构建富集分析背景库)")
        
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

def run_extraction():
    base_dir = get_base_dir()
    check_gffread()
    
    print("\n" + "="*50)
    print(" 🧬 目标蛋白序列提取系统 (支持全量背景提取)")
    print("="*50)

    # 1. 交互式选择文件
    ref_list = glob.glob(os.path.join(base_dir, "**", "*.fasta"), recursive=True) + \
               glob.glob(os.path.join(base_dir, "**", "*.fna"), recursive=True)
    ref_fasta = interactive_select(ref_list, "参考基因组 (.fasta)")

    gff_list = glob.glob(os.path.join(base_dir, "**", "*.gff*"), recursive=True)
    gff_file = interactive_select(gff_list, "基因注释文件 (.gff3)")

    txt_list = glob.glob(os.path.join(base_dir, "**", "*.txt"), recursive=True) + \
               glob.glob(os.path.join(base_dir, "**", "*.tsv"), recursive=True)
    # 允许选择全量提取
    target_txt = interactive_select(txt_list, "目标基因列表 (.txt / .tsv)", allow_all=True)

    # 2. 创建输出目录
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    out_dir = os.path.join(base_dir, f"11_IPRSCAN_Input_{timestamp}")
    os.makedirs(out_dir, exist_ok=True)
    
    temp_all_faa = os.path.join(out_dir, "temp_all_proteins.faa")

    # 3. 调用 gffread 提取全量蛋白
    print(f"\n⚙️ 正在调用 gffread 提取全基因组蛋白序列...")
    cmd = f"gffread '{gff_file}' -g '{ref_fasta}' -y '{temp_all_faa}'"
    subprocess.run(cmd, shell=True, check=True)

    # 4. 根据用户选择决定是否进行过滤
    if target_txt == 'ALL':
        print(f"🔍 检测到 [全基因组模式]，正在跳过列表过滤...")
        final_faa = os.path.join(out_dir, "All_Proteins_for_IPRSCAN.faa")
        os.rename(temp_all_faa, final_faa) # 直接重命名保留
        clean_fasta(final_faa)
        
        # 统计序列总数
        count = sum(1 for line in open(final_faa) if line.startswith('>'))
        
        print("\n" + "="*50)
        print(f"🎉 全基因组提取完成！共获得 {count} 条蛋白质序列。")
        print(f"📂 最终蛋白序列文件: {final_faa}")
        print("💡 下一步：请将此文件投喂给 InterProScan。注意，这可能需要运行数十个小时！")
        
    else:
        final_faa = os.path.join(out_dir, "Target_Proteins_for_IPRSCAN.faa")
        print(f"🔍 正在匹配目标列表...")
        target_ids = parse_target_list(target_txt)
        
        count = 0
        with open(temp_all_faa, 'r') as fin, open(final_faa, 'w') as fout:
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
                        clean_seq = "".join(filter(str.isalpha, line.strip()))
                        fout.write(clean_seq + "\n")
                        count += 1
                    else:
                        write_flag = False
                else:
                    if write_flag:
                        fout.write(line)

        if os.path.exists(temp_all_faa):
            os.remove(temp_all_faa)

        print("\n" + "="*50)
        print(f"🎉 提取完成！共匹配到 {count} 条目标序列。")
        print(f"📂 最终蛋白序列文件: {final_faa}")

def clean_fasta(file_path):
    """移除 FASTA 序列中的非法字符（. 和 *）"""
    print(f"🧹 正在清洗序列文件中的非法字符...")
    temp_path = file_path + ".tmp"
    with open(file_path, 'r') as f_in, open(temp_path, 'w') as f_out:
        for line in f_in:
            if line.startswith(">"):
                f_out.write(line)
            else:
                # 只保留字母，移除点、星号、空格等
                clean_seq = "".join(filter(str.isalpha, line.strip()))
                f_out.write(clean_seq + "\n")
    os.replace(temp_path, file_path)

if __name__ == "__main__":
    try:
        run_extraction()
    except KeyboardInterrupt:
        print("\n🛑 用户强制退出。")