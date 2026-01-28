#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import sys
import glob
import subprocess
import pandas as pd
import time

# ==============================================================================
# 1. Cite2 交互逻辑模块 (源自您提供的 cite_v2.py)
# ==============================================================================
def get_base_dir():
    """获取脚本所在目录，作为搜索的起点"""
    if getattr(sys, 'frozen', False):
        return os.path.dirname(sys.executable)
    return os.path.dirname(os.path.abspath(__file__))

def find_files(ext, path=None):
    """查找指定扩展名文件"""
    if not ext.startswith('.'):
        ext = '.' + ext
    if path is None:
        path = get_base_dir()
    search_pattern = os.path.join(path, '**', f'*{ext}')
    # 排除之前的输出结果，防止混淆
    all_files = sorted(glob.glob(search_pattern, recursive=True))
    return [f for f in all_files if "cleaned" not in f]

def choose_file(files, desc="文件"):
    """支持多种选取模式 (单选/多选/范围/全部)"""
    if not files:
        print(f"提示：在当前目录及其子目录下未找到任何 {desc}")
        return []

    print(f"\n2. 找到 {len(files)} 个 {desc}:")
    for i, f in enumerate(files, 1):
        print(f"  [{i}] {os.path.basename(f)}  (路径: {os.path.relpath(f, get_base_dir())})")

    while True:
        try:
            prompt = (
                f"\n3. 请输入欲操作的 {desc} 编号\n"
                f" (通常请输入单个编号，如 1): \n"
            )
            user_input = input(prompt).strip().lower()

            if not user_input: continue

            selected_indices = set()
            if user_input in ['all', 'a']:
                selected_indices = set(range(len(files)))
            else:
                parts = user_input.split(',')
                for part in parts:
                    part = part.strip()
                    if '-' in part:
                        start_str, end_str = part.split('-')
                        for idx in range(int(start_str) - 1, int(end_str)):
                            if 0 <= idx < len(files): selected_indices.add(idx)
                    else:
                        idx = int(part) - 1
                        if 0 <= idx < len(files): selected_indices.add(idx)
            
            if not selected_indices:
                print("未匹配到有效编号，请重新输入。")
                continue

            selected_paths = [files[i] for i in sorted(list(selected_indices))]
            print(f"\n4. 已选择: {os.path.basename(selected_paths[0])}")
            if len(selected_paths) > 1:
                print("注意：本脚本一次仅处理一个表型文件，将只使用第一个选择的文件。")
            
            return selected_paths

        except ValueError:
            print("输入错误：请输入数字编号。")

def make_sure(action_name="执行操作"):
    response = input(f"\n5. 确认{action_name}? (y/n): ").strip().lower()
    return response in ['y', 'yes']

# ==============================================================================
# 2. 核心业务逻辑
# ==============================================================================

def get_vcf_samples(vcf_path):
    """调用 bcftools 获取 VCF/BCF 中的样本 ID 列表"""
    print(f"正在从基因型文件中读取样本列表: {os.path.basename(vcf_path)} ...")
    try:
        # bcftools query -l <file>
        result = subprocess.run(
            ["bcftools", "query", "-l", vcf_path],
            capture_output=True, text=True, check=True
        )
        samples = result.stdout.strip().split('\n')
        # 过滤空行
        samples = [s for s in samples if s]
        print(f"-> 基因型文件中共有 {len(samples)} 个样本。")
        return set(samples)
    except FileNotFoundError:
        print("[错误] 未找到 bcftools，请确保已安装 (conda install bcftools)。")
        sys.exit(1)
    except subprocess.CalledProcessError as e:
        print(f"[错误] 读取 VCF 样本失败: {e}")
        sys.exit(1)

def process_phenotype(tsv_path, vcf_samples_set):
    """处理 TSV 表型文件并进行过滤"""
    print(f"\n正在读取表型文件: {os.path.basename(tsv_path)} ...")
    
    try:
        df = pd.read_csv(tsv_path, sep='\t')
    except Exception as e:
        print(f"[错误] 读取 TSV 文件失败: {e}")
        return

    # --- 交互式选择 ID 列 ---
    columns = list(df.columns)
    id_col = columns[0] # 默认第一列
    
    if len(columns) > 1:
        print("\n检测到表型文件包含以下列:")
        for i, col in enumerate(columns, 1):
            print(f"  [{i}] {col}")
        
        try:
            col_idx = input(f"\n请选择代表 '样本ID' 的列编号 (默认 1: {columns[0]}): ").strip()
            if col_idx:
                idx = int(col_idx) - 1
                if 0 <= idx < len(columns):
                    id_col = columns[idx]
        except ValueError:
            print("输入无效，将使用第一列作为 ID 列。")
    
    print(f"-> 使用 '{id_col}' 作为样本 ID 进行对齐。")
    
    # 确保 ID 列是字符串类型，防止数字/字符串匹配错误
    df[id_col] = df[id_col].astype(str)

    # --- 核心对齐逻辑 ---
    # 找出 TSV 中有，但 VCF 中没有的样本
    pheno_samples_set = set(df[id_col])
    
    # 交集 (保留的)
    keep_ids = pheno_samples_set.intersection(vcf_samples_set)
    # 差集 (剔除的)
    remove_ids = pheno_samples_set - vcf_samples_set
    
    # 过滤 DataFrame
    clean_df = df[df[id_col].isin(keep_ids)]
    
    # --- 打印报告 ---
    print("\n" + "="*50)
    print("【对齐报告】")
    print(f"表型文件原始样本数 : {len(df)}")
    print(f"基因型文件样本数   : {len(vcf_samples_set)}")
    print("-" * 50)
    print(f"对齐成功 (保留)    : {len(clean_df)}")
    print(f"对齐失败 (剔除)    : {len(remove_ids)}")
    print("="*50)
    
    if remove_ids:
        print("\n[以下样本因无 VCF 数据被剔除]:")
        for i, rid in enumerate(sorted(list(remove_ids)), 1):
            print(f"  {i}. {rid}")
            if i >= 20:
                print(f"  ... 以及其他 {len(remove_ids)-20} 个")
                break
    else:
        print("\n完美匹配！没有任何样本被剔除。")

    # --- 保存文件 ---
    if make_sure("保存清洗后的 TSV 文件"):
        base_name = os.path.basename(tsv_path)
        file_name_no_ext = os.path.splitext(base_name)[0]
        output_filename = os.path.join(get_base_dir(), f"{file_name_no_ext}_cleaned.tsv")
        
        # 保存为 TSV，保留表头，不保留索引
        clean_df.to_csv(output_filename, sep='\t', index=False)
        print(f"\n[成功] 文件已保存至: {output_filename}")

# ==============================================================================
# 主函数
# ==============================================================================
def main():
    print("==============================================")
    print("   Step 5: 表型文件对齐与清洗 (TSV 工具)")
    print("==============================================")
    
    # 1. 选择 VCF/BCF 文件 (基准)
    print("\n>>> 第一步: 选择基准基因型文件 (VCF/BCF)")
    vcf_files = find_files('.gz', get_base_dir()) + find_files('.bcf', get_base_dir())
    # 去重
    vcf_files = sorted(list(set(vcf_files)))
    
    selected_vcfs = choose_file(vcf_files, "基因型文件")
    if not selected_vcfs: return
    vcf_file = selected_vcfs[0]

    # 2. 获取 VCF 样本
    vcf_samples = get_vcf_samples(vcf_file)

    # 3. 选择 TSV 表型文件
    print("\n>>> 第二步: 选择待清洗的表型文件 (TSV)")
    tsv_files = find_files('.tsv')
    selected_tsvs = choose_file(tsv_files, "TSV 表型文件")
    if not selected_tsvs: return
    tsv_file = selected_tsvs[0]

    # 4. 执行清洗
    process_phenotype(tsv_file, vcf_samples)

if __name__ == "__main__":
    main()