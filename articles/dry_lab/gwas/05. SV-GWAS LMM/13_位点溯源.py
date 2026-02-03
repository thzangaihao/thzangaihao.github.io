#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import sys
import glob
import subprocess
import pandas as pd
import datetime
import time

# ================= 0. 核心配置 =================
OUTPUT_DIR_NAME = "13_Variant_Details"
BCFTOOLS_EXEC = "bcftools"

# ================= 1. 交互逻辑 (支持多选) =================
def get_base_dir():
    return os.path.dirname(os.path.abspath(__file__))

def find_files(exts, path=None):
    if path is None: path = get_base_dir()
    all_files = []
    if isinstance(exts, str): exts = [exts]
    for ext in exts:
        if not ext.startswith('.'): ext = '.' + ext
        search_pattern = os.path.join(path, '**', f'*{ext}')
        all_files.extend(glob.glob(search_pattern, recursive=True))
    # 排除输出目录，防止循环读取自己生成的结果
    return sorted(list(set([f for f in all_files if OUTPUT_DIR_NAME not in f])))

def choose_files_multi(files, desc="文件"):
    """多选模式：支持 1,2 | 1-5 | all"""
    if not files:
        print(f"[提示] 未找到任何 {desc}")
        return []

    print(f"\n2. 找到 {len(files)} 个 {desc}:")
    limit = 15
    for i, f in enumerate(files[:limit], 1):
        print(f"  [{i}] {os.path.basename(f)}  ({os.path.relpath(f, get_base_dir())})")
    if len(files) > limit:
        print(f"  ... (共 {len(files)} 个)")

    while True:
        try:
            prompt = (
                f"\n3. 请输入欲处理的 {desc} 编号\n"
                f" (支持格式: 1,2 | 1-4 | all): \n"
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
                        s, e = map(int, part.split('-'))
                        selected_indices.update(range(s-1, e))
                    else:
                        selected_indices.add(int(part)-1)
            
            selected_files = [files[i] for i in sorted(selected_indices) if 0 <= i < len(files)]
            if selected_files:
                print(f"\n4. 已选择 {len(selected_files)} 个文件。")
                return selected_files
        except ValueError:
            print("输入错误，请重试。")

def choose_file_single(files, desc="文件"):
    """单选模式：用于选择 VCF"""
    if not files: return None
    print(f"\n>>> 请选择作为基准的 {desc} (单选):")
    for i, f in enumerate(files, 1):
        print(f"  [{i}] {os.path.basename(f)}")
    while True:
        try:
            u_input = input(f"编号: ").strip()
            idx = int(u_input) - 1
            if 0 <= idx < len(files): return files[idx]
        except: pass

# ================= 2. 核心处理逻辑 =================

def get_ids_from_file(file_path):
    """从显著结果文件中智能提取 ID"""
    try:
        # 兼容 csv/tsv/txt
        if file_path.endswith('.csv'):
            df = pd.read_csv(file_path)
        else:
            # 尝试 Tab 分隔，如果只有一列则兼容
            df = pd.read_csv(file_path, sep='\t')
            if len(df.columns) == 1 and ',' in df.iloc[0,0]:
                 df = pd.read_csv(file_path, sep=',')

        # 智能匹配 ID 列名
        possible_cols = ['rs', 'RS', 'snp', 'SNP', 'id', 'ID', 'Variant', 'Target_ID']
        target_col = None
        for col in df.columns:
            if col in possible_cols:
                target_col = col
                break
        
        if target_col:
            # 转字符串并去重
            return df[target_col].dropna().astype(str).unique().tolist()
        else:
            print(f"[跳过] {os.path.basename(file_path)}: 未找到 ID 列 (列名: {list(df.columns)})")
            return []
    except Exception as e:
        print(f"[错误] 读取 {os.path.basename(file_path)} 失败: {e}")
        return []

def trace_worker(vcf_path, ids, source_name):
    """具体的溯源工兵函数"""
    if not ids: return

    # 1. 构建输出文件名 (时间戳 + 原文件名)
    timestamp = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
    clean_name = os.path.splitext(source_name)[0].replace("_sig_sites", "").replace("_filtered", "")
    out_filename = f"{timestamp}_{clean_name}_Details.tsv"
    
    out_dir = os.path.join(get_base_dir(), OUTPUT_DIR_NAME)
    if not os.path.exists(out_dir): os.makedirs(out_dir)
    out_path = os.path.join(out_dir, out_filename)

    print(f"  -> 目标 ID 数: {len(ids)}")
    print(f"  -> 输出文件名: {out_filename}")

    # 2. 准备集合用于精确匹配
    target_set = set(ids)
    results = []

    # 3. 流式读取 VCF (高效)
    cmd = [
        BCFTOOLS_EXEC, "query", 
        "-f", '%CHROM\t%POS\t%INFO/END\t%ID\t%INFO/SVLEN\t%INFO/SVTYPE\t%REF\t%ALT\n',
        vcf_path
    ]
    
    try:
        process = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        
        # 进度条逻辑
        matched_count = 0
        
        for line in process.stdout:
            cols = line.strip().split('\t')
            if len(cols) < 8: continue
            
            # 精确匹配：VCF 里的这一行 ID 必须在我们的清单里
            # (针对 Merged VCF 里的 "ID1;ID2" 这种情况，必须完全一致才能匹配)
            # 如果您希望 "ID1;ID2" 中的 "ID1" 也能匹配，请将下方改为:
            # if not target_set.isdisjoint(cols[3].split(';')): 
            
            # 目前保持精确全字符匹配 (Robust Exact Match)
            if cols[3] in target_set:
                results.append({
                    "Target_ID": cols[3],
                    "Chr": cols[0],
                    "Start_Pos": cols[1],
                    "End_Pos": cols[2] if cols[2] != '.' else "NA",
                    "SV_Length": cols[4] if cols[4] != '.' else "NA",
                    "SV_Type": cols[5],
                    "Ref": cols[6],
                    "Alt": cols[7]
                })
                matched_count += 1
                
        process.wait()
        
        # 4. 保存结果
        if results:
            df = pd.DataFrame(results)
            # 优选列顺序
            cols_seq = ["Target_ID", "Chr", "Start_Pos", "End_Pos", "SV_Length", "SV_Type", "Ref", "Alt"]
            df = df[cols_seq]
            
            df.to_csv(out_path, sep='\t', index=False)
            print(f"  [成功] 提取到 {len(df)} 条记录")
        else:
            print(f"  [警告] 0 条匹配！请检查 VCF 是否包含这些 ID。")

    except Exception as e:
        print(f"  [异常] {e}")

# ================= 主函数 =================
def main():
    print("==============================================")
    print("   Step 13: 批量自动化变异溯源 (Batch Auto-Trace)")
    print("==============================================")
    
    # 1. 多选显著位点文件
    print("\n>>> 第一步: 选择显著位点结果文件 (支持多选/全选)")
    # 自动搜索所有可能的 txt/tsv/csv 结果文件
    sig_files = find_files(['.txt', '.tsv', '.csv'])
    # 简单过滤，优先展示看起来像结果的文件
    sig_candidates = [f for f in sig_files if "sig" in f or "assoc" in f or "08_" in f]
    
    selected_sig_files = choose_files_multi(sig_candidates, "显著位点文件")
    if not selected_sig_files: return

    # 2. 单选 VCF 文件
    print("\n>>> 第二步: 选择基准 VCF 文件 (用于查户口)")
    vcf_files = find_files(['.vcf.gz'])
    # 优先推荐 Merged 文件
    vcf_files = sorted(vcf_files, key=lambda x: "Merged_Population.vcf.gz" not in x)
    
    target_vcf = choose_file_single(vcf_files, "VCF 文件")
    if not target_vcf: return

    # 3. 批量循环处理
    print(f"\n" + "="*50)
    print(f"开始批量处理 {len(selected_sig_files)} 个任务...")
    print("="*50)

    for i, sig_file in enumerate(selected_sig_files, 1):
        fname = os.path.basename(sig_file)
        print(f"\n>>> 任务 [{i}/{len(selected_sig_files)}]: {fname}")
        
        # 提取 ID
        ids = get_ids_from_file(sig_file)
        if ids:
            # 执行溯源
            trace_worker(target_vcf, ids, fname)
        else:
            print("  [跳过] 空 ID 列表或读取失败")
            
    print("\n" + "="*50)
    print("所有任务执行完毕！结果已存入 13_Variant_Details 文件夹。")
    print("="*50)

if __name__ == "__main__":
    main()