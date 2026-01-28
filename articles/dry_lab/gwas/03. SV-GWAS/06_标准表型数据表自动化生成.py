#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import sys
import glob
import pandas as pd

# ==============================================================================
# 0. 配置参数
# ==============================================================================
OUTPUT_DIR_NAME = "06_GEMMA_Phenotypes"
MISSING_VALUE = "-9"  # GEMMA 默认缺失值

# ==============================================================================
# 1. Cite2 交互逻辑模块
# ==============================================================================
def get_base_dir():
    if getattr(sys, 'frozen', False):
        return os.path.dirname(sys.executable)
    return os.path.dirname(os.path.abspath(__file__))

def find_files(ext, path=None):
    if not ext.startswith('.'): ext = '.' + ext
    if path is None: path = get_base_dir()
    search_pattern = os.path.join(path, '**', f'*{ext}')
    all_files = sorted(glob.glob(search_pattern, recursive=True))
    # 排除输出目录
    return [f for f in all_files if OUTPUT_DIR_NAME not in f]

def choose_file(files, desc="文件"):
    if not files:
        print(f"提示：未找到任何 {desc}")
        return []

    print(f"\n找到 {len(files)} 个 {desc}:")
    for i, f in enumerate(files, 1):
        print(f"  [{i}] {os.path.basename(f)}  (路径: {os.path.relpath(f, get_base_dir())})")

    while True:
        try:
            user_input = input(f"\n请输入 {desc} 编号 (单选): ").strip()
            if not user_input: continue
            
            idx = int(user_input) - 1
            if 0 <= idx < len(files):
                return files[idx]
            print("编号无效，请重试。")
        except ValueError:
            print("请输入数字。")

def make_sure(action_name="执行操作"):
    response = input(f"\n确认{action_name}? (y/n): ").strip().lower()
    return response in ['y', 'yes']

# ==============================================================================
# 2. 核心处理模块
# ==============================================================================
def process_splitting(tsv_path, fam_path):
    print(f"\n--- 正在读取数据 ---")
    
    # 1. 读取 FAM 文件 (获取标准顺序)
    try:
        # fam 文件前两列是 FID, IID
        fam_df = pd.read_csv(fam_path, sep='\s+', header=None, usecols=[1], names=['IID'], dtype=str)
        print(f"FAM 文件样本数: {len(fam_df)}")
    except Exception as e:
        print(f"[错误] 读取 FAM 文件失败: {e}")
        return

    # 2. 读取 TSV 表型文件
    try:
        pheno_df = pd.read_csv(tsv_path, sep='\t', dtype=str)
        print(f"TSV 文件样本数: {len(pheno_df)}")
        print(f"包含列: {list(pheno_df.columns)}")
    except Exception as e:
        print(f"[错误] 读取 TSV 文件失败: {e}")
        return

    # 3. 确定 ID 列
    columns = list(pheno_df.columns)
    id_col = columns[0]
    
    print("\n请确认 TSV 文件中的【样本 ID】列:")
    for i, col in enumerate(columns, 1):
        print(f"  [{i}] {col}")
    
    try:
        idx = input(f"请输入列编号 (默认 1: {id_col}): ").strip()
        if idx:
            id_col = columns[int(idx) - 1]
    except:
        print(f"使用默认列: {id_col}")
    
    # 重命名 ID 列以便合并
    pheno_df = pheno_df.rename(columns={id_col: 'IID'})
    
    # 4. 确定性状列 (所有非 ID 列)
    trait_cols = [c for c in pheno_df.columns if c != 'IID']
    print(f"\n检测到 {len(trait_cols)} 个性状: {trait_cols}")
    
    if not make_sure("开始拆分"): return

    # 5. 对齐数据 (Left Join)
    # 这一步至关重要，确保输出的行数和顺序与 FAM 完全一致
    merged_df = pd.merge(fam_df, pheno_df, on='IID', how='left')
    
    # 6. 创建输出目录
    out_dir = os.path.join(get_base_dir(), OUTPUT_DIR_NAME)
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
    
    print(f"\n--- 开始输出 (保存至 {OUTPUT_DIR_NAME}) ---")
    
    count = 0
    for trait in trait_cols:
        # 提取单列
        series = merged_df[trait]
        
        # 填充缺失值
        series = series.fillna(MISSING_VALUE)
        
        # 替换可能的空字符串
        series = series.replace(r'^\s*$', MISSING_VALUE, regex=True)
        
        # 构造文件名 (使用列名)
        # 建议剔除特殊字符
        safe_name = "".join([c if c.isalnum() or c in ['_','-'] else '_' for c in trait])
        out_path = os.path.join(out_dir, f"{safe_name}.txt")
        
        # 保存：无表头，无索引，纯数值
        series.to_csv(out_path, index=False, header=False)
        print(f"  -> 生成: {safe_name}.txt")
        count += 1
        
    print(f"\n[成功] 已拆分 {count} 个性状文件。")
    print("提示：这些文件可以直接用于 GEMMA 的 -p 参数。")

# ==============================================================================
# 主函数
# ==============================================================================
def main():
    print("==============================================")
    print("   Step 6: 表型拆分与格式化 (GEMMA Ready)")
    print("==============================================")
    
    # 1. 选择 TSV (包含多个性状)
    print("\n>>> 选择 1: 包含多个性状的 TSV 文件")
    tsvs = find_files('.tsv')
    tsv_file = choose_file(tsvs, "表型总表")
    if not tsv_file: return

    # 2. 选择 FAM (定义样本顺序)
    print("\n>>> 选择 2: PLINK .fam 文件 (用于对齐顺序)")
    fams = find_files('.fam')
    fam_file = choose_file(fams, "FAM 文件")
    if not fam_file: return

    # 3. 执行
    process_splitting(tsv_file, fam_file)

if __name__ == "__main__":
    main()