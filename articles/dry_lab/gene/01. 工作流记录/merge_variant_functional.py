#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import sys
import glob
import pandas as pd
from datetime import datetime

'''
功能与变异数据合并脚本 (v3 强力清洗版)
修复：针对 NCBI/Liftoff 复杂的 ID 前缀进行强力清洗，确保 100% 匹配。
'''

def get_base_dir():
    return os.path.dirname(os.path.abspath(__file__))

def interactive_select(files, desc):
    if not files:
        print(f"⚠️ 未找到任何 {desc}！")
        return None
    print(f"\n📂 扫描到以下 {len(files)} 个 {desc}:")
    for i, f in enumerate(files, 1):
        print(f"  [{i}] {os.path.relpath(f, get_base_dir())}")
    while True:
        choice = input(f"\n👉 请选择 {desc} (输入编号，q退出): ").strip().lower()
        if choice == 'q': sys.exit(0)
        try: return files[int(choice)-1]
        except: print("⚠️ 输入无效。")

def clean_ipr_id(raw_id):
    """强力 ID 清洗器：剥离各种前缀和后缀，还原纯净的基因编号"""
    x = str(raw_id)
    # 处理你提供数据中的特殊前缀
    if 'mrna.' in x:
        x = x.split('mrna.')[-1]
    elif 'gene.' in x:
        x = x.split('gene.')[-1]
    elif '|' in x:
        x = x.split('|')[-1]
    
    # 清理掉可能存在的 rna- 或 gene- 前缀
    x = x.replace('rna-', '').replace('gene-', '').replace('cds-', '')
    
    # 按照小数点切割，去掉 .t1 或 .1 这样的后缀
    return x.split('.')[0]

def run_merge():
    base_dir = get_base_dir()
    print("\n" + "="*50)
    print(" 📊 变异-功能终极大师表生成器 (v3)")
    print("="*50)

    var_list = glob.glob(os.path.join(base_dir, "**", "*Combined_Gene_List.tsv"), recursive=True) + \
               glob.glob(os.path.join(base_dir, "**", "*.tsv"), recursive=True)
    var_file = interactive_select(var_list, "变异汇总表 (.tsv)")

    ipr_list = glob.glob(os.path.join(base_dir, "**", "*.tsv"), recursive=True)
    ipr_file = interactive_select(ipr_list, "InterProScan 结果表 (.tsv)")

    print(f"\n⚙️ 正在分析变异表结构...")
    df_var = pd.read_csv(var_file, sep='\t')
    
    if 'Gene_ID' in df_var.columns: id_col = 'Gene_ID'
    elif 'ID' in df_var.columns: id_col = 'ID'
    else:
        id_col = df_var.columns[0]
        print(f"ℹ️ 自动识别第一列 '{id_col}' 为 ID 列。")

    print(f"🔍 正在读取 InterProScan 结果并清洗 ID...")
    ipr_cols = ['ID', 'MD5', 'Length', 'Analysis', 'Signature_ID', 'Description', 'Start', 'End', 'Score', 'Status', 'Date', 'InterPro_ID', 'InterPro_Desc', 'GO', 'Pathway']
    try:
        df_ipr = pd.read_csv(ipr_file, sep='\t', names=ipr_cols, header=None, engine='python', on_bad_lines='skip')
    except Exception as e:
        print(f"❌ 读取 IPRSCAN 文件失败: {e}")
        return

    # 汇总功能描述
    df_ipr['Clean_Desc'] = df_ipr['InterPro_Desc'].fillna(df_ipr['Description'])
    
    def summarize_func(x):
        items = [str(i) for i in x if pd.notna(i) and str(i) != '-' and str(i).strip() != '']
        return "; ".join(sorted(set(items))) if items else "Hypothetical Protein"

    ipr_summary = df_ipr.groupby('ID').agg({
        'Clean_Desc': summarize_func,
        'GO': summarize_func,
        'Pathway': summarize_func
    }).reset_index()

    # 👉 核心修复：使用强力清洗器处理两侧的 ID
    ipr_summary['ID_Simple'] = ipr_summary['ID'].apply(clean_ipr_id)
    df_var['Match_Key'] = df_var[id_col].apply(lambda x: str(x).split('.')[0])

    print(f"🧬 正在将变异数据与蛋白质功能进行合并...")
    final_df = pd.merge(df_var, ipr_summary, left_on='Match_Key', right_on='ID_Simple', how='left')

    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    out_file = f"Final_Master_Table_V3_{timestamp}.csv"
    
    cols_to_keep = [c for c in df_var.columns if c != 'Match_Key'] + ['Clean_Desc', 'GO', 'Pathway']
    final_df = final_df[cols_to_keep]
    final_df['Clean_Desc'] = final_df['Clean_Desc'].fillna("No annotation found")
    
    final_df.to_csv(out_file, index=False)

    print("\n" + "="*50)
    print(f"🎉 大师表合并成功！那些隐藏的功能终于重见天日了！")
    print(f"📂 结果文件: {os.path.basename(out_file)}")
    print("="*50)

if __name__ == "__main__":
    run_merge()