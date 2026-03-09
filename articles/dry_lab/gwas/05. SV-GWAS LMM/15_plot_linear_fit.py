#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import sys
import glob
import subprocess
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats

# ================= 0. 核心配置 =================
OUTPUT_DIR_NAME = "15_Boxplots_From_Table"
BCFTOOLS_EXEC = "bcftools" 

# ================= 1. 基础交互逻辑 =================
def get_base_dir():
    return os.path.dirname(os.path.abspath(__file__))

def find_files(exts, path=None):
    if path is None: path = get_base_dir()
    if isinstance(exts, str): exts = [exts]
    all_files = []
    for ext in exts:
        if not ext.startswith('.'): ext = '.' + ext
        search_pattern = os.path.join(path, '**', f'*{ext}')
        all_files.extend(glob.glob(search_pattern, recursive=True))
    return sorted(list(set([f for f in all_files if OUTPUT_DIR_NAME not in f and "output" not in f])))

def choose_file_unlimited(files, desc="文件"):
    if not files:
        print(f"[提示] 未找到任何 {desc}")
        return None

    print(f"\n--- 找到 {len(files)} 个 {desc}: ---")
    for i, f in enumerate(files, 1):
        print(f"  [{i}] {os.path.basename(f)}")

    while True:
        try:
            u_input = input(f"\n>>> 请输入 {desc} 编号 (单选): ").strip()
            idx = int(u_input) - 1
            if 0 <= idx < len(files): return files[idx]
        except: pass

# ================= 2. 表格解析与选择 =================
def select_variants_from_tsv(tsv_path):
    print(f"\n正在读取表格: {os.path.basename(tsv_path)} ...")
    try:
        df = pd.read_csv(tsv_path, sep='\t')
        if 'Target_ID' not in df.columns:
            print("[错误] 未找到 'Target_ID' 列。")
            return []
        
        ids = df['Target_ID'].astype(str).tolist()
        types = df['SV_Type'].astype(str).tolist() if 'SV_Type' in df.columns else ['NA'] * len(df)
        counts = df['Carriers_Count'].astype(str).tolist() if 'Carriers_Count' in df.columns else ['?'] * len(df)
        
        print(f"\n--- 表中包含 {len(df)} 个变异位点 ---")
        print(f"{'编号':<6} | {'变异ID (预览)':<20} | {'类型':<6} | {'携带样本数'}")
        print("-" * 60)
        
        for i, (vid, vtype, vcount) in enumerate(zip(ids, types, counts), 1):
            # ID 预览截断，防止刷屏
            short_id = (vid[:18] + '..') if len(vid) > 18 else vid
            print(f"[{i:<4}] | {short_id:<20} | {vtype:<6} | {vcount}")
            
        print("-" * 60)
        
        user_input = input(">>> 请选择欲绘图的变异编号 (如 1,3,5-8): ").strip().lower()
        if not user_input: return []
        
        selected_indices = set()
        if user_input in ['all', 'a']:
            selected_indices = set(range(len(ids)))
        else:
            parts = user_input.split(',')
            for part in parts:
                part = part.strip()
                if '-' in part:
                    s, e = map(int, part.split('-'))
                    selected_indices.update(range(s-1, e))
                else:
                    selected_indices.add(int(part)-1)
        
        final_ids = [ids[i] for i in sorted(selected_indices) if 0 <= i < len(ids)]
        print(f"-> 已选中 {len(final_ids)} 个变异。")
        return final_ids

    except Exception as e:
        print(f"[读取失败] {e}")
        return []

# ================= 3. 数据提取与绘图 =================
def extract_genotype_data(vcf_path, variant_id):
    # 如果 ID 是一长串 "ID1;ID2;...", 我们只用第一个去 query 基因型
    # 因为它们在 VCF 里是同一行
    query_id = variant_id
    if ';' in variant_id:
        # 使用精确匹配查找整个字符串比较困难，不如直接根据行号提取？
        # 这里为了稳妥，我们尝试只用第一个 ID 去查，看看能不能查到
        # 或者使用 -i 'ID="ID1;ID2"' 的精确匹配
        pass
    
    # 尝试精确匹配整个长 ID
    cmd = [
        BCFTOOLS_EXEC, "query",
        "-i", f'ID="{variant_id}"',
        "-f", '[%SAMPLE\t%GT\n]',
        vcf_path
    ]
    
    # 如果失败，再尝试拆分策略 (此处简化，直接跑)
    try:
        res = subprocess.run(cmd, capture_output=True, text=True, check=True)
        # 如果结果为空，尝试用第一个 ID 查 (针对 Delly 合并 ID 的情况)
        if not res.stdout.strip() and ';' in variant_id:
             first_id = variant_id.split(';')[0]
             # 模糊匹配第一个 ID
             cmd = [BCFTOOLS_EXEC, "query", "-i", f'ID~"{first_id}"', "-f", '[%SAMPLE\t%GT\n]', vcf_path]
             res = subprocess.run(cmd, capture_output=True, text=True, check=True)

        lines = res.stdout.strip().split('\n')
        data = []
        for line in lines:
            if not line: continue
            parts = line.split('\t')
            if len(parts) < 2: continue
            sample, gt = parts[0], parts[1]
            if '.' in gt: dosage = np.nan
            else: dosage = gt.count('1')
            data.append({'SampleID': sample, 'Genotype': dosage})
        return pd.DataFrame(data)
    except: return pd.DataFrame()

def load_phenotype_data(trait_path, fam_path):
    try:
        fam_df = pd.read_csv(fam_path, sep=r'\s+', header=None, usecols=[1], names=['SampleID'])
        pheno_df = pd.read_csv(trait_path, sep=r'\s+', header=None, names=['Phenotype'])
        if len(fam_df) != len(pheno_df): return None
        return pd.concat([fam_df, pheno_df], axis=1)
    except: return None

def plot_boxplot(merged_df, variant_id, trait_name):
    df = merged_df.dropna(subset=['Genotype', 'Phenotype'])
    if len(df) < 5:
        print(f"  [跳过] 有效样本不足")
        return

    slope, intercept, r_value, p_value, std_err = stats.linregress(df['Genotype'], df['Phenotype'])

    plt.figure(figsize=(8, 6))
    sns.set_style("ticks")
    
    # [修复] hue=x 参数，解决 FutureWarning
    sns.boxplot(x='Genotype', y='Phenotype', data=df, 
                hue='Genotype', legend=False,  # 显式指定 hue 并关闭 legend
                width=0.5, palette="Pastel1", showfliers=False)
    
    sns.stripplot(x='Genotype', y='Phenotype', data=df, 
                  color='#2c3e50', size=4, alpha=0.6, jitter=0.2)
    
    x_range = np.array([-0.2, 2.2])
    y_pred = intercept + slope * x_range
    plt.plot(x_range, y_pred, color='#e74c3c', lw=2, ls='--',
             label=f'$R^2$={r_value**2:.3f}, $P$={p_value:.2e}')

    # [修复] 标题截断：防止 ID 太长导致 Raster Overflow
    display_title = variant_id
    if len(display_title) > 40:
        display_title = display_title[:40] + "..."
        
    plt.title(f"{display_title}\nvs {trait_name}", fontsize=14)
    plt.xlabel("Genotype (0:Ref, 1:Het, 2:Alt)", fontsize=12)
    plt.ylabel("Phenotype", fontsize=12)
    plt.legend(frameon=False)
    sns.despine()

    base_dir = get_base_dir()
    out_dir = os.path.join(base_dir, OUTPUT_DIR_NAME)
    if not os.path.exists(out_dir): os.makedirs(out_dir)
    
    # 文件名也要截断一下，或者只用哈希/前缀
    safe_fname_id = variant_id.replace(':', '_').replace(';', '_')
    if len(safe_fname_id) > 50:
        safe_fname_id = safe_fname_id[:50] + "_etc"
        
    out_name = f"{safe_fname_id}_{trait_name}_boxplot.png"
    out_path = os.path.join(out_dir, out_name)
    plt.savefig(out_path, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"  [成功] 绘图完成: {out_name}")

# ================= 主函数 =================
def main():
    print("==============================================")
    print("   Step 15: 交互式变异效应可视化")
    print("==============================================")
    
    tsv_files = find_files('.tsv')
    target_tsvs = [f for f in tsv_files if "Details" in f]
    target_tsv = choose_file_unlimited(target_tsvs, "溯源表格")
    if not target_tsv: return

    selected_ids = select_variants_from_tsv(target_tsv)
    if not selected_ids: return

    print("\n>>> 选择 VCF 文件")
    vcf_files = sorted(find_files('.vcf.gz'), key=lambda x: "Merged" not in x)
    target_vcf = choose_file_unlimited(vcf_files, "VCF")
    if not target_vcf: return

    print("\n>>> 选择表型文件 (.txt)")
    trait_files = [f for f in find_files('.txt') if "covar" not in f and "kinship" not in f]
    target_trait = choose_file_unlimited(trait_files, "表型")
    if not target_trait: return

    print("\n>>> 选择 FAM 文件")
    target_fam = choose_file_unlimited(find_files('.fam'), "FAM")
    if not target_fam: return

    trait_name = os.path.basename(target_trait).replace('.txt', '')
    print(f"\n--- 开始批量绘图 (共 {len(selected_ids)} 个) ---")
    
    pheno_df = load_phenotype_data(target_trait, target_fam)
    if pheno_df is None: return

    for i, vid in enumerate(selected_ids, 1):
        # 打印时也截断一下 ID，防止刷屏太乱
        short_id = (vid[:30] + '..') if len(vid) > 30 else vid
        print(f"[{i}/{len(selected_ids)}] 处理: {short_id} ...")
        
        geno_df = extract_genotype_data(target_vcf, vid)
        if not geno_df.empty:
            merged = pd.merge(geno_df, pheno_df, on='SampleID')
            plot_boxplot(merged, vid, trait_name)
            
    print(f"\n所有图片已保存在: {OUTPUT_DIR_NAME}")

if __name__ == "__main__":
    main()