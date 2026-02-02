#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import sys
import glob
import subprocess
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# ================= 0. 核心配置 =================
OUTPUT_DIR_NAME = "15_PCA"
PLINK_EXEC = "plink"  # 请确保 plink 已安装并在环境变量中
PCA_COMPONENTS = 10   # 计算多少个 PC
COVAR_PCS = 3         # 选多少个 PC 放入 GEMMA 协变量文件

# ================= 1. 交互逻辑 (复用) =================
def get_base_dir():
    return os.path.dirname(os.path.abspath(__file__))

def find_files(ext, path=None):
    if not ext.startswith('.'): ext = '.' + ext
    if path is None: path = get_base_dir()
    search_pattern = os.path.join(path, '**', f'*{ext}')
    # 排除输出目录
    files = sorted(glob.glob(search_pattern, recursive=True))
    return [f for f in files if OUTPUT_DIR_NAME not in f]

def choose_file(files, desc="文件"):
    if not files:
        print(f"[提示] 未找到任何 {desc}")
        return None

    print(f"\n2. 找到 {len(files)} 个 {desc}:")
    for i, f in enumerate(files, 1):
        print(f"  [{i}] {os.path.basename(f)}  ({os.path.relpath(f, get_base_dir())})")

    while True:
        try:
            u_input = input(f"\n3. 请选择 {desc} 编号: ").strip()
            idx = int(u_input) - 1
            if 0 <= idx < len(files): return files[idx]
        except: pass

def make_sure(action):
    return input(f"\n4. 确认{action}? (y/n): ").strip().lower() in ['y', 'yes']

# ================= 2. 核心处理逻辑 =================
def run_pca_pipeline(fam_path):
    base_dir = get_base_dir()
    out_dir = os.path.join(base_dir, OUTPUT_DIR_NAME)
    if not os.path.exists(out_dir): os.makedirs(out_dir)

    # 1. 准备文件前缀
    bfile_prefix = os.path.splitext(fam_path)[0] # 去掉 .fam
    file_name = os.path.basename(bfile_prefix)
    out_prefix = os.path.join(out_dir, f"{file_name}_pca")

    # 2. 调用 PLINK 计算 PCA
    print(f"\n--- [1/3] 正在调用 PLINK 计算 PCA (Top {PCA_COMPONENTS}) ---")
    cmd = [
        PLINK_EXEC, 
        "--bfile", bfile_prefix,
        "--pca", str(PCA_COMPONENTS),
        "--out", out_prefix,
        "--allow-extra-chr"
    ]
    
    try:
        subprocess.run(cmd, check=True)
        print(f"  [成功] PCA 计算完成。")
    except FileNotFoundError:
        print(f"  [错误] 找不到 '{PLINK_EXEC}' 命令。请确保 PLINK 已安装。")
        return
    except subprocess.CalledProcessError as e:
        print(f"  [失败] PLINK 运行出错: {e}")
        return

    eigenvec_file = out_prefix + ".eigenvec"
    eigenval_file = out_prefix + ".eigenval"

    if not os.path.exists(eigenvec_file):
        print("  [错误] 未生成 .eigenvec 文件，流程终止。")
        return

    # 3. 绘制 PCA 散点图
    print(f"\n--- [2/3] 正在绘制 PCA 散点图 (PC1 vs PC2) ---")
    try:
        # PLINK 的 eigenvec 文件没有表头，前两列是 FID, IID，后面是 PC1, PC2...
        # 也就是: FID, IID, PC1, PC2, ...
        col_names = ["FID", "IID"] + [f"PC{i}" for i in range(1, PCA_COMPONENTS + 1)]
        df = pd.read_csv(eigenvec_file, sep='\s+', header=None, names=col_names)
        
        plt.figure(figsize=(10, 8))
        sns.scatterplot(x="PC1", y="PC2", data=df, alpha=0.7, edgecolor=None)
        
        plt.title(f"PCA Plot: {file_name}", fontsize=15)
        plt.xlabel(f"PC1", fontsize=12)
        plt.ylabel(f"PC2", fontsize=12)
        plt.grid(True, linestyle='--', alpha=0.3)
        
        plot_out = out_prefix + "_plot.png"
        plt.savefig(plot_out, dpi=300)
        plt.close()
        print(f"  [成功] 散点图已保存: {os.path.basename(plot_out)}")
        
    except Exception as e:
        print(f"  [绘图失败] {e}")

    # 4. 生成 GEMMA 协变量文件
    print(f"\n--- [3/3] 正在生成 GEMMA 协变量文件 (Top {COVAR_PCS} PCs) ---")
    try:
        # GEMMA 要求协变量文件:
        # 1. 纯数字矩阵，无表头
        # 2. 第一列通常是 1 (Intercept)
        # 3. 顺序必须与 .fam 文件中的样本顺序严格一致 (PLINK 输出通常是一致的，但为了保险我们这里直接用 eigenvec 的行)
        
        # 提取前 N 个 PC
        # df 已经在上面读取了
        
        # 创建 Intercept 列 (全为 1)
        covar_df = pd.DataFrame()
        covar_df['Intercept'] = [1] * len(df)
        
        # 加入 PC 列
        for i in range(1, COVAR_PCS + 1):
            covar_df[f'PC{i}'] = df[f'PC{i}']
            
        # 保存为 .txt (Tab 分隔，无索引，无表头)
        covar_out = os.path.join(out_dir, f"{file_name}_gemma_covar.txt")
        covar_df.to_csv(covar_out, sep='\t', index=False, header=False)
        
        print(f"  [成功] 协变量文件已生成: {os.path.basename(covar_out)}")
        print(f"      -> 包含列: Intercept + PC1 ~ PC{COVAR_PCS}")
        print(f"      -> 可在 GEMMA 中使用参数: -c {os.path.basename(covar_out)}")
        
    except Exception as e:
        print(f"  [生成失败] {e}")

# ================= 主函数 =================
def main():
    print("==============================================")
    print("   Step 15: 主成分分析 (PCA) & 协变量准备")
    print("==============================================")
    
    # 1. 选择 PLINK 文件 (.fam)
    # 因为 PCA 是基于基因型的，所以要选 .fam (对应 .bed/.bim)
    fam_files = find_files('.fam')
    target_fam = choose_file(fam_files, "PLINK 基因型文件 (.fam)")
    
    if not target_fam: return

    # 2. 执行
    if make_sure("开始计算 PCA"):
        run_pca_pipeline(target_fam)
        
    print(f"\n" + "="*50)
    print(f"任务完成！结果存放在: {OUTPUT_DIR_NAME}")
    print("提示: 生成的 _gemma_covar.txt 可用于后续 LMM 分析以校正群体结构。")
    print("="*50)

if __name__ == "__main__":
    main()