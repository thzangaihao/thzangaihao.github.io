#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import sys
import glob
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# ==============================================================================
# 0. 核心配置
# ==============================================================================
OUTPUT_DIR_NAME = "09_Manhattan_Plots"
FIG_DPI = 300

# ==============================================================================
# 1. Cite2 交互逻辑 (升级为多选版)
# ==============================================================================
def get_base_dir():
    return os.path.dirname(os.path.abspath(__file__))

def find_files(exts, path=None):
    if path is None: path = get_base_dir()
    all_found = []
    for ext in exts:
        if not ext.startswith('.'): ext = '.' + ext
        search_pattern = os.path.join(path, '**', f'*{ext}')
        files = glob.glob(search_pattern, recursive=True)
        # 排除之前的输出目录，防止重复
        files = [f for f in files if OUTPUT_DIR_NAME not in f and "png" not in f]
        all_found.extend(files)
    return sorted(list(set(all_found)))

def choose_files(files, desc="文件"):
    if not files:
        print(f"[提示] 未找到任何 {desc}。")
        return []

    print(f"\n2. 找到 {len(files)} 个 {desc}:")
    limit = 15
    if len(files) > limit:
        for i, f in enumerate(files[:limit], 1):
            print(f"  [{i}] {os.path.basename(f)}")
        print(f"  ... (共 {len(files)} 个文件) ...")
    else:
        for i, f in enumerate(files, 1):
            print(f"  [{i}] {os.path.basename(f)}")

    while True:
        try:
            prompt = (
                f"\n3. 请输入欲绘图的 {desc} 编号\n"
                f" (输入 'all' 全选，或格式: 1,2 | 1-4): \n"
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
            print("选择无效，请重试。")
        except ValueError:
            print("输入错误。")

def make_sure(action_name="执行操作"):
    response = input(f"\n5. 确认{action_name}? (y/n): ").strip().lower()
    return response in ['y', 'yes']

# ==============================================================================
# 2. 核心绘图逻辑
# ==============================================================================

def draw_qq(p_values, base_name, output_dir):
    """绘制 QQ 图"""
    print("   -> 正在绘制 QQ 图...")
    p_observed = -np.log10(np.sort(p_values))
    n = len(p_values)
    p_expected = -np.log10(np.arange(1, n + 1) / (n + 1))
    
    fig, ax = plt.subplots(figsize=(6, 6))
    ax.scatter(p_expected, p_observed, s=10, color='blue', alpha=0.5)
    
    max_val = max(np.max(p_expected), np.max(p_observed))
    if pd.isna(max_val) or np.isinf(max_val): max_val = 10
        
    ax.plot([0, max_val], [0, max_val], color='red', linestyle='--')
    
    ax.set_xlabel('Expected $-log_{10}(P)$')
    ax.set_ylabel('Observed $-log_{10}(P)$')
    ax.set_title(f'QQ Plot: {base_name}')
    
    output_path = os.path.join(output_dir, f"{base_name}.qq.png")
    plt.tight_layout()
    plt.savefig(output_path, dpi=FIG_DPI)
    plt.close(fig) # 关闭画布释放内存
    print(f"      已保存: {os.path.basename(output_path)}")

def draw_manhattan(file_path, output_dir):
    base_name = os.path.splitext(os.path.basename(file_path))[0]
    # 如果文件名是 .assoc.txt，去掉双重后缀
    if base_name.endswith('.assoc'): base_name = base_name[:-6]
    
    print(f"\n--- 正在处理: {base_name} ---")
    print(f"   [1/3] 读取数据...")
    
    try:
        df = pd.read_csv(file_path, sep='\t')
        if len(df.columns) < 3:
            df = pd.read_csv(file_path, delim_whitespace=True)
    except Exception as e:
        print(f"   [错误] 读取失败: {e}")
        return

    # 数据预处理 (防报错)
    df = df.loc[:, ~df.columns.duplicated()]
    
    required_cols = ['chr', 'ps', 'p_wald']
    for col in required_cols:
        if col not in df.columns:
            print(f"   [跳过] 文件缺少列 {col}，可能不是结果文件。")
            return

    # 清洗 P 值
    df = df.dropna(subset=['p_wald'])
    df = df[df['p_wald'] > 0] # 避免 log(0)
    df['minuslog10p'] = -np.log10(df['p_wald'])

    if df.empty:
        print("   [错误] 有效数据为空。")
        return

    # 清洗染色体 (兼容 Chr1 和 1)
    df['chr_raw'] = df['chr'].astype(str).str.extract(r'(\d+)')
    df['chr_no'] = pd.to_numeric(df['chr_raw'], errors='coerce')
    df = df.dropna(subset=['chr_no', 'ps'])
    df['chr_no'] = df['chr_no'].astype(int)
    
    # 排序
    df = df.sort_values(['chr_no', 'ps'])
    
    # 计算坐标
    print("   [2/3] 计算曼哈顿坐标...")
    chr_len = df.groupby('chr_no')['ps'].max()
    chr_offset = chr_len.cumsum().shift(1).fillna(0)
    offset_map = chr_offset.to_dict()
    df['plot_pos'] = df['ps'] + df['chr_no'].map(offset_map)
    
    max_pos = df['plot_pos'].max()
    if pd.isna(max_pos) or np.isinf(max_pos):
        print("   [错误] 坐标计算异常。")
        return

    # 绘图
    print("   [3/3] 绘制曼哈顿图...")
    fig, ax = plt.subplots(figsize=(14, 6))
    
    colors = ['#4A4A4A', '#87CEFA'] 
    x_labels = []
    x_labels_pos = []
    
    grouped = df.groupby('chr_no')
    for num, (name, group) in enumerate(grouped):
        ax.scatter(group['plot_pos'], group['minuslog10p'], 
                   color=colors[num % len(colors)], s=10, linewidth=0)
        x_labels.append(f"Chr{int(name)}")
        x_labels_pos.append(group['plot_pos'].mean())
        
    ax.set_xticks(x_labels_pos)
    ax.set_xticklabels(x_labels, fontsize=9, rotation=0)
    ax.set_xlim([0, max_pos])
    ax.set_ylim(bottom=0)
    ax.set_xlabel('Chromosome')
    ax.set_ylabel(r'$-log_{10}(P)$')
    
    # 阈值线 (Bonferroni)
    threshold = -np.log10(0.05 / len(df))
    ax.axhline(threshold, color='red', linestyle='--', linewidth=1)
    ax.text(0, threshold + 0.2, f'Threshold', color='red', fontsize=8)
    
    ax.set_title(f'Manhattan Plot: {base_name}', fontsize=14)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    
    output_png = os.path.join(output_dir, f"{base_name}.manhattan.png")
    plt.tight_layout()
    plt.savefig(output_png, dpi=FIG_DPI)
    plt.close(fig)
    print(f"      已保存: {os.path.basename(output_png)}")
    
    # 顺便画 QQ 图
    draw_qq(df['p_wald'], base_name, output_dir)

# ==============================================================================
# 主函数
# ==============================================================================
def main():
    print("==============================================")
    print("   Step 9: 批量曼哈顿图绘制工具")
    print("==============================================")
    
    # 1. 查找结果文件
    print("\n>>> 第一步: 搜索 GWAS 结果文件 (.assoc.txt)")
    files = find_files(['assoc.txt', 'assoc'])
    
    # 2. 交互选择
    selected_files = choose_files(files, "GWAS 结果文件")
    if not selected_files: return

    # 3. 准备目录
    base_dir = get_base_dir()
    output_dir = os.path.join(base_dir, OUTPUT_DIR_NAME)
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
        print(f"\n[系统] 创建图片保存目录: {output_dir}")

    # 4. 执行绘图
    if make_sure("开始批量绘图"):
        print(f"\n--- 开始处理 {len(selected_files)} 个文件 ---")
        for i, f in enumerate(selected_files, 1):
            print(f"\n>>> 任务 [{i}/{len(selected_files)}]")
            draw_manhattan(f, output_dir)
            
        print("\n" + "="*50)
        print("所有图片已生成完毕！")
        print(f"存放目录: {output_dir}")
        print("="*50)

if __name__ == "__main__":
    main()