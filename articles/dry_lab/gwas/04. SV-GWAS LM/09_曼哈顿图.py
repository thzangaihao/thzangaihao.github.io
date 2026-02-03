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
# 1. 辅助与交互逻辑
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

def get_threshold_config():
    """
    交互式获取阈值配置
    """
    print("\n>>> 设置阈值线 (Threshold Line)")
    print(" 1. 自动 Bonferroni 矫正 (0.05 / SNP数量) [默认]")
    print(" 2. 自定义 P-value (例如 1e-5)")
    print(" 3. 不显示阈值线")
    
    choice = input(" 请输入选项 (1/2/3): ").strip()
    
    if choice == '2':
        while True:
            try:
                val_str = input(" 请输入 P-value (如 1e-5): ").strip()
                val = float(val_str)
                if 0 < val < 1:
                    return {'type': 'fixed', 'val': val}
                print(" 数值必须在 0 到 1 之间。")
            except ValueError:
                print(" 输入格式错误，请输入浮点数。")
    elif choice == '3':
        return {'type': 'none'}
    else:
        # 默认为 Bonferroni
        return {'type': 'bonferroni', 'val': 0.05}

def make_sure(action_name="执行操作"):
    response = input(f"\n6. 确认{action_name}? (y/n): ").strip().lower()
    return response in ['y', 'yes']

# ==============================================================================
# 2. 核心绘图逻辑
# ==============================================================================

def draw_qq(p_values, base_name, output_dir):
    """绘制 QQ 图"""
    print("   -> 正在绘制 QQ 图...")
    # 移除无效值
    p_values = p_values.dropna()
    p_values = p_values[p_values > 0]
    
    p_observed = -np.log10(np.sort(p_values))
    n = len(p_values)
    p_expected = -np.log10(np.arange(1, n + 1) / (n + 1))
    
    fig, ax = plt.subplots(figsize=(6, 6))
    ax.scatter(p_expected, p_observed, s=10, color='blue', alpha=0.5)
    
    # 动态调整坐标轴范围
    max_val = max(np.max(p_expected), np.max(p_observed))
    if pd.isna(max_val) or np.isinf(max_val): max_val = 10
    limit = max_val * 1.05
        
    ax.plot([0, limit], [0, limit], color='red', linestyle='--')
    
    ax.set_xlabel('Expected $-log_{10}(P)$')
    ax.set_ylabel('Observed $-log_{10}(P)$')
    ax.set_title(f'QQ Plot: {base_name}')
    
    output_path = os.path.join(output_dir, f"{base_name}.qq.png")
    plt.tight_layout()
    plt.savefig(output_path, dpi=FIG_DPI)
    plt.close(fig)
    print(f"      已保存: {os.path.basename(output_path)}")

def draw_manhattan(file_path, output_dir, threshold_config):
    base_name = os.path.splitext(os.path.basename(file_path))[0]
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

    df = df.loc[:, ~df.columns.duplicated()]
    required_cols = ['chr', 'ps', 'p_wald']
    for col in required_cols:
        if col not in df.columns:
            print(f"   [跳过] 文件缺少列 {col}，可能不是结果文件。")
            return

    df = df.dropna(subset=['p_wald'])
    df = df[df['p_wald'] > 0]
    df['minuslog10p'] = -np.log10(df['p_wald'])

    if df.empty:
        print("   [错误] 有效数据为空。")
        return

    df['chr_raw'] = df['chr'].astype(str).str.extract(r'(\d+)')
    df['chr_no'] = pd.to_numeric(df['chr_raw'], errors='coerce')
    df = df.dropna(subset=['chr_no', 'ps'])
    df['chr_no'] = df['chr_no'].astype(int)
    
    df = df.sort_values(['chr_no', 'ps'])
    
    print("   [2/3] 计算曼哈顿坐标...")
    chr_len = df.groupby('chr_no')['ps'].max()
    chr_offset = chr_len.cumsum().shift(1).fillna(0)
    offset_map = chr_offset.to_dict()
    df['plot_pos'] = df['ps'] + df['chr_no'].map(offset_map)
    
    max_pos = df['plot_pos'].max()
    if pd.isna(max_pos) or np.isinf(max_pos):
        print("   [错误] 坐标计算异常。")
        return

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
    
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.set_title(f'Manhattan Plot: {base_name}', fontsize=14)

    # =======================================================
    # 新增逻辑：根据配置和实际数据决定是否画线
    # =======================================================
    threshold_val = None
    threshold_label = ""
    
    if threshold_config['type'] == 'bonferroni':
        # 计算动态 Bonferroni
        threshold_val = -np.log10(threshold_config['val'] / len(df))
        threshold_label = 'Bonferroni'
    elif threshold_config['type'] == 'fixed':
        # 固定值
        threshold_val = -np.log10(threshold_config['val'])
        threshold_label = f'P={threshold_config["val"]}'
        
    if threshold_val is not None:
        max_log_p = df['minuslog10p'].max()
        # 只有当最大值超过阈值时才画线
        if max_log_p >= threshold_val:
            ax.axhline(threshold_val, color='red', linestyle='--', linewidth=1)
            # 标签位置微调
            ax.text(0, threshold_val, f' {threshold_label}', 
                    color='red', fontsize=8, va='bottom', ha='left')
            print(f"      [阈值] 绘制线条于 -log10(p) = {threshold_val:.2f}")
        else:
            print(f"      [阈值] 数据最大值 ({max_log_p:.2f}) 未超过阈值 ({threshold_val:.2f})，隐藏线条。")
    # =======================================================
    
    output_png = os.path.join(output_dir, f"{base_name}.manhattan.png")
    plt.tight_layout()
    plt.savefig(output_png, dpi=FIG_DPI)
    plt.close(fig)
    print(f"      已保存: {os.path.basename(output_png)}")
    
    draw_qq(df['p_wald'], base_name, output_dir)

# ==============================================================================
# 主函数
# ==============================================================================
def main():
    print("==============================================")
    print("   Step 9: 批量曼哈顿图绘制工具 (升级版)")
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

    # 4. 获取阈值设置 (新增步骤)
    threshold_config = get_threshold_config()

    # 5. 执行绘图
    if make_sure("开始批量绘图"):
        print(f"\n--- 开始处理 {len(selected_files)} 个文件 ---")
        for i, f in enumerate(selected_files, 1):
            print(f"\n>>> 任务 [{i}/{len(selected_files)}]")
            # 将配置传入绘图函数
            draw_manhattan(f, output_dir, threshold_config)
            
        print("\n" + "="*50)
        print("所有图片已生成完毕！")
        print(f"存放目录: {output_dir}")
        print("="*50)

if __name__ == "__main__":
    main()