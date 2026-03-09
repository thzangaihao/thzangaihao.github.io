#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import sys
import glob
import pandas as pd
import matplotlib.pyplot as plt

# ============= 图表字体设置 (支持中文) =============
plt.rcParams['font.sans-serif'] = ['SimHei', 'Arial Unicode MS', 'DejaVu Sans'] 
plt.rcParams['axes.unicode_minus'] = False

# ============= 基础路径获取 =============
def get_base_dir():
    if getattr(sys, 'frozen', False):
        return os.path.dirname(sys.executable)
    return os.path.dirname(os.path.abspath(__file__))

# ============= 高级文件查找 =============
def find_files(ext, path=None):
    if not ext.startswith('.'):
        ext = '.' + ext
    if path is None:
        path = get_base_dir()
    search_pattern = os.path.join(path, '**', f'*{ext}')
    return sorted(glob.glob(search_pattern, recursive=True))

# ============= 交互式文件选择 =============
def choose_file(files, desc="文件"):
    if not files:
        print(f"当前目录下未找到 {desc}！")
        return []
    
    print(f"\n找到以下 {len(files)} 个 {desc}:")
    for i, f in enumerate(files):
        print(f"  [{i+1}] {os.path.basename(f)}")
        
    print("\n💡 选择模式: 单选(1), 多选(1,3), 范围(1-3), 全部(all)")
    
    while True:
        try:
            choice = input(f"👉 请选择要可视化的 {desc} 编号 (输入 q 退出): ").strip().lower()
            if choice == 'q':
                sys.exit()
            if choice == 'all':
                return files
                
            selected_paths = []
            parts = choice.replace(' ', '').split(',')
            indices = set()
            for part in parts:
                if '-' in part:
                    start, end = map(int, part.split('-'))
                    indices.update(range(start, end + 1))
                else:
                    indices.add(int(part))
            
            for idx in sorted(list(indices)):
                if 1 <= idx <= len(files):
                    selected_paths.append(files[idx-1])
                else:
                    print(f"⚠️ 警告：编号 {idx} 超出范围，已忽略。")
            
            if selected_paths:
                return selected_paths
        except ValueError:
            print("输入错误：请确保输入的是数字编号、范围（如1-5）或 'all'。")
        except KeyboardInterrupt:
            print("\n用户取消操作。")
            sys.exit()

# ============= 确认函数 =============
def make_sure(selected_paths):
    print("\n" + "="*40)
    print("准备对以下文件进行可视化分析:")
    for p in selected_paths:
        print(f"  - {os.path.basename(p)}")
    print("="*40)
    response = input("确认执行? (y/n): ").strip().lower()
    return response in ['y', 'yes']

# ============= 解析与绘图核心 =============
def plot_kraken_report(report_path):
    data = []
    unclassified_pct = 0
    with open(report_path, 'r', encoding='utf-8') as f:
        for line in f:
            parts = line.strip('\n').split('\t')
            if len(parts) < 6:
                continue
            
            # 关键修复：加入 .strip() 清除 Kraken2 排版用的空格
            pct = float(parts[0].strip())
            rank = parts[3].strip() 
            name = parts[5].strip() 
            
            if rank == 'U':
                unclassified_pct = pct
            elif rank.startswith('S'): 
                data.append((name, pct))
    
    if not data:
        print(f"⚠️ {os.path.basename(report_path)} 中未找到物种级别(S)的数据，跳过绘制。")
        return

    # 转换为 DataFrame 并提取 Top 10
    df = pd.DataFrame(data, columns=['Species', 'Percentage'])
    # 由于存在同种下的亚种 (比如 S 和 S1)，可能出现重复名字或父子包含，这里直接按丰度排
    df = df.sort_values(by='Percentage', ascending=False).head(10)
    df = df.sort_values(by='Percentage', ascending=True) 

    # 绘制条形图
    plt.figure(figsize=(10, 6))
    bars = plt.barh(df['Species'], df['Percentage'], color='#4C72B0', edgecolor='black')
    
    # 添加数据标签
    for bar in bars:
        width = bar.get_width()
        plt.text(width + 0.3, bar.get_y() + bar.get_height()/2, 
                 f'{width:.2f}%', va='center', fontsize=10)

    # 图表美化
    sample_name = os.path.basename(report_path).replace('_report.txt', '').replace('.txt', '')
    plt.title(f'物种丰度分析 - 样本: {sample_name}', fontsize=14, pad=20)
    plt.xlabel('相对丰度 (%)', fontsize=12)
    plt.ylabel('物种名称', fontsize=12)
    
    # 动态调整 X 轴，防止标签超出边界
    max_pct = max(df['Percentage']) if not df.empty else 0
    plt.xlim(0, max_pct * 1.25) 
    
    # 在右上角标注未分类比例
    plt.text(0.95, 0.05, f'Unclassified: {unclassified_pct:.2f}%', 
             transform=plt.gca().transAxes, ha='right', va='bottom',
             bbox=dict(boxstyle='round,pad=0.5', facecolor='lightyellow', alpha=0.8))
    
    plt.tight_layout()
    
    # 保存图片
    out_name = os.path.join(os.path.dirname(report_path), f"{sample_name}_abundance_barplot.png")
    plt.savefig(out_name, dpi=300)
    plt.close()
    print(f"✅ 生成可视化结果: {os.path.basename(out_name)}")

# ============= 主函数 =============
def main():
    print("--- Kraken2 鉴定报告自动可视化工具 ---")
    
    all_txt = find_files('report.txt') 
    if not all_txt:
        all_txt = find_files('.txt')
        
    selected_files = choose_file(all_txt, "Kraken2 报告 (.txt)")
    
    if selected_files and make_sure(selected_files):
        print("\n🚀 开始生成图表...")
        for report in selected_files:
            plot_kraken_report(report)
        print("\n🎉 所有选中样本的可视化已完成！")

if __name__ == "__main__":
    main()