#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import sys
import glob
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

'''
GO 功能分类交互式绘图工具
功能：读取 GO_Classification_Stats.csv，支持自定义 Top N 和预览/保存逻辑。
'''

# ============= 基础路径获取 =============
def get_base_dir():
    if getattr(sys, 'frozen', False):
        return os.path.dirname(sys.executable)
    return os.path.dirname(os.path.abspath(__file__))

# ============= 搜索数据文件 =============
def find_files(ext):
    path = get_base_dir()
    search_pattern = os.path.join(path, '**', f'*{ext}')
    return sorted(glob.glob(search_pattern, recursive=True))

# ============= 交互式单选文件 =============
def choose_file(files):
    if not files:
        print("⚠️ 未找到任何相关文件，请检查目录或后缀名！")
        return None
    print("\n📂 找到以下数据文件:")
    for i, f in enumerate(files):
        print(f"[{i+1}] {os.path.relpath(f, get_base_dir())}")
    
    while True:
        choice = input("\n👉 请选择要绘图的数据文件编号 (例如 '1', 或 'q' 退出): ").strip()
        if choice.lower() == 'q':
            sys.exit(0)
        try:
            return files[int(choice)-1]
        except (ValueError, IndexError):
            print("⚠️ 输入无效，请重新输入正确的编号。")

# ============= 绘图主逻辑 =============
def plot_go_stats(file_path, top_n, save_prefix):
    print(f"\n⚙️ 正在加载数据: {os.path.basename(file_path)}")
    df = pd.read_csv(file_path)
    
    # 检查必要的列
    required_cols = {'Category', 'Description', 'Count'}
    if not required_cols.issubset(df.columns):
        print(f"❌ 错误：数据文件缺失必要的列，需包含 {required_cols}")
        return

    # 清洗异常值并提取 Top N
    df = df.dropna(subset=['Category'])
    df_sorted = df.sort_values(by=['Category', 'Count'], ascending=[True, False])
    df_top = df_sorted.groupby('Category').head(top_n)
    
    # 为了让图中每个类别内数值最大的排在最上面，需要调整 DataFrame 的排序
    df_top = df_top.sort_values(by=['Category', 'Count'], ascending=[True, False])

    # 设置画布大小（根据条目多少自适应高度）
    fig_height = max(6, int(len(df_top) * 0.4))
    plt.figure(figsize=(12, fig_height))
    sns.set_theme(style="whitegrid")

    # 定义三大类的经典配色
    palette = {
        'Biological Process (BP)': '#e74c3c',  # 红色
        'Cellular Component (CC)': '#2ecc71',  # 绿色
        'Molecular Function (MF)': '#3498db'   # 蓝色
    }

    # 绘制水平柱状图
    ax = sns.barplot(
        x='Count', 
        y='Description', 
        hue='Category', 
        data=df_top, 
        dodge=False, 
        palette=palette
    )

    # 细节调控：标题与标签
    plt.title(f'GO Classification of Mutated Target Genes (Top {top_n} per category)', 
              fontsize=16, fontweight='bold', pad=20)
    plt.xlabel('Number of Mutated Genes', fontsize=12, fontweight='bold')
    plt.ylabel('Gene Ontology (GO) Term', fontsize=12, fontweight='bold')
    
    # 细节调控：图例位置
    plt.legend(title='GO Category', loc='lower right', bbox_to_anchor=(1, 0.1))
    
    plt.tight_layout()

    # 预览与保存逻辑 (参考 CR_v4)
    if save_prefix:
        pdf_path = f"{save_prefix}_GO_Barplot.pdf"
        png_path = f"{save_prefix}_GO_Barplot.png"
        plt.savefig(pdf_path, dpi=300, bbox_inches='tight')
        plt.savefig(png_path, dpi=300, bbox_inches='tight')
        print(f"\n🎉 图表已成功保存至:")
        print(f"  - {pdf_path}")
        print(f"  - {png_path}")
    else:
        print("\n👀 正在打开预览窗口 (关闭窗口后程序结束)...")
        plt.show()

# ============= 主函数 =============
def main():
    print("="*50)
    print("   📊 交互式 GO 功能分类绘图工具")
    print("="*50)
    
    # 1. 自动搜索 CSV 文件并选择
    files = find_files('.csv')
    selected_file = choose_file(files)
    if not selected_file:
        return
        
    # 2. 交互式选择 Top N
    print("\n" + "="*50)
    top_n_str = input("👉 1. 请输入每个大类要显示的最多 GO 词条数量 (直接回车默认 10): ").strip()
    top_n = int(top_n_str) if top_n_str.isdigit() else 10

    # 3. 预览或保存逻辑
    print("\n" + "="*50)
    save_prefix = input("👉 2. 请输入保存文件的前缀名 (例如 'Salt_Tolerance_Target')\n   [直接回车则仅弹出窗口预览，不保存文件]: ").strip()
    
    plot_go_stats(selected_file, top_n, save_prefix)

if __name__ == "__main__":
    try:
        main()
    except KeyboardInterrupt:
        print("\n🛑 用户强制退出。")