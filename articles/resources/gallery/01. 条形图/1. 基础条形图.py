#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import os
import glob
import sys
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import scipy.stats

# ====================================文件查找
def get_base_dir():
    """获取脚本所在目录，兼容 PyInstaller 打包"""
    if getattr(sys, 'frozen', False):
        return os.path.dirname(sys.executable)
    return os.path.dirname(os.path.abspath(__file__))
def find_files(ext, path=None):
    """
    查找指定扩展名文件，如 ext='.tsv'
    自动在当前目录 + 子目录递归查找
    """
    if not ext.startswith('.'):
        ext = '.' + ext
    if path is None:
        path = get_base_dir()
    return glob.glob(os.path.join(path, '**', f'*{ext}'), recursive=True)
def choose_file(files, desc="文件"):
    """
    自动选择文件：
      - 如果 0 个 → 返回 None
      - 如果 1 个 → 自动选择
      - 如果多个 → 用户交互选择
    """
    if not files:
        print(f"未找到任何 {desc}")
        return None

    if len(files) == 1:
        print(f"自动选择 {desc}: {files[0]}")
        return files[0]

    print(f"找到 {len(files)} 个 {desc}:")
    for i, f in enumerate(files, 1):
        print(f"  [{i}] {f}")

    while True:
        try:
            c = input(f"请选择 {desc} 编号（1-{len(files)}）: ").strip()
            return files[int(c) - 1]
        except:
            print("无效输入，请重新输入。")
def make_sure():
    print("=" * 50)
    response = input("是否继续? (y/n): ").strip().lower()
    if response not in ['y', 'yes']:
        print("操作已取消, 程序已退出")
        sys.exit()
    return True
def save_dataframe(select_df):
    if make_sure():
        save_name = input(f'请输入保存文件名: \n').strip()
        current_path = get_base_dir()
        output_filename = os.path.join(current_path, f'{save_name}_matrix.tsv')
        select_df.to_csv(output_filename, sep='\t', index=True)
        print("=" * 50)
        print(f"矩阵已保存到脚本同级目录下: {output_filename}")
# ====================================文件查找
def draw(df):
    # 设置中文字体
    plt.rcParams['font.sans-serif'] = ['SimHei', 'Arial Unicode MS', 'DejaVu Sans']
    plt.rcParams['axes.unicode_minus'] = False

    # 1. 数据准备
    x_labels = df.iloc[:, 0].tolist()
    # 1.1 柱状图数据
    bar_data = []
    for index in df.iloc[:, 1:4].values.tolist():
        bar_data.append(np.mean(index))
    # 1.2 散点图数据
    scatter_data = []
    max_scatter_data = []
    for index in df.iloc[:, 1:4].values.tolist():
        scatter_data.append(index)
        max_scatter_data.append(max(index))
        
    # 1.3 标准差数据
    err_data = []
    for index in df.iloc[:, 1:4].values.tolist():
        err_data.append(np.std(index))
    # 1.4 显著性数据
    pvalues_data = []
    for index in df.iloc[1:, 4].values:
        pvalues_data.append(index)

    # 2. 创建画布并分离对象
    fig, ax = plt.subplots(figsize=(16, 10), dpi=80)
    x = np.arange(len(x_labels)) # 设置x轴初始刻度

    # 3. 图像绘制
    # 3.1 柱状图绘制
    bar_width = 0.5
    ax.bar(x=x, height=bar_data, width=bar_width,
           yerr = err_data, capsize=5, error_kw={'elinewidth': 1, 'markeredgewidth': 1},
           color = '#f5785f')
    # 3.2 散点图绘制
    for i, me_list in enumerate(scatter_data):
        # 为每个柱子的三个重复值生成x坐标
        x_positions = [x[i] - bar_width/4, x[i], x[i] + bar_width/4]
        # 绘制三个散点
        ax.scatter(x_positions, me_list, color="#000000", s=10, zorder=5)
    # 3.3 显著性绘制
    pliter = 3 # 上下偏移量
    start_pliter = pliter # 起始偏移量
    bar_plot_pliter = 2 # 最大点与显著性标记偏移量
    max_y2 = max(max_scatter_data) + bar_plot_pliter
    left_y1 = max(scatter_data[0]) + bar_plot_pliter
    left_x = 0
    for right_x, p_sign, right_y in zip(x[1:], pvalues_data, max_scatter_data[1:]):
        
        ax.plot([left_x, left_x, left_x, right_x, right_x],
                [left_y1, max_y2, max_y2+start_pliter, max_y2+start_pliter, right_y+bar_plot_pliter],
                lw=1.5, c='black')

        ax.text((left_x+right_x)*0.5, max_y2+start_pliter, p_sign,
                ha='center', va='bottom', color='black', size=15)
        
        start_pliter += pliter

    # 4. 图像设置
    ax.set_ylabel('生长直径 (mm)', size=15)
    ax.set_xlabel('菌株编号', size=15)

    ax.set_xticks(x, x_labels, rotation=25, size=15)
    plt.subplots_adjust(bottom=0.15)  # 0.1~0.3 之间自行调节

    ax.set_title(input('请输入图表标题：'), loc='left')  # 设置图表标题

    # 5. 显示图像
    plt.show()

# ============= 主函数（保持简洁，可被其它脚本调用） =============
def main():
    # 选择 .tsv 文件
    files = find_files('.tsv')
    chosen = choose_file(files, desc=".tsv 文件")
    if chosen is None:
        print("未找到 .tsv 文件，程序退出。")
        return

    df = pd.read_csv(chosen, sep='\t')
    draw(df=df)

if __name__ == "__main__":
    main()
