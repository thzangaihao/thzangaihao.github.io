#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import sys
import glob
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

# ==============================================================================
# 0. 核心配置
# ==============================================================================
OUTPUT_DIR_NAME = "12_Kinship_Plots"
FIG_SIZE = (10, 10) # 图片尺寸
COLOR_MAP = "YlGnBu"  # 颜色方案：黄色-绿色-蓝色

# ==============================================================================
# 1. Cite2 交互逻辑模块 (复用)
# ==============================================================================
def get_base_dir():
    return os.path.dirname(os.path.abspath(__file__))

def find_files(ext, path=None):
    if not ext.startswith('.'): ext = '.' + ext
    if path is None: path = get_base_dir()
    search_pattern = os.path.join(path, '**', f'*{ext}')
    all_files = sorted(glob.glob(search_pattern, recursive=True))
    return [f for f in all_files if OUTPUT_DIR_NAME not in f]

def choose_file(files, desc="文件"):
    if not files:
        print(f"[提示] 未找到任何 {desc}")
        return None

    print(f"\n2. 找到 {len(files)} 个 {desc}:")
    for i, f in enumerate(files, 1):
        print(f"  [{i}] {os.path.basename(f)}  (路径: {os.path.relpath(f, get_base_dir())})")

    while True:
        try:
            u_input = input(f"\n3. 请选择 {desc} 编号: ").strip()
            idx = int(u_input) - 1
            if 0 <= idx < len(files): return files[idx]
        except ValueError: print("请输入有效数字。")

# ==============================================================================
# 2. 核心绘图模块 (优化版)
# ==============================================================================
def plot_kinship(kinship_path, fam_path):
    print(f"\n--- 正在准备绘图数据 ---")
    
    # 1. 读取样本 ID
    try:
        fam_df = pd.read_csv(fam_path, sep='\s+', header=None)
        sample_ids = fam_df[1].astype(str).tolist()
        print(f"-> 成功获取 {len(sample_ids)} 个样本 ID")
    except Exception as e:
        print(f"[错误] 读取 FAM 文件失败: {e}")
        return

    # 2. 读取 Kinship 矩阵
    try:
        kin_matrix = pd.read_csv(kinship_path, sep='\t', header=None)
        
        # 维度校验与截断
        if kin_matrix.shape[0] != len(sample_ids):
            print(f"[警告] 矩阵维度 ({kin_matrix.shape[0]}) 与样本数 ({len(sample_ids)}) 不匹配！")
            size = min(kin_matrix.shape[0], len(sample_ids))
            kin_matrix = kin_matrix.iloc[:size, :size]
            sample_ids = sample_ids[:size]
            
        kin_matrix.index = sample_ids
        kin_matrix.columns = sample_ids
    except Exception as e:
        print(f"[错误] 读取 Kinship 矩阵失败: {e}")
        return

    # 3. 绘图 (优化重点)
    print(f"   正在生成热图 (含聚类)... 这可能需要一点时间...")
    
    try:
        # 使用 clustermap

        # 定义 Colorbar 的位置 (left, bottom, width, height)
        # 例如：放在画布右侧外面，高度适中
        # left=1.02 (刚好在画布右边缘外)
        # bottom=0.3 (距离底部30%的位置)
        # width=0.03 (宽度为画布宽度的3%)
        # height=0.4 (高度为画布高度的40%)
        my_cbar_pos = (1.02, 0.7, 0.03, 0.2)
        
        g = sns.clustermap(
            kin_matrix, 
            cmap=COLOR_MAP, 
            figsize=FIG_SIZE,
            
            # [优化1] 移除坐标轴标签 (xticklabels, yticklabels)
            xticklabels=False,
            yticklabels=False,

            cbar_pos=my_cbar_pos,
            
            # [优化2] 调整树状图比例，留出更多空间
            dendrogram_ratio=(0.05, 0.05),
            
            # [优化3] 移除手动 cbar_pos，使用默认布局防止重叠
            # 或者将其放在更安全的位置 (如右上角)
            # 这里我们使用默认值，Seaborn 会自动安排不重叠的位置
            cbar_kws={"label": "Kinship Coefficient"}
        )
        
        # [优化4] 调整标题位置，防止与顶部聚类树重叠
        # y=1.02 表示在图片最顶端再往上一点点
        title_text = f'Kinship Matrix: {os.path.basename(kinship_path)}'
        g.fig.suptitle(title_text, fontsize=16, y=1.02)
        
        # 准备输出
        base_dir = get_base_dir()
        out_dir = os.path.join(base_dir, OUTPUT_DIR_NAME)
        if not os.path.exists(out_dir): os.makedirs(out_dir)
        
        prefix = os.path.basename(kinship_path).replace(".cXX.txt", "")
        output_path = os.path.join(out_dir, f"{prefix}_heatmap_opt.png")
        
        # 保存时使用 bbox_inches='tight' 自动裁剪空白
        plt.savefig(output_path, dpi=300, bbox_inches='tight')
        plt.close()
        print(f"\n[成功] 优化后的亲缘关系热图已保存至: {output_path}")

    except Exception as e:
        print(f"[绘图失败] {e}")
        import traceback
        traceback.print_exc()

# ==============================================================================
# 主函数
# ==============================================================================
def main():
    print("==============================================")
    print("   Step 12: 亲缘关系矩阵可视化 (优化版)")
    print("==============================================")
    
    # 1. 选择 Kinship 文件
    kin_files = find_files('.cXX.txt')
    target_kin = choose_file(kin_files, "Kinship 矩阵 (.cXX.txt)")
    if not target_kin: return

    # 2. 选择对应的 FAM 文件
    fam_files = find_files('.fam')
    target_fam = choose_file(fam_files, "对应的 FAM 文件 (.fam)")
    if not target_fam: return

    # 3. 绘图
    plot_kinship(target_kin, target_fam)

if __name__ == "__main__":
    main()