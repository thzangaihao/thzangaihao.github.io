import os
import glob
import sys
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

def current_path_function():
    # 获取基础路径，兼容打包模式
    if getattr(sys, 'frozen', False):
        # 打包后的情况 - 使用可执行文件所在目录
        current_path = os.path.dirname(sys.executable)
    else:
        # 普通Python脚本运行的情况
        current_path = os.path.dirname(os.path.abspath(__file__))

    return current_path

def find_file(file_type):
    current_path = current_path_function()

    # 在当前及其子目录下搜索
    current_file_path = os.path.join(current_path, '**', f'*{file_type}')
    current_file_list = glob.glob(current_file_path, recursive=True)

    print('='*50)
    if len(current_file_list) >= 1:
        return current_file_list

    else:
        response = input(f'在当前目录及子目录下没有发现{file_type}文件, 是否在父级目录搜索(y/n): \n').strip().lower()
        if response not in ['y', 'yes']:
            print("操作已取消, 程序已退出")
            exit()
        else:
            parent_path = os.path.dirname(current_path)

            # 在父级目录下搜索
            parent_file_path = os.path.join(parent_path, '**', f'*{file_type}')
            parent_file_list = glob.glob(parent_file_path, recursive=True)

            if len(parent_file_list) >= 1:
                return parent_file_list
            else:
                print('未发现任何文件, 程序已退出')
                exit()

def make_sure():
    print('='*50)
    response = input("是否继续? (y/n): \n").strip().lower()
    if response not in ['y', 'yes']:
        print("操作已取消, 程序已退出")
    else:
        return True

def choose_file(file_list):
    number = 1
    print('='*50)
    print(f'共发现{len(file_list)}个文件: ')
    for file in file_list:
        print(f'  [{number}]{file}')
        number += 1
    print('='*50)
    response = int(input('请选择目标文件: '))
    choosing_file = file_list[response - 1]

    return choosing_file

def save_dataframe(select_df, save_name='default'):
    # 获取当前路径
    current_path = current_path_function()
    
    # 构建完整的文件路径
    output_filename = os.path.join(current_path, save_name)
    
    # 保存到TSV文件
    select_df.to_csv(output_filename, sep='\t', index=True)
    print('='*50)
    print(f"矩阵已保存到脚本同级目录下: {output_filename}")

def create_volcano_plot(deseq_df, padj_threshold=0.05, log2fc_threshold=1, 
                       top_genes=10, title=None):

    
    # 创建图形
    fig, ax = plt.subplots(figsize=(12, 8), dpi=100)
    
    # 准备数据
    df_plot = deseq_df.copy()
    df_plot['-log10_padj'] = -np.log10(df_plot['padj'])
    
    # 分类基因
    conditions = [
        (df_plot['padj'] < padj_threshold) & (df_plot['log2FoldChange'] > log2fc_threshold),
        (df_plot['padj'] < padj_threshold) & (df_plot['log2FoldChange'] < -log2fc_threshold),
        (df_plot['padj'] < padj_threshold) & (abs(df_plot['log2FoldChange']) <= log2fc_threshold),
        (df_plot['padj'] >= padj_threshold) & (abs(df_plot['log2FoldChange']) > log2fc_threshold)
    ]
    choices = ['Up-regulated', 'Down-regulated', 'FDR significant only', 'FC significant only']
    df_plot['category'] = np.select(conditions, choices, default='Not significant')
    
    # 颜色映射
    color_map = {
        'Up-regulated': 'red',
        'Down-regulated': 'blue', 
        'FDR significant only': 'purple',
        'FC significant only': 'orange',
        'Not significant': 'gray'
    }
    
    # 绘制散点图
    for category, color in color_map.items():
        mask = df_plot['category'] == category
        if mask.sum() > 0:  # 只有当该类别有点时才绘制
            ax.scatter(df_plot.loc[mask, 'log2FoldChange'],
                      df_plot.loc[mask, '-log10_padj'],
                      c=color, label=category, alpha=0.6, s=1)
    
    # 添加阈值线
    ax.axhline(y=-np.log10(padj_threshold), color='black', linestyle='--', alpha=0.8, linewidth=1)
    ax.axvline(x=log2fc_threshold, color='black', linestyle='--', alpha=0.8, linewidth=1)
    ax.axvline(x=-log2fc_threshold, color='black', linestyle='--', alpha=0.8, linewidth=1)
    
    # 标注top基因
    if top_genes > 0:
        # 选择最显著的差异表达基因
        significant_genes = df_plot[
            (df_plot['padj'] < padj_threshold) & 
            (abs(df_plot['log2FoldChange']) > log2fc_threshold)
        ].copy()
        
        if len(significant_genes) > 0:
            # 按显著性排序并取前top_genes个
            significant_genes = significant_genes.nsmallest(min(top_genes, len(significant_genes)), 'padj')
            
            for idx, gene in significant_genes.iterrows():
                ax.annotate(idx,  # 使用索引作为基因名
                           xy=(gene['log2FoldChange'], gene['-log10_padj']),
                           xytext=(5, 5), textcoords='offset points',
                           fontsize=8, alpha=0.8,
                           arrowprops=dict(arrowstyle='->', color='black', alpha=0.6, lw=0.5))
    
    # 设置标签和标题
    ax.set_xlabel('log2 Fold Change', fontsize=12)
    ax.set_ylabel('-log10(adjusted p-value)', fontsize=12)
    
    if title is None:
        title = f'Volcano Plot (FDR < {padj_threshold}, |log2FC| > {log2fc_threshold})'
    ax.set_title(title, fontsize=14)
    
    # 添加图例
    ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
    
    # 统计信息
    up_regulated = len(df_plot[df_plot['category'] == 'Up-regulated'])
    down_regulated = len(df_plot[df_plot['category'] == 'Down-regulated'])
    fdr_only = len(df_plot[df_plot['category'] == 'FDR significant only'])
    fc_only = len(df_plot[df_plot['category'] == 'FC significant only'])
    
    # 在图上添加统计信息
    stats_text = f'Up-regulated: {up_regulated}\nDown-regulated: {down_regulated}\nFDR only: {fdr_only}\nFC only: {fc_only}'
    ax.text(0.02, 0.98, stats_text, transform=ax.transAxes, 
            verticalalignment='top', fontsize=10,
            bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
    
    # 设置网格
    ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    return fig, ax

def get_deg_stats(deseq_df, padj_threshold=0.05, log2fc_threshold=1):
    """获取差异表达基因的统计信息"""
    up_regulated = len(deseq_df[
        (deseq_df['padj'] < padj_threshold) & 
        (deseq_df['log2FoldChange'] > log2fc_threshold)
    ])
    
    down_regulated = len(deseq_df[
        (deseq_df['padj'] < padj_threshold) & 
        (deseq_df['log2FoldChange'] < -log2fc_threshold)
    ])
    
    total_deg = up_regulated + down_regulated
    
    return {
        'up_regulated': up_regulated,
        'down_regulated': down_regulated,
        'total_deg': total_deg,
        'thresholds': f'FDR < {padj_threshold}, |log2FC| > {log2fc_threshold}'
    }

def main():
    choosing_file = choose_file(file_list=find_file('.tsv'))
    df = pd.DataFrame(pd.read_csv(choosing_file, sep='\t'))

    """
    创建DESeq2结果的火山图
    
    参数:
    deseq_df: 包含DESeq2结果的DataFrame, 必须包含'log2FoldChange'和'padj'列
    padj_threshold: 显著性阈值 (默认: 0.05)
    log2fc_threshold: log2 fold change(LFC)阈值 (默认: 1)
    top_genes: 要标注的top基因数量 (默认: 10)
    figsize: 图形大小 (默认: (12, 8))
    title: 图表标题 (默认: 自动生成)
    
    返回:
    fig, ax: matplotlib的图形和坐标轴对象
    """
    create_volcano_plot(
        deseq_df=df, 
        padj_threshold=0.05,
        log2fc_threshold=1,
        top_genes=15,
        title='Differential Expression Analysis - Volcano Plot'
    )
    
    # 获取统计信息
    stats = get_deg_stats(df)
    print('='*50)
    print("差异表达基因统计:")
    for key, value in stats.items():
        print(f"{key}: {value}")

    # 显示图形
    plt.show()

main()