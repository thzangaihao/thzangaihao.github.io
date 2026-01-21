#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import os
import sys
import glob
import json
import pandas as pd

# ============= 基础路径获取（复用 cite_v1 逻辑） =============
def get_base_dir():
    """获取脚本所在目录，兼容 PyInstaller 打包"""
    if getattr(sys, 'frozen', False):
        return os.path.dirname(sys.executable)
    return os.path.dirname(os.path.abspath(__file__))

# ============= 高级文件查找（复用 cite_v1 逻辑） =============
def find_files(filename_pattern, path=None):
    """
    查找指定模式的文件
    自动在当前目录 + 子目录递归查找
    """
    if path is None:
        path = get_base_dir()
    # 使用 recursive=True 进行递归查找
    return glob.glob(os.path.join(path, '**', filename_pattern), recursive=True)

# ============= 核心处理逻辑 =============
def process_json_files(file_list):
    """
    读取文件列表，提取 ID 和 ptm 分数
    """
    data = []
    print(f"正在处理 {len(file_list)} 个文件...")

    for i, file_path in enumerate(file_list, 1):
        try:
            # 1. 获取文件名
            basename = os.path.basename(file_path)
            
            # 2. 清洗文件名以提取 ID
            # 假设文件名格式为: fold_fs000286_t01_summary_confidences_0.json
            # 去除前缀 'fold_'
            if basename.startswith('fold_'):
                gene_id = basename[5:] 
            else:
                gene_id = basename
            
            # 去除后缀 '_summary_confidences_0.json'
            suffix = '_summary_confidences_0.json'
            if gene_id.endswith(suffix):
                gene_id = gene_id[:-len(suffix)]
            
            # 3. 读取 JSON 获取 ptm
            with open(file_path, 'r', encoding='utf-8') as f:
                content = json.load(f)
                # 提取 ptm，如果没有则设为 None
                ptm_score = content.get('ptm', None)
                
                # 如果你想顺便提取 ranking_score，可以把下面这行注释打开
                # ranking_score = content.get('ranking_score', None)

            # 4. 存入列表
            data.append({
                'Gene_ID': gene_id,
                'pTM_Score': ptm_score
            })

            # 简单的进度打印
            if i % 100 == 0:
                print(f"已处理 {i}/{len(file_list)}...")

        except Exception as e:
            print(f"[Warning] 处理文件 {basename} 时出错: {e}")

    return pd.DataFrame(data)

# ============= 主函数 =============
def main():
    # 1. 确定搜索模式
    target_pattern = '*summary_confidences_0.json'
    base_dir = get_base_dir()
    
    print(f"正在目录 {base_dir} 下搜索 {target_pattern} ...")
    
    # 2. 搜索文件
    json_files = find_files(target_pattern, base_dir)
    
    if not json_files:
        print(f"未找到任何 {target_pattern} 文件，程序退出。")
        return

    print(f"共找到 {len(json_files)} 个相关 JSON 文件。")

    # 3. 提取数据
    df = process_json_files(json_files)

    # 4. 预览数据
    print("\n生成的矩阵预览:")
    print(df.head())
    
    # 5. 保存文件
    output_filename = os.path.join(base_dir, 'ptm_matrix.tsv')
    try:
        df.to_csv(output_filename, sep='\t', index=False)
        print("=" * 50)
        print(f"成功！矩阵已保存至: {output_filename}")
        print("=" * 50)
    except Exception as e:
        print(f"保存文件失败: {e}")

if __name__ == "__main__":
    main()