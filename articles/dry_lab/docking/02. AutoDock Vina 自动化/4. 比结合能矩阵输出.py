#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import os
import sys
import glob
import pandas as pd

# ============= 1. 基础路径 =============
def get_base_dir():
    if getattr(sys, 'frozen', False):
        return os.path.dirname(sys.executable)
    return os.path.dirname(os.path.abspath(__file__))

# ============= 2. 获取完整名单 (确保矩阵完整性) =============
def get_all_names(directory, ext):
    """扫描原始文件夹，获取所有 受体/配体 名称"""
    if not ext.startswith('.'): ext = '.' + ext
    files = glob.glob(os.path.join(directory, f"*{ext}"))
    names = [os.path.splitext(os.path.basename(f))[0] for f in files]
    return sorted(list(set(names)))

# ============= 3. 解析 Log 文件 =============
def parse_vina_log(log_path):
    """提取 Best Mode 结合能"""
    try:
        with open(log_path, 'r') as f:
            start_reading = False
            for line in f:
                if "-----+------------" in line:
                    start_reading = True
                    continue
                if start_reading:
                    parts = line.strip().split()
                    if len(parts) >= 2 and parts[0] == "1":
                        return float(parts[1])
    except:
        pass
    return None

# ============= 4. 主流程 =============
def main():
    base_path = get_base_dir()
    
    # 路径定义
    lig_dir = os.path.join(base_path, "1. Ligand_molecular_pdbqt")
    rec_dir = os.path.join(base_path, "1. Receptor_protein_pdbqt")
    result_dir = os.path.join(base_path, "3. Docking_Results")
    output_file = os.path.join(base_path, "Final_Matrix_ReceptorRow.tsv")

    print("=== Step 4 (v3): 生成结合能矩阵 (行=受体, 列=配体) ===")

    # 1. 扫描原始文件，建立理论全名单
    print("正在扫描原始文件列表...")
    all_ligands = get_all_names(lig_dir, '.pdbqt')
    all_receptors = get_all_names(rec_dir, '.pdbqt')

    if not all_ligands or not all_receptors:
        print("错误: 在 1. 文件夹下未找到 PDBQT 文件，无法建立列表。")
        return

    print(f"  - 配体 (列): {len(all_ligands)} 个")
    print(f"  - 受体 (行): {len(all_receptors)} 个")

    # 2. 提取对接结果
    print("\n正在解析对接日志...")
    log_files = glob.glob(os.path.join(result_dir, "**", "*.log"), recursive=True)
    
    data_list = []
    
    for log_path in log_files:
        # 提取受体名
        parent_dir = os.path.dirname(log_path)
        rec_name = os.path.basename(parent_dir)
        
        # 提取配体名
        filename = os.path.basename(log_path)
        if "_vs_" in filename:
            lig_name = filename.split("_vs_")[0]
        else:
            lig_name = filename.replace(".log", "")

        # 提取结合能
        score = parse_vina_log(log_path)
        if score is not None:
            data_list.append({
                "Ligand": lig_name,
                "Receptor": rec_name,
                "Affinity": score
            })

    # 3. 构建矩阵
    print(f"提取到 {len(data_list)} 条有效数据，正在构建矩阵...")
    
    if data_list:
        df = pd.DataFrame(data_list)
        
        # 【关键改动】透视表：行(Index)=受体, 列(Columns)=配体
        matrix = df.pivot(index="Receptor", columns="Ligand", values="Affinity")
        
        # 强制重索引：确保所有受体和配体都在表中，哪怕全是空的
        matrix = matrix.reindex(index=all_receptors, columns=all_ligands)
        
        # 【关键改动】将 NaN (空值) 填充为 '*'
        matrix = matrix.fillna('*')
        
        # 保存
        matrix.to_csv(output_file, sep='\t')
        
        print("\n" + "="*50)
        print(f"矩阵已保存: {output_file}")
        print(f"维度: {matrix.shape[0]} 行 (受体) x {matrix.shape[1]} 列 (配体)")
        print("说明: 表格中标记为 '*' 的位置表示对接失败或数据缺失。")
    else:
        # 如果完全没有数据，也生成一个全是 * 的空表
        print("警告: 未提取到有效数据，正在生成全空表...")
        empty_df = pd.DataFrame('*', index=all_receptors, columns=all_ligands)
        empty_df.to_csv(output_file, sep='\t')
        print(f"空表已保存: {output_file}")

    input("按回车键退出...")

if __name__ == "__main__":
    main()