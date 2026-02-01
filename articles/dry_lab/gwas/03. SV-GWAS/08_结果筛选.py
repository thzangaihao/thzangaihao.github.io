#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import sys
import glob
import pandas as pd

# ==============================================================================
# 0. 核心配置参数
# ==============================================================================
# 显著性阈值 (P-value Threshold)
# 常用标准:
#   - 严谨 (Bonferroni): 0.05 / 变异总数 (例如 1e-6)
#   - 宽松 (初筛): 1e-4 或 1e-5
P_VALUE_THRESHOLD = 1e-5

# 输出目录名称
OUTPUT_DIR_NAME = "08_Significant_Results"

# ==============================================================================
# 1. Cite2 交互逻辑模块 (内置)
# ==============================================================================
def get_base_dir():
    if getattr(sys, 'frozen', False):
        return os.path.dirname(sys.executable)
    return os.path.dirname(os.path.abspath(__file__))

def find_files(ext, path=None):
    if not ext.startswith('.'): ext = '.' + ext
    if path is None: path = get_base_dir()
    search_pattern = os.path.join(path, '**', f'*{ext}')
    all_files = sorted(glob.glob(search_pattern, recursive=True))
    # 排除之前的输出目录，避免重复处理结果文件
    return [f for f in all_files if OUTPUT_DIR_NAME not in f and "plot" not in f]

def choose_files(files, desc="文件"):
    if not files:
        print(f"[提示] 未找到任何 {desc}")
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
                f"\n3. 请输入欲筛选的 {desc} 编号\n"
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
# 2. 核心处理模块
# ==============================================================================
def extract_significant_sites(file_list):
    # 1. 准备输出目录
    base_dir = get_base_dir()
    out_dir = os.path.join(base_dir, OUTPUT_DIR_NAME)
    
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
        print(f"\n[系统] 创建结果存储目录: {out_dir}")

    print(f"\n--- 开始筛选显著位点 (阈值 P < {P_VALUE_THRESHOLD}) ---")
    
    summary_list = [] # 用于最后汇总统计

    for idx, file_path in enumerate(file_list, 1):
        filename = os.path.basename(file_path)
        print(f"\n>>> [{idx}/{len(file_list)}] 正在处理: {filename}")
        
        try:
            # 读取数据 (自动处理空格或Tab分隔)
            try:
                df = pd.read_csv(file_path, sep='\t')
            except:
                df = pd.read_csv(file_path, delim_whitespace=True)
            
            # 检查关键列名
            if 'p_wald' not in df.columns:
                print(f"  [跳过] 文件缺少 'p_wald' 列，可能不是 GEMMA 结果。")
                continue
                
            # 核心筛选逻辑
            sig_df = df[df['p_wald'] < P_VALUE_THRESHOLD].copy()
            
            # 按照 P 值从小到大排序 (最显著的在前面)
            sig_df = sig_df.sort_values(by='p_wald', ascending=True)
            
            num_sig = len(sig_df)
            print(f"  -> 发现 {num_sig} 个显著位点")
            
            summary_list.append({'File': filename, 'Significant_Hits': num_sig})

            if num_sig > 0:
                # 构造输出文件名
                name_prefix = os.path.splitext(filename)[0]
                # 移除 .assoc 后缀如果存在
                if name_prefix.endswith('.assoc'): name_prefix = name_prefix[:-6]
                
                out_name = f"{name_prefix}_sig_sites.txt"
                out_path = os.path.join(out_dir, out_name)
                
                # 保存结果
                sig_df.to_csv(out_path, sep='\t', index=False)
                print(f"  [保存] 已导出至: {out_name}")
                
                # 可选：导出 Excel 版本方便查看 (需要 openpyxl 库，如果没有则跳过)
                try:
                    excel_path = os.path.join(out_dir, f"{name_prefix}_sig_sites.xlsx")
                    sig_df.to_excel(excel_path, index=False)
                    print(f"  [保存] Excel版本: {os.path.basename(excel_path)}")
                except:
                    pass
            else:
                print("  [提示] 未找到满足阈值的位点，不生成文件。")

        except Exception as e:
            print(f"  [错误] 处理失败: {e}")

    # --- 最终汇总报告 ---
    print("\n" + "="*50)
    print("【筛选工作完成】")
    print(f"设定阈值: P < {P_VALUE_THRESHOLD}")
    print(f"结果目录: {out_dir}")
    print("-" * 50)
    print(f"{'Trait File':<40} | {'Hits':<10}")
    print("-" * 50)
    for item in summary_list:
        print(f"{item['File']:<40} | {item['Significant_Hits']:<10}")
    print("="*50)

# ==============================================================================
# 主函数
# ==============================================================================
def main():
    print("==============================================")
    print("   Step 8: 显著变异位点筛选工具 (GWAS Miner)")
    print("==============================================")
    print(f"当前筛选阈值: P < {P_VALUE_THRESHOLD}")
    print("提示: 如需修改阈值，请用文本编辑器打开本脚本修改第 16 行。")
    
    # 1. 搜索 GEMMA 结果文件 (.assoc.txt)
    print("\n>>> 第一步: 搜索 GWAS 结果文件 (.assoc.txt)")
    files = find_files('.assoc.txt')
    
    # 2. 交互式选择
    selected = choose_files(files, "GWAS 结果文件")
    if not selected: return

    # 3. 执行筛选
    if make_sure("开始筛选"):
        extract_significant_sites(selected)

if __name__ == "__main__":
    main()