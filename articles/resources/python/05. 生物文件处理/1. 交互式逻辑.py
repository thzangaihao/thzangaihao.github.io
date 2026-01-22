#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import os
import sys
import glob
import pandas as pd
import time

# ============= 基础路径获取 =============
def get_base_dir():
    """获取脚本所在目录，作为搜索的起点"""
    if getattr(sys, 'frozen', False):
        return os.path.dirname(sys.executable)
    return os.path.dirname(os.path.abspath(__file__))

# ============= 高级文件查找（优化：仅限同级及子级） =============
def find_files(ext, path=None):
    """
    查找指定扩展名文件。
    - 仅在当前路径(path)及其子目录下递归查找。
    """
    if not ext.startswith('.'):
        ext = '.' + ext
    if path is None:
        path = get_base_dir()
    search_pattern = os.path.join(path, '**', f'*{ext}')
    return sorted(glob.glob(search_pattern, recursive=True))

# ============= 文件选择（支持多种选取模式） =============
def choose_file(files, desc="文件"):
    """
    支持多种选取模式：
    1. 单选：1
    2. 多选：1,2,3
    3. 范围：1-4
    4. 组合：1,2,5-10
    5. 全部：all
    """
    if not files:
        print(f"提示：在当前目录及其子目录下未找到任何 {desc}")
        return []

    print(f"\n2. 找到 {len(files)} 个 {desc}:")
    for i, f in enumerate(files, 1):
        print(f"  [{i}] {f}")

    while True:
        try:
            prompt = (
                f"\n3. 请输入欲操作的 {desc} 编号\n"
                f" (支持格式: 1,2 | 1-4 | 1,3,5-8 | all): \n"
            )
            user_input = input(prompt).strip().lower()

            if not user_input:
                continue

            selected_indices = set()

            # 模式 1: 全部选择
            if user_input == 'all' or 'ALL' or 'a' or 'A':
                selected_indices = set(range(len(files)))
            
            # 模式 2: 解析组合输入
            else:
                parts = user_input.split(',')
                for part in parts:
                    part = part.strip()
                    if '-' in part:
                        # 处理 1-4 格式
                        start_str, end_str = part.split('-')
                        start, end = int(start_str), int(end_str)
                        # 转换为 0-based 索引并加入集合
                        for idx in range(start - 1, end):
                            if 0 <= idx < len(files):
                                selected_indices.add(idx)
                    else:
                        # 处理单个数字
                        idx = int(part) - 1
                        if 0 <= idx < len(files):
                            selected_indices.add(idx)
            
            if not selected_indices:
                print("未匹配到有效编号，请重新输入。")
                continue

            # 转换为路径列表并按原始顺序排序
            selected_paths = [files[i] for i in sorted(list(selected_indices))]
            
            print(f"\n4. 已选择 {len(selected_paths)} 个文件:")
            for p in selected_paths:
                print(f"  - {os.path.basename(p)}")
            
            return selected_paths

        except ValueError:
            print("输入错误：请确保输入的是数字编号、范围（如1-5）或 'all'。")
        except KeyboardInterrupt:
            print("\n用户取消操作。")
            sys.exit()

# ============= 确认函数 =============
def make_sure():
    response = input(f"\n5. 确认对上述选中的文件执行后续操作? (y/n): ").strip().lower()
    if response not in ['y', 'yes']:
        print("操作已取消。")
        return False
    return True

# ============= DataFrame 保存 =============
def save_dataframe(select_df):
    if make_sure():
        save_name = input(f'请输入保存文件名（无需后缀）: \n').strip()
        current_path = get_base_dir()
        output_filename = os.path.join(current_path, f'{save_name}_matrix.tsv')
        select_df.to_csv(output_filename, sep='\t', index=True)
        print("-" * 50)
        print(f"矩阵已成功保存: {output_filename}")

# ============= 主函数 =============
def main():
    print("--- 批量文件搜索与多选工具 ---")
    
    target_ext = input("1. 请输入要搜索的文件后缀名: ").strip()
    if not target_ext:
        print("未输入后缀，程序退出")
        time.sleep(1)
        exit()

    # 1. 查找文件
    files = find_files(target_ext)
    
    # 2. 获取选中的路径列表
    chosen_list = choose_file(files, desc=f"{target_ext} 文件")
    
    # 3. 后续处理示例
    if chosen_list and make_sure():
        print(f"\n程序开始处理这 {len(chosen_list)} 个文件...")
        # 此处可接入你的批量处理逻辑

if __name__ == "__main__":
    main()