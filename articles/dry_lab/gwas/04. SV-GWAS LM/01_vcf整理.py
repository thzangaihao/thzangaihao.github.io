#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import os
import sys
import glob
import shutil
import time

# ============= 基础路径获取 =============
def get_base_dir():
    """获取脚本所在目录，作为搜索的起点"""
    if getattr(sys, 'frozen', False):
        return os.path.dirname(sys.executable)
    return os.path.dirname(os.path.abspath(__file__))

# ============= 高级文件查找 =============
def find_files(ext, path=None):
    """查找指定扩展名文件，递归子目录"""
    if not ext.startswith('.'):
        ext = '.' + ext
    if path is None:
        path = get_base_dir()
    # 搜索当前目录及所有子目录
    search_pattern = os.path.join(path, '**', f'*{ext}')
    return sorted(glob.glob(search_pattern, recursive=True))

# ============= 文件选择 =============
def choose_file(files, desc="文件"):
    """交互式文件选择逻辑"""
    if not files:
        print(f"提示：在当前目录及其子目录下未找到任何 {desc}")
        return []

    print(f"\n2. 找到 {len(files)} 个 {desc}:")
    if len(files) > 20:
         print(f"  ... (共 {len(files)} 个文件，列表过长仅显示部分示例) ...")
         for i, f in enumerate(files[:5], 1):
             print(f"  [{i}] {os.path.basename(f)}")
         print(f"  ...")
         for i, f in enumerate(files[-5:], len(files)-4):
             print(f"  [{i}] {os.path.basename(f)}")
    else:
        for i, f in enumerate(files, 1):
            print(f"  [{i}] {os.path.basename(f)}")

    while True:
        try:
            prompt = (
                f"\n3. 请输入欲处理的 {desc} 编号\n"
                f" (支持格式: 1,2 | 1-4 | 1,3,5-8 | all): \n"
            )
            user_input = input(prompt).strip().lower()

            if not user_input:
                continue

            selected_indices = set()

            if user_input in ['all', 'a']:
                selected_indices = set(range(len(files)))
            else:
                parts = user_input.split(',')
                for part in parts:
                    part = part.strip()
                    if '-' in part:
                        start_str, end_str = part.split('-')
                        start, end = int(start_str), int(end_str)
                        for idx in range(start - 1, end):
                            if 0 <= idx < len(files):
                                selected_indices.add(idx)
                    else:
                        idx = int(part) - 1
                        if 0 <= idx < len(files):
                            selected_indices.add(idx)
            
            if not selected_indices:
                print("未匹配到有效编号，请重新输入。")
                continue

            selected_paths = [files[i] for i in sorted(list(selected_indices))]
            print(f"\n4. 已选择 {len(selected_paths)} 个文件准备复制。")
            return selected_paths

        except ValueError:
            print("输入错误：请确保输入的是数字编号、范围（如1-5）或 'all'。")
        except KeyboardInterrupt:
            print("\n用户取消操作。")
            sys.exit()

# ============= 确认函数 =============
def make_sure():
    response = input(f"是否确认开始执行复制操作? (y/n): ").strip().lower()
    if response not in ['y', 'yes']:
        print("操作已取消。")
        return False
    return True

# ============= 核心逻辑：复制处理 (已修改) =============
def copy_vcf_gz(file_list, output_dir_name="01_Raw_VCFs"):
    """
    将选中的 .vcf.gz 复制到指定目录 (保持压缩状态)
    """
    base_dir = get_base_dir()
    output_dir = os.path.join(base_dir, output_dir_name)
    
    # 创建输出目录
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
        print(f"\n已创建输出目录: {output_dir}")
    else:
        print(f"\n输出目录已存在: {output_dir}")

    success_count = 0
    total = len(file_list)
    
    print("-" * 50)
    for idx, gz_path in enumerate(file_list, 1):
        try:
            # 获取文件名（保持原文件名，不去除 .gz）
            file_name = os.path.basename(gz_path)
            
            # 目标路径
            out_path = os.path.join(output_dir, file_name)
            
            # 显示进度
            print(f"[{idx}/{total}] 正在复制: {file_name} ...", end="", flush=True)
            
            # 使用 shutil.copy2 进行复制 (保留文件元数据如修改时间)
            shutil.copy2(gz_path, out_path)
            
            print(" 完成")
            success_count += 1
            
        except Exception as e:
            print(f" 失败! \n错误信息: {e}")
    
    print("-" * 50)
    print(f"处理结束。成功复制 {success_count} / {total} 个文件。")
    print(f"文件保存在: {output_dir}")

# ============= 主函数 =============
def main():
    print("--- 基因型文件(.vcf.gz) 批量搜集工具 ---")
    print("该脚本将扫描子目录，将选中的 VCF.GZ 文件复制到统一文件夹 (不解压)。")
    
    # 1. 默认搜索 .vcf.gz
    target_ext = "vcf.gz"
    print(f"\n1. 正在搜索后缀为 .{target_ext} 的文件...")
    
    # 2. 查找文件
    files = find_files(target_ext)
    
    # 3. 交互选择
    chosen_list = choose_file(files, desc=f".{target_ext} 文件")
    
    # 4. 执行复制
    if chosen_list and make_sure():
        # 这里调用新的复制函数
        copy_vcf_gz(chosen_list, output_dir_name="01_Raw_VCFs")

if __name__ == "__main__":
    main()