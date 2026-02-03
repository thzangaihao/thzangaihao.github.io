#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import sys
import glob
import subprocess
import shutil

# ==============================================================================
# 0. 核心配置参数
# ==============================================================================
# 输出目录名称
OUTPUT_DIR_NAME = "03_QC_Filtered"

# 过滤条件 (使用 bcftools 表达式)
# 保留 (FILTER列为PASS) 并且 (次等位基因频率 > 0.05) 的位点
# FILTER_EXPRESSION = 'FILTER="PASS" && MAF > 0.05'
FILTER_EXPRESSION = 'MAF > 0.05'

# 线程数
THREADS = 16

# 输出格式: 'z' (vcf.gz), 'b' (bcf), 'v' (vcf)
OUTPUT_FORMAT = "z" 

# 可执行程序路径
BCFTOOLS_EXEC = "bcftools"

# ==============================================================================
# 1. Cite2 交互逻辑模块 (内置)
# ==============================================================================
def get_base_dir():
    return os.path.dirname(os.path.abspath(__file__))

def find_files(exts, path=None):
    if path is None: path = get_base_dir()
    all_found = []
    for ext in exts:
        if not ext.startswith('.'): ext = '.' + ext
        search_pattern = os.path.join(path, '**', f'*{ext}')
        # 排除之前的输出目录，防止循环处理
        files = glob.glob(search_pattern, recursive=True)
        files = [f for f in files if OUTPUT_DIR_NAME not in f and "04_" not in f]
        all_found.extend(files)
    return sorted(list(set(all_found)))

def choose_file(files, desc="文件"):
    if not files:
        print(f"[提示] 未找到任何 {desc}")
        return []

    print(f"\n2. 找到 {len(files)} 个 {desc}:")
    limit = 10
    if len(files) > limit * 2:
        for i, f in enumerate(files[:limit], 1):
            print(f"  [{i}] {os.path.basename(f)} ({os.path.relpath(f, get_base_dir())})")
        print(f"  ... (中间省略 {len(files) - limit*2} 个) ...")
        for i, f in enumerate(files[-limit:], len(files)-limit+1):
            print(f"  [{i}] {os.path.basename(f)} ({os.path.relpath(f, get_base_dir())})")
    else:
        for i, f in enumerate(files, 1):
            print(f"  [{i}] {os.path.basename(f)} ({os.path.relpath(f, get_base_dir())})")

    while True:
        try:
            prompt = f"\n3. 请输入编号 (格式: 1,2 | 1-4 | all): "
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
            print("输入格式错误。")

def make_sure():
    response = input(f"\n5. 确认执行过滤操作? (y/n): ").strip().lower()
    return response in ['y', 'yes']

# ==============================================================================
# 2. 核心处理模块
# ==============================================================================
def run_filtering(file_list):
    # 创建输出目录
    base_dir = get_base_dir()
    out_dir = os.path.join(base_dir, OUTPUT_DIR_NAME)
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
        print(f"\n[系统] 创建输出目录: {out_dir}")

    for idx, input_file in enumerate(file_list, 1):
        filename = os.path.basename(input_file)
        # 构造输出文件名: OriginalName_filtered.vcf.gz
        prefix = filename.split('.')[0] # 简单取前缀
        if OUTPUT_FORMAT == 'z':
            suffix = ".vcf.gz"
        elif OUTPUT_FORMAT == 'b':
            suffix = ".bcf"
        else:
            suffix = ".vcf"
            
        output_filename = f"{prefix}_filtered{suffix}"
        output_path = os.path.join(out_dir, output_filename)
        
        print(f"\n--- [{idx}/{len(file_list)}] 正在处理: {filename} ---")
        print(f"过滤条件: {FILTER_EXPRESSION}")
        
        cmd = [
            BCFTOOLS_EXEC, "view",
            "-i", FILTER_EXPRESSION,
            input_file,
            "-O", OUTPUT_FORMAT,
            "-o", output_path,
            "--threads", str(THREADS)
        ]
        
        print(f"执行命令: {' '.join(cmd)}")
        try:
            subprocess.run(cmd, check=True)
            print(f"[成功] 输出至: {output_path}")
            
            # 顺手建立索引，方便下一步使用
            if OUTPUT_FORMAT in ['z', 'b']:
                print("正在建立索引...")
                subprocess.run([BCFTOOLS_EXEC, "index", "-f", output_path], check=False)
                
        except subprocess.CalledProcessError as e:
            print(f"[失败] 处理文件 {filename} 时出错: {e}")

# ==============================================================================
# 主函数
# ==============================================================================
def main():
    print("==============================================")
    print("   Step 3: VCF/BCF 质控与过滤工具 ")
    print("==============================================")
    
    # 1. 查找
    files = find_files(['vcf.gz', 'bcf'])
    
    # 2. 选择
    selected = choose_file(files, "待过滤文件 (建议选择 02 文件夹中的合并结果)")
    if not selected: return

    # 3. 确认与执行
    if make_sure():
        run_filtering(selected)

if __name__ == "__main__":
    main()