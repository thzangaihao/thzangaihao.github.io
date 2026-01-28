#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import sys
import glob
import subprocess

# ==============================================================================
# 0. 核心配置参数
# ==============================================================================
# 输出目录名称
OUTPUT_DIR_NAME = "04_PLINK_Binary"

# PLINK 额外参数
# --allow-extra-chr: 必须开启，因为大多数生物染色体命名通常不是 1-22
# --const-fid: 如果 VCF 没有 Family ID，自动将 Family ID 设为 0
EXTRA_FLAGS = ["--allow-extra-chr", "--const-fid"]

# 可执行程序路径
PLINK_EXEC = "plink"  # 如果需要 PLINK2，请改为 "plink2" 并相应调整参数

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
        files = glob.glob(search_pattern, recursive=True)
        # 排除输出目录，但为了方便，优先推荐 03 目录下的文件
        files = [f for f in files if OUTPUT_DIR_NAME not in f]
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
    response = input(f"\n5. 确认执行转换操作? (y/n): ").strip().lower()
    return response in ['y', 'yes']

# ==============================================================================
# 2. 核心处理模块
# ==============================================================================
def run_conversion(file_list):
    base_dir = get_base_dir()
    out_dir = os.path.join(base_dir, OUTPUT_DIR_NAME)
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
        print(f"\n[系统] 创建输出目录: {out_dir}")

    for idx, input_file in enumerate(file_list, 1):
        filename = os.path.basename(input_file)
        # 去除扩展名作为前缀 (处理 .vcf.gz 和 .bcf)
        if filename.endswith(".vcf.gz"):
            prefix = filename[:-7]
        elif filename.endswith(".bcf"):
            prefix = filename[:-4]
        else:
            prefix = os.path.splitext(filename)[0]
            
        # 结果文件前缀
        out_prefix = os.path.join(out_dir, prefix + "_plink")
        
        print(f"\n--- [{idx}/{len(file_list)}] 正在转换: {filename} ---")
        
        # 判断输入格式
        input_arg = "--bcf" if input_file.endswith(".bcf") else "--vcf"
        
        cmd = [
            PLINK_EXEC,
            input_arg, input_file,
            "--make-bed",
            "--out", out_prefix
        ] + EXTRA_FLAGS
        
        print(f"执行命令: {' '.join(cmd)}")
        try:
            # PLINK 的日志通常输出到标准输出，这里我们让它显示出来以便查错
            subprocess.run(cmd, check=True)
            print(f"[成功] 生成: {out_prefix}.bed / .bim / .fam")
        except subprocess.CalledProcessError as e:
            print(f"[失败] 转换 {filename} 时出错。请检查文件格式或 PLINK 版本。")

# ==============================================================================
# 主函数
# ==============================================================================
def main():
    print("==============================================")
    print("   Step 4: 格式转换工具 (VCF -> PLINK Binary)")
    print("==============================================")
    
    # 1. 查找
    # 优先找过滤后的文件，通常是 .vcf.gz 或 .bcf
    files = find_files(['vcf.gz', 'bcf'])
    
    # 2. 选择
    selected = choose_file(files, "待转换文件 (VCF/BCF)")
    if not selected: return

    # 3. 确认与执行
    if make_sure():
        run_conversion(selected)

if __name__ == "__main__":
    main()