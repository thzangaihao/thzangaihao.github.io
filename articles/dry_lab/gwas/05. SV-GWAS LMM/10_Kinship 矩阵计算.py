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
OUTPUT_DIR_NAME = "10_Kinship"
GEMMA_EXEC = "gemma"
GEMMA_ARGS = ["-gk", "1"] # -gk 1: Centered relatedness matrix

# ==============================================================================
# 1. Cite2 交互逻辑模块
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
    return [f for f in all_files if OUTPUT_DIR_NAME not in f and "output" not in f]

def choose_file(files, desc="文件"):
    if not files:
        print(f"[提示] 未找到任何 {desc}")
        return []

    print(f"\n2. 找到 {len(files)} 个 {desc}:")
    for i, f in enumerate(files, 1):
        print(f"  [{i}] {os.path.basename(f)}  (路径: {os.path.relpath(f, get_base_dir())})")

    while True:
        try:
            prompt = (
                f"\n3. 请输入欲操作的 {desc} 编号\n"
                f" (支持格式: 1,2 | 1-4 | all): \n"
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
                        for idx in range(s - 1, e):
                            selected_indices.add(idx)
                    else:
                        selected_indices.add(int(part)-1)
            
            selected_files = [files[i] for i in sorted(selected_indices) if 0 <= i < len(files)]
            if selected_files:
                print(f"\n4. 已选择 {len(selected_files)} 个文件。")
                return selected_files
        except ValueError:
            print("输入错误。")

def make_sure(action_name="执行操作"):
    response = input(f"\n5. 确认{action_name}? (y/n): ").strip().lower()
    return response in ['y', 'yes']

# ==============================================================================
# 2. 核心处理模块 (含自动修复功能)
# ==============================================================================
def patch_fam_phenotype(fam_path):
    """
    检查 .fam 文件的表型列（第6列）。
    如果是 -9 (缺失)，则将其修改为 1 (Control)，防止 GEMMA 剔除所有样本。
    """
    print(f"正在检查 .fam 文件格式: {os.path.basename(fam_path)} ...")
    
    needs_patch = False
    with open(fam_path, 'r') as f:
        first_line = f.readline()
        if not first_line: return # 空文件
        parts = first_line.split()
        if len(parts) >= 6 and parts[5] == '-9':
            needs_patch = True
            
    if needs_patch:
        print("  [警告] 检测到表型列全为 -9 (缺失)。")
        print("  [自动修复] 正在将表型重置为 1，以允许 GEMMA 进行计算...")
        
        # 备份原文件
        bak_path = fam_path + ".bak"
        if not os.path.exists(bak_path):
            shutil.copy(fam_path, bak_path)
            print(f"  [备份] 原文件已备份至: {os.path.basename(bak_path)}")
        
        # 读取并修改
        with open(fam_path, 'r') as f:
            lines = f.readlines()
            
        with open(fam_path, 'w') as f:
            for line in lines:
                parts = line.strip().split()
                if len(parts) >= 6:
                    parts[5] = '1' # 强制设为 1
                f.write('\t'.join(parts) + '\n')
        print("  [成功] .fam 文件修复完成。")
    else:
        print("  [状态] .fam 文件格式正常。")

def run_kinship(bfile_path):
    # 1. 自动修复 .fam 文件
    patch_fam_phenotype(bfile_path)

    # 2. 准备目录
    base_dir = get_base_dir()
    final_dir = os.path.join(base_dir, OUTPUT_DIR_NAME)
    gemma_default_out_dir = os.path.join(base_dir, "output")

    if not os.path.exists(gemma_default_out_dir): os.makedirs(gemma_default_out_dir)
    if not os.path.exists(final_dir): os.makedirs(final_dir)

    prefix = os.path.basename(bfile_path).replace(".fam", "")
    output_name = f"{prefix}_kinship"
    bfile_arg = os.path.splitext(bfile_path)[0] # 去掉 .fam

    print(f"\n--- 开始计算 Kinship 矩阵 ---")
    
    cmd = [
        GEMMA_EXEC,
        "-bfile", bfile_arg,
        "-o", output_name
    ] + GEMMA_ARGS
    
    print(f"执行命令: {' '.join(cmd)}")
    
    try:
        subprocess.run(cmd, check=True)
        
        src_cxx = os.path.join(gemma_default_out_dir, f"{output_name}.cXX.txt")
        src_log = os.path.join(gemma_default_out_dir, f"{output_name}.log.txt")
        
        dst_cxx = os.path.join(final_dir, f"{output_name}.cXX.txt")
        dst_log = os.path.join(final_dir, f"{output_name}.log.txt")
        
        if os.path.exists(src_cxx):
            shutil.move(src_cxx, dst_cxx)
            print(f"\n[成功] Kinship 矩阵已保存至: {dst_cxx}")
        else:
            print("\n[警告] 未找到 .cXX.txt 输出，GEMMA 可能运行失败。")
            
        if os.path.exists(src_log):
            shutil.move(src_log, dst_log)

    except subprocess.CalledProcessError as e:
        print(f"\n[失败] 运行出错: {e}")

# ==============================================================================
# 主函数
# ==============================================================================
def main():
    print("==============================================")
    print("   Step 10: 计算亲缘关系矩阵 (Kinship) - 修复版")
    print("==============================================")
    
    fam_files = find_files('.fam')
    selected = choose_file(fam_files, "PLINK 基因型文件 (.fam)")
    
    if not selected: return
    
    if make_sure("开始计算"):
        run_kinship(selected[0])

if __name__ == "__main__":
    main()