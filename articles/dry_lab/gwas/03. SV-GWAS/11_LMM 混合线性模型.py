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
# 最终结果存放目录
FINAL_OUTPUT_DIR = "11_LMM_GWAS"

# GEMMA 可执行文件路径
GEMMA_EXEC = "gemma"

# GEMMA 模型参数
# -lmm 1: Wald Test (混合线性模型)
# -lmm 2: Likelihood Ratio Test
# -lmm 4: Score Test
GEMMA_ARGS = ["-lmm", "1"]

# ==============================================================================
# 1. Cite2 交互逻辑模块 (复用版)
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
    # 排除输出目录
    return [f for f in all_files if FINAL_OUTPUT_DIR not in f and "output" not in f]

def choose_file(files, desc="文件"):
    if not files:
        print(f"[提示] 未找到任何 {desc}")
        return []

    print(f"\n--- 找到 {len(files)} 个 {desc}: ---")
    for i, f in enumerate(files, 1):
        print(f"  [{i}] {os.path.basename(f)}  (路径: {os.path.relpath(f, get_base_dir())})")

    while True:
        try:
            prompt = (
                f"\n请输入 {desc} 编号\n"
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
                print(f"已选择 {len(selected_files)} 个文件。")
                return selected_files
            print("选择无效，请重试。")
        except ValueError:
            print("输入错误。")

def make_sure(action_name="执行操作"):
    response = input(f"\n确认{action_name}? (y/n): ").strip().lower()
    return response in ['y', 'yes']

# ==============================================================================
# 2. 核心处理模块
# ==============================================================================
def run_lmm_batch(bfile_path, kinship_path, trait_files):
    base_dir = get_base_dir()
    final_dir = os.path.join(base_dir, FINAL_OUTPUT_DIR)
    gemma_default_out_dir = os.path.join(base_dir, "output")

    if not os.path.exists(gemma_default_out_dir): os.makedirs(gemma_default_out_dir)
    if not os.path.exists(final_dir): os.makedirs(final_dir)

    # PLINK 前缀
    bfile_prefix = os.path.splitext(bfile_path)[0]

    print(f"\n--- 开始批量 LMM 分析 (共 {len(trait_files)} 个性状) ---")
    print(f"基因型: {os.path.basename(bfile_path)}")
    print(f"Kinship: {os.path.basename(kinship_path)}")

    for idx, trait_file in enumerate(trait_files, 1):
        trait_name = os.path.splitext(os.path.basename(trait_file))[0]
        output_name = f"{trait_name}_lmm"
        
        print(f"\n>>> [{idx}/{len(trait_files)}] 正在分析性状: {trait_name}")
        
        # 构造命令
        # gemma -bfile ... -k ... -p ... -lmm 1 -o ...
        cmd = [
            GEMMA_EXEC,
            "-bfile", bfile_prefix,
            "-k", kinship_path,
            "-p", trait_file,
            "-o", output_name
        ] + GEMMA_ARGS
        
        print(f"执行命令: {' '.join(cmd)}")
        
        try:
            # 运行 (capture_output=True 保持界面整洁，只看进度)
            subprocess.run(cmd, check=True)
            
            # 移动结果
            src_assoc = os.path.join(gemma_default_out_dir, f"{output_name}.assoc.txt")
            src_log   = os.path.join(gemma_default_out_dir, f"{output_name}.log.txt")
            
            dst_assoc = os.path.join(final_dir, f"{output_name}.assoc.txt")
            dst_log   = os.path.join(final_dir, f"{output_name}.log.txt")
            
            if os.path.exists(src_assoc):
                shutil.move(src_assoc, dst_assoc)
                print(f"  [成功] 结果已保存: {os.path.basename(dst_assoc)}")
            else:
                print(f"  [警告] 结果文件缺失，GEMMA 可能运行失败。")
                
            if os.path.exists(src_log):
                shutil.move(src_log, dst_log)

        except subprocess.CalledProcessError as e:
            print(f"  [失败] 运行出错: {e}")

    print("\n" + "="*50)
    print(f"批量 LMM 分析完成！")
    print(f"结果存放于: {final_dir}")
    print("="*50)

# ==============================================================================
# 主函数
# ==============================================================================
def main():
    print("==============================================")
    print("   Step 11: 批量 LMM (混合线性模型) 分析工具")
    print("==============================================")

    # 1. 选择基因型文件 (.fam)
    print("\n>>> 第一步: 选择基因型数据 (.fam)")
    fam_files = find_files('.fam')
    sel_fam = choose_file(fam_files, "基因型文件")
    if not sel_fam: return
    bfile_path = sel_fam[0] # 单选

    # 2. 选择 Kinship 矩阵 (.cXX.txt)
    print("\n>>> 第二步: 选择 Kinship 矩阵 (.cXX.txt)")
    # 优先推荐 10_Kinship 文件夹
    kinship_files = find_files('.cXX.txt')
    sel_kin = choose_file(kinship_files, "Kinship 矩阵")
    if not sel_kin: return
    kinship_path = sel_kin[0] # 单选

    # 3. 选择性状文件 (.txt)
    print("\n>>> 第三步: 选择性状文件 (支持多选/全选)")
    # 查找所有 .txt，并过滤掉明显不是性状的文件
    txt_files = find_files('.txt')
    trait_files = [f for f in txt_files if "log" not in f and "kinship" not in f and "cXX" not in f]
    
    sel_traits = choose_file(trait_files, "性状文件")
    if not sel_traits: return

    # 4. 执行
    if make_sure("开始批量 LMM 分析"):
        run_lmm_batch(bfile_path, kinship_path, sel_traits)

if __name__ == "__main__":
    main()