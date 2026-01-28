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
FINAL_OUTPUT_DIR = "07_GWAS_Results"

# GEMMA 可执行文件路径
GEMMA_EXEC = "gemma"

# GEMMA 模型参数
# -lm 1: Wald Test (线性模型)，速度快，适合初筛
# -lmm 1: 如果你有 Kinship 矩阵并想跑混合线性模型，请改为 "-lmm 1" 并添加 -k 参数逻辑
GEMMA_ARGS = ["-lm", "1"]

# ==============================================================================
# 1. Cite2 交互逻辑模块 (复用)
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
    # 排除之前的输出目录
    return [f for f in all_files if FINAL_OUTPUT_DIR not in f and "output" not in f]

def choose_files(files, desc="文件"):
    if not files:
        print(f"提示：未找到任何 {desc}")
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
                f"\n3. 请输入欲分析的 {desc} 编号\n"
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
                print(f"\n4. 已选择 {len(selected_files)} 个文件待分析。")
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
def run_gemma_batch(bfile_path, trait_files):
    # 1. 准备目录
    base_dir = get_base_dir()
    final_dir = os.path.join(base_dir, FINAL_OUTPUT_DIR)
    
    # GEMMA 强制输出到当前目录下的 'output' 文件夹
    # 我们需要确保这个文件夹存在，否则 GEMMA 会报错
    gemma_default_out_dir = os.path.join(base_dir, "output")
    if not os.path.exists(gemma_default_out_dir):
        os.makedirs(gemma_default_out_dir)

    if not os.path.exists(final_dir):
        os.makedirs(final_dir)
        print(f"\n[系统] 创建最终结果目录: {final_dir}")

    # PLINK 前缀 (去除 .fam后缀)
    bfile_prefix = bfile_path.replace(".fam", "").replace(".bed", "").replace(".bim", "")

    print(f"\n--- 开始批量 GWAS 分析 (共 {len(trait_files)} 个任务) ---")

    for idx, trait_file in enumerate(trait_files, 1):
        # 获取性状名称 (文件名无后缀)
        trait_name = os.path.splitext(os.path.basename(trait_file))[0]
        
        print(f"\n>>> [{idx}/{len(trait_files)}] 正在分析性状: {trait_name}")
        
        # 构造命令
        # gemma -bfile <plink> -p <trait> -lm 1 -o <trait_name>
        cmd = [
            GEMMA_EXEC,
            "-bfile", bfile_prefix,
            "-p", trait_file,
            "-o", trait_name  # 输出文件前缀 (会自动保存到 output/trait_name.assoc.txt)
        ] + GEMMA_ARGS
        
        print(f"执行命令: {' '.join(cmd)}")
        
        try:
            # 运行 GEMMA
            # capture_output=True 可以避免屏幕被 GEMMA 的刷屏日志淹没，只打印错误
            # 如果想看实时日志，去掉 capture_output=True
            subprocess.run(cmd, check=True)
            
            # 移动结果文件
            # GEMMA 输出通常是: output/trait_name.assoc.txt 和 output/trait_name.log.txt
            src_assoc = os.path.join(gemma_default_out_dir, f"{trait_name}.assoc.txt")
            src_log   = os.path.join(gemma_default_out_dir, f"{trait_name}.log.txt")
            
            dst_assoc = os.path.join(final_dir, f"{trait_name}.assoc.txt")
            dst_log   = os.path.join(final_dir, f"{trait_name}.log.txt")
            
            if os.path.exists(src_assoc):
                shutil.move(src_assoc, dst_assoc)
                print(f"  [成功] 结果已保存: {os.path.basename(dst_assoc)}")
            else:
                print(f"  [警告] 未找到结果文件，GEMMA 可能运行失败。")
                
            if os.path.exists(src_log):
                shutil.move(src_log, dst_log)

        except subprocess.CalledProcessError as e:
            print(f"  [失败] 运行出错: {e}")

    # 善后：如果 output 文件夹空了，可以删掉，保持整洁 (可选)
    # try:
    #     os.rmdir(gemma_default_out_dir)
    # except:
    #     pass

    print("\n" + "="*50)
    print(f"批量分析完成！")
    print(f"所有结果文件已存放于: {final_dir}")
    print("="*50)

# ==============================================================================
# 主函数
# ==============================================================================
def main():
    print("==============================================")
    print("   Step 7: 批量 GWAS 分析运行器 (GEMMA)")
    print("==============================================")

    # 1. 选择 PLINK 文件 (基因型)
    print("\n>>> 第一步: 选择基因型数据 (.fam/.bed)")
    # 搜索 .fam 最准确，因为 bed/bim 肯定和它在一起
    fam_files = find_files('.fam')
    if not fam_files:
        print("未找到 .fam 文件，请先完成 Step 4。")
        return
    
    # 这里我们只选一个基准基因型文件
    print(f"\n2. 找到 {len(fam_files)} 个基因型文件:")
    for i, f in enumerate(fam_files, 1):
        print(f"  [{i}] {os.path.basename(f)}")
        
    try:
        idx = int(input("\n3. 请输入基因型文件编号 (单选): ").strip()) - 1
        bfile_path = fam_files[idx]
    except:
        print("输入无效。")
        return

    # 2. 选择性状文件 (支持多选/全选)
    print("\n>>> 第二步: 选择要分析的性状文件 (.txt)")
    # 优先推荐 Step 6 的输出目录，但这里通用搜索 txt
    trait_files = find_files('.txt')
    # 排除掉 log 文件或其他无关文本
    trait_files = [f for f in trait_files if "gemma_input" not in f and "log" not in f]
    
    selected_traits = choose_files(trait_files, "性状文件")
    if not selected_traits: return

    # 3. 执行
    if make_sure("开始批量运行"):
        run_gemma_batch(bfile_path, selected_traits)

if __name__ == "__main__":
    main()