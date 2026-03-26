#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import os
import sys
import glob
import subprocess
import time

# ============= 基础路径获取 =============
def get_base_dir():
    if getattr(sys, 'frozen', False):
        return os.path.dirname(sys.executable)
    return os.path.dirname(os.path.abspath(__file__))

# ============= 文件查找 =============
def find_files(exts, path=None):
    if path is None:
        path = get_base_dir()
    all_files = []
    for ext in exts:
        if not ext.startswith('.'):
            ext = '.' + ext
        search_pattern = os.path.join(path, '**', f'*{ext}')
        all_files.extend(glob.glob(search_pattern, recursive=True))
    return sorted(list(set(all_files)))

# ============= 核心选择逻辑（沿用你的代码） =============
def choose_file(files, desc="文件"):
    if not files:
        print(f"提示：未找到任何 {desc}")
        return []

    print(f"\n找到 {len(files)} 个 {desc}:")
    for i, f in enumerate(files, 1):
        print(f"  [{i}] {os.path.relpath(f)}")

    while True:
        try:
            prompt = f"\n请输入欲操作的 {desc} 编号 (如: 1 | 1,2 | 1-4 | all): "
            user_input = input(prompt).strip().lower()
            if not user_input: continue
            
            selected_indices = set()
            if user_input in ['all', 'a']:
                selected_indices = set(range(len(files)))
            else:
                for part in user_input.split(','):
                    part = part.strip()
                    if '-' in part:
                        start, end = map(int, part.split('-'))
                        for idx in range(start - 1, end):
                            if 0 <= idx < len(files): selected_indices.add(idx)
                    else:
                        idx = int(part) - 1
                        if 0 <= idx < len(files): selected_indices.add(idx)
            
            if not selected_indices:
                print("未匹配到有效编号。")
                continue
            return [files[i] for i in sorted(list(selected_indices))]
        except ValueError:
            print("输入错误：请输入数字、范围或 'all'。")

# ============= IQ-TREE 参数配置 =============
def get_iqtree_params():
    print("\n--- IQ-TREE 参数配置 ---")
    model = input("1. 进化模型 (默认 MFP): ").strip() or "MFP"
    bootstrap = input("2. UFBoot 自展值 (默认 1000): ").strip() or "1000"
    threads = input("3. 线程数 (默认 AUTO): ").strip() or "AUTO"
    
    use_partition = input("4. 是否使用分区文件 (.nex/.txt)? (y/n, 默认 n): ").strip().lower()
    partition_file = None
    if use_partition == 'y':
        p_files = find_files(['.txt', '.nex', '.partition'])
        if p_files:
            partition_file = choose_file(p_files, "分区文件")[0]
        else:
            print("未在当前目录找到分区文件，将跳过分区设置。")
            
    return model, bootstrap, threads, partition_file

# ============= 执行函数 =============
def run_iqtree(fasta_path, model, bb, nt, part_file):
    # 构建基础命令
    cmd = [
        "iqtree", 
        "-s", fasta_path, 
        "-m", model, 
        "-bb", bb, 
        "-nt", nt
    ]
    
    # 如果有分区文件，使用 -spp 参数（允许分区速率缩放）
    if part_file:
        cmd.extend(["-spp", part_file])
    
    print(f"\n[执行命令]: {' '.join(cmd)}")
    try:
        # 实时输出 IQ-TREE 的日志
        subprocess.run(cmd, check=True)
    except FileNotFoundError:
        print("错误：未找到 iqtree 程序，请确保它已添加到系统环境变量 (PATH) 中。")
        sys.exit(1)
    except subprocess.CalledProcessError as e:
        print(f"运行出错：{e}")

# ============= 主函数 =============
def main():
    print("="*40)
    print("  IQ-TREE 交互式建树工具 (基于拼接蛋白序列)")
    print("="*40)
    
    # 1. 查找 MSA 文件
    msa_files = find_files(['.fasta', '.fa', '.phy', '.aln'])
    # selected_msas = choose_file(msa_files, "MSA 比对文件")
    selected_msas = choose_file(msa_files, "序列文件")
    
    if not selected_msas:
        return

    # 2. 获取参数
    model, bb, nt, part_file = get_iqtree_params()
    
    # 3. 确认并执行
    print(f"\n选中的文件将依次使用模型 {model} 进行 {bb} 次自展检验。")
    confirm = input("确认开始运行? (y/n): ").strip().lower()
    
    if confirm in ['y', 'yes']:
        for i, msa in enumerate(selected_msas, 1):
            print(f"\n>>> 正在处理第 {i}/{len(selected_msas)} 个文件: {os.path.basename(msa)}")
            run_iqtree(msa, model, bb, nt, part_file)
        print("\n所有任务已完成！")
    else:
        print("操作已取消。")

if __name__ == "__main__":
    main()