#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import os
import sys
import glob
import subprocess
import time
import re
from concurrent.futures import ProcessPoolExecutor

# ============= 基础路径与文件查找 =============
def get_base_dir():
    return os.path.dirname(os.path.abspath(__file__))

def find_fasta_files(path=None):
    if path is None: path = get_base_dir()
    extensions = ['*.fasta', '*.fa', '*.fna']
    files = []
    for ext in extensions:
        files.extend(glob.glob(os.path.join(path, '**', ext), recursive=True))
    return sorted(list(set(files)))

# ============= 交互式选择逻辑 (保持你的风格) =============
def choose_files(files):
    if not files:
        print("提示：当前目录下未找到基因组 FASTA 文件。")
        return []
    print(f"\n[1] 找到 {len(files)} 个待预测基因组:")
    for i, f in enumerate(files, 1):
        print(f"  [{i}] {os.path.basename(f)}")
    while True:
        try:
            prompt = "\n[2] 请输入编号选择 (如: 1,2 | 1-3 | all): "
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
            if not selected_indices: continue
            selected_paths = [files[i] for i in sorted(list(selected_indices))]
            return selected_paths
        except Exception as e:
            print(f"输入错误: {e}")

# ============= 蛋白提取函数 (独立出来供并行调用) =============
def extract_proteins(gff_file, faa_file):
    try:
        with open(gff_file, 'r') as f, open(faa_file, 'w') as out:
            current_id = ""
            is_protein = False
            for line in f:
                if line.startswith("# start gene"):
                    current_id = line.strip().split()[-1]
                if line.startswith("# protein sequence = ["):
                    is_protein = True
                    out.write(f">{current_id}\n")
                    seq = re.search(r'\[(.*)\]', line)
                    if seq:
                        out.write(seq.group(1).replace("\n", "").replace("# ", "") + "\n")
                    continue
                if is_protein and line.startswith("# "):
                    seq_part = line.replace("# ", "").strip()
                    if "]" in seq_part:
                        out.write(seq_part.split("]")[0] + "\n")
                        is_protein = False
                    else:
                        out.write(seq_part + "\n")
    except Exception as e:
        print(f"提取蛋白序列出错: {e}")

# ============= 单个任务的封装 =============
def worker(task_info):
    file_path, species, out_dir = task_info
    base_name = os.path.splitext(os.path.basename(file_path))[0]
    gff_out = os.path.join(out_dir, f"{base_name}.gff")
    faa_out = os.path.join(out_dir, f"{base_name}.faa")
    
    cmd = ["augustus", f"--species={species}", "--gff3=on", "--protein=on", file_path]
    
    print(f"开始预测: {base_name} ...")
    try:
        with open(gff_out, 'w') as out_f:
            process = subprocess.run(cmd, stdout=out_f, stderr=subprocess.PIPE, text=True)
        
        if process.returncode == 0:
            extract_proteins(gff_out, faa_out)
            return f"√ {base_name} 完成"
        else:
            return f"× {base_name} 失败: {process.stderr[:100]}"
    except Exception as e:
        return f"× {base_name} 崩溃: {str(e)}"

# ============= 主程序 =============
def main():
    print("=== Augustus 基因预测 过易并行 ===")
    
    files = find_fasta_files()
    selected_files = choose_files(files)
    if not selected_files: return

    # 参数设置
    species = input("\n请输入训练物种 (默认: neurospora_crassa): ").strip() or "neurospora_crassa"
    
    # 【新增核心数选择】
    max_workers = input("请输入并行任务数: ").strip()
    max_workers = int(max_workers) if max_workers else 4
    
    out_dir = "_Augustus_Out"
    if not os.path.exists(out_dir): os.makedirs(out_dir)
    
    # 准备任务列表
    tasks = [(f, species, out_dir) for f in selected_files]
    
    print(f"\n[即将开始] 并行核心数: {max_workers} | 任务总数: {len(tasks)}")
    confirm = input("确认执行? (y/n): ").strip().lower()
    
    if confirm in ['y', 'yes']:
        start_time = time.time()
        
        # 使用进程池执行
        with ProcessPoolExecutor(max_workers=max_workers) as executor:
            results = list(executor.map(worker, tasks))
        
        for res in results:
            print(res)
            
        end_time = time.time()
        print(f"\n全部任务完成！总耗时: {(end_time - start_time)/60:.2f} 分钟。")

if __name__ == "__main__":
    main()