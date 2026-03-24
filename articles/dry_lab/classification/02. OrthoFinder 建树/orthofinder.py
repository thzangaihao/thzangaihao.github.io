#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import os
import sys
import glob
import shutil
import subprocess
import time

def get_base_dir():
    return os.path.dirname(os.path.abspath(__file__))

def find_protein_files(path=None):
    if path is None: path = get_base_dir()
    extensions = ['*.faa', '*.pep', '*.fasta']
    files = []
    for ext in extensions:
        files.extend(glob.glob(os.path.join(path, '**', ext), recursive=True))
    return sorted(list(set(files)))

def choose_files(files):
    if not files:
        print("提示：未找到序列文件。")
        return []
    print(f"\n[1] 找到 {len(files)} 个序列文件:")
    for i, f in enumerate(files, 1):
        print(f"  [{i}] {os.path.basename(f)}")
    while True:
        try:
            user_input = input("\n[2] 请输入编号 (如: 1,2 | 1-5 | all): ").strip().lower()
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
            return [files[i] for i in sorted(list(selected_indices))]
        except: print("输入错误。")

def prepare_input_dir(selected_files, project_name):
    input_dir = os.path.join(get_base_dir(), f"{project_name}_Input")
    if os.path.exists(input_dir): shutil.rmtree(input_dir)
    os.makedirs(input_dir)
    
    print(f"\n>>> 正在导入并预处理序列 (清洗非法字符)...")
    for f in selected_files:
        dest_path = os.path.join(input_dir, os.path.basename(f))
        # 【核心修复】：读取并清洗序列中的 * 号，确保 Diamond 不报错
        with open(f, 'r') as fin, open(dest_path, 'w') as fout:
            for line in fin:
                if line.startswith(">"):
                    fout.write(line)
                else:
                    fout.write(line.replace("*", "").strip() + "\n")
        print(f"  - 已导入: {os.path.basename(f)}")
    return input_dir

def run_orthofinder(input_dir, cpu, use_msa):
    cmd = ["orthofinder", "-f", input_dir, "-t", str(cpu), "-a", str(cpu)]
    if use_msa: cmd.extend(["-M", "msa"])
    
    print("\n🚀 正在启动 OrthoFinder ...\n" + "="*40)
    try:
        # 使用 check=True，如果报错会直接抛出异常
        process = subprocess.run(cmd, text=True)
        if process.returncode == 0:
            print("\n" + "="*40 + "\n🎉 OrthoFinder 分析圆满完成！")
        else:
            print("\n" + "="*40 + "\n❌ OrthoFinder 运行中途夭折，请检查上方日志。")
            sys.exit(1)
    except Exception as e:
        print(f"\n❌ 脚本执行崩溃: {e}")
        sys.exit(1)

def main():
    print("=== OrthoFinder 自动化建树脚本 (修复版) ===")
    all_files = find_protein_files()
    selected_files = choose_files(all_files)
    if not selected_files: return
    
    project_name = input("项目名称: ").strip() or "Ortho_Project"
    cpu = input("CPU核心数 (默认 4): ").strip() or "4"
    
    input_dir = prepare_input_dir(selected_files, project_name)
    run_orthofinder(input_dir, cpu, True)

if __name__ == "__main__":
    main()