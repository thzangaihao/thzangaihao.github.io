#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import os
import sys
import glob
import subprocess
import time

# ============= 基础路径与文件查找 =============
def get_base_dir():
    return os.path.dirname(os.path.abspath(__file__))

def find_fasta_files(path=None):
    """查找常见的基因组后缀文件"""
    if path is None:
        path = get_base_dir()
    extensions = ['*.fasta', '*.fa', '*.fna', '*.gz']
    files = []
    for ext in extensions:
        files.extend(glob.glob(os.path.join(path, '**', ext), recursive=True))
    return sorted(list(set(files)))

# ============= 交互式选择逻辑 (继承自你的 cite_v2.py) =============
def choose_files(files):
    if not files:
        print("提示：当前目录下未找到 .fasta 或 .fa 文件。")
        return []

    print(f"\n[1] 找到 {len(files)} 个待测基因组文件:")
    for i, f in enumerate(files, 1):
        print(f"  [{i}] {os.path.basename(f)}")

    while True:
        try:
            prompt = "\n[2] 请输入编号进行选择 (如: 1,2 | 1-3 | all): "
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
                print("输入无效，请重新选择。")
                continue
                
            selected_paths = [files[i] for i in sorted(list(selected_indices))]
            print(f"\n已选择: {[os.path.basename(p) for p in selected_paths]}")
            return selected_paths
        except Exception as e:
            print(f"输入错误: {e}")

# ============= BUSCO 参数设置 =============
def get_busco_config():
    print("\n" + "="*30)
    print("      BUSCO 参数设置")
    print("="*30)
    
    # 默认针对木霉菌推荐 hypocreales_odb10
    lineage = input("请输入 Lineage 数据库 (默认: hypocreales_odb10): ").strip()
    if not lineage: lineage = "hypocreales_odb10"
    
    cpu = input("请输入使用的 CPU 核心数 (默认: 8): ").strip()
    if not cpu: cpu = "8"
    
    mode = input("请输入运行模式 (genome/proteins, 默认: genome): ").strip()
    if not mode: mode = "genome"
    
    return lineage, cpu, mode

# ============= 执行核心 =============
def run_busco(file_path, lineage, cpu, mode):
    file_name = os.path.basename(file_path)
    # 自动生成输出文件夹名，例如 BUSCO_CK-3-5
    out_name = f"BUSCO_{os.path.splitext(file_name)[0]}"
    
    cmd = [
        "busco",
        "-i", file_path,
        "-o", out_name,
        "-l", lineage,
        "-m", mode,
        "-c", cpu#,
        #"--offline"  # 如果你已下好数据库，建议开启；若未下载，请去掉此行
    ]
    
    print(f"\n>>> 正在处理: {file_name}")
    print(f"执行命令: {' '.join(cmd)}")
    
    try:
        # 使用 subprocess 运行并实时打印输出
        process = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
        for line in process.stdout:
            print(line, end="")
        process.wait()
        
        if process.returncode == 0:
            print(f"√ {file_name} 分析完成。")
        else:
            print(f"× {file_name} 分析失败，请检查输出日志。")
            
    except FileNotFoundError:
        print("错误：未在系统 PATH 中找到 'busco' 命令，请确认 conda 环境已激活。")
        sys.exit(1)

# ============= 主程序 =============
def main():
    print("=== BUSCO 自动化分析交互脚本 ===")
    
    # 1. 查找并选择文件
    all_files = find_fasta_files()
    selected_files = choose_files(all_files)
    
    if not selected_files:
        print("未选择文件，退出。")
        return

    # 2. 设置参数
    lineage, cpu, mode = get_busco_config()
    
    # 3. 最后的确认
    confirm = input(f"\n准备对 {len(selected_files)} 个文件进行分析，确认开始? (y/n): ").strip().lower()
    if confirm not in ['y', 'yes']:
        print("已取消。")
        return

    start_time = time.time()
    for f in selected_files:
        run_busco(f, lineage, cpu, mode)
    
    end_time = time.time()
    duration = (end_time - start_time) / 60
    print(f"\n所有任务已完成！总耗时: {duration:.2f} 分钟。")

if __name__ == "__main__":
    main()