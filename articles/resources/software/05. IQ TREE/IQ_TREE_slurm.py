#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import os
import sys
import glob
import subprocess

# ============= 基础文件选择逻辑 (沿用你的 cite_v2.py) =============
def find_files(exts):
    path = os.getcwd()
    all_files = []
    for ext in exts:
        search_pattern = os.path.join(path, '**', f'*{ext}')
        all_files.extend(glob.glob(search_pattern, recursive=True))
    return sorted(list(set(all_files)))

def choose_file(files, desc="文件"):
    if not files:
        print(f"提示：未找到 {desc}"); return []
    print(f"\n[1] 找到 {len(files)} 个 {desc}:")
    for i, f in enumerate(files, 1):
        print(f"  [{i}] {os.path.relpath(f)}")
    while True:
        try:
            user_input = input(f"\n请输入欲操作的编号 (如: 1,2 | 1-4 | all): ").strip().lower()
            if not user_input: continue
            idx_list = []
            if user_input in ['all', 'a']:
                idx_list = list(range(len(files)))
            else:
                for part in user_input.split(','):
                    if '-' in part:
                        s, e = map(int, part.split('-'))
                        idx_list.extend(range(s-1, e))
                    else:
                        idx_list.append(int(part)-1)
            return [files[i] for i in idx_list if 0 <= i < len(files)]
        except: print("输入格式错误，请重新输入。")

# ============= Slurm 资源参数获取 =============
def get_slurm_config():
    print("\n" + "="*30 + "\n[2] 配置 Slurm 集群资源\n" + "="*30)
    config = {
        "partition": input("  队列/分区名 (Partition, 默认 batch): ") or "batch",
        "cpus": input("  每项任务 CPU 核心数 (默认 16): ") or "16",
        "mem": input("  内存需求 (如 64G, 默认 32G): ") or "32G",
        "time": input("  运行限时 (格式 D-HH:MM:SS, 默认 7-00:00:00): ") or "7-00:00:00",
        "iqtree_cmd": input("  IQ-TREE 参数 (默认 -m MFP -bb 1000): ") or "-m MFP -bb 1000"
    }
    return config

# ============= 生成并提交脚本 =============
def submit_job(msa_path, cfg):
    job_name = os.path.basename(msa_path).split('.')[0]
    sh_name = f"run_{job_name}.sh"
    
    # 构建 Slurm 脚本内容
    slurm_content = f"""#!/bin/bash
#SBATCH --job-name={job_name}
#SBATCH --partition={cfg['partition']}
#SBATCH --nodes=1
#SBATCH --cpus-per-task={cfg['cpus']}
#SBATCH --mem={cfg['mem']}
#SBATCH --time={cfg['time']}
#SBATCH --output=%j_{job_name}.log
#SBATCH --error=%j_{job_name}.err

# 自动加载模块 (如果你的集群需要)
# module load iqtree/2.2.0 

echo "Job started at: `date`"
echo "Running on node: `hostname`"

# 执行 IQ-TREE
# -nt $SLURM_CPUS_PER_TASK 确保 IQ-TREE 使用你申请的核心数
iqtree -s {msa_path} {cfg['iqtree_cmd']} -nt $SLURM_CPUS_PER_TASK

echo "Job finished at: `date`"
"""
    
    # 写入文件
    with open(sh_name, 'w') as f:
        f.write(slurm_content)
    
    # 提交任务
    try:
        res = subprocess.run(['sbatch', sh_name], capture_output=True, text=True, check=True)
        print(f"  成功投递: {sh_name} -> {res.stdout.strip()}")
    except Exception as e:
        print(f"  投递失败: {sh_name}, 错误: {e}")

# ============= 主流程 =============
def main():
    print("--- IQ-TREE Slurm 批量提交工具 ---")
    
    # 1. 选择比对文件
    msas = find_files(['.fasta', '.fa', '.phy'])
    selected_msas = choose_file(msas, "MSA 文件")
    if not selected_msas: return

    # 2. 配置资源
    cfg = get_slurm_config()
    
    # 3. 确认提交
    confirm = input(f"\n确认提交 {len(selected_msas)} 个任务到集群? (y/n): ").lower()
    if confirm in ['y', 'yes']:
        for msa in selected_msas:
            submit_job(msa, cfg)
        print("\n所有任务已尝试提交。你可以使用 'squeue -u $USER' 查看进度。")
    else:
        print("已取消。")

if __name__ == "__main__":
    main()