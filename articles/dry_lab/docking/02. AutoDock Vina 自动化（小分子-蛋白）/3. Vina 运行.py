#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# 正式对接，不知道vina为什么只能跑在一个核心上，即使按照vina官方manual写也不行，因此采用过易并行策略 #
import os
import sys
import glob
import subprocess
import time
import shutil
from multiprocessing import Pool, cpu_count

# ============= 1. 基础配置 =============
MAX_PARALLEL_TASKS = 4  # 同时运行的最大 Vina 进程数，如果出问题那就由过易并行兜底
CPU_PER_VINA = 16         # 每个 Vina 进程占用的核心数 (建议设为 1 以实现最大并行)

def get_base_dir():
    if getattr(sys, 'frozen', False):
        return os.path.dirname(sys.executable)
    return os.path.dirname(os.path.abspath(__file__))

def find_files_recursive(directory, ext):
    if not ext.startswith('.'): ext = '.' + ext
    search_pattern = os.path.join(directory, '**', f'*{ext}')
    return sorted(glob.glob(search_pattern, recursive=True))

# ============= 2. 单个对接任务函数 (Worker) =============
def run_single_docking(task_info):
    """
    该函数会被多个进程同时调用
    task_info: (lig_path, rec_path, conf_file, out_pdbqt, log_file, lig_base, rec_base)
    """
    lig_path, rec_path, conf_file, out_pdbqt, log_file, lig_base, rec_base = task_info
    
    # 构建 Vina 命令
    # 强制设置 --cpu 为指定的 CPU_PER_VINA
    cmd = [
        'vina',
        '--config', conf_file,
        '--ligand', lig_path,
        '--out', out_pdbqt,
        '--log', log_file,
        '--cpu', str(CPU_PER_VINA)
    ]

    try:
        # 执行对接
        subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.PIPE)
        return True, f"完成: {lig_base} vs {rec_base}"
    except subprocess.CalledProcessError as e:
        return False, f"出错: {lig_base} vs {rec_base} | 错误信息: {e.stderr.decode()[:100]}"
    except Exception as e:
        return False, f"异常: {lig_base} vs {rec_base} | {str(e)}"

# ============= 3. 主控制程序 =============
def main():
    print("=== Step 3 (Parallel): 多进程批量对接程序 ===")
    
    # 环境检查
    if not shutil.which("vina"):
        print("错误: 未找到 'vina' 命令。请确保 conda 环境已激活。")
        return

    base_path = get_base_dir()
    lig_dir = os.path.join(base_path, "1. Ligand_molecular_pdbqt")
    rec_dir = os.path.join(base_path, "1. Receptor_protein_pdbqt")
    conf_dir = os.path.join(base_path, "2. Docking_Config")
    result_base_dir = os.path.join(base_path, "3. Docking_Results")

    # 获取输入列表
    ligands = find_files_recursive(lig_dir, '.pdbqt')
    receptors = find_files_recursive(rec_dir, '.pdbqt')

    if not ligands or not receptors:
        print("错误: 未找到配体或受体文件。")
        return

    # 准备任务列表
    all_tasks = []
    for rec_path in receptors:
        rec_name = os.path.basename(rec_path)
        rec_base = os.path.splitext(rec_name)[0]
        conf_file = os.path.join(conf_dir, f"conf_{rec_base}.txt")
        
        if not os.path.exists(conf_file):
            continue

        # 创建受体结果文件夹
        rec_result_dir = os.path.join(result_base_dir, rec_base)
        os.makedirs(rec_result_dir, exist_ok=True)

        for lig_path in ligands:
            lig_base = os.path.splitext(os.path.basename(lig_path))[0]
            out_pdbqt = os.path.join(rec_result_dir, f"{lig_base}_vs_{rec_base}.pdbqt")
            log_file = os.path.join(rec_result_dir, f"{lig_base}_vs_{rec_base}.log")
            
            # 跳过已存在的结果 (可选)
            if os.path.exists(out_pdbqt) and os.path.exists(log_file):
                continue

            all_tasks.append((lig_path, rec_path, conf_file, out_pdbqt, log_file, lig_base, rec_base))

    total_tasks = len(all_tasks)
    if total_tasks == 0:
        print("所有对接任务均已完成，无需重复运行。")
        return

    print(f"待处理总任务数: {total_tasks}")
    print(f"并行策略: 同时启动 {MAX_PARALLEL_TASKS} 个进程，每个进程使用 {CPU_PER_VINA} 核。")
    
    if input("是否开始大规模并行对接? (y/n): ").lower() not in ['y', 'yes']:
        return

    # 启动进程池
    start_time = time.time()
    success_count = 0
    
    with Pool(processes=MAX_PARALLEL_TASKS) as pool:
        # 使用 imap_unordered 可以实时获取完成的结果
        for i, (success, msg) in enumerate(pool.imap_unordered(run_single_docking, all_tasks), 1):
            if success:
                success_count += 1
            
            # 每完成 1% 或每 10 个打印一次进度
            if i % 10 == 0 or i == total_tasks:
                elapsed = time.time() - start_time
                avg_time = elapsed / i
                eta = avg_time * (total_tasks - i)
                sys.stdout.write(f"\r  进度: {i}/{total_tasks} | 成功: {success_count} | 预计剩余时间: {eta/60:.1f} 分钟")
                sys.stdout.flush()

    end_time = time.time()
    print(f"\n\n" + "="*50)
    print(f"任务结束！成功完成: {success_count}/{total_tasks}")
    print(f"总耗时: {(end_time - start_time)/306:.2f} 小时")
    print("="*50)

if __name__ == "__main__":
    main()