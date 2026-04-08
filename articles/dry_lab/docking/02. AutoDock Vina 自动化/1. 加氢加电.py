#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import os
import sys
import glob
import subprocess
import time

# ============= 0. 环境检查 =============
try:
    import pymol
    from pymol import cmd
except ImportError:
    print("错误: 未找到 'pymol' 模块。请先安装: sudo apt-get install python3-pymol")
    sys.exit(1)

# ============= 1. 基础工具函数 =============
def get_base_dir():
    if getattr(sys, 'frozen', False):
        return os.path.dirname(sys.executable)
    return os.path.dirname(os.path.abspath(__file__))

def find_files_recursive(directory, ext):
    if not ext.startswith('.'): ext = '.' + ext
    search_pattern = os.path.join(directory, '**', f'*{ext}')
    return sorted(glob.glob(search_pattern, recursive=True))

def choose_items(files, base_path, desc="文件"):
    if not files:
        print(f"提示：在 {base_path} 下未找到任何 {desc}")
        return []

    print(f"\n--- 找到以下 {desc} ---")
    for i, f in enumerate(files, 1):
        rel_path = os.path.relpath(f, base_path)
        print(f"  [{i}] {rel_path}")

    while True:
        try:
            prompt = (f"\n请输入欲操作的 {desc} 编号\n"
                      f" (例如: 1 | 1,2,5 | 1-10 | all): ")
            user_input = input(prompt).strip().lower()

            if not user_input: continue

            selected_indices = set()
            if user_input == 'all':
                selected_indices = set(range(len(files)))
            else:
                parts = user_input.split(',')
                for part in parts:
                    part = part.strip()
                    if '-' in part:
                        start_str, end_str = part.split('-')
                        for idx in range(int(start_str) - 1, int(end_str)):
                            if 0 <= idx < len(files): selected_indices.add(idx)
                    else:
                        idx = int(part) - 1
                        if 0 <= idx < len(files): selected_indices.add(idx)
            
            if not selected_indices:
                print("无效编号，请重试。")
                continue
            return [files[i] for i in sorted(list(selected_indices))]

        except Exception as e:
            print(f"输入错误: {e}")

# ============= 2. 核心处理 (PyMOL加氢转pdb + OpenBabel加电荷转pdbqt) =============
def process_one_file(input_file, output_dir, mode='ligand'):
    """
    单文件处理流：CIF -> PyMOL(加氢) -> TMP.PDB -> OpenBabel(电荷) -> PDBQT
    """
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    file_name = os.path.basename(input_file)
    name_no_ext = os.path.splitext(file_name)[0]
    
    # 定义路径
    temp_pdb = os.path.join(output_dir, f"temp_{name_no_ext}.pdb")
    final_pdbqt = os.path.join(output_dir, f"{name_no_ext}.pdbqt")

    try:
        # --- 阶段 1: PyMOL 加氢修复 ---
        # 启动 PyMOL (Quiet mode)
        pymol.finish_launching(['pymol', '-c', '-q'])
        cmd.reinitialize()
        cmd.load(input_file)
        cmd.h_add()            # 智能加氢
        cmd.save(temp_pdb)     # 保存为临时 PDB
        cmd.delete('all')
        
        # 检查中间文件是否生成
        if not os.path.exists(temp_pdb):
            print(f"  [失败] PyMOL 未能生成中间文件: {file_name}")
            return False

        # --- 阶段 2: Open Babel 格式转换与电荷计算 ---
        # 构造命令
        ob_cmd = ['obabel', '-ipdb', temp_pdb, '-opdbqt', '-O', final_pdbqt]
        
        # 添加参数
        if mode == 'ligand':
            # 配体: 允许扭转，计算电荷
            ob_cmd.extend(['--partialcharge', 'gasteiger'])
        else:
            # 受体: 刚性分子(-xr)，计算电荷
            ob_cmd.extend(['-xr', '--partialcharge', 'gasteiger'])
        
        # 执行 Open Babel
        subprocess.run(ob_cmd, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.PIPE)

        # --- 阶段 3: 清理 ---
        if os.path.exists(temp_pdb):
            os.remove(temp_pdb)

        return True

    except subprocess.CalledProcessError as e:
        print(f"\n  [OpenBabel 报错] {file_name}: {e.stderr.decode()}")
        if os.path.exists(temp_pdb): os.remove(temp_pdb)
        return False
    except Exception as e:
        print(f"\n  [PyMOL 报错] {file_name}: {e}")
        return False

# ============= 3. 主流程 =============
def main():
    base_path = get_base_dir()
    
    # 目录配置
    lig_in = os.path.join(base_path, "0. Ligand_molecular_cif")
    rec_in = os.path.join(base_path, "0. Receptor_protein_cif")
    lig_out = os.path.join(base_path, "1. Ligand_molecular_pdbqt")
    rec_out = os.path.join(base_path, "1. Receptor_protein_pdbqt")

    print("=== 流程1: CIF -> PDBQT (PyMOL内核 + OpenBabel) ===")

    # 1. 选择配体
    lig_files = find_files_recursive(lig_in, '.cif')
    chosen_ligs = choose_items(lig_files, lig_in, desc="配体")
    if not chosen_ligs: return

    # 2. 选择受体
    rec_files = find_files_recursive(rec_in, '.cif')
    chosen_recs = choose_items(rec_files, rec_in, desc="受体")
    if not chosen_recs: return

    # 3. 确认
    print("-" * 50)
    print(f"即将处理:\n - 配体: {len(chosen_ligs)} 个 -> 存入 {os.path.basename(lig_out)}\n - 受体: {len(chosen_recs)} 个 -> 存入 {os.path.basename(rec_out)}")
    if input("确认开始? (y/n): ").lower() not in ['y', 'yes']:
        return

    # 4. 执行配体
    print("\n>>> 正在处理配体...")
    for f in chosen_ligs:
        print(f"  处理: {os.path.basename(f)} ...")
        if process_one_file(f, lig_out, mode='ligand'):
            print("    -> 完成")

    # 5. 执行受体
    print("\n>>> 正在处理受体...")
    success_count = 0
    total_recs = len(chosen_recs)
    
    for i, f in enumerate(chosen_recs, 1):
        # 打印进度条
        sys.stdout.write(f"\r  进度: {i}/{total_recs} | 当前: {os.path.basename(f)}                    ")
        sys.stdout.flush()
        
        if process_one_file(f, rec_out, mode='receptor'):
            success_count += 1
            
    print(f"\n\n" + "="*50)
    print(f"全部任务结束。")
    print(f"受体成功率: {success_count}/{total_recs}")
    print(f"请检查 '1. ...' 文件夹下的最终 PDBQT 文件。")
    input("按回车键退出...")

if __name__ == "__main__":
    main()