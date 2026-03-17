#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import os
import sys
import glob
import json

# ============= 基础功能函数 (保持原汁原味) =============
def get_base_dir():
    if getattr(sys, 'frozen', False):
        return os.path.dirname(sys.executable)
    return os.path.dirname(os.path.abspath(__file__))

def find_files(exts, path=None):
    """支持查找多种后缀"""
    if path is None:
        path = get_base_dir()
    
    files = []
    if isinstance(exts, str): exts = [exts]
    
    for ext in exts:
        if not ext.startswith('.'): ext = '.' + ext
        search_pattern = os.path.join(path, '**', f'*{ext}')
        files.extend(glob.glob(search_pattern, recursive=True))
    return sorted(list(set(files)))

def choose_file(files, desc="文件"):
    if not files:
        print(f"提示：未找到任何 {desc}")
        return []
    print(f"\n找到 {len(files)} 个 {desc}:")
    for i, f in enumerate(files, 1):
        print(f"  [{i}] {os.path.basename(f)}") # 只显示文件名，简洁
    
    while True:
        try:
            prompt = f"\n请选择 {desc} 编号 (如: 1,2 | all): "
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
                        start_str, end_str = part.split('-')
                        for idx in range(int(start_str)-1, int(end_str)):
                            if 0 <= idx < len(files): selected_indices.add(idx)
                    else:
                        idx = int(part) - 1
                        if 0 <= idx < len(files): selected_indices.add(idx)
            
            if not selected_indices: continue
            return [files[i] for i in sorted(list(selected_indices))]
        except ValueError:
            print("输入错误，请输入有效数字。")

def make_sure(text="\n确认执行操作? (y/n): "):
    response = input(text).strip().lower()
    return response in ['y', 'yes']

# ============= 核心解析逻辑 =============

def parse_fasta(fasta_path):
    """解析 FASTA 文件，返回 [(header, sequence), ...]"""
    sequences = []
    current_header = None
    current_seq = []
    
    with open(fasta_path, 'r', encoding='utf-8') as f:
        for line in f:
            line = line.strip()
            if not line: continue
            if line.startswith(">"):
                if current_header:
                    sequences.append((current_header, "".join(current_seq)))
                current_header = line[1:].split()[0] # 取>后的第一个单词作为ID
                current_seq = []
            else:
                current_seq.append(line)
        # 添加最后一个
        if current_header:
            sequences.append((current_header, "".join(current_seq)))
    return sequences

def parse_custom_smiles(file_path):
    """
    解析自定义 SMILES 格式:
    >CAS_ID
    SMILES_STRING
    """
    ligands = [] # 格式: [(cas_id, smiles), ...]
    current_id = None
    
    with open(file_path, 'r', encoding='utf-8') as f:
        for line in f:
            line = line.strip()
            if not line: continue
            
            if line.startswith(">"):
                current_id = line[1:].strip()
            else:
                # 如果有ID，这行就是SMILES
                if current_id:
                    ligands.append((current_id, line))
                    current_id = None # 重置，等待下一个 >
    return ligands

def generate_af3_local_json(name, protein_seq, ligand_smiles):
    """
    生成符合 AlphaFold 3 本地版 (Local Dialect) 的 JSON 结构
    包含必须的 'dialect' 和 'version' 字段
    """
    return {
        "name": name,
        "dialect": "alphafold3",
        "version": 1,
        "sequences": [
            {
                "protein": {
                    "id": "A",
                    "sequence": protein_seq
                }
            },
            {
                "ligand": {
                    "id": "B",
                    "smiles": ligand_smiles
                }
            }
        ],
        "modelSeeds": [1]
    }

# ============= 主程序 =============
def main():
    print("="*60)
    print("      AF3 大规模对接任务生成器 (Local Dialect 版)")
    print("      功能：FASTA蛋白 x SMILES小分子 -> 批量JSON")
    print("="*60)

    # 1. 选择 FASTA 文件 (蛋白)
    fasta_files = find_files(["fasta", "fa"])
    selected_fastas = choose_file(fasta_files, "FASTA 蛋白序列文件")
    if not selected_fastas: return

    # 2. 选择 SMILES 文件 (配体)
    # 支持 txt, smi 等格式
    smiles_files = find_files(["smiles", "SMILES", "txt"])
    selected_smiles_files = choose_file(smiles_files, "SMILES 配体列表文件")
    if not selected_smiles_files: return

    # 3. 准备数据
    all_proteins = []
    print("\n[1/3] 正在读取蛋白序列...")
    for f in selected_fastas:
        prots = parse_fasta(f)
        all_proteins.extend(prots)
        print(f"  - {os.path.basename(f)}: 读取到 {len(prots)} 条序列")

    all_ligands = []
    print("\n[2/3] 正在读取配体 SMILES...")
    for f in selected_smiles_files:
        ligs = parse_custom_smiles(f)
        all_ligands.extend(ligs)
        print(f"  - {os.path.basename(f)}: 读取到 {len(ligs)} 个配体")

    if not all_proteins or not all_ligands:
        print("\n错误：蛋白或配体数据为空，无法继续。")
        return

    total_tasks = len(all_proteins) * len(all_ligands)
    print(f"\n[统计] 蛋白: {len(all_proteins)} x 配体: {len(all_ligands)} = 总任务数: {total_tasks}")

    # 4. 执行生成
    if not make_sure(f"\n即将生成 {total_tasks} 个 JSON 文件，确认执行? (y/n): "):
        print("操作已取消。")
        return

    base_output_dir = os.path.join(get_base_dir(), "JSON")
    print(f"\n[3/3] 正在写入文件到 {base_output_dir} ...")

    count = 0
    for cas_id, smiles in all_ligands:
        # 1. 创建 CAS 号文件夹
        # 替换非法字符，防止文件名报错
        safe_cas = cas_id.replace("/", "_").replace("\\", "_")
        cas_dir = os.path.join(base_output_dir, safe_cas)
        
        if not os.path.exists(cas_dir):
            os.makedirs(cas_dir)
        
        # 2. 遍历所有蛋白
        for prot_id, prot_seq in all_proteins:
            # 组合任务名: CAS号_蛋白ID
            # 例如: 328968-36-1_Fs000286.t01
            job_name = f"{safe_cas}_{prot_id}"
            
            json_content = generate_af3_local_json(job_name, prot_seq, smiles)
            
            # 写入文件
            file_name = f"{job_name}.json"
            file_path = os.path.join(cas_dir, file_name)
            
            with open(file_path, 'w', encoding='utf-8') as jf:
                json.dump(json_content, jf, indent=2, ensure_ascii=False)
            
            count += 1
            if count % 100 == 0:
                print(f"  已生成 {count}/{total_tasks} ...")

    print("-" * 50)
    print(f"全部完成！共生成 {count} 个 JSON 任务文件。")
    print(f"请将 '{base_output_dir}' 文件夹上传至集群进行计算。")

if __name__ == "__main__":
    main()