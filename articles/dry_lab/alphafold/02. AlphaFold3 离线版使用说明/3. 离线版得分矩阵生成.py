import os
import json
import csv
import glob

# ================= 配置区域 =================
RESULTS_DIR = "Results"           # 结果根目录
OUTPUT_IPTM = "matrix_iptm.tsv"    # ipTM 矩阵
OUTPUT_PTM  = "matrix_ptm.tsv"     # pTM 矩阵
OUTPUT_RANK = "matrix_ranking.tsv" # Ranking Score 矩阵

def get_scores_from_file(json_path):
    """从 summary_confidences JSON 文件中提取三个核心分数"""
    try:
        with open(json_path, 'r') as f:
            data = json.load(f)
            # 提取三个分数，如果不存在则返回 None 或 0.0
            iptm = data.get('iptm', 0.0)
            ptm = data.get('ptm', 0.0)
            ranking = data.get('ranking_score', 0.0)
            return iptm, ptm, ranking
    except Exception as e:
        print(f"读取警告: 无法解析 {json_path} - {e}")
        return None, None, None

def main():
    print(f"正在深度扫描 {RESULTS_DIR} 目录下的最佳结果...")
    
    # data[protein_id][ligand_id] = (iptm, ptm, ranking)
    matrix_data = {}
    all_ligands = set()
    all_proteins = set()
    
    # 1. 递归搜索所有的 summary_confidences.json
    search_pattern = os.path.join(RESULTS_DIR, "**", "*_summary_confidences.json")
    all_files = glob.glob(search_pattern, recursive=True)
    
    if not all_files:
        print(f"错误: 在 {RESULTS_DIR} 下没找到任何 summary JSON 文件。")
        return

    print(f"初步找到 {len(all_files)} 个文件，正在处理...")

    valid_count = 0
    for file_path in all_files:
        # 2. 过滤掉 seed-X 文件夹内部的文件 (只留外面的最佳汇总)
        if "seed-" in file_path or "sample-" in file_path:
            continue
            
        # ---------------- 解析 ID 逻辑 (保持 V2 版本的鲁棒性) ----------------
        rel_path = os.path.relpath(file_path, RESULTS_DIR)
        path_parts = rel_path.split(os.sep)
        
        if len(path_parts) < 2:
            continue
            
        ligand_id = path_parts[0]  # CAS号
        
        # 解析蛋白 ID
        filename = os.path.basename(file_path)
        base_name = filename.replace("_summary_confidences.json", "")
        
        # 尝试去掉 CAS 前缀以获得纯蛋白名
        if "_" in base_name:
            if base_name.lower().startswith(ligand_id.lower()):
                temp_name = base_name[len(ligand_id):]
                if temp_name.startswith('_'):
                    temp_name = temp_name[1:]
                protein_id = temp_name
            else:
                protein_id = base_name
        else:
            protein_id = base_name

        protein_id = protein_id.split('.')[0] # 去后缀
        
        # ---------------- 提取分数 ----------------
        iptm, ptm, ranking = get_scores_from_file(file_path)
        
        if iptm is not None:
            if protein_id not in matrix_data:
                matrix_data[protein_id] = {}
            
            # 存储三个分数
            matrix_data[protein_id][ligand_id] = (iptm, ptm, ranking)
            all_ligands.add(ligand_id)
            all_proteins.add(protein_id)
            valid_count += 1

    print(f"提取完成！共找到 {valid_count} 个有效数据。")
    print(f"统计范围: {len(all_proteins)} 个蛋白 x {len(all_ligands)} 个配体")

    # 3. 排序
    sorted_ligands = sorted(list(all_ligands))
    sorted_proteins = sorted(list(all_proteins))
    
    # 4. 通用写入函数 (避免写三遍重复代码)
    def write_matrix(filename, index):
        print(f"正在生成 {filename} ...")
        with open(filename, 'w', newline='', encoding='utf-8') as f:
            writer = csv.writer(f, delimiter='\t')
            writer.writerow(["Protein"] + sorted_ligands)
            for prot in sorted_proteins:
                row = [prot]
                for lig in sorted_ligands:
                    val = matrix_data.get(prot, {}).get(lig)
                    # index: 0=iptm, 1=ptm, 2=ranking
                    row.append(f"{val[index]:.4f}" if val else "N/A")
                writer.writerow(row)

    # 分别生成三个矩阵
    write_matrix(OUTPUT_IPTM, 0)
    write_matrix(OUTPUT_PTM, 1)
    write_matrix(OUTPUT_RANK, 2)

    print("-" * 50)
    print("全部完成！请查看生成的三个 .tsv 文件。")

if __name__ == "__main__":
    main()