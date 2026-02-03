#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import sys
import glob
import pandas as pd
import datetime

# ================= 0. 核心配置 =================
OUTPUT_DIR_NAME = "14_IGV_Tracks"

# ================= 1. 交互逻辑 =================
def get_base_dir():
    return os.path.dirname(os.path.abspath(__file__))

def find_files(ext, path=None):
    if not ext.startswith('.'): ext = '.' + ext
    if path is None: path = get_base_dir()
    search_pattern = os.path.join(path, '**', f'*{ext}')
    # 排除之前的输出目录
    files = sorted(glob.glob(search_pattern, recursive=True))
    return [f for f in files if OUTPUT_DIR_NAME not in f and "13_Variant_Details" in f]

def choose_files_multi(files, desc="文件"):
    if not files:
        print(f"[提示] 未找到任何 {desc} (请先运行 Step 13)")
        return []

    print(f"\n2. 找到 {len(files)} 个 {desc}:")
    for i, f in enumerate(files, 1):
        print(f"  [{i}] {os.path.basename(f)}")

    while True:
        try:
            user_input = input(f"\n3. 请输入编号 (all / 1,2): ").strip().lower()
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
            if selected_files: return selected_files
        except: pass

# ================= 2. 核心转换逻辑 =================
def convert_to_gff3(tsv_path):
    print(f"\n--- 正在转换: {os.path.basename(tsv_path)} ---")
    
    try:
        df = pd.read_csv(tsv_path, sep='\t')
        
        # 检查必要的列
        required_cols = ['Chr', 'Start_Pos', 'Target_ID']
        if not all(col in df.columns for col in required_cols):
            print(f"  [跳过] 文件缺少必要列 (需要 Step 13 的输出格式)")
            return

        # 准备 GFF3 数据列表
        gff_lines = []
        gff_lines.append("##gff-version 3")
        
        for _, row in df.iterrows():
            chrom = str(row['Chr'])
            start = int(row['Start_Pos'])
            sv_type = str(row.get('SV_Type', 'Sequence_Variant'))
            
            # 计算终止位置 (End)
            end = start # 默认点突变
            
            # 逻辑1: 优先使用 End_Pos
            if 'End_Pos' in row and pd.notna(row['End_Pos']) and str(row['End_Pos']) != 'NA':
                end = int(row['End_Pos'])
            
            # 逻辑2: 如果没有 End_Pos，尝试用 Length 推算
            elif 'SV_Length' in row and pd.notna(row['SV_Length']) and str(row['SV_Length']) != 'NA':
                try:
                    length = abs(int(row['SV_Length']))
                    end = start + length
                except: pass
            
            # 确保 end >= start
            if end < start: end = start
            
            # 构建属性栏 (Attributes)
            # ID=...;Name=...;Note=...
            target_id = str(row['Target_ID'])
            sv_len = str(row.get('SV_Length', 'NA'))
            
            attributes = f"ID={target_id};Name={target_id};Type={sv_type};Length={sv_len}"
            
            # GFF3 9列标准
            # seqid source type start end score strand phase attributes
            line_cols = [
                chrom,
                "GWAS_Significant", # source
                sv_type,            # type (e.g., DEL, INS)
                str(start),
                str(end),
                ".",                # score
                ".",                # strand
                ".",                # phase
                attributes
            ]
            gff_lines.append("\t".join(line_cols))

        # 保存文件
        base_dir = get_base_dir()
        out_dir = os.path.join(base_dir, OUTPUT_DIR_NAME)
        if not os.path.exists(out_dir): os.makedirs(out_dir)
        
        timestamp = datetime.datetime.now().strftime("%Y%m%d")
        original_name = os.path.splitext(os.path.basename(tsv_path))[0]
        # 清理文件名
        original_name = original_name.replace("_Details", "").replace("Significant_", "")
        
        out_filename = f"{timestamp}_{original_name}_IGV.gff3"
        out_path = os.path.join(out_dir, out_filename)
        
        with open(out_path, 'w') as f:
            f.write("\n".join(gff_lines))
            
        print(f"  [成功] GFF3 已生成: {out_filename}")
        print(f"  -> 包含 {len(gff_lines)-1} 个变异区间")

    except Exception as e:
        print(f"  [错误] {e}")

# ================= 主函数 =================
def main():
    print("==============================================")
    print("   Step 14: IGV 专用轨道生成器 (GFF3)")
    print("==============================================")
    
    # 1. 寻找 Step 13 的结果文件
    tsv_files = find_files('.tsv')
    if not tsv_files:
        print("未找到 Step 13 的结果文件。请先运行溯源脚本。")
        return

    # 2. 多选
    selected = choose_files_multi(tsv_files, "溯源结果文件 (TSV)")
    if not selected: return

    # 3. 批量转换
    for f in selected:
        convert_to_gff3(f)
        
    print(f"\n" + "="*50)
    print(f"转换完成！文件存放在: {OUTPUT_DIR_NAME}")
    print("请将 .gff3 文件拖入 IGV 即可查看。")
    print("="*50)

if __name__ == "__main__":
    main()