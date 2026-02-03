#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import sys
import glob
import subprocess
import pandas as pd
import datetime

# ================= 0. 核心配置 =================
OUTPUT_DIR_NAME = "10_Variant_Details_With_Samples"
BCFTOOLS_EXEC = "bcftools"

# ================= 1. 交互逻辑 (复用) =================
def get_base_dir():
    return os.path.dirname(os.path.abspath(__file__))

def find_files(exts, path=None):
    if path is None: path = get_base_dir()
    all_files = []
    if isinstance(exts, str): exts = [exts]
    for ext in exts:
        if not ext.startswith('.'): ext = '.' + ext
        search_pattern = os.path.join(path, '**', f'*{ext}')
        all_files.extend(glob.glob(search_pattern, recursive=True))
    return sorted(list(set([f for f in all_files if OUTPUT_DIR_NAME not in f])))

def choose_files_multi(files, desc="文件"):
    if not files:
        print(f"[提示] 未找到任何 {desc}")
        return []
    print(f"\n2. 找到 {len(files)} 个 {desc}:")
    limit = 15
    for i, f in enumerate(files[:limit], 1):
        print(f"  [{i}] {os.path.basename(f)} ({os.path.relpath(f, get_base_dir())})")
    if len(files) > limit: print(f"  ... (共 {len(files)} 个)")
    while True:
        try:
            prompt = f"\n3. 请输入编号 (支持 1,2 | all): "
            u = input(prompt).strip().lower()
            if not u: continue
            sel = set()
            if u in ['all', 'a']: sel = set(range(len(files)))
            else:
                for p in u.split(','):
                    if '-' in p: s, e = map(int, p.split('-')); sel.update(range(s-1, e))
                    else: sel.add(int(p)-1)
            res = [files[i] for i in sorted(sel) if 0 <= i < len(files)]
            if res: return res
        except: pass

def choose_file_single(files, desc="文件"):
    if not files: return None
    print(f"\n>>> 请选择作为基准的 {desc} (单选):")
    for i, f in enumerate(files, 1): print(f"  [{i}] {os.path.basename(f)}")
    while True:
        try:
            idx = int(input("编号: ").strip()) - 1
            if 0 <= idx < len(files): return files[idx]
        except: pass

# ================= 2. 核心处理 (含样本提取) =================

def get_ids_from_file(file_path):
    try:
        if file_path.endswith('.csv'): df = pd.read_csv(file_path)
        else:
            df = pd.read_csv(file_path, sep='\t')
            if len(df.columns) == 1 and ',' in df.iloc[0,0]: df = pd.read_csv(file_path, sep=',')
        
        possible_cols = ['rs', 'RS', 'snp', 'SNP', 'id', 'ID', 'Variant', 'Target_ID']
        for col in df.columns:
            if col in possible_cols:
                return df[col].dropna().astype(str).unique().tolist()
        return []
    except: return []

def trace_worker(vcf_path, ids, source_name):
    if not ids: return

    # 输出设置
    timestamp = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
    clean_name = os.path.splitext(source_name)[0].replace("_sig_sites", "").replace("_filtered", "")
    out_filename = f"{timestamp}_{clean_name}_Sample_Details.tsv"
    
    out_dir = os.path.join(get_base_dir(), OUTPUT_DIR_NAME)
    if not os.path.exists(out_dir): os.makedirs(out_dir)
    out_path = os.path.join(out_dir, out_filename)

    print(f"  -> 正在提取 {len(ids)} 个变异的携带者信息...")
    target_set = set(ids)
    results = []

    # 构造 bcftools 命令
    # 关键修改: 增加 [%SAMPLE:%GT;] 来获取所有样本的基因型
    # 格式: CHROM, POS, END, ID, LEN, TYPE, REF, ALT, [SAMPLE:GT;SAMPLE:GT...]
    cmd = [
        BCFTOOLS_EXEC, "query", 
        "-f", '%CHROM\t%POS\t%INFO/END\t%ID\t%INFO/SVLEN\t%INFO/SVTYPE\t%REF\t%ALT\t[%SAMPLE:%GT;]\n',
        vcf_path
    ]
    
    try:
        process = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        
        for line in process.stdout:
            cols = line.strip().split('\t')
            if len(cols) < 9: continue
            
            # 精确匹配 ID
            if cols[3] in target_set:
                # 解析样本信息
                # cols[8] 是一长串 "SF1:0/0;SF2:0/1;..."
                raw_samples = cols[8].split(';')
                carriers = []
                
                for s_info in raw_samples:
                    if ':' not in s_info: continue
                    sample_id, gt = s_info.split(':')
                    
                    # 筛选逻辑: 基因型包含 '1' (即 0/1, 1/0, 1/1) 且不是缺失
                    if '1' in gt: 
                        carriers.append(sample_id)
                
                results.append({
                    "Target_ID": cols[3],
                    "Chr": cols[0],
                    "Start_Pos": cols[1],
                    "End_Pos": cols[2] if cols[2] != '.' else "NA",
                    "SV_Type": cols[5],
                    "SV_Length": cols[4] if cols[4] != '.' else "NA",
                    "Carriers_Count": len(carriers),
                    "Carriers_List": ",".join(carriers) # 用逗号连接所有样本ID
                })

        process.wait()
        
        if results:
            df = pd.DataFrame(results)
            # 调整顺序
            cols_seq = ["Target_ID", "Chr", "Start_Pos", "End_Pos", "SV_Type", "SV_Length", "Carriers_Count", "Carriers_List"]
            df = df[cols_seq]
            df.to_csv(out_path, sep='\t', index=False)
            print(f"  [成功] 结果已保存 (含样本ID): {out_filename}")
        else:
            print(f"  [警告] 未找到匹配变异。")

    except Exception as e:
        print(f"  [异常] {e}")

# ================= 主函数 =================
def main():
    print("==============================================")
    print("   Step 10: 变异溯源 + 样本锁定 (Sample Tracer)")
    print("==============================================")
    
    # 1. 选显著文件
    sig_files = find_files(['.txt', '.tsv', '.csv'])
    sig_candidates = [f for f in sig_files if "sig" in f or "assoc" in f or "08_" in f]
    selected_sig_files = choose_files_multi(sig_candidates, "显著位点文件")
    if not selected_sig_files: return

    # 2. 选 VCF
    print("\n>>> 第二步: 选择基准 VCF 文件")
    vcf_files = find_files(['.vcf.gz'])
    vcf_files = sorted(vcf_files, key=lambda x: "Merged_Population.vcf.gz" not in x)
    target_vcf = choose_file_single(vcf_files, "VCF 文件")
    if not target_vcf: return

    # 3. 运行
    print(f"\n开始处理...")
    for i, f in enumerate(selected_sig_files, 1):
        print(f"\n>>> 任务 [{i}/{len(selected_sig_files)}]: {os.path.basename(f)}")
        ids = get_ids_from_file(f)
        trace_worker(target_vcf, ids, os.path.basename(f))
    
    print(f"\n所有任务完成！结果存放在: {OUTPUT_DIR_NAME}")

if __name__ == "__main__":
    main()