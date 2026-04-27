#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import os
import sys
import pandas as pd
import subprocess
import shutil
import math
import json

# ================= 配置与常量 =================
CONFIG_FILE = "saf_config.json"
SUMMARY_FILE = "assembly_summary_genbank.txt"
TASKS_DIR = "tasks_chunks"
DOWNLOADS_DIR = "downloads"
EXTRACT_DIR = "extracted_data"

FILE_TYPES = {
    "1": {"name": "Genome", "suffix": "_genomic.fna.gz"},
    "2": {"name": "Protein", "suffix": "_protein.faa.gz"},
    "3": {"name": "CDS", "suffix": "_cds_from_genomic.fna.gz"},
    "4": {"name": "GFF3", "suffix": "_genomic.gff.gz"}
}

def load_config():
    if os.path.exists(CONFIG_FILE):
        with open(CONFIG_FILE, 'r') as f:
            return json.load(f)
    return {"api_key": "", "chunk_size": 50, "timeout": 600}

def save_config(config):
    with open(CONFIG_FILE, 'w') as f:
        json.dump(config, f, indent=4)

def print_header():
    print("\n" + "=" * 70)
    print(" 🌱 SAF 数据库: 模块化拉取框架 (Task Split -> Download -> Extract)")
    print("=" * 70)

# ================= 模块 1: 任务拆封 =================
def step1_split_tasks(config):
    print("\n--- [模块 1: 任务拆封] ---")
    file_path = input("请输入包含 Assembly ID 的表格路径 (如 input.csv): ").strip()
    if not os.path.exists(file_path):
        print("❌ 找不到文件！")
        return

    # 读取文件逻辑
    try:
        with open(file_path, 'r', encoding='utf-8-sig') as f:
            first_line = f.readline()
        if '\t' in first_line: df = pd.read_csv(file_path, sep='\t', dtype=str)
        elif ',' in first_line: df = pd.read_csv(file_path, sep=',', dtype=str)
        else: df = pd.read_csv(file_path, sep=r'\s+', engine='python', dtype=str)
        
        columns = df.columns.tolist()
        print("找到以下列:")
        for i, col in enumerate(columns, 1): print(f"  [{i}] {col}")
        col_idx = input(f"请输入包含【Assembly ID】的列编号 (1-{len(columns)}): ").strip()
        col_name = columns[int(col_idx) - 1]
        
        ids = df[col_name].dropna().astype(str).str.strip().unique().tolist()
        valid_ids = [x for x in ids if x.startswith('GCA_') or x.startswith('GCF_')]
        print(f"✅ 成功提取到 {len(valid_ids)} 个有效的 Assembly ID。")
    except Exception as e:
        print(f"❌ 读取失败: {e}")
        return

    # 下载/读取总目录
    if not os.path.exists(SUMMARY_FILE):
        print(f"\n📥 正在下载 NCBI GenBank 总目录...")
        url = "https://ftp.ncbi.nlm.nih.gov/genomes/ASSEMBLY_REPORTS/assembly_summary_genbank.txt"
        subprocess.run(["wget", "-4", "-c", "-O", SUMMARY_FILE, url], check=True)
    
    print("🔍 正在与总目录匹配物种 ID...")
    df_summary = pd.read_csv(SUMMARY_FILE, sep='\t', skiprows=1, low_memory=False)
    df_summary.columns = [c.replace('#', '').strip() for c in df_summary.columns]
    df_matched = df_summary[df_summary['assembly_accession'].isin(valid_ids)].copy()
    
    print(f"🎯 成功匹配到 {len(df_matched)} 个基因组记录！")
    if len(df_matched) == 0: return

    # 处理 HTTPS 和 Basename
    df_matched['https_dir'] = df_matched['ftp_path'].astype(str).apply(lambda x: x.replace("ftp://", "https://").rstrip("/") if x != 'nan' else "")
    df_matched = df_matched[df_matched['https_dir'] != ""]
    df_matched['basename'] = df_matched['https_dir'].apply(lambda x: x.split("/")[-1])
    
    # 分块并保存
    os.makedirs(TASKS_DIR, exist_ok=True)
    chunk_size = config.get("chunk_size", 50)
    total_chunks = math.ceil(len(df_matched) / chunk_size)
    
    # 清空旧任务
    for f in os.listdir(TASKS_DIR): os.remove(os.path.join(TASKS_DIR, f))
        
    for i in range(total_chunks):
        chunk_df = df_matched.iloc[i*chunk_size : (i+1)*chunk_size]
        chunk_file = os.path.join(TASKS_DIR, f"chunk_{i+1:03d}.csv")
        chunk_df[['assembly_accession', 'https_dir', 'basename']].to_csv(chunk_file, index=False)
        
    print(f"✅ 任务已拆分为 {total_chunks} 个文件，保存在 ./{TASKS_DIR}/ 目录下。")

# ================= 模块 2: 选择性下载 =================
def step2_download(config):
    print("\n--- [模块 2: 分块下载] ---")
    if not os.path.exists(TASKS_DIR) or not os.listdir(TASKS_DIR):
        print("❌ 没有找到拆分的任务文件，请先执行模块 1。")
        return
        
    chunks = sorted(os.listdir(TASKS_DIR))
    print("可用任务块:")
    for i, chunk in enumerate(chunks, 1): print(f"  [{i}] {chunk}")
    print("  [0] 下载所有块")
    
    chunk_choice = input("\n请选择要下载的任务编号 (例如 1，或 0 下载全部): ").strip()
    target_chunks = chunks if chunk_choice == "0" else [chunks[int(chunk_choice)-1]]
    
    print("\n请选择需要下载的文件类型 (可多选，用逗号分隔，如 1,3):")
    for k, v in FILE_TYPES.items(): print(f"  [{k}] {v['name']}")
    type_choices = input("选择: ").strip().split(',')
    selected_suffixes = [FILE_TYPES[c.strip()]['suffix'] for c in type_choices if c.strip() in FILE_TYPES]
    
    if not selected_suffixes:
        print("❌ 未选择有效的文件类型！")
        return

    api_key_header = f"--header=api-key: {config['api_key']}" if config.get("api_key") else ""
    timeout_val = str(config.get("timeout", 600))
    
    for chunk in target_chunks:
        chunk_path = os.path.join(TASKS_DIR, chunk)
        chunk_name = chunk.replace('.csv', '')
        out_dir = os.path.join(DOWNLOADS_DIR, chunk_name)
        os.makedirs(out_dir, exist_ok=True)
        
        df_chunk = pd.read_csv(chunk_path)
        urls = []
        for _, row in df_chunk.iterrows():
            for suffix in selected_suffixes:
                urls.append(f"{row['https_dir']}/{row['basename']}{suffix}")
                
        url_list_file = os.path.join(out_dir, "temp_urls.txt")
        with open(url_list_file, "w") as f:
            f.write("\n".join(urls))
            
        print(f"\n🚀 开始下载 {chunk_name} (共 {len(urls)} 个文件) ...")
        # 优化 wget 参数：设置超时，精简输出
        cmd = ["wget", "-4", "-c", "-nc", f"--timeout={timeout_val}", "--tries=3"]
        if api_key_header: cmd.append(api_key_header)
        cmd.extend(["-q", "--show-progress", "-P", out_dir, "-i", url_list_file])
        
        try:
            # 增加子进程整体超时机制
            subprocess.run(cmd, timeout=len(urls) * int(timeout_val))
            print(f"✅ {chunk_name} 下载进程结束。")
        except subprocess.TimeoutExpired:
            print(f"⚠️ {chunk_name} 下载超时被强制中断！下次运行可断点续传。")
        except KeyboardInterrupt:
            print(f"\n⚠️ 用户手动中断。")
            break
            
        if os.path.exists(url_list_file): os.remove(url_list_file)

# ================= 模块 3: 文件提取与报告 =================
def step3_extract():
    print("\n--- [模块 3: 提取与核验] ---")
    if not os.path.exists(TASKS_DIR):
        print("❌ 找不到任务清单，无法核验。")
        return
        
    chunks = sorted(os.listdir(TASKS_DIR))
    print("选择要提取的任务块对应的数据:")
    for i, chunk in enumerate(chunks, 1): print(f"  [{i}] {chunk}")
    print("  [0] 提取所有块")
    
    chunk_choice = input("\n选择编号: ").strip()
    target_chunks = chunks if chunk_choice == "0" else [chunks[int(chunk_choice)-1]]
    
    print("\n请选择要提取的文件类型 (单选):")
    for k, v in FILE_TYPES.items(): print(f"  [{k}] {v['name']}")
    type_choice = input("选择: ").strip()
    if type_choice not in FILE_TYPES: return
    
    target_suffix = FILE_TYPES[type_choice]['suffix']
    target_type_name = FILE_TYPES[type_choice]['name']
    
    out_dir = os.path.join(EXTRACT_DIR, target_type_name)
    os.makedirs(out_dir, exist_ok=True)
    
    report_lines = [f"Extraction Report: {target_type_name}", "="*40]
    success_cnt, fail_cnt, corrupt_cnt = 0, 0, 0
    
    for chunk in target_chunks:
        chunk_path = os.path.join(TASKS_DIR, chunk)
        chunk_name = chunk.replace('.csv', '')
        source_dir = os.path.join(DOWNLOADS_DIR, chunk_name)
        
        df_chunk = pd.read_csv(chunk_path)
        for _, row in df_chunk.iterrows():
            acc_id = row['assembly_accession']
            expected_file = f"{row['basename']}{target_suffix}"
            gz_path = os.path.join(source_dir, expected_file)
            
            if os.path.exists(gz_path):
                # 【新增逻辑】：静默测试 gzip 文件的完整性
                test_result = subprocess.run(["gzip", "-t", gz_path], stderr=subprocess.PIPE)
                
                if test_result.returncode != 0:
                    # 返回码非0，说明文件不完整/损坏（下载了一半）
                    report_lines.append(f"{acc_id}\tCORRUPTED (Incomplete Download)")
                    corrupt_cnt += 1
                else:
                    # 文件完整，进行解压操作
                    unzipped_name = expected_file[:-3]
                    unzipped_path = os.path.join(out_dir, unzipped_name)
                    with open(unzipped_path, "w") as outfile:
                        subprocess.run(["gzip", "-d", "-c", gz_path], stdout=outfile)
                    
                    report_lines.append(f"{acc_id}\tSUCCESS")
                    success_cnt += 1
            else:
                report_lines.append(f"{acc_id}\tNOT FOUND")
                fail_cnt += 1
                
    report_file = os.path.join(out_dir, "extraction_report.txt")
    with open(report_file, "w") as f:
        f.write("\n".join(report_lines))
        
    print(f"\n🎉 提取完成！成功: {success_cnt} | 缺失: {fail_cnt} | 损坏/未下完: {corrupt_cnt}")
    print(f"📂 文件保存在: {out_dir}")
    print(f"📄 报告保存在: {report_file}")

# ================= 系统设置 =================
def settings(config):
    print("\n--- [设置] ---")
    print(f"1. NCBI API Key (当前: {config['api_key'] or '未设置'})")
    print(f"2. 单个分块大小 (当前: {config['chunk_size']})")
    print(f"3. wget 下载单次超时秒数 (当前: {config['timeout']}s)")
    
    c = input("选择要修改的项 (回车返回): ").strip()
    if c == "1": config['api_key'] = input("输入 API Key: ").strip()
    elif c == "2": config['chunk_size'] = int(input("输入分块大小: ").strip())
    elif c == "3": config['timeout'] = int(input("输入超时时间(秒): ").strip())
    save_config(config)

def main():
    config = load_config()
    while True:
        print_header()
        print(" [1] 任务拆封 (生成下载清单与分块)")
        print(" [2] 执行下载 (支持指定格式、断点续传)")
        print(" [3] 文件提取与核验报告生成")
        print(" [4] 全局设置 (API Key, 分块大小, 超时控制)")
        print(" [0] 退出")
        
        choice = input("\n请输入操作编号: ").strip()
        if choice == "1": step1_split_tasks(config)
        elif choice == "2": step2_download(config)
        elif choice == "3": step3_extract()
        elif choice == "4": settings(config)
        elif choice == "0": break
        else: print("无效输入。")

if __name__ == "__main__":
    main()