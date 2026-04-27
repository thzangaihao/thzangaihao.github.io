#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import os
import sys
import time
import pandas as pd
from Bio import Entrez
import urllib.request
import urllib.error
import socket

# 设置全局默认超时时间
socket.setdefaulttimeout(15)

# ================= 配置区 =================
# 已根据你的偏好设置默认邮箱
Entrez.email = "thzangaihao@outlook.com" 
# ==========================================

def get_base_dir():
    if getattr(sys, 'frozen', False):
        return os.path.dirname(sys.executable)
    return os.path.dirname(os.path.abspath(__file__))

def check_url_exists(url):
    """
    轻量级验证链接真实性
    """
    try:
        req = urllib.request.Request(url, method='HEAD', headers={'User-Agent': 'Mozilla/5.0'})
        with urllib.request.urlopen(req, timeout=5) as response:
            return response.status == 200
    except:
        return False

def fetch_assembly_links(accession_id):
    """
    逻辑改动：
    1. 爬取 Assembly Level。
    2. 如果 Level 是 Contig，则链接全部强制为 NA（即使物理链接存在也不返回）。
    """
    if pd.isna(accession_id) or not str(accession_id).strip():
        return []

    accession_id = str(accession_id).strip()
    term = f'"{accession_id}"[Assembly Accession]'
    
    try:
        # 1. 检索 UID
        search_handle = Entrez.esearch(db="assembly", term=term, retmax=1)
        search_record = Entrez.read(search_handle)
        search_handle.close()

        id_list = search_record.get("IdList", [])
        if not id_list:
            return {"Accession": accession_id, "Query_Status": "Not Found"}

        uid = id_list[0]

        # 2. 获取元数据摘要
        summary_handle = Entrez.esummary(db="assembly", id=uid)
        summaries = Entrez.read(summary_handle)
        summary_handle.close()

        doc = summaries['DocumentSummarySet']['DocumentSummary'][0]

        # 3. 提取官方组装层级 (AssemblyStatus)
        # 可能的值: 'Contig', 'Scaffold', 'Chromosome', 'Complete Genome'
        assembly_level = doc.get("AssemblyStatus", "Unknown")
        
        # 核心逻辑：判断是否达到 Scaffold 或更高层级
        # 如果只有 Contig 层级，即便有下载链接也将其标记为 NA
        has_scaffold_or_above = assembly_level in ['Scaffold', 'Chromosome', 'Complete Genome']

        # 4. 获取 FTP 根路径
        ftp_path = doc.get("FtpPath_RefSeq") or doc.get("FtpPath_GenBank") or ""

        if not ftp_path:
            return {
                "Accession": accession_id,
                "Assembly_Level": assembly_level,
                "Query_Status": "No FTP Path"
            }

        # 5. 生成结果
        if not has_scaffold_or_above:
            # 如果是 Contig 层级，直接置为 NA，不再进行后续验证
            return {
                "Accession": accession_id,
                "Assembly_Level": assembly_level,
                "Genome_FNA_Link": "NA",
                "Annotation_GFF_Link": "NA",
                "Protein_FAA_Link": "NA",
                "Query_Status": f"Level is {assembly_level} (Filtered)"
            }
        else:
            # 达到 Scaffold 以上层级，拼装并验证链接
            base_url = ftp_path.replace("ftp://", "https://")
            basename = base_url.split("/")[-1]

            cand_fna = f"{base_url}/{basename}_genomic.fna.gz"
            cand_gff = f"{base_url}/{basename}_genomic.gff.gz"
            cand_faa = f"{base_url}/{basename}_protein.faa.gz"

            return {
                "Accession": accession_id,
                "Assembly_Level": assembly_level,
                "Genome_FNA_Link": cand_fna if check_url_exists(cand_fna) else "Not Found",
                "Annotation_GFF_Link": cand_gff if check_url_exists(cand_gff) else "Not Found",
                "Protein_FAA_Link": cand_faa if check_url_exists(cand_faa) else "Not Found",
                "Query_Status": "Success"
            }

    except Exception as e:
        return {"Accession": accession_id, "Query_Status": f"Error: {str(e)}"}

def load_dataframe(filepath):
    ext = os.path.splitext(filepath)[1].lower()
    try:
        if ext == '.csv': return pd.read_csv(filepath)
        elif ext in ['.tsv', '.txt']: return pd.read_csv(filepath, sep='\t')
        elif ext in ['.xlsx', '.xls']: return pd.read_excel(filepath)
        else: return pd.read_csv(filepath)
    except Exception as e:
        print(f"读取表格失败: {e}")
        sys.exit(1)

def main():
    print("=" * 75)
    print(" " * 5 + "NCBI Assembly 链接爬取工具 (严格层级校验: Contig -> NA)")
    print("=" * 75)
    
    mode = input("\n[1] 手动输入 ID\n[2] 批量导入表格\n请选择: ").strip()
    
    accessions_to_search = []
    df = None
    query_col = ""

    if mode == '1':
        acc = input("\n请输入 Assembly ID: ").strip()
        if acc: accessions_to_search.append(acc)
    elif mode == '2':
        file_path = input("\n请输入/拖拽表格路径: ").strip().strip('\'"')
        df = load_dataframe(file_path)
        cols = df.columns.tolist()
        for i, c in enumerate(cols, 1): print(f"  [{i}] {c}")
        col_idx = int(input(f"\n请选择包含 ID 的列编号 (1-{len(cols)}): ")) - 1
        query_col = cols[col_idx]
        accessions_to_search = df[query_col].dropna().astype(str).str.strip().unique().tolist()
    else:
        sys.exit()

    print(f"\n🚀 开始处理 {len(accessions_to_search)} 个 ID...")
    
    results = []
    for idx, acc in enumerate(accessions_to_search, 1):
        print(f"[{idx}/{len(accessions_to_search)}] {acc} ... ", end="", flush=True)
        info = fetch_assembly_links(acc)
        results.append(info)
        
        level = info.get('Assembly_Level', 'N/A')
        status = info.get('Query_Status', '')
        
        if "Success" in status:
            print(f"✅ {level}")
        elif "Filtered" in status:
            print(f"⚠️ {level} (已置 NA)")
        else:
            print(f"❌ {status}")
            
        time.sleep(0.3)

    # 结果导出
    results_df = pd.DataFrame(results)
    output_name = "NCBI_Verified_Links_Filtered.csv"
    
    if mode == '2':
        final_df = pd.merge(df, results_df, left_on=query_col, right_on="Accession", how="left")
        final_df.to_csv(output_name, index=False, encoding='utf-8-sig')
    else:
        results_df.to_csv(output_name, index=False, encoding='utf-8-sig')

    print(f"\n🎉 完成！结果已存入: {output_name}")

if __name__ == "__main__":
    main()