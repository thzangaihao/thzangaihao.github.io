#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import os
import sys
import time
import pandas as pd
from Bio import Entrez
import urllib.error
import socket

# 设置全局默认超时时间（防止 NCBI 响应卡死）
socket.setdefaulttimeout(15)

# ================= 配置区 =================
# 【重要】请替换为你自己的邮箱，防止被 NCBI 封禁 IP
Entrez.email = "your_email@example.com" 
# ==========================================

def get_base_dir():
    if getattr(sys, 'frozen', False):
        return os.path.dirname(sys.executable)
    return os.path.dirname(os.path.abspath(__file__))

def fetch_genomes_by_genus(genus_name):
    """
    通过属名查询 NCBI Assembly 数据库，获取属下所有带有 Reference/Representative 基因组的种
    """
    if pd.isna(genus_name) or not str(genus_name).strip():
        return []

    genus_name = str(genus_name).strip()
    
    # 【修复Bug的关键点】：为带空格的 filter 加上了严格的双引号
    term = f'{genus_name}[Organism] AND ("reference genome"[Filter] OR "representative genome"[Filter]) AND latest[Filter]'
    
    results = []
    
    try:
        # 1. 搜索符合条件的 Assembly ID
        search_handle = Entrez.esearch(db="assembly", term=term, retmax=1000)
        search_record = Entrez.read(search_handle)
        search_handle.close()

        id_list = search_record.get("IdList", [])
        if not id_list:
            return [{"Genus(属)": genus_name, "Query_Status": "No Reference Genomes Found"}]

        # 2. 批量获取这些 ID 的详细摘要信息
        ids_str = ",".join(id_list)
        summary_handle = Entrez.esummary(db="assembly", id=ids_str)
        summaries = Entrez.read(summary_handle)
        summary_handle.close()

        # 提取数据 (Entrez 返回的 XML 结构为 DocumentSummarySet -> DocumentSummary)
        doc_summaries = summaries['DocumentSummarySet']['DocumentSummary']
        
        # 如果只有一个结果，Entrez 可能会将其解析为字典而非列表，这里做个兼容处理
        if isinstance(doc_summaries, dict):
            doc_summaries = [doc_summaries]

        for doc in doc_summaries:
            # 提取物种与分类 ID
            species_name = doc.get("SpeciesName", "")
            species_taxid = doc.get("SpeciesTaxid", "")
            strain_taxid = doc.get("Taxid", "") # 具体测序菌株的 TaxID
            
            # 提取 Accession ID
            synonyms = doc.get("Synonym", {})
            genbank_id = synonyms.get("Genbank", "")
            refseq_id = synonyms.get("RefSeq", "")
            
            category = doc.get("RefSeq_category", "")
            assembly_name = doc.get("AssemblyName", "")

            results.append({
                "Genus(属)": genus_name,
                "Species(种)": species_name,
                "Species_TaxID(种分类ID)": species_taxid,
                "GenBank_ID": genbank_id,
                "RefSeq_ID": refseq_id,
                "Category(级别)": category,
                "AssemblyName(组装名称)": assembly_name,
                "Strain_TaxID(菌株分类ID)": strain_taxid,
                "Query_Status": "Success"
            })

        return results

    except urllib.error.URLError as e:
        # 捕获网络异常
        return [{"Genus(属)": genus_name, "Query_Status": f"Network Error: {e}"}]
    except Exception as e:
        # 捕获其他未知异常
        return [{"Genus(属)": genus_name, "Query_Status": f"Error: {str(e)}"}]

def load_dataframe(filepath):
    ext = os.path.splitext(filepath)[1].lower()
    try:
        try:
            df = pd.read_csv(filepath, sep='\t')
            if len(df.columns) > 1: return df
        except:
            pass
        
        if ext == '.csv': return pd.read_csv(filepath)
        elif ext in ['.tsv', '.txt']: return pd.read_csv(filepath, sep='\t')
        elif ext in ['.xlsx', '.xls']: return pd.read_excel(filepath)
        else: return pd.read_csv(filepath)
    except Exception as e:
        print(f"读取数据表失败: {e}")
        sys.exit(1)

def main():
    print("=" * 60)
    print(" " * 8 + "NCBI 属下参考/代表性基因组 自动爬取工具 (修正版)")
    print("=" * 60)
    
    if Entrez.email == "your_email@example.com":
        print("⚠️ 警告：请在脚本代码第16行修改为你自己的邮箱地址！\n")

    print("请选择操作模式:")
    print("  [1] 手动输入单个属名 (例如: Alternaria)")
    print("  [2] 导入包含属名的表格进行批量查询")
    
    mode = input("\n你的选择 (1 或 2): ").strip()
    
    genera_to_search = []
    file_path = ""

    if mode == '1':
        single_genus = input("\n请输入要查询的属名 (英文): ").strip()
        if single_genus:
            genera_to_search.append(single_genus)
        else:
            print("未输入有效属名，退出程序。")
            sys.exit()
    elif mode == '2':
        file_path = input("\n请拖拽或输入表格文件路径 (.csv/.tsv/.xlsx): ").strip().strip('\'"')
        if not os.path.isfile(file_path):
            print("文件不存在，请检查路径。")
            sys.exit(1)
            
        df = load_dataframe(file_path)
        columns = df.columns.tolist()
        
        print(f"\n找到以下 {len(columns)} 列:")
        for i, col in enumerate(columns, 1):
            print(f"  [{i}] {col}")

        while True:
            try:
                col_idx = int(input(f"\n请输入包含【属名】的列编号 (1-{len(columns)}): ").strip()) - 1
                if 0 <= col_idx < len(columns):
                    query_col = columns[col_idx]
                    # 获取去重且非空的属名列表
                    genera_to_search = df[query_col].dropna().astype(str).str.strip().unique().tolist()
                    break
                print("编号超出范围。")
            except ValueError:
                print("请输入有效的数字。")
    else:
        print("无效的选择，退出程序。")
        sys.exit()

    print(f"\n🚀 开始查询，共 {len(genera_to_search)} 个属需要检索...")
    
    all_results = []
    
    for idx, genus in enumerate(genera_to_search, 1):
        print(f"[{idx}/{len(genera_to_search)}] 正在检索: {genus} ... ", end="", flush=True)
        
        results = fetch_genomes_by_genus(genus)
        all_results.extend(results)
        
        if results and results[0].get("Query_Status") == "Success":
            print(f"✅ 找到 {len(results)} 个参考/代表性基因组")
        else:
            status = results[0].get("Query_Status") if results else "Unknown Error"
            print(f"❌ 未找到或失败 ({status})")
            
        time.sleep(0.5) # 稍微延长了暂停时间，让请求更稳

    # 汇总并保存数据
    final_df = pd.DataFrame(all_results)
    
    # 构建输出文件名
    if mode == '1':
        save_name = f"{genera_to_search[0]}_genomes.csv"
    else:
        base_name = os.path.splitext(os.path.basename(file_path))[0]
        save_name = f"{base_name}_genomes_crawled.csv"
        
    output_path = os.path.join(get_base_dir(), save_name)
    
    # 保存为标准 CSV
    final_df.to_csv(output_path, sep=',', index=False, encoding='utf-8-sig')
    
    print("-" * 60)
    print(f"🎉 全部处理完成！共收集到 {len(all_results)} 条记录。")
    print(f"📄 结果已保存至:\n   {output_path}")

if __name__ == "__main__":
    main()