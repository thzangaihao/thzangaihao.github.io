#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import os
import sys
import time
import pandas as pd
from Bio import Entrez
import urllib.error
import socket
import xml.etree.ElementTree as ET

# 设置全局默认超时时间
socket.setdefaulttimeout(15)

# ================= 配置区 =================
# 【重要】请替换为你自己的邮箱
Entrez.email = "thzangaihao@outlook.com" 
# ==========================================

def get_base_dir():
    if getattr(sys, 'frozen', False):
        return os.path.dirname(sys.executable)
    return os.path.dirname(os.path.abspath(__file__))

def parse_assembly_meta(meta_xml_str):
    """
    解析 NCBI Assembly 返回的 Meta XML 字符串，提取组装统计信息
    """
    stats = {
        "Total_Length": "",
        "Scaffold_Count": "",
        "Scaffold_N50": "",
        "Scaffold_L50": "",
        "Contig_Count": "",
        "Contig_N50": "",
        "GC_Percent": ""
    }
    
    if not meta_xml_str:
        return stats

    try:
        # 将 Meta 字符串包裹在 root 节点中以确保是合法的 XML
        root = ET.fromstring(f"<root>{meta_xml_str}</root>")
        
        # 遍历所有的 Stat 节点
        for stat in root.findall(".//Stat"):
            category = stat.attrib.get('category', '')
            seq_tag = stat.attrib.get('sequence_tag', '')
            
            # 我们通常提取 "all" 序列层级的统计数据
            if seq_tag == "all":
                if category == "total_length":
                    stats["Total_Length"] = stat.text
                elif category == "scaffold_count":
                    stats["Scaffold_Count"] = stat.text
                elif category == "scaffold_n50":
                    stats["Scaffold_N50"] = stat.text
                elif category == "scaffold_l50":
                    stats["Scaffold_L50"] = stat.text
                elif category == "contig_count":
                    stats["Contig_Count"] = stat.text
                elif category == "contig_n50":
                    stats["Contig_N50"] = stat.text
                # GC 含量可能不带 seq_tag="all"，所以单独处理
            
            if category == "gc_percent" or category == "gc_count":
                 stats["GC_Percent"] = stat.text
                 
    except Exception as e:
        pass # 解析错误时静默，返回空值
        
    return stats

def fetch_assembly_details(accession_id):
    """
    通过 Assembly Accession ID (如 GCA_055944425.1) 抓取统计和方法信息
    """
    if pd.isna(accession_id) or not str(accession_id).strip():
        return []

    accession_id = str(accession_id).strip()
    
    # 严格限定查询字段为 Assembly Accession
    term = f'"{accession_id}"[Assembly Accession]'
    
    try:
        # 1. 获取内部 UID
        search_handle = Entrez.esearch(db="assembly", term=term, retmax=1)
        search_record = Entrez.read(search_handle)
        search_handle.close()

        id_list = search_record.get("IdList", [])
        if not id_list:
            return {"Accession": accession_id, "Query_Status": "Not Found"}

        uid = id_list[0]

        # 2. 获取详细摘要
        summary_handle = Entrez.esummary(db="assembly", id=uid)
        summaries = Entrez.read(summary_handle)
        summary_handle.close()

        doc_summaries = summaries['DocumentSummarySet']['DocumentSummary']
        if isinstance(doc_summaries, list):
            doc = doc_summaries[0]
        else:
            doc = doc_summaries

        # 3. 提取常规信息（方法、测序技术、覆盖度）
        assembly_name = doc.get("AssemblyName", "")
        coverage = doc.get("Coverage", "")
        seq_tech = doc.get("SequencingTech", "")
        assembly_method = doc.get("AssemblyMethod", "")

        # 4. 解析 Meta 字段获取详细统计数据 (N50等)
        meta_str = doc.get("Meta", "")
        stats_info = parse_assembly_meta(meta_str)

        result = {
            "Accession": accession_id,
            "AssemblyName": assembly_name,
            "Assembly_Method": assembly_method,
            "Sequencing_Tech": seq_tech,
            "Coverage": coverage,
            "Total_Length": stats_info.get("Total_Length"),
            "Scaffold_Count": stats_info.get("Scaffold_Count"),
            "Scaffold_N50": stats_info.get("Scaffold_N50"),
            "Scaffold_L50": stats_info.get("Scaffold_L50"),
            "Contig_Count": stats_info.get("Contig_Count"),
            "Contig_N50": stats_info.get("Contig_N50"),
            "GC_Percent": stats_info.get("GC_Percent"),
            "Query_Status": "Success"
        }

        return result

    except urllib.error.URLError as e:
        return {"Accession": accession_id, "Query_Status": f"Network Error: {e}"}
    except Exception as e:
        return {"Accession": accession_id, "Query_Status": f"Error: {str(e)}"}

def load_dataframe(filepath):
    ext = os.path.splitext(filepath)[1].lower()
    try:
        try:
            df = pd.read_csv(filepath, sep='\t')
            if len(df.columns) > 1: return df
        except: pass
        if ext == '.csv': return pd.read_csv(filepath)
        elif ext in ['.tsv', '.txt']: return pd.read_csv(filepath, sep='\t')
        elif ext in ['.xlsx', '.xls']: return pd.read_excel(filepath)
        else: return pd.read_csv(filepath)
    except Exception as e:
        print(f"读取数据表失败: {e}")
        sys.exit(1)

def main():
    print("=" * 65)
    print(" " * 8 + "NCBI Assembly 统计信息与组装方法自动爬取工具")
    print("=" * 65)
    
    if Entrez.email == "your_email@example.com":
        print("⚠️ 警告：请在脚本代码第17行修改为你自己的邮箱地址！\n")

    print("请选择操作模式:")
    print("  [1] 手动输入单个 Assembly ID (例如: GCA_055944425.1)")
    print("  [2] 导入包含 Assembly ID 的表格进行批量查询")
    
    mode = input("\n你的选择 (1 或 2): ").strip()
    
    accessions_to_search = []
    file_path = ""
    df = None
    query_col = ""

    if mode == '1':
        single_acc = input("\n请输入要查询的 Assembly ID: ").strip()
        if single_acc:
            accessions_to_search.append(single_acc)
        else:
            print("未输入有效 ID，退出程序。")
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
                col_idx = int(input(f"\n请输入包含【Assembly ID】的列编号 (1-{len(columns)}): ").strip()) - 1
                if 0 <= col_idx < len(columns):
                    query_col = columns[col_idx]
                    # 获取去重且非空的 ID 列表
                    accessions_to_search = df[query_col].dropna().astype(str).str.strip().unique().tolist()
                    break
                print("编号超出范围。")
            except ValueError:
                print("请输入有效的数字。")
    else:
        print("无效的选择，退出程序。")
        sys.exit()

    print(f"\n🚀 开始查询，共 {len(accessions_to_search)} 个 Assembly ID 需要检索...")
    
    results_list = []
    
    for idx, acc in enumerate(accessions_to_search, 1):
        print(f"[{idx}/{len(accessions_to_search)}] 正在检索: {acc} ... ", end="", flush=True)
        
        info = fetch_assembly_details(acc)
        results_list.append(info)
        
        if info.get("Query_Status") == "Success":
            print(f"✅ 成功 (N50: {info.get('Scaffold_N50', '未知')})")
        else:
            print(f"❌ 失败 ({info.get('Query_Status')})")
            
        time.sleep(0.5) # 防止请求过快

    # 保存数据
    results_df = pd.DataFrame(results_list)
    
    if mode == '1':
        save_name = f"{accessions_to_search[0]}_assembly_stats.csv"
        results_df.to_csv(os.path.join(get_base_dir(), save_name), sep=',', index=False, encoding='utf-8-sig')
        print("-" * 65)
        print(f"🎉 处理完成！结果已保存至: {os.path.join(get_base_dir(), save_name)}")
    else:
        # 如果是表格模式，将爬取结果与原表格合并 (左连接)
        final_df = pd.merge(df, results_df, left_on=query_col, right_on="Accession", how="left")
        # 删掉重复的 Accession 列以保持表格整洁
        if "Accession" in final_df.columns and query_col != "Accession":
            final_df = final_df.drop(columns=["Accession"])
            
        base_name = os.path.splitext(os.path.basename(file_path))[0]
        output_path = os.path.join(get_base_dir(), f"{base_name}_with_stats.csv")
        final_df.to_csv(output_path, sep=',', index=False, encoding='utf-8-sig')
        
        print("-" * 65)
        print(f"🎉 全部处理完成！原表格已附带上统计信息。")
        print(f"📄 结果已保存至:\n   {output_path}")

if __name__ == "__main__":
    main()