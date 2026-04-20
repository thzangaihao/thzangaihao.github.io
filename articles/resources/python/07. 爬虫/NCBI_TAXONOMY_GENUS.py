#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import os
import sys
import time
import pandas as pd
from Bio import Entrez
import urllib.error

# ================= 配置区 =================
# 【重要】请将这里替换为你自己的邮箱！NCBI 要求提供邮箱，否则高频访问可能导致IP被封禁。
Entrez.email = "thzangaihao@outlook.com" 
# ==========================================

def get_base_dir():
    if getattr(sys, 'frozen', False):
        return os.path.dirname(sys.executable)
    return os.path.dirname(os.path.abspath(__file__))

def load_dataframe(filepath):
    ext = os.path.splitext(filepath)[1].lower()
    try:
        # 兼容读取可能用制表符分割的"伪csv"
        try:
            df = pd.read_csv(filepath, sep='\t')
            if len(df.columns) > 1: return df
        except:
            pass
        
        # 常规读取
        if ext == '.csv':
            return pd.read_csv(filepath)
        elif ext in ['.tsv', '.txt']:
            return pd.read_csv(filepath, sep='\t')
        elif ext in ['.xlsx', '.xls']:
            return pd.read_excel(filepath)
        else:
            return pd.read_csv(filepath)
    except Exception as e:
        print(f"读取数据表失败: {e}")
        sys.exit(1)

def fetch_taxonomy(query_name):
    """
    通过分类单元名查询 NCBI Taxonomy，获取属向上的各层级及其 TaxID
    """
    # 初始化目标层级结构
    target_ranks = {
        "superkingdom": "Kingdom(界)", 
        "phylum": "Phylum(门)", 
        "class": "Class(纲)", 
        "order": "Order(目)", 
        "family": "Family(科)", 
        "genus": "Genus(属)"
    }
    
    # 初始化空结果字典
    result = {"Query_Status": "Not Found"}
    for rank_cn in target_ranks.values():
        result[rank_cn] = ""
        result[f"{rank_cn}_TaxID"] = ""

    if pd.isna(query_name) or not str(query_name).strip():
        return result

    query_name = str(query_name).strip()
    
    try:
        # 1. 查询 TaxID
        search_handle = Entrez.esearch(db="taxonomy", term=query_name)
        search_record = Entrez.read(search_handle)
        search_handle.close()

        if not search_record["IdList"]:
            return result

        tax_id = search_record["IdList"][0]

        # 2. 获取详细的 XML 信息
        fetch_handle = Entrez.efetch(db="taxonomy", id=tax_id, retmode="xml")
        fetch_record = Entrez.read(fetch_handle)
        fetch_handle.close()
        
        tax_data = fetch_record[0]
        lineage_ex = tax_data.get("LineageEx", [])
        
        result["Query_Status"] = "Success"
        
        # 3. 提取祖先节点 (Ancestors) 的名称和 TaxID
        for node in lineage_ex:
            rank = node.get("Rank")
            if rank in target_ranks:
                col_name = target_ranks[rank]
                result[col_name] = node.get("ScientificName")
                result[f"{col_name}_TaxID"] = node.get("TaxId")
                
        # 4. 提取当前查询节点本身的信息（如果它是属、科等）
        main_rank = tax_data.get("Rank")
        if main_rank in target_ranks:
            col_name = target_ranks[main_rank]
            result[col_name] = tax_data.get("ScientificName")
            result[f"{col_name}_TaxID"] = tax_data.get("TaxId")

        return result

    except urllib.error.URLError as e:
        print(f"\n⚠️ 网络连接失败: {e}，请检查网络环境。")
        result["Query_Status"] = "Network Error"
        return result
    except Exception as e:
        result["Query_Status"] = f"Error: {str(e)}"
        return result

def main():
    print("=" * 60)
    print(" " * 10 + "NCBI 分类与 TaxID 批量补齐工具 (属及以上)")
    print("=" * 60)
    
    if Entrez.email == "your_email@example.com":
        print("⚠️ 警告：请在脚本代码的第12行修改为你自己的邮箱地址！\n")

    file_path = input("\n1. 请拖拽或输入包含【属名】等分类单元的数据表格路径 (.csv/.tsv/.xlsx): \n").strip().strip('\'"')
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
            col_idx = int(input(f"\n2. 请输入包含【查询名称】的列编号 (1-{len(columns)}): ").strip()) - 1
            if 0 <= col_idx < len(columns):
                query_col = columns[col_idx]
                break
            print("编号超出范围。")
        except ValueError:
            print("请输入有效的数字。")

    print(f"\n3. 开始联网向 NCBI 查询 [{query_col}] 列的详细分类信息与 TaxID...")
    
    taxonomy_results = []
    total_rows = len(df)
    
    for idx, row in df.iterrows():
        query_val = row[query_col]
        print(f"[{idx+1}/{total_rows}] 正在查询: {query_val} ... ", end="", flush=True)
        
        tax_info = fetch_taxonomy(query_val)
        taxonomy_results.append(tax_info)
        
        status = tax_info.get("Query_Status")
        if status == "Success":
            print("✅ 成功")
        else:
            print(f"❌ 失败 ({status})")
            
        time.sleep(0.4) 

    tax_df = pd.DataFrame(taxonomy_results)
    
    # 定义期望的列输出顺序，剔除种(Species)层级
    ordered_cols = [
        "Kingdom(界)", "Kingdom(界)_TaxID", 
        "Phylum(门)", "Phylum(门)_TaxID", 
        "Class(纲)", "Class(纲)_TaxID", 
        "Order(目)", "Order(目)_TaxID", 
        "Family(科)", "Family(科)_TaxID", 
        "Genus(属)", "Genus(属)_TaxID",
        "Query_Status"
    ]
    
    existing_cols = [c for c in ordered_cols if c in tax_df.columns]
    tax_df = tax_df.reindex(columns=existing_cols)
    
    final_df = pd.concat([df, tax_df], axis=1)

    # 强制保存为标准 CSV 格式
    base_name = os.path.splitext(os.path.basename(file_path))[0]
    output_path = os.path.join(get_base_dir(), f"{base_name}_with_taxids.csv")
    
    final_df.to_csv(output_path, sep=',', index=False, encoding='utf-8-sig')
    print("-" * 60)
    print(f"🎉 全部处理完成！包含全层级 TaxID 的 CSV 文件已保存至:\n   {output_path}")

if __name__ == "__main__":
    main()