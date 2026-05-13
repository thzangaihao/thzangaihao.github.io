#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import sys
import glob
import pandas as pd
import urllib.request
from goatools.obo_parser import GODag
from collections import defaultdict

'''
GO 功能追溯数据生成工具
功能：解析大师表中的 GO 编号，关联具体基因 ID，生成用于后续绘图的追溯统计表。
'''

def get_base_dir():
    return os.path.dirname(os.path.abspath(__file__))

def check_and_get_obo():
    """检查并获取 GO 字典"""
    obo_file = os.path.join(get_base_dir(), "go-basic.obo")
    if not os.path.exists(obo_file):
        print("\n⏳ 正在下载 GO 词典 (go-basic.obo)...")
        url = "http://purl.obolibrary.org/obo/go/go-basic.obo"
        urllib.request.urlretrieve(url, obo_file)
    return obo_file

def run_data_generation():
    base_dir = get_base_dir()
    print("\n" + "="*50)
    print(" 📑 GO 功能-基因追溯数据生成器")
    print("="*50)

    # 1. 自动定位大师表 (Master Table V3)
    csv_list = glob.glob(os.path.join(base_dir, "**", "*Master_Table_V3*.csv"), recursive=True)
    if not csv_list:
        print("⚠️ 未找到大师表！请确保已运行 merge_variant_functional_v3.py")
        return
    
    print("\n📂 找到以下数据表:")
    for i, f in enumerate(csv_list, 1):
        print(f"  [{i}] {os.path.relpath(f, base_dir)}")
    
    while True:
        choice = input("\n👉 请选择要处理的数据表编号 (q退出): ").strip().lower()
        if choice == 'q': sys.exit(0)
        try:
            data_file = csv_list[int(choice)-1]
            break
        except:
            print("⚠️ 输入无效。")

    # 2. 加载本体数据库
    obo_file = check_and_get_obo()
    print("\n📚 正在加载 GO 字典结构...")
    go_dag = GODag(obo_file)
    
    # 3. 读取数据
    print(f"📖 正在读取: {os.path.basename(data_file)}")
    df = pd.read_csv(data_file)
    
    # 自动定位 ID 列
    if 'Gene_ID' in df.columns:
        id_col = 'Gene_ID'
    elif 'ID' in df.columns:
        id_col = 'ID'
    else:
        id_col = df.columns[0]

    # 4. 核心映射逻辑：GO -> 基因 ID 集合
    print("🔍 正在建立功能与基因的对应关系...")
    go_to_genes = defaultdict(set)
    
    for _, row in df.iterrows():
        gene_id = str(row[id_col])
        go_str = str(row.get('GO', ''))
        
        if pd.isna(go_str) or go_str == '-' or "No annotation" in go_str:
            continue
            
        # 兼容多种分隔符
        parts = go_str.replace('|', ';').replace(',', ';').split(';')
        for p in parts:
            clean_go = p.split('(')[0].strip() # 移除类似 (InterPro) 的后缀
            if clean_go.startswith("GO:"):
                go_to_genes[clean_go].add(gene_id)

    # 5. 构建统计数据框
    results = []
    for go_id, genes in go_to_genes.items():
        if go_id in go_dag:
            term = go_dag[go_id]
            results.append({
                'GO_ID': go_id,
                'Description': term.name,
                'Category': term.namespace,
                'Count': len(genes),
                'Gene_List': ";".join(sorted(list(genes))) # 使用分号分隔，方便 Excel 进一步处理
            })

    df_stats = pd.DataFrame(results)
    
    # 类别美化命名
    category_map = {
        'biological_process': 'Biological Process (BP)',
        'cellular_component': 'Cellular Component (CC)',
        'molecular_function': 'Molecular Function (MF)'
    }
    df_stats['Category'] = df_stats['Category'].map(category_map)
    df_stats = df_stats.dropna(subset=['Category'])

    # 6. 保存结果
    out_file = "GO_Classification_Traceable_Stats.csv"
    df_stats.sort_values(['Category', 'Count'], ascending=[True, False]).to_csv(out_file, index=False)

    print("\n" + "="*50)
    print(f"🎉 统计表生成成功！")
    print(f"📄 文件名: {out_file}")
    print(f"💡 提示：现在你可以运行绘图脚本来可视化此文件了。")
    print("="*50)

if __name__ == "__main__":
    try:
        run_data_generation()
    except KeyboardInterrupt:
        print("\n🛑 用户强制退出。")