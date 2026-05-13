#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import sys
import glob
import pandas as pd
import urllib.request
import matplotlib.pyplot as plt
import seaborn as sns
from goatools.obo_parser import GODag
from collections import Counter

def check_and_get_obo():
    """检查当前目录下是否有 GO 词典，没有则自动从官网拉取最新版"""
    obo_file = "go-basic.obo"
    if not os.path.exists(obo_file):
        print("\n⏳ 未检测到 GO 词典，正在自动下载最新的 go-basic.obo ...")
        url = "http://purl.obolibrary.org/obo/go/go-basic.obo"
        try:
            urllib.request.urlretrieve(url, obo_file)
            print("✅ 下载完成！")
        except Exception as e:
            print(f"❌ 下载失败: {e}")
            print("请手动在浏览器输入 http://purl.obolibrary.org/obo/go/go-basic.obo 下载并放入该目录。")
            sys.exit(1)
    return obo_file

def select_master_table():
    """交互式选择刚才生成的变异-功能合并大师表"""
    csv_list = glob.glob("**/*Master_Table_V3*.csv", recursive=True)
    if not csv_list:
        print("⚠️ 未找到大师表！请确保目录下有 Final_Master_Table_V3...csv 文件。")
        sys.exit(1)
        
    print("\n📂 扫描到以下大师表:")
    for i, f in enumerate(csv_list, 1):
        print(f"  [{i}] {f}")
        
    while True:
        choice = input(f"\n👉 请选择你的数据表 (输入编号，q退出): ").strip().lower()
        if choice == 'q': sys.exit(0)
        try: 
            return csv_list[int(choice)-1]
        except: 
            print("⚠️ 输入无效。")

def run_analysis():
    print("\n" + "="*50)
    print(" 📊 独立版 GO 功能分类与绘图系统 (兼容版)")
    print("="*50)

    data_file = select_master_table()
    obo_file = check_and_get_obo()

    print("\n📚 正在解析 GO 本体词典结构...")
    go_dag = GODag(obo_file)

    print("🔍 正在从数据表中清洗和提取 GO 编号...")
    df = pd.read_csv(data_file)
    if 'GO' not in df.columns:
        print("❌ 错误：表格中缺失 'GO' 这一列！")
        sys.exit(1)

    all_go_ids = []
    for go_str in df['GO'].dropna():
        if go_str and str(go_str) != '-' and "No annotation" not in str(go_str):
            parts = str(go_str).replace('|', ';').replace(',', ';').split(';')
            for p in parts:
                clean_go = p.split('(')[0].strip() 
                if clean_go.startswith("GO:"):
                    all_go_ids.append(clean_go)

    if not all_go_ids:
        print("❌ 警告：未提取到有效的 GO 编号。")
        sys.exit(1)

    print(f"✅ 共提取到 {len(all_go_ids)} 个 GO 标签，正在映射分类...")
    go_counts = Counter(all_go_ids)
    
    plot_data = []
    for go_id, count in go_counts.items():
        if go_id in go_dag:
            term = go_dag[go_id]
            plot_data.append({
                'GO_ID': go_id,
                'Description': term.name,
                'Category': term.namespace,
                'Count': count
            })

    df_plot = pd.DataFrame(plot_data)
    
    category_map = {
        'biological_process': 'Biological Process (BP)',
        'cellular_component': 'Cellular Component (CC)',
        'molecular_function': 'Molecular Function (MF)'
    }
    df_plot['Category'] = df_plot['Category'].map(category_map)

    # 🚀 核心修复区：完美兼容全版本 Pandas 的 Top10 提取逻辑
    df_plot = df_plot.dropna(subset=['Category']) # 剔除没匹配上三大类的异常值
    # 1. 先全局按 Count 降序排
    df_sorted = df_plot.sort_values('Count', ascending=False)
    # 2. 对 Category 分组，直接取前 10 行
    df_top = df_sorted.groupby('Category').head(10)
    # 3. 再次排序，让画图时的同类别聚在一起，且条带长度按从大到小排列
    df_top = df_top.sort_values(by=['Category', 'Count'], ascending=[True, False])

    print("🎨 正在渲染科研级柱状图...")
    sns.set_theme(style="whitegrid")
    plt.figure(figsize=(12, 8)) # 稍微加宽一点防止英文名字被截断

    palette = {
        'Biological Process (BP)': '#e74c3c', 
        'Cellular Component (CC)': '#2ecc71', 
        'Molecular Function (MF)': '#3498db'  
    }

    sns.barplot(
        x='Count', 
        y='Description', 
        hue='Category', 
        data=df_top, 
        dodge=False, 
        palette=palette
    )

    plt.title('GO Classification of Mutated Target Genes', fontsize=16, fontweight='bold', pad=20)
    plt.xlabel('Number of Mutated Genes', fontsize=12, fontweight='bold')
    plt.ylabel('Gene Ontology (GO) Term', fontsize=12, fontweight='bold')
    plt.legend(title='GO Category', loc='lower right', bbox_to_anchor=(1, 0.1))
    plt.tight_layout()

    plt.savefig("GO_Classification_Barplot.pdf", dpi=300)
    plt.savefig("GO_Classification_Barplot.png", dpi=300)
    df_plot.sort_values(['Category', 'Count'], ascending=[True, False]).to_csv("GO_Classification_Stats.csv", index=False)

    print("\n" + "="*50)
    print("🎉 分析与绘图圆满完成！")
    print(f"📄 统计明细表: GO_Classification_Stats.csv")
    print(f"📊 论文用 PDF: GO_Classification_Barplot.pdf")
    print("="*50)

if __name__ == "__main__":
    run_analysis()