#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import glob
import sys
import pandas as pd
from goatools.obo_parser import GODag
from goatools.go_enrichment import GOEnrichmentStudy
import urllib.request

# ==========================================================
#  cite.py 功能整合：文件查找 + 文件选择
# ==========================================================

def current_path_function():
    """返回脚本运行目录"""
    if getattr(sys, 'frozen', False):
        return os.path.dirname(sys.executable)
    return os.path.dirname(os.path.abspath(__file__))

def find_file(file_type):
    """在当前目录和子目录搜索文件"""
    current_path = current_path_function()
    search_pattern = os.path.join(current_path, '**', f'*{file_type}')
    file_list = glob.glob(search_pattern, recursive=True)

    print("="*60)
    if len(file_list) >= 1:
        print(f"在当前目录找到 {len(file_list)} 个 {file_type} 文件：")
        for f in file_list:
            print(" -", f)
        return file_list

    # 询问是否搜索父目录
    response = input(f"未找到 {file_type}，是否搜索父目录？ (y/n): ").strip().lower()
    if response not in ["y", "yes"]:
        print("用户取消，程序退出。")
        sys.exit()

    parent_path = os.path.dirname(current_path)
    search_pattern = os.path.join(parent_path, '**', f'*{file_type}')
    file_list = glob.glob(search_pattern, recursive=True)

    if len(file_list) >= 1:
        print(f"在父目录找到 {len(file_list)} 个 {file_type} 文件：")
        for f in file_list:
            print(" -", f)
        return file_list

    print("仍未找到文件，程序退出。")
    sys.exit()

def choose_file(file_list):
    print("="*60)
    print("请选择文件：")
    for i, f in enumerate(file_list, 1):
        print(f"[{i}] {f}")
    idx = int(input("输入编号："))
    return file_list[idx - 1]


# ==========================================================
#  Step 1. 选择 DESeq2 结果 + InterProScan 注释文件
# ==========================================================

def load_files():
    print("\n请选择差异基因文件（包含 gene, padj）：")
    deg_file = choose_file(find_file(".tsv"))

    print("\n请选择 InterProScan 注释文件（包含 GO:xxxx）：")
    ipr_file = choose_file(find_file(".tsv"))

    deg = pd.read_csv(deg_file, sep="\t")
    ipr = pd.read_csv(ipr_file, sep="\t")

    return deg, ipr


# ==========================================================
#  Step 2. 从 InterProScan 解析 GO 注释
# ==========================================================

def parse_go(ipr):
    go_col = None
    for col in ipr.columns:
        if ipr[col].astype(str).str.contains("GO:").any():
            go_col = col
            break

    if go_col is None:
        print("❌ 未在注释文件中找到 GO 信息！")
        sys.exit()

    print(f"✔ 识别到 GO 注释列：{go_col}")

    def extract_go(x):
        if pd.isna(x):
            return []
        return [i.split("(")[0] for i in str(x).split("|") if "GO:" in i]

    ipr["GO_terms"] = ipr[go_col].apply(extract_go)
    return ipr[["gene", "GO_terms"]]


# ==========================================================
#  Step 3. 构建 gene2go 映射
# ==========================================================

def build_gene2go(ipr):
    rows = []
    for gene, gos in zip(ipr["gene"], ipr["GO_terms"]):
        for go in gos:
            rows.append([gene, go])

    g2g = pd.DataFrame(rows, columns=["GeneID", "GO"])
    g2g_file = os.path.join(current_path_function(), "gene2go.tsv")
    g2g.to_csv(g2g_file, sep="\t", index=False)

    print(f"✔ gene2go 已输出：{g2g_file}")

    # 转为 GOATOOLS 格式
    mapping = {}
    for gene, go in zip(g2g["GeneID"], g2g["GO"]):
        mapping.setdefault(gene, set()).add(go)

    return mapping


# ==========================================================
#  Step 4. 提取 DEGs（padj < 0.05）
# ==========================================================

def extract_deg(deg):
    sig = set(deg[deg["padj"] < 0.05]["gene"])
    bg = set(deg["gene"])
    print(f"显著基因数：{len(sig)}，背景基因数：{len(bg)}")
    return sig, bg


# ==========================================================
#  Step 5. GO 富集分析
# ==========================================================

def go_enrich(sig, bg, gene2go):

    # 自动下载 go-basic.obo（若不存在）
    if not os.path.exists("go-basic.obo"):
        print("⏳ 正在下载 go-basic.obo ...")
        url = "http://purl.obolibrary.org/obo/go/go-basic.obo"
        try:
            urllib.request.urlretrieve(url, "go-basic.obo")
            print("✔ 已成功下载 go-basic.obo")
        except Exception as e:
            print("❌ go-basic.obo 下载失败，请手动下载：")
            print("   http://purl.obolibrary.org/obo/go/go-basic.obo")
            print("错误信息：", e)
            sys.exit()

    go_dag = GODag("go-basic.obo")

    goea = GOEnrichmentStudy(
        bg,
        gene2go,
        go_dag,
        propagate_counts=True,
        alpha=0.05,
        methods=['fdr_bh']
    )
    results = goea.run_study(sig)

    df = pd.DataFrame([r.__dict__ for r in results])
    out = os.path.join(current_path_function(), "GO_enrichment_results.tsv")
    df.to_csv(out, sep="\t", index=False)

    print(f"\n🎉 GO 富集完成！结果输出：{out}")


# ==========================================================
#  主程序入口
# ==========================================================

def main():
    print("\n===== 🚀 启动 GO 富集分析 =====\n")

    deg, ipr = load_files()
    ipr = parse_go(ipr)
    gene2go = build_gene2go(ipr)
    sig, bg = extract_deg(deg)
    go_enrich(sig, bg, gene2go)

    print("\n🎉 分析流程全部结束！")


if __name__ == "__main__":
    main()
