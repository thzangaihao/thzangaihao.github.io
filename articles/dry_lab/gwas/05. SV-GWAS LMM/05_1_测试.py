#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import sys
import glob
import subprocess
import pandas as pd
import re

# ==============================================================================
# 1. Cite2 交互逻辑模块 (完全复用您的代码)
# ==============================================================================
def get_base_dir():
    """获取脚本所在目录，作为搜索的起点"""
    if getattr(sys, 'frozen', False):
        return os.path.dirname(sys.executable)
    return os.path.dirname(os.path.abspath(__file__))

def find_files(ext, path=None):
    """查找指定扩展名文件"""
    if not ext.startswith('.'):
        ext = '.' + ext
    if path is None:
        path = get_base_dir()
    search_pattern = os.path.join(path, '**', f'*{ext}')
    # 排除之前的输出目录，防止混淆
    all_files = sorted(glob.glob(search_pattern, recursive=True))
    # 简单的过滤，排除非原始数据
    return [f for f in all_files if "cleaned" not in f and "08_" not in f]

def choose_file(files, desc="文件"):
    """支持多种选取模式"""
    if not files:
        print(f"提示：在当前目录及其子目录下未找到任何 {desc}")
        return []

    print(f"\n2. 找到 {len(files)} 个 {desc}:")
    for i, f in enumerate(files, 1):
        print(f"  [{i}] {os.path.basename(f)}  (路径: {os.path.relpath(f, get_base_dir())})")

    while True:
        try:
            prompt = (
                f"\n3. 请输入欲操作的 {desc} 编号\n"
                f" (通常请输入单个编号，如 1): \n"
            )
            user_input = input(prompt).strip().lower()

            if not user_input: continue

            selected_indices = set()
            if user_input in ['all', 'a']:
                selected_indices = set(range(len(files)))
            else:
                parts = user_input.split(',')
                for part in parts:
                    part = part.strip()
                    if '-' in part:
                        start_str, end_str = part.split('-')
                        for idx in range(int(start_str) - 1, int(end_str)):
                            if 0 <= idx < len(files): selected_indices.add(idx)
                    else:
                        idx = int(part) - 1
                        if 0 <= idx < len(files): selected_indices.add(idx)
            
            if not selected_indices:
                print("未匹配到有效编号，请重新输入。")
                continue

            selected_paths = [files[i] for i in sorted(list(selected_indices))]
            print(f"\n4. 已选择: {os.path.basename(selected_paths[0])}")
            if len(selected_paths) > 1:
                print("注意：本脚本用于诊断，将只使用第一个选择的文件。")
            
            return selected_paths

        except ValueError:
            print("输入错误：请输入数字编号。")

def make_sure(action_name="执行操作"):
    response = input(f"\n5. 确认{action_name}? (y/n): ").strip().lower()
    return response in ['y', 'yes']

# ==============================================================================
# 2. 核心诊断逻辑
# ==============================================================================

def get_vcf_samples(vcf_path):
    """读取 VCF 中的 ID"""
    print(f"\n正在读取基因型 ID: {os.path.basename(vcf_path)} ...")
    try:
        cmd = ["bcftools", "query", "-l", vcf_path]
        res = subprocess.run(cmd, capture_output=True, text=True, check=True)
        ids = res.stdout.strip().split()
        print(f"-> 成功读取 {len(ids)} 个 VCF 样本 ID。")
        return set(ids)
    except Exception as e:
        print(f"[错误] bcftools 读取失败: {e}")
        sys.exit(1)

def diagnose_mismatch(vcf_path, tsv_path):
    """对比 ID 并给出修改建议"""
    # 1. 获取 ID 集合
    vcf_ids = get_vcf_samples(vcf_path)
    
    print(f"正在读取表型 ID: {os.path.basename(tsv_path)} ...")
    # 强制读取为字符串，避免 SF01 变成 1
    df = pd.read_csv(tsv_path, sep='\t', dtype=str)
    
    # 简单的列选择交互
    cols = list(df.columns)
    id_col = cols[0]
    if len(cols) > 1:
        print("\n表型文件包含以下列，请确认 ID 列:")
        for i, c in enumerate(cols, 1):
            print(f"  [{i}] {c}")
        try:
            sel = input(f"请输入列编号 (默认 1): ").strip()
            if sel: id_col = cols[int(sel)-1]
        except:
            pass
            
    # 提取表型 ID (原始)
    raw_pheno_ids = df[id_col].dropna().tolist()
    pheno_ids_set = set(raw_pheno_ids)
    
    # 2. 计算交集与差集
    common = vcf_ids.intersection(pheno_ids_set)
    missing_in_vcf = pheno_ids_set - vcf_ids  # 这就是那 88 个被剔除的
    
    # 3. 打印基础报告
    print("\n" + "="*50)
    print(f"【ID 匹配诊断报告】")
    print(f"VCF 样本总数    : {len(vcf_ids)}")
    print(f"表型 样本总数   : {len(pheno_ids_set)}")
    print("-" * 50)
    print(f"完全匹配 (保留) : {len(common)}")
    print(f"匹配失败 (剔除) : {len(missing_in_vcf)}  <-- 重点分析这里")
    print("="*50)
    
    if len(missing_in_vcf) == 0:
        print("所有表型样本都匹配上了！无需修复。")
        return

    # 4. 深度侦探模式：分析为什么匹配不上
    print("\n【深度侦探：为何这 88 个样本对不上？】")
    
    # 取几个典型的失败案例
    sample_examples = sorted(list(missing_in_vcf))[:5]
    vcf_examples = sorted(list(vcf_ids))[:5]
    
    print(f"表型 ID 示例: {sample_examples}")
    print(f"VCF  ID 示例: {vcf_examples}")
    print("-" * 30)
    
    # --- 假设检验 1: 空格问题 ---
    stripped_ids = {s.strip() for s in raw_pheno_ids}
    if len(stripped_ids.intersection(vcf_ids)) > len(common):
        print("[发现] 表型 ID 中可能包含看不见的空格！")
        print("  -> 建议：在读取时去除首尾空格。")
        
    # --- 假设检验 2: 大小写问题 ---
    lower_pheno = {s.lower(): s for s in missing_in_vcf}
    lower_vcf = {s.lower(): s for s in vcf_ids}
    common_lower = set(lower_pheno.keys()).intersection(lower_vcf.keys())
    
    if common_lower:
        print(f"[发现] 存在大小写不一致的情况！(共 {len(common_lower)} 个)")
        ex = list(common_lower)[0]
        print(f"  -> 例子: 表型 '{lower_pheno[ex]}' vs VCF '{lower_vcf[ex]}'")
    
    # --- 假设检验 3: 前缀/后缀问题 (SF1 vs SF01) ---
    # 尝试提取数字部分进行比对
    def extract_num(s):
        nums = re.findall(r'\d+', s)
        return int(nums[0]) if nums else None
        
    vcf_num_map = {extract_num(s): s for s in vcf_ids if extract_num(s) is not None}
    
    recoverable_count = 0
    print("\n[尝试智能推断对应关系]:")
    for pid in sample_examples:
        p_num = extract_num(pid)
        if p_num in vcf_num_map:
            match_vid = vcf_num_map[p_num]
            print(f"  表型 '{pid}' 可能对应 VCF '{match_vid}'")
            recoverable_count += 1
            
    if recoverable_count > 0:
        print(f"\n  -> 看来主要是 ID 格式 (如前导零、前缀) 不一致。")
        print("  -> 建议：修改 Excel/TSV 中的 ID，使其与 VCF 完全一致。")
    else:
        print("\n  -> 未发现明显规律，可能是完全不同的样本命名，或者 VCF 里确实没有这些样本。")

# ==============================================================================
# 主函数
# ==============================================================================
def main():
    print("==============================================")
    print("   ID 匹配侦探工具 (集成 cite2 交互)")
    print("==============================================")
    
    # 1. 选择 VCF
    print("\n>>> 第一步: 选择基因型文件 (VCF/BCF)")
    vcfs = find_files('.gz') + find_files('.bcf')
    # 优先找合并后的文件
    vcfs = sorted(list(set(vcfs)))
    sel_vcf = choose_file(vcfs, "基因型文件")
    if not sel_vcf: return
    
    # 2. 选择 TSV
    print("\n>>> 第二步: 选择表型文件 (TSV)")
    tsvs = find_files('.tsv')
    sel_tsv = choose_file(tsvs, "表型文件")
    if not sel_tsv: return
    
    # 3. 执行诊断
    if make_sure("开始诊断"):
        diagnose_mismatch(sel_vcf[0], sel_tsv[0])

if __name__ == "__main__":
    main()