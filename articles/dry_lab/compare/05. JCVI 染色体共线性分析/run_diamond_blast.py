#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import sys
import glob
import subprocess
import time
import multiprocessing

'''
💎 全基因组蛋白序列双向比对自动流程 (Diamond 终极防撞车版)
功能：
1. 强制类型校验：严禁使用 DNA 序列建库，防止报错。
2. 算力动态分配：交互式设置 CPU 线程数，榨干服务器性能。
3. 🌟 新增核心功能：自动为 BLAST 结果的第一列和第二列注入物种前缀，彻底解决 Liftoff 同名 ID 冲突！
4. 自动合并比对结果，无缝衔接 JCVI/MCScanX。
'''

def get_base_dir():
    return os.path.dirname(os.path.abspath(__file__))

def check_dependency():
    """检查 diamond 环境"""
    if subprocess.call(['which', 'diamond'], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL) != 0:
        print("❌ 错误：未找到 diamond！")
        print("💡 请在终端执行: conda install -c bioconda diamond")
        sys.exit(1)

def is_protein_fasta(filepath):
    """快速检查文件是否包含蛋白序列，防止用户误传基因组 DNA"""
    try:
        with open(filepath, 'r') as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith('>'):
                    continue
                protein_chars = set('DEFHIKLMNPQRSTVWY')
                if any(char in protein_chars for char in line.upper()):
                    return True
                return False
    except Exception:
        return False
    return False

def interactive_select(files, desc):
    """交互式单选逻辑"""
    if not files:
        print(f"⚠️ 未在当前目录下找到任何候选 {desc}！")
        sys.exit(0)
    
    valid_files = [f for f in files if f.endswith(('.faa', '.pep')) or is_protein_fasta(f)]
    
    if not valid_files:
        print(f"❌ 严重错误：你提供的文件中似乎没有【蛋白质序列】！")
        print("💡 共线性分析必须使用翻译后的氨基酸序列（.faa / .pep）。")
        sys.exit(1)

    print(f"\n📂 扫描到以下 {len(valid_files)} 个 {desc}:")
    for i, f in enumerate(valid_files, 1):
        print(f"  [{i}] {os.path.relpath(f, get_base_dir())}")
        
    while True:
        choice = input(f"\n👉 请选择 {desc} (输入编号，或 q 退出): ").strip().lower()
        if choice == 'q': sys.exit(0)
        try:
            return valid_files[int(choice)-1]
        except Exception:
            print("⚠️ 输入无效，请重新输入。")

def get_prefix(filepath):
    """智能猜测物种前缀"""
    basename = os.path.basename(filepath)
    guess = basename.split('.')[0].split('_')[0]
    
    print(f"\n💡 猜测该物种的前缀代号为: \033[92m{guess}\033[0m")
    print("⚠️ 警告：此代号必须与你转换生成的 .bed 文件前缀【完全一致】！")
    user_input = input(f"👉 请输入确认的前缀 (直接回车使用 '{guess}'): ").strip()
    return user_input if user_input else guess

def set_threads():
    """交互式设置运行线程数"""
    max_cores = multiprocessing.cpu_count()
    default_cores = max(1, int(max_cores * 0.8))

    print("\n" + "-"*40)
    print(" ⚡ 算力分配配置")
    print(f"💡 检测到系统共有 \033[93m{max_cores}\033[0m 个 CPU 核心可用。")
    
    while True:
        user_input = input(f"👉 请输入 Diamond 调用的线程数 (直接回车使用推荐值 {default_cores}): ").strip()
        if not user_input:
            return default_cores
        try:
            threads = int(user_input)
            if 0 < threads <= max_cores: return threads
            if threads > max_cores:
                confirm = input(f"⚠️ 警告：设定的线程数 ({threads}) 超过系统核心数 ({max_cores})！坚持使用？(y/n): ").strip().lower()
                if confirm == 'y': return threads
        except ValueError:
            print("⚠️ 请输入一个有效的整数！")

def run_cmd(cmd, desc):
    """执行系统命令并处理异常"""
    print(f"\n▶️ 正在执行: {desc}")
    print(f"   命令: {cmd}")
    start_time = time.time()
    try:
        subprocess.run(cmd, shell=True, check=True)
        elapsed = time.time() - start_time
        print(f"✅ {desc} 完成！耗时: {elapsed:.2f} 秒")
    except subprocess.CalledProcessError:
        print(f"\n❌ 错误：{desc} 失败，请检查上方报错信息。")
        sys.exit(1)

def main():
    print("=" * 60)
    print(" 💎 全基因组蛋白双向比对自动流程 (终极前缀注入版)")
    print("=" * 60)

    check_dependency()
    base_dir = get_base_dir()

    all_files = glob.glob(os.path.join(base_dir, "**", "*.faa"), recursive=True) + \
                glob.glob(os.path.join(base_dir, "**", "*.pep"), recursive=True) + \
                glob.glob(os.path.join(base_dir, "**", "*.fasta"), recursive=True)

    # 1. 选择文件与确认前缀
    print("\n" + "-"*40)
    print(" 第一步：配置【物种 A】(参考库构建方)")
    file_a = interactive_select(all_files, "物种 A 的蛋白序列文件")
    prefix_a = get_prefix(file_a)

    print("\n" + "-"*40)
    print(" 第二步：配置【物种 B】(查询比对方)")
    file_b = interactive_select([f for f in all_files if f != file_a], "物种 B 的蛋白序列文件")
    prefix_b = get_prefix(file_b)

    # 2. 线程设置与前缀注入确认
    threads = set_threads()
    
    print("\n" + "-"*40)
    print(" 🛡️ 格式安全配置")
    print("强烈建议为比对结果注入前缀，以匹配 BED 文件，防止 Liftoff 同名冲突。")
    inject_prefix = input(f"👉 是否为输出的基因 ID 注入 '{prefix_a}_' 和 '{prefix_b}_' 前缀？(y/n, 默认 y): ").strip().lower()
    do_inject = False if inject_prefix == 'n' else True

    print("\n" + "="*60)
    print(f"🚀 即将启动: [{prefix_a}] ⚔️ [{prefix_b}] (使用 {threads} 线程)")
    print("="*60)

    # 3. 运行构建与比对
    run_cmd(f"diamond makedb --threads {threads} --in '{file_a}' -d {prefix_a}", f"构建物种 A ({prefix_a}) 数据库")
    run_cmd(f"diamond makedb --threads {threads} --in '{file_b}' -d {prefix_b}", f"构建物种 B ({prefix_b}) 数据库")

    blast_a_vs_b = os.path.join(base_dir, f"{prefix_a}_{prefix_b}.blast")
    blast_b_vs_a = os.path.join(base_dir, f"{prefix_b}_{prefix_a}.blast")
    
    diamond_params = f"--threads {threads} --evalue 1e-5 --outfmt 6 --max-target-seqs 5"
    
    run_cmd(f"diamond blastp -q '{file_a}' -d '{prefix_b}' -o '{blast_a_vs_b}' {diamond_params}", 
            f"并发比对: {prefix_a} ➔ {prefix_b}")
    run_cmd(f"diamond blastp -q '{file_b}' -d '{prefix_a}' -o '{blast_b_vs_a}' {diamond_params}", 
            f"并发比对: {prefix_b} ➔ {prefix_a}")

    # 4. 合并与智能前缀注入
    merged_blast = os.path.join(base_dir, f"{prefix_a}.{prefix_b}.blast")
    print(f"\n▶️ 正在执行: 合并双向比对结果 & 数据清洗")
    
    try:
        with open(merged_blast, 'w') as outfile:
            # 处理 A vs B (Query=A, Subject=B)
            with open(blast_a_vs_b) as infile:
                for line in infile:
                    if not line.strip(): continue
                    parts = line.strip().split('\t')
                    if do_inject and len(parts) >= 2:
                        if not parts[0].startswith(f"{prefix_a}_"): parts[0] = f"{prefix_a}_{parts[0]}"
                        if not parts[1].startswith(f"{prefix_b}_"): parts[1] = f"{prefix_b}_{parts[1]}"
                    outfile.write('\t'.join(parts) + '\n')
            
            # 处理 B vs A (Query=B, Subject=A)
            with open(blast_b_vs_a) as infile:
                for line in infile:
                    if not line.strip(): continue
                    parts = line.strip().split('\t')
                    if do_inject and len(parts) >= 2:
                        if not parts[0].startswith(f"{prefix_b}_"): parts[0] = f"{prefix_b}_{parts[0]}"
                        if not parts[1].startswith(f"{prefix_a}_"): parts[1] = f"{prefix_a}_{parts[1]}"
                    outfile.write('\t'.join(parts) + '\n')
                    
        print(f"✅ 数据处理成功！生成了带物种标签的最终比对文件: {merged_blast}")
    except Exception as e:
        print(f"❌ 合并失败: {e}")
        sys.exit(1)

    print("\n🧹 正在清理单向比对临时文件...")
    os.remove(blast_a_vs_b)
    os.remove(blast_b_vs_a)

    print("\n" + "="*60)
    print("🎉 Diamond 全基因组双向比对火力全开完成！")
    print(f"📂 核心产出文件：")
    print(f"   💎 {merged_blast}")
    print("="*60 + "\n")

if __name__ == "__main__":
    main()