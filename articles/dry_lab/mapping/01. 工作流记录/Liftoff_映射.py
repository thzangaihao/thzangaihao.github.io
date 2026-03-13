#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import sys
import glob
import subprocess
from datetime import datetime

'''
Liftoff 基因组注释转移流水线 (多任务交互版)
功能：以某个标准基因组及其 GFF 为模板，将注释精确投射到多个自建草图上。
工具：Liftoff
'''

def get_base_dir():
    return os.path.dirname(os.path.abspath(__file__))

def check_env():
    """检查环境中是否安装了 liftoff"""
    if subprocess.call(['which', 'liftoff'], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL) != 0:
        print("❌ 错误: 未找到 liftoff 工具！请确保已激活 annotation_env 环境。")
        sys.exit(1)

def interactive_select(files, desc, single_choice=False):
    """通用的交互式选择逻辑 (支持单选或多选)"""
    if not files:
        print(f"⚠️ 未找到任何 {desc}！请检查文件是否放在了当前目录或子目录中。")
        return []
    
    print(f"\n📂 扫描到以下 {len(files)} 个 {desc}:")
    for i, f in enumerate(files, 1):
        print(f"  [{i}] {os.path.relpath(f, get_base_dir())}")
        
    while True:
        prompt = f"\n👉 请选择 {desc} (输入编号"
        if single_choice:
            prompt += ", 仅限单选): "
        else:
            prompt += ", 多选如1,3 或 1-3, 全部选all): "
            
        choice = input(prompt).strip().lower()
        if choice == 'q': sys.exit(0)
        
        if not single_choice and choice == 'all': return files
        
        try:
            selected = []
            parts = choice.replace(' ', '').split(',')
            for part in parts:
                if '-' in part:
                    start, end = map(int, part.split('-'))
                    selected.extend(files[start-1:end])
                else:
                    selected.append(files[int(part)-1])
                    
            if single_choice and len(selected) > 1:
                print("⚠️ 此项只能选择一个文件，请重新输入！")
                continue
                
            return list(set(selected))
        except:
            print("⚠️ 输入格式错误，请重新选择。")

def run_pipeline():
    base_dir = get_base_dir()
    
    print("\n" + "="*50)
    print("🧬 欢迎使用 Liftoff 基因组注释转移系统")
    print("="*50)

    # 1. 选择【模板基因组】和【模板GFF3】 (即你的 Gamsii 参考)
    print("\n--- 第一步：选择模板 (Reference) ---")
    all_fastas = glob.glob(os.path.join(base_dir, "**", "*.fasta"), recursive=True) + \
                 glob.glob(os.path.join(base_dir, "**", "*.fna"), recursive=True)
                 
    all_gffs = glob.glob(os.path.join(base_dir, "**", "*.gff3"), recursive=True) + \
               glob.glob(os.path.join(base_dir, "**", "*.gff"), recursive=True)

    ref_fasta = interactive_select(all_fastas, "【模板】基因组序列 (.fasta / .fna)", single_choice=True)
    if not ref_fasta: return
    ref_fasta = ref_fasta[0]

    ref_gff = interactive_select(all_gffs, "【模板】注释文件 (.gff / .gff3)", single_choice=True)
    if not ref_gff: return
    ref_gff = ref_gff[0]

    # 2. 选择【目标靶基因组】 (即你需要注释的 CK-3-5, CK-5-4, CK-7-9)
    print("\n--- 第二步：选择目标靶序列 (Target) ---")
    # 过滤掉刚刚选为模板的那个文件，免得自己映射自己产生混乱
    target_candidates = [f for f in all_fastas if f != ref_fasta]
    target_fastas = interactive_select(target_candidates, "【待注释】的自建基因组序列 (.fasta)")
    if not target_fastas: return

    # 3. 设置多线程参数 (针对 N100 优化)
    print("\n--- 第三步：参数配置 ---")
    threads = input("👉 请输入分配的 CPU 核心数 (默认 4): ") or "4"

    # 4. 创建输出目录
    timestamp = datetime.now().strftime("%Y%m%d_%H%M")
    out_root = os.path.join(base_dir, f"07_Annotation_Liftoff_{timestamp}")
    os.makedirs(out_root, exist_ok=True)
    print(f"\n📂 所有的注释结果将保存在: {out_root}")

    # 5. 开始批量映射
    total = len(target_fastas)
    for i, target in enumerate(target_fastas, 1):
        target_name = os.path.basename(target).replace('.fasta', '').replace('.fna', '')
        output_gff = os.path.join(out_root, f"{target_name}_annotated.gff3")
        
        print(f"\n🚀 [任务 {i}/{total}] 正在注释: {target_name}")
        print(f"  - 正在使用 Liftoff 进行同源序列比对与坐标投射 (约需 5-15 分钟)...")
        
        # 核心命令：liftoff -g 模板GFF -o 输出GFF 目标序列 模板序列 -p 线程数
        cmd = f"liftoff -g {ref_gff} -o {output_gff} {target} {ref_fasta} -p {threads}"
        
        try:
            subprocess.run(cmd, shell=True, check=True)
            print(f"✅ {target_name} 注释成功！生成文件: {output_gff}")
        except subprocess.CalledProcessError:
            print(f"❌ {target_name} 注释失败，请检查终端报错信息。")

    print(f"\n🎉 所有基因组注释均已完成。")

if __name__ == "__main__":
    try:
        check_env()
        run_pipeline()
    except KeyboardInterrupt:
        print("\n🛑 用户强制退出。")