#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import sys
import glob
import subprocess
from datetime import datetime

'''
FASTA 序列深度清洗工具 (SeqKit 版)
功能：从组装结果中筛选高质量 Scaffold，剔除短片段。
特性：支持多任务交互、时间戳文件夹输出、SeqKit 强力驱动。
'''

def get_base_dir():
    return os.path.dirname(os.path.abspath(__file__))

def check_env():
    """检查 SeqKit 是否安装"""
    if subprocess.call(['which', 'seqkit'], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL) != 0:
        print("❌ 错误: 未找到 seqkit！请先执行: conda install -c bioconda seqkit")
        sys.exit(1)

def choose_fasta_tasks():
    """动态扫描并选择组装好的 FASTA 文件"""
    base_dir = get_base_dir()
    # 扫描之前组装脚本产生的 ref_draft.fasta
    search_pattern = os.path.join(base_dir, '04_Assembly_Result_*', '*_ref_draft.fasta')
    fasta_files = sorted(glob.glob(search_pattern, recursive=True))
    
    if not fasta_files:
        print("⚠️ 提示：未找到待清洗的组装结果文件 (*_ref_draft.fasta)！")
        return []

    print(f"\n📂 扫描到以下 {len(fasta_files)} 个组装草图:")
    for i, f in enumerate(fasta_files, 1):
        rel_path = os.path.relpath(f, base_dir)
        print(f"  [{i}] {rel_path}")

    while True:
        try:
            choice = input("\n👉 请选择要清洗的样本编号 (单选:1, 多选:1,3, 范围:1-3, 全部:all, 退出:q): ").strip().lower()
            if choice == 'q': sys.exit(0)
            
            selected_files = []
            if choice == 'all':
                selected_files = fasta_files
            else:
                parts = choice.replace(' ', '').split(',')
                indices = set()
                for part in parts:
                    if '-' in part:
                        start, end = map(int, part.split('-'))
                        indices.update(range(start, end + 1))
                    else:
                        indices.add(int(part))
                for idx in sorted(list(indices)):
                    if 1 <= idx <= len(fasta_files):
                        selected_files.append(fasta_files[idx-1])
            
            if selected_files: return selected_files
            print("⚠️ 未选中有效样本。")
        except ValueError:
            print("⚠️ 输入格式有误。")

def run_cleaning(fasta_path, output_dir):
    """使用 SeqKit 执行清洗"""
    sample_name = os.path.basename(fasta_path).replace('_ref_draft.fasta', '')
    cleaned_fasta = os.path.join(output_dir, f"{sample_name}_cleaned_v1.fasta")
    
    # --- 核心逻辑：使用 seqkit 过滤 ---
    # -m 1000: 仅保留长度 >= 1000 bp 的序列
    # -g: 压缩输出 (可选，这里我们输出解压的以便后续建索引)
    print(f"\n✨ 正在清洗样本: {sample_name}")
    
    # 统计清洗前的数据
    print("  📊 清洗前统计:")
    subprocess.run(f"seqkit stats {fasta_path}", shell=True)
    
    # 执行过滤
    cmd = f"seqkit seq -g -m 1000 {fasta_path} -o {cleaned_fasta}"
    subprocess.run(cmd, shell=True, check=True)
    
    # 统计清洗后的数据
    print("  ✅ 清洗后统计:")
    subprocess.run(f"seqkit stats {cleaned_fasta}", shell=True)
    return cleaned_fasta

if __name__ == "__main__":
    try:
        check_env()
        tasks = choose_fasta_tasks()
        
        if tasks:
            # 1. 创建带有当前时间戳的输出目录 (如 202603111923)
            timestamp = datetime.now().strftime("%Y%m%d%H%M")
            output_dir = os.path.join(get_base_dir(), f"05_Cleaned_Ref_{timestamp}")
            os.makedirs(output_dir, exist_ok=True)
            
            print(f"\n📂 所有清洗后的结果将保存在: {output_dir}")
            
            # 2. 依次执行清洗任务
            for fasta in tasks:
                run_cleaning(fasta, output_dir)
                
            print(f"\n🎉 恭喜！{len(tasks)} 个参考基因组已精炼完毕，准备进行变异比对！\n")
            
    except KeyboardInterrupt:
        print("\n🛑 用户中断。")
        sys.exit(0)