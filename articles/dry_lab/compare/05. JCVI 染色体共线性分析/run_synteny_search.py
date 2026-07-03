#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import sys
import glob
import subprocess
from datetime import datetime

'''
🧬 共线性区块识别系统 (修复版：高级工作区与软链接伪装术)
功能：
1. 交互式选择 BED 和 BLAST 文件。
2. 自动创建独立分析工作区 (Workspace)。
3. 🔥 神级操作：强制将 blast 文件伪装成 .last 后缀，完美绕过 JCVI 的序列检查机制。
'''

def get_base_dir():
    return os.path.dirname(os.path.abspath(__file__))

def check_dependency():
    try:
        subprocess.run(["python", "-m", "jcvi"], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    except Exception:
        print("❌ 错误：未检测到 jcvi 环境！请确认已激活 conda 环境。")
        sys.exit(1)

def interactive_select(files, desc):
    if not files:
        print(f"⚠️ 未扫描到任何 {desc}！")
        sys.exit(0)
    
    print(f"\n📂 扫描到以下 {len(files)} 个 {desc}:")
    for i, f in enumerate(files, 1):
        print(f"  [{i}] {os.path.relpath(f, get_base_dir())}")
        
    while True:
        choice = input(f"\n👉 请选择 {desc} (输入编号，或 q 退出): ").strip()
        if choice.lower() == 'q': sys.exit(0)
        try:
            return files[int(choice)-1]
        except Exception:
            print("⚠️ 输入无效，请重新输入。")

def main():
    print("=" * 60)
    print(" 🧬 共线性区块识别系统 (完美伪装版)")
    print("=" * 60)

    check_dependency()
    base_dir = get_base_dir()

    # 1. 扫描可用文件
    bed_files = glob.glob(os.path.join(base_dir, "**", "*.bed"), recursive=True)
    blast_files = glob.glob(os.path.join(base_dir, "**", "*.blast"), recursive=True) + \
                  glob.glob(os.path.join(base_dir, "**", "*.tsv"), recursive=True)

    if len(bed_files) < 2:
        print("❌ 错误：当前目录下找不到足够的 .bed 文件！至少需要 2 个。")
        sys.exit(1)
    if not blast_files:
        print("❌ 错误：找不到全基因组比对结果 (.blast/.tsv) 文件！")
        sys.exit(1)

    # 2. 交互式选择文件
    print("\n" + "-"*40)
    print(" 第一步：配置【物种 A】")
    bed_a = interactive_select(bed_files, "物种 A 的位置文件 (.bed)")
    prefix_a = os.path.basename(bed_a).replace('.bed', '')

    print("\n" + "-"*40)
    print(" 第二步：配置【物种 B】")
    bed_files_remaining = [f for f in bed_files if f != bed_a]
    bed_b = interactive_select(bed_files_remaining, "物种 B 的位置文件 (.bed)")
    prefix_b = os.path.basename(bed_b).replace('.bed', '')

    print("\n" + "-"*40)
    print(" 第三步：配置【全基因组双向比对结果】")
    blast_file = interactive_select(blast_files, "双向比对结果文件 (.blast / .tsv)")

    print("\n" + "="*60)
    print(f" 🚀 即将进行共线性分析: \033[96m{prefix_a}\033[0m ⚔️ \033[96m{prefix_b}\033[0m")
    print("="*60)

    # 3. 创建独立的工作目录
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    work_dir = os.path.join(base_dir, f"Synteny_Result_{prefix_a}_vs_{prefix_b}_{timestamp}")
    os.makedirs(work_dir, exist_ok=True)
    print(f"\n📁 已创建独立分析工作区: {os.path.relpath(work_dir, base_dir)}")

    # 4. 创建软链接
    print("🔗 正在将文件软链接至工作区并强制规范命名...")
    def make_symlink(src_path, target_name):
        abs_src = os.path.abspath(src_path)
        dest_path = os.path.join(work_dir, target_name)
        os.symlink(abs_src, dest_path)

    try:
        make_symlink(bed_a, f"{prefix_a}.bed")
        make_symlink(bed_b, f"{prefix_b}.bed")
        
        # 🔥🔥🔥 核心破局点：不管原文件叫啥，强制伪装成 .last 文件！
        # 这样 JCVI 就会把它当做现成的比对结果直接读取，绝不会去搜寻 .cds 文件
        make_symlink(blast_file, f"{prefix_a}.{prefix_b}.last")
        
    except Exception as e:
        print(f"❌ 创建软链接失败: {e}")
        sys.exit(1)

    # 5. 运行 JCVI 核心命令
    print("\n" + "="*60)
    print("🧠 正在调用 JCVI MCscan 算法寻找 Anchors...")
    print("="*60)
    
    cmd = f"python -m jcvi.compara.catalog ortholog {prefix_a} {prefix_b} --no_strip_names"
    print(f"▶️ 执行命令: {cmd}")
    
    try:
        subprocess.run(cmd, shell=True, check=True, cwd=work_dir)
        
        print("\n" + "="*60)
        print("🎉 恭喜！共线性区块识别成功！(傲娇的 JCVI 被成功制服)")
        print(f"📂 所有的结果文件均已保存在独立文件夹下: \n   \033[92m{os.path.relpath(work_dir, base_dir)}/\033[0m")
        print(f"   💎 核心产物: \033[93m{prefix_a}.{prefix_b}.anchors\033[0m")
        print("="*60 + "\n")
    except subprocess.CalledProcessError:
        print("\n❌ 运行失败，请检查上方 JCVI 的报错信息。")

if __name__ == "__main__":
    main()