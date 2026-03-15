#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import sys
import glob
import subprocess
from datetime import datetime

'''
InterProScan 交互式运行工具
'''

# === 配置区：请确保路径正确 ===
IPRSCAN_PATH = "/home/thz/software/interproscan-5.76-107.0/interproscan.sh"

def get_base_dir():
    return os.path.dirname(os.path.abspath(__file__))

def interactive_select(files, desc):
    """交互式单选逻辑"""
    if not files:
        print(f"⚠️ 未找到任何 {desc}！")
        return None
    print(f"\n📂 扫描到以下 {len(files)} 个 {desc}:")
    for i, f in enumerate(files, 1):
        print(f"  [{i}] {os.path.relpath(f, get_base_dir())}")
    while True:
        choice = input(f"\n👉 请选择输入文件 (输入编号，q退出): ").strip().lower()
        if choice == 'q': sys.exit(0)
        try: return files[int(choice)-1]
        except: print("⚠️ 输入无效。")

def run_iprscan():
    base_dir = get_base_dir()
    print("\n" + "="*50)
    print("InterProScan")
    print("="*50)

    # 1. 自动寻找蛋白序列文件
    faa_list = glob.glob(os.path.join(base_dir, "**", "*.faa"), recursive=True) + \
               glob.glob(os.path.join(base_dir, "**", "*.fasta"), recursive=True)
    input_faa = interactive_select(faa_list, "待分析蛋白序列 (.faa)")
    if not input_faa: return

    # 2. 配置分析资源
    print("\n💻 --- 运算资源配置 ---")
    cpu_cores = input("👉 使用多少个 CPU 核心？(默认 2): ") or "2"
    
    # 3. 配置输出格式
    print("\n📝 --- 输出格式选择 ---")
    print("  [1] TSV (表格, 最常用)")
    print("  [2] XML (后续脚本提取数据必选)")
    print("  [3] GFF3 (查看基因组分布)")
    formats_choice = input("👉 请选择输出格式 (默认 1,2): ") or "1,2"
    format_map = {"1": "tsv", "2": "xml", "3": "gff3"}
    selected_formats_list = [format_map[c] for c in formats_choice.replace(' ', '').split(',') if c in format_map]
    formats_str = ",".join(selected_formats_list)

    # 4. 创建输出目录
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    out_dir = os.path.abspath(os.path.join(base_dir, f"_IPRSCAN_Results_{timestamp}"))
    os.makedirs(out_dir, exist_ok=True)

    # 5. 构建命令 (修复核心逻辑)
    # 如果只有一种格式，可以使用 -o 指定文件名
    # 如果有多种格式，必须使用 -d 指定目录，文件会自动按输入文件名命名
    cmd = [
        IPRSCAN_PATH,
        "-i", f"'{input_faa}'",
        "-f", formats_str,
        "-cpu", cpu_cores,
        "-goterms",
        "-pa",
        "-dp"
    ]

    if len(selected_formats_list) == 1:
        # 单一格式，指定具体输出文件
        out_file = os.path.join(out_dir, f"Functional_Annotation.{selected_formats_list[0]}")
        cmd.extend(["-o", f"'{out_file}'"])
    else:
        # 多种格式，指定输出目录
        cmd.extend(["-d", f"'{out_dir}'"])

    full_cmd = " ".join(cmd)
    
    print("\n" + "="*50)
    print(f"🚀 准备启动分析！")
    print(f"📂 输入文件: {input_faa}")
    print(f"🧵 使用核心: {cpu_cores}")
    print(f"📊 输出格式: {formats_str}")
    print(f"🛠️  执行方式: {'目录输出模式' if len(selected_formats_list) > 1 else '单文件输出模式'}")
    print("="*50)
    
    confirm = input("\n⚠️ 是否开始？(y/n): ").strip().lower()
    if confirm != 'y':
        print("🛑 已取消。")
        return

    print(f"\n🔥 正在运行，请勿关闭终端...\n")
    try:
        process = subprocess.Popen(full_cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
        for line in process.stdout:
            print(line, end='')
        process.wait()
        
        if process.returncode == 0:
            print(f"\n🎉 分析成功完成！")
            print(f"📂 结果文件存放在：{out_dir}")
        else:
            print(f"\n❌ 分析失败，错误代码：{process.returncode}")
            print("💡 建议：检查 Java 环境或内存占用。")
    except Exception as e:
        print(f"\n❌ 运行出错：{str(e)}")

if __name__ == "__main__":
    try:
        run_iprscan()
    except KeyboardInterrupt:
        print("\n🛑 用户强制退出。")