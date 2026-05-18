#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import sys
import glob
import subprocess

'''
📄 JCVI .simple 块文件交互式生成工具
功能：
1. 自动扫描当前目录及子目录下的所有 .anchors / .lifted.anchors 文件。
2. 交互式选择数据源，并动态微调 --minspan 过滤阈值。
3. 调用 JCVI 核心过滤算法，生成绘制染色体连线图（Karyotype）必需的 .simple 文件。
'''

def get_base_dir():
    return os.getcwd()

def main():
    print("=" * 60)
    print(" 📄 JCVI .simple 块文件交互式生成工具")
    print("=" * 60)

    # 1. 检查 JCVI 环境
    try:
        subprocess.run(["python", "-m", "jcvi"], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    except Exception:
        print("❌ 错误：未检测到 jcvi 环境！请确认已激活 conda 环境。")
        sys.exit(1)

    base_dir = get_base_dir()

    # 2. 智能扫描所有的 anchors 文件
    anchor_files = glob.glob(os.path.join(base_dir, "*.anchors")) + \
                   glob.glob(os.path.join(base_dir, "**", "*.anchors"), recursive=True)
    
    # 过滤掉已经生成的 .simple 文件本身（防止干扰）
    anchor_files = [f for f in anchor_files if not f.endswith('.simple')]
    # 去重
    anchor_files = list(set([os.path.abspath(f) for f in anchor_files]))

    if not anchor_files:
        print("❌ 错误：当前目录及其子目录下找不到任何 .anchors 或 .lifted.anchors 文件！")
        print("💡 提示：请确保你处于刚刚运行完共线性搜索的工作区文件夹中。")
        sys.exit(1)

    print(f"🔍 扫描到以下 {len(anchor_files)} 个 Anchors 结果文件:")
    for i, f in enumerate(anchor_files, 1):
        print(f"  [{i}] {os.path.relpath(f, base_dir)}")

    while True:
        choice = input("\n👉 请选择要简化的 Anchors 文件编号 (或输入 q 退出): ").strip()
        if choice.lower() == 'q': sys.exit(0)
        try:
            selected_anchor = anchor_files[int(choice) - 1]
            break
        except:
            print("⚠️ 输入无效，请重新选择。")

    # 3. 交互式配置 minspan 阈值
    print("\n" + "-" * 40)
    print("💡 关键步骤：设置最小区块跨度过滤阈值 (--minspan)")
    print("  - 代表一个共线性区块中至少要包含多少个核心基因对。")
    print("  - 推荐值: 30 (适合宏观染色体大图，过滤掉小的碎片噪声，画面更干净) 。")
    print("  - 如果你的物种基因组非常小（如真菌），可以尝试降低到 10 或 20。")
    
    while True:
        minspan_input = input("\n👉 请输入 minspan 值 (直接回车使用推荐值 30): ").strip()
        if not minspan_input:
            minspan = 30
            break
        try:
            minspan = int(minspan_input)
            if minspan > 0:
                break
            print("⚠️ 阈值必须大于 0！")
        except ValueError:
            print("❌ 输入错误，请输入一个有效的正整数。")

    # 4. 构建输出路径与执行命令
    anchor_dir, anchor_name = os.path.split(selected_anchor)
    output_anchor = selected_anchor + ".new"

    print("\n" + "=" * 60)
    print("🚀 正在运行 JCVI 宏观区块聚合滤波算法...")
    print("=" * 60)

    # 核心命令：python -m jcvi.compara.synteny screen --minspan=X --simple A B
    cmd = f"python -m jcvi.compara.synteny screen --minspan={minspan} --simple '{selected_anchor}' '{output_anchor}'"
    print(f"▶️ 执行命令: {cmd}\n")

    try:
        # 在文件所在子目录下安全执行，防止路径错乱
        subprocess.run(cmd, shell=True, check=True, cwd=anchor_dir)

        # JCVI 算法在使用 --simple 时，会自动在输入文件同目录下生成一个以 .simple 结尾的文件
        expected_simple = selected_anchor.replace(".anchors", ".simple")

        print("\n" + "=" * 60)
        print("🎉 宏观共线性块过滤与简化成功！")
        print(f"📂 核心大产物已生成:")
        if os.path.exists(expected_simple):
            print(f"   💎 \033[92m{os.path.relpath(expected_simple, base_dir)}\033[0m")
        else:
            print(f"   💎 \033[92m{anchor_name.replace('.anchors', '.simple')} (保存在子目录下)\033[0m")
        
        print("\n💡 下一步全套作图冲刺指南:")
        print("  1. 运行 `python generate_seqids.py` 选择并重排你要上图的染色体。")
        print("  2. 在当前文件夹下创建或编辑 `layout` 文本文件。")
        print("  3. 运行终极绘图命令：")
        print("     \033[96mpython -m jcvi.graphics.karyotype seqids layout\033[0m")
        print("=" * 60 + "\n")

    except subprocess.CalledProcessError:
        print("\n❌ 运行失败，请检查上方 JCVI 的报错信息。")

if __name__ == "__main__":
    main()