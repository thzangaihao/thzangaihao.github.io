#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import sys
import glob
import re

'''
✨ SVG 矢量图描边一键清洗系统
功能：
1. 自动扫描当前目录及其子目录下的所有 .svg 文件。
2. 交互式选择需要处理的 SVG 图形。
3. 利用正则表达式，直接修改底层的 XML/SVG 代码，无损抹除所有描边属性。
4. 自动生成无描边的纯净版文件，保留原文件。
'''

def get_base_dir():
    return os.path.dirname(os.path.abspath(__file__))

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

def clean_svg_strokes(input_svg, output_svg):
    """
    核心清洗逻辑：使用正则暴力替换 SVG 中的描边属性
    """
    try:
        with open(input_svg, 'r', encoding='utf-8') as f:
            svg_content = f.read()

        # 替换 1：处理 style 属性中的 stroke-width (例如 style="...;stroke-width:1.0;...")
        svg_content = re.sub(r'stroke-width\s*:\s*[\d.]+', 'stroke-width:0', svg_content)
        
        # 替换 2：处理直接的 stroke-width 属性 (例如 stroke-width="1.0")
        svg_content = re.sub(r'stroke-width\s*=\s*"[\d.]+"', 'stroke-width="0"', svg_content)
        
        # 替换 3：处理 style 属性中的 stroke 颜色，将其置为空或 none (针对 matplotlib 的常见黑色描边)
        svg_content = re.sub(r'stroke\s*:\s*(#[0-9a-fA-F]{6}|#000000|black)', 'stroke:none', svg_content, flags=re.IGNORECASE)
        
        # 替换 4：处理直接的 stroke 属性
        svg_content = re.sub(r'stroke\s*=\s*"(#[0-9a-fA-F]{6}|#000000|black)"', 'stroke="none"', svg_content, flags=re.IGNORECASE)

        with open(output_svg, 'w', encoding='utf-8') as f:
            f.write(svg_content)
            
        return True
    except Exception as e:
        print(f"❌ 处理文件时发生错误: {e}")
        return False

def main():
    print("=" * 60)
    print(" ✨ SVG 矢量图描边一键清洗系统")
    print("=" * 60)

    base_dir = get_base_dir()

    # 1. 扫描所有 SVG 文件
    # 排除掉之前可能已经清洗过的文件，避免重复套娃
    all_svg_files = glob.glob(os.path.join(base_dir, "**", "*.svg"), recursive=True)
    svg_files = [f for f in all_svg_files if not f.endswith("_clean.svg")]

    if not svg_files:
        print("❌ 错误：当前目录及子目录下找不到任何 .svg 文件！")
        sys.exit(1)

    # 2. 交互式选择文件
    selected_svg = interactive_select(svg_files, "需要去除描边的 SVG 文件")
    
    # 3. 确定输出文件名
    dir_name = os.path.dirname(selected_svg)
    base_name = os.path.basename(selected_svg)
    name_part, ext_part = os.path.splitext(base_name)
    output_svg = os.path.join(dir_name, f"{name_part}_clean{ext_part}")

    # 4. 执行清洗
    print("\n" + "="*60)
    print(f"🔍 正在读取文件: {base_name}")
    print("🧹 正在通过正则匹配抹除所有 stroke 属性...")
    
    success = clean_svg_strokes(selected_svg, output_svg)
    
    if success:
        print("\n🎉 恭喜！SVG 描边清洗完成！原图已脱胎换骨。")
        print(f"📁 原始文件: {os.path.relpath(selected_svg, base_dir)}")
        print(f"💎 清洗结果: \033[93m{os.path.relpath(output_svg, base_dir)}\033[0m")
    print("=" * 60 + "\n")

if __name__ == "__main__":
    main()