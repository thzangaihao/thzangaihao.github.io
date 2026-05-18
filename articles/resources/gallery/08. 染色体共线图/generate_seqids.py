#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import sys
import glob

'''
📄 JCVI seqids 配置文件自动生成工具
功能：
1. 自动解析当前目录下的 .bed 文件，提取所有独一无二的染色体/Scaffold ID。
2. 交互式选择哪些染色体需要上图，并支持手动排序。
3. 严格按照 JCVI 规范输出无尾部空行的 seqids 文件。
'''

def parse_chromosomes_from_bed(bed_path):
    """从 BED 文件中提取所有染色体名，并保持其在基因组中的先后出现顺序"""
    chroms = []
    seen = set()
    with open(bed_path, 'r', encoding='utf-8') as f:
        for line in f:
            parts = line.strip().split('\t')
            if parts and parts[0]:
                chrom_name = parts[0]
                if chrom_name not in seen:
                    seen.add(chrom_name)
                    chroms.append(chrom_name)
    return chroms

def interactive_select_chroms(chroms, species_name):
    """交互式筛选和排序染色体"""
    print(f"\n📊 在 【{species_name}】 中共检测到 {len(chroms)} 条染色体/Scaffolds:")
    for i, chr_id in enumerate(chroms, 1):
        print(f"  [{i}] {chr_id}")
        
    print(f"\n💡 请选择要在图中绘制的染色体（依照你希望展示的从左到右的顺序输入）：")
    print("  - 直接回车（按 Enter）：默认选中上面【所有】染色体。")
    print("  - 自定义输入编号：例如输入 '1,2,3,5'（用英文逗号隔开）将只画这 4 条。")
    print("  - 改变显示顺序：例如输入 '3,2,1' 会在画布上颠倒它们的物理排列。")
    
    while True:
        user_input = input("👉 请输入你的选择: ").strip()
        if not user_input:
            return chroms
            
        try:
            # 解析用户输入的数字编号
            indices = [int(x.strip()) for x in user_input.split(',')]
            selected_chroms = [chroms[idx - 1] for idx in indices if 0 < idx <= len(chroms)]
            if selected_chroms:
                print(f"✅ 已确认选择的染色体顺序: {', '.join(selected_chroms)}")
                return selected_chroms
            else:
                print("⚠️ 未选中任何有效的染色体编号，请重新输入。")
        except (ValueError, IndexError):
            print("❌ 输入格式错误！请使用英文逗号分隔纯数字编号（如: 1,2,3,4）。")

def main():
    print("=" * 60)
    print(" 📄 JCVI seqids 配置文件交互式生成器")
    print("=" * 60)

    base_dir = os.getcwd()
    bed_files = glob.glob(os.path.join(base_dir, "*.bed"))

    if len(bed_files) < 2:
        print("❌ 错误：当前目录下找不到足够的 .bed 文件！至少需要 2 个。")
        sys.exit(1)

    print(f"🔍 扫描到以下 {len(bed_files)} 个位置文件 (.bed):")
    for i, f in enumerate(bed_files, 1):
        print(f"  [{i}] {os.path.basename(f)}")

    # 1. 选择放在第一行的物种 (通常在图的上方轨)
    while True:
        try:
            choice_a = int(input("\n👉 请选择放在【第一轨 (图上方)】的物种编号: ").strip())
            bed_a = bed_files[choice_a - 1]
            prefix_a = os.path.basename(bed_a).replace('.bed', '')
            break
        except:
            print("⚠️ 输入无效，请重新选择。")

    # 2. 选择放在第二行的物种 (通常在图的下方轨)
    while True:
        try:
            choice_b = int(input("👉 请选择放在【第二轨 (图下方)】的物种编号: ").strip())
            if choice_b == choice_a:
                print("⚠️ 不能选择同一个物种作为比较对象！")
                continue
            bed_b = bed_files[choice_b - 1]
            prefix_b = os.path.basename(bed_b).replace('.bed', '')
            break
        except:
            print("⚠️ 输入无效，请重新选择。")

    # 3. 分别提取并筛选两个物种的染色体
    chroms_a = parse_chromosomes_from_bed(bed_a)
    selected_a = interactive_select_chroms(chroms_a, prefix_a)

    chroms_b = parse_chromosomes_from_bed(bed_b)
    selected_b = interactive_select_chroms(chroms_b, prefix_b)

    # 4. 严格按照 JCVI 规范写入 seqids 文件
    output_path = os.path.join(base_dir, "seqids")
    
    line1 = ",".join(selected_a)
    line2 = ",".join(selected_b)
    
    # ⚠️ 核心防爆点：两行数据之间用 \n 隔开，但第二行末尾绝对不加 \n，确保没有尾部空行
    with open(output_path, 'w', encoding='utf-8') as f:
        f.write(line1 + "\n" + line2)

    print("\n" + "=" * 60)
    print("🎉 配置文件 `seqids` 生成成功！")
    print(f"📂 保存路径: {output_path}")
    print("\n📝 文件内容预览:")
    print(f"   第 1 行 ({prefix_a}): {line1}")
    print(f"   第 2 行 ({prefix_b}): {line2}")

if __name__ == "__main__":
    main()