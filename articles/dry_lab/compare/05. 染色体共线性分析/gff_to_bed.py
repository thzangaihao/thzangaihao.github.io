#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import sys
import glob
import subprocess

'''
📂 JCVI 简化位置文件 (BED) 交互式转换工具 (终极防撞车优化版)
功能：
1. 自动扫描当前目录及子目录下的 .gff / .gff3 文件。
2. 交互式选择文件与特征类型 (mRNA / transcript / gene)。
3. 🌟 新增：强制前缀注入功能！在生成 BED 后，自动为所有 ID 添加物种前缀，彻底根除 Liftoff 同名陷阱。
'''

def get_base_dir():
    return os.path.dirname(os.path.abspath(__file__))

def main():
    print("=" * 60)
    print(" 📂 JCVI 简化位置文件 (BED) 交互式转换工具 (终极版)")
    print("=" * 60)

    # 1. 检查当前 Conda 环境
    try:
        subprocess.run(
            ["python", "-m", "jcvi"], 
            stdout=subprocess.DEVNULL, 
            stderr=subprocess.DEVNULL
        )
    except Exception:
        print("❌ 错误：未在当前环境中检测到 jcvi！")
        print("💡 请确保你已经执行了 `conda activate synteny`。")
        sys.exit(1)

    base_dir = get_base_dir()
    
    # 2. 交互式循环
    while True:
        gff_files = glob.glob(os.path.join(base_dir, "**", "*.gff*"), recursive=True)
        
        if not gff_files:
            print("\n⚠️ 未在当前目录及子目录下找到任何 .gff 或 .gff3 文件！")
            sys.exit(0)
            
        print(f"\n🔍 扫描到以下 {len(gff_files)} 个 GFF 注释文件:")
        for i, f in enumerate(gff_files, 1):
            print(f"  [{i}] {os.path.relpath(f, base_dir)}")
            
        choice = input("\n👉 请选择要转换的 GFF 文件编号 (或输入 q 退出): ").strip()
        if choice.lower() == 'q':
            print("[INFO] 程序已退出。")
            break
            
        try:
            selected_gff = gff_files[int(choice) - 1]
        except Exception:
            print("❌ 输入无效，请重新选择。")
            continue
            
        # 3. 选择特征类型
        print("\n" + "-" * 40)
        print("💡 关键步骤 1：选择提取的特征类型 (Feature Type)")
        print("  [1] mRNA (默认推荐)")
        print("  [2] transcript")
        print("  [3] gene")
        print("  [4] 自定义")
        
        type_choice = input("👉 请选择 (直接回车默认选 1): ").strip()
        
        feature_type = "mRNA"
        if type_choice == '2': feature_type = "transcript"
        elif type_choice == '3': feature_type = "gene"
        elif type_choice == '4': feature_type = input("👉 请输入自定义类型: ").strip()

        # 4. 提取文件前缀
        gff_dir, gff_name = os.path.split(selected_gff)
        name_root = gff_name.split('.')[0]
        output_bed = os.path.join(gff_dir, f"{name_root}.bed")
        
        print(f"\n🚀 正在转换: {gff_name} ➔ {name_root}.bed (提取特征: {feature_type}) ...")
        
        # 5. 调用 JCVI 转换
        cmd = f"python -m jcvi.formats.gff bed '{selected_gff}' -o '{output_bed}' --type={feature_type}"
        
        try:
            subprocess.run(cmd, shell=True, check=True)
            print(f"🎉 基础位置文件生成成功！")
            
            # 🌟 6. 核心优化：Liftoff 陷阱防护 (前缀注入)
            print("\n" + "-" * 40)
            print("🛡️ 关键步骤 2：防止 Liftoff 同名 ID 冲突 (0 Anchors 陷阱)")
            print("建议为 ID 加上当前物种的前缀，以此区分不同物种的同名基因。")
            prefix_input = input(f"👉 请输入物种前缀 (直接回车默认使用 '{name_root}'; 输入 'n' 跳过): ").strip()

            if prefix_input.lower() != 'n':
                prefix = prefix_input if prefix_input else name_root
                print(f"🔧 正在将 '{prefix}_' 注入到 {name_root}.bed 的所有基因 ID 中...")
                
                lines = []
                with open(output_bed, 'r') as f:
                    for line in f:
                        parts = line.strip().split('\t')
                        if len(parts) >= 4:
                            # 避免重复添加前缀
                            if not parts[3].startswith(f"{prefix}_"):
                                parts[3] = f"{prefix}_{parts[3]}"
                        lines.append('\t'.join(parts))
                        
                with open(output_bed, 'w') as f:
                    f.write('\n'.join(lines) + '\n')
                    
                print(f"✅ 完美！该物种的基因 ID 已全部变更为 '{prefix}_XXX' 格式。")
                print("⚠️ 【重要提示】：请确保你在比对前，给 .faa 蛋白序列的 ID 也加上了相同的前缀。")
                print("   （或者跑完 Diamond 后，用我们之前的 fix_liftoff_ids.py 脚本给 blast 文件加上前缀）")
            
        except subprocess.CalledProcessError:
            print(f"❌ 转换失败！jcvi 报错。请检查该 GFF 文件的格式是否规范。")
            continue
            
        # 7. 询问是否继续
        cont = input("\n[?] 是否继续转换另一个物种的 GFF 文件？(y/n): ").strip().lower()
        if cont != 'y':
            print("\n" + "=" * 60)
            print("[INFO] 所有位置文件转换结束！")
            print("=" * 60 + "\n")
            break

if __name__ == "__main__":
    main()