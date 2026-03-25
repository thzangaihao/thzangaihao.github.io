#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import os
import sys

def simplify_fasta():
    print("=== FASTA 标签简化与重命名工具 ===")
    
    # 1. 自动定位 OrthoFinder 生成的文件
    input_file = "SpeciesTreeAlignment.fa"
    if not os.path.exists(input_file):
        input_file = input("未在当前目录找到 SpeciesTreeAlignment.fa，请输入完整路径: ").strip()
    
    if not os.path.exists(input_file):
        print("错误：找不到输入文件！")
        return

    # 2. 提取所有原始 ID
    original_ids = []
    with open(input_file, 'r') as f:
        for line in f:
            if line.startswith(">"):
                original_ids.append(line.strip()[1:])
    
    print(f"\n找到 {len(original_ids)} 个物种。现在开始设置简短名：")
    print("提示：输入为空则保留原名，输入 's' 则仅保留 GCA 编号。")
    
    # 3. 交互式建立映射字典
    name_map = {}
    for old_id in original_ids:
        # 如果是你的实验样本，默认给个好听的名字
        default_suggestion = old_id.replace("_cleaned_v1", "")
        
        new_name = input(f"原名: {old_id}\n改名 [默认 {default_suggestion}]: ").strip()
        
        if not new_name:
            name_map[old_id] = default_suggestion
        elif new_name.lower() == 's':
            # 自动提取 GCA 编号部分
            gca_match = old_id.split('_')[0:2]
            name_map[old_id] = "_".join(gca_match)
        else:
            name_map[old_id] = new_name
            
    # 4. 执行替换并写出新文件
    output_file = "Cleaned_Alignment.fa"
    with open(input_file, 'r') as fin, open(output_file, 'w') as fout:
        for line in fin:
            if line.startswith(">"):
                old_id = line.strip()[1:]
                fout.write(f">{name_map.get(old_id, old_id)}\n")
            else:
                fout.write(line)
                
    print(f"\n" + "="*40)
    print(f"🎉 转换完成！新文件已生成: {output_file}")
    print(f"下一步指令：\niqtree -s {output_file} -m MFP -bb 1000 -nt AUTO")

if __name__ == "__main__":
    simplify_fasta()