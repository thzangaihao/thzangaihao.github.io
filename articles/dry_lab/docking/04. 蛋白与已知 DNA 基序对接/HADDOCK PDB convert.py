#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import os
import sys
import glob
import time

# ============= 基础路径与文件查找 (保留你的逻辑) =============
def get_base_dir():
    if getattr(sys, 'frozen', False):
        return os.path.dirname(sys.executable)
    return os.path.dirname(os.path.abspath(__file__))

def find_files(ext='.pdb'):
    path = get_base_dir()
    search_pattern = os.path.join(path, '**', f'*{ext}')
    return sorted(glob.glob(search_pattern, recursive=True))

# ============= 核心修复逻辑 (集成 DNA 改名 + 全局重编号) =============
def process_pdb_line(line, res_map, global_info):
    """处理单行 PDB 数据，确保列宽绝对固定"""
    if line.startswith(("ATOM", "HETATM")):
        # 1. 提取原始字段
        res_name = line[17:20].strip()
        chain_id = line[21]
        res_seq_raw = line[22:26].strip()
        
        # 2. 判断是否是新残基，更新全局编号
        current_res_uid = f"{chain_id}_{res_seq_raw}"
        if current_res_uid != global_info["last_res_uid"]:
            global_info["count"] += 1
            global_info["last_res_uid"] = current_res_uid
        
        # 3. 修复残基名 (A -> DA 等)
        new_res_name = res_map.get(res_name, res_name)
        
        # 4. 组装新行 (严格遵守 PDB Column 格式)
        # ATOM/HETATM(6) + ID(5) + Space(1) + Name(4) + Alt(1) + ResName(3) + Chain(1) + ResSeq(4) ...
        # 我们采用切片替换法最安全
        line_list = list(line.ljust(80)) 
        
        # 替换 ATOM (强制转换 HETATM 为 ATOM)
        line_list[0:6] = list("ATOM  ")
        # 替换残基名 (17-20)
        line_list[17:20] = list(f"{new_res_name:>3}")
        # 替换残基编号 (22-26)
        line_list[22:26] = list(f"{global_info['count']:>4}")
        
        return "".join(line_list).rstrip() + "\n"
    
    elif line.startswith("TER"):
        # TER 也要对应最后的编号
        return f"TER   {line[6:11]}  {line[17:21]}{global_info['count']:>4}\n"
    
    return line

def fix_file(file_path):
    res_map = {'G': 'DG', 'A': 'DA', 'T': 'DT', 'C': 'DC', 'U': 'U '}
    global_info = {"count": 0, "last_res_uid": None}
    
    new_content = []
    with open(file_path, 'r') as f:
        for line in f:
            new_content.append(process_pdb_line(line, res_map, global_info))
    
    output_path = os.path.join(os.path.dirname(file_path), f"fixed_{os.path.basename(file_path)}")
    with open(output_path, 'w') as f:
        f.writelines(new_content)
    return output_path

# ============= 交互界面 (保留 cite_v2.py 的风格) =============
def main():
    print("="*50)
    print("   HADDOCK PDB 格式交互式修复工具 v3.0")
    print("   功能：1.DNA单字母转标准双字母  2.全局重编号防止冲突")
    print("="*50)

    files = find_files()
    if not files:
        print("❌ 未能在当前目录及子目录下找到 .pdb 文件")
        return

    print(f"\n1. 找到以下 {len(files)} 个文件:")
    for i, f in enumerate(files, 1):
        print(f"  [{i}] {os.path.relpath(f, get_base_dir())}")

    user_input = input("\n2. 请输入编号进行处理 (例如: 1,2 或 1-3 或 all): ").strip().lower()
    
    selected_indices = []
    if user_input in ['all', 'a']:
        selected_indices = range(len(files))
    else:
        try:
            for part in user_input.split(','):
                if '-' in part:
                    start, end = map(int, part.split('-'))
                    selected_indices.extend(range(start-1, end))
                else:
                    selected_indices.append(int(part)-1)
        except:
            print("❌ 输入格式错误。")
            return

    print("\n3. 正在处理...")
    for idx in selected_indices:
        if 0 <= idx < len(files):
            p = files[idx]
            out = fix_file(p)
            print(f"  ✅ 已修复: {os.path.basename(out)}")
    
    print("\n所有任务完成！现在可以使用 fixed_ 开头的文件进行对接了。")

if __name__ == "__main__":
    main()