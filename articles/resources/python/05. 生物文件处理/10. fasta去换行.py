#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import os
import sys
import glob

# ============= 基础路径获取 =============
def get_base_dir():
    """获取脚本所在目录"""
    if getattr(sys, 'frozen', False):
        return os.path.dirname(sys.executable)
    return os.path.dirname(os.path.abspath(__file__))

# ============= 文件查找 =============
def find_fasta_files(path=None):
    """搜索常见的FASTA后缀文件"""
    if path is None:
        path = get_base_dir()
    
    # 定义匹配的后缀名
    extensions = ['*.fasta', '*.fna', '*.faa', '*.fa']
    files = []
    for ext in extensions:
        # 递归查找子目录下的文件
        files.extend(glob.glob(os.path.join(path, '**', ext), recursive=True))
    return sorted(list(set(files)))

# ============= 多选逻辑 (参考 cite_v2.py) =============
def choose_multiple_files(files):
    if not files:
        print("未找到任何序列文件。")
        return []

    print("\n找到以下文件:")
    for i, f in enumerate(files, 1):
        print(f"  {i}. {os.path.basename(f)}")
    
    print("\n[选择方式说明]:")
    print("  - 输入单个数字: 1")
    print("  - 输入多个数字: 1,2,5")
    print("  - 输入范围: 1-4")
    print("  - 输入全部: all")
    
    while True:
        try:
            choice = input("\n请输入要处理的文件编号: ").strip().lower()
            if not choice: continue
            
            selected_indices = []
            if choice == 'all':
                selected_indices = list(range(len(files)))
            else:
                parts = choice.split(',')
                for part in parts:
                    if '-' in part:
                        start, end = map(int, part.split('-'))
                        selected_indices.extend(range(start-1, end))
                    else:
                        selected_indices.append(int(part)-1)
            
            # 过滤掉无效索引并去重
            selected_paths = [files[i] for i in sorted(set(selected_indices)) if 0 <= i < len(files)]
            
            if not selected_paths:
                print("未选中任何有效文件，请重新输入。")
                continue
                
            print(f"\n已选中 {len(selected_paths)} 个文件:")
            for p in selected_paths:
                print(f"  - {os.path.basename(p)}")
            return selected_paths

        except ValueError:
            print("输入错误：请确保输入的是数字、范围（如1-5）或 'all'。")

# ============= 核心转换函数 =============
def process_single_file(input_path):
    """将多行FASTA转换为单行"""
    base, ext = os.path.splitext(input_path)
    output_path = f"{base}_oneline{ext}"
    
    with open(input_path, 'r', encoding='utf-8') as f_in, \
         open(output_path, 'w', encoding='utf-8') as f_out:
        
        header = None
        seq_parts = []
        
        for line in f_in:
            line = line.strip()
            if not line: continue
            
            if line.startswith('>'):
                if header:
                    f_out.write(f"{header}\n{''.join(seq_parts)}\n")
                header = line
                seq_parts = []
            else:
                seq_parts.append(line)
        
        # 写入最后一个条目
        if header:
            f_out.write(f"{header}\n{''.join(seq_parts)}\n")
    
    return output_path

# ============= 主程序 =============
def main():
    print("="*50)
    print("      批量 FASTA 序列单行化转换工具      ")
    print("="*50)

    # 1. 自动搜索
    all_files = find_fasta_files()
    
    # 2. 多项选择
    selected_files = choose_multiple_files(all_files)
    
    if not selected_files:
        return

    # 3. 确认执行
    confirm = input(f"\n确认要处理这 {len(selected_files)} 个文件吗? (y/n): ").strip().lower()
    if confirm not in ['y', 'yes']:
        print("操作已取消。")
        return

    # 4. 循环处理
    print("\n正在开始批量处理...")
    success_count = 0
    for f_path in selected_files:
        try:
            out_path = process_single_file(f_path)
            print(f"  [成功] {os.path.basename(f_path)} -> {os.path.basename(out_path)}")
            success_count += 1
        except Exception as e:
            print(f"  [失败] {os.path.basename(f_path)}: {e}")

    print("-" * 50)
    print(f"处理结束！共成功处理 {success_count} 个文件。")

if __name__ == '__main__':
    main()