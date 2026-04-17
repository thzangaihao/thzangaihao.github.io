#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import os
import sys
import glob
import re
import time

# ============= 基础路径获取 =============
def get_base_dir():
    """获取脚本所在目录，作为搜索的起点"""
    if getattr(sys, 'frozen', False):
        return os.path.dirname(sys.executable)
    return os.path.dirname(os.path.abspath(__file__))

# ============= 高级文件查找 =============
def find_files(ext, path=None):
    """查找指定扩展名文件"""
    if not ext.startswith('.'):
        ext = '.' + ext
    if path is None:
        path = get_base_dir()
    search_pattern = os.path.join(path, '**', f'*{ext}')
    return sorted(glob.glob(search_pattern, recursive=True))

# ============= 文件选择（支持多种选取模式） =============
def choose_file(files, desc="文件"):
    if not files:
        print(f"提示：在当前目录及其子目录下未找到任何 {desc}")
        return []

    print(f"\n2. 找到 {len(files)} 个 {desc}:")
    for i, f in enumerate(files, 1):
        print(f"  [{i}] {f}")

    while True:
        try:
            prompt = (
                f"\n3. 请输入欲操作的 {desc} 编号\n"
                f" (支持格式: 1,2 | 1-4 | 1,3,5-8 | all): \n"
            )
            user_input = input(prompt).strip().lower()

            if not user_input:
                continue

            selected_indices = set()

            # [修复 Bug] 修改了原代码中错误的 or 逻辑判断
            if user_input in ['all', 'a']:
                selected_indices = set(range(len(files)))
            else:
                parts = user_input.split(',')
                for part in parts:
                    part = part.strip()
                    if '-' in part:
                        start_str, end_str = part.split('-')
                        start, end = int(start_str), int(end_str)
                        for idx in range(start - 1, end):
                            if 0 <= idx < len(files):
                                selected_indices.add(idx)
                    else:
                        idx = int(part) - 1
                        if 0 <= idx < len(files):
                            selected_indices.add(idx)
            
            if not selected_indices:
                print("未匹配到有效编号，请重新输入。")
                continue

            selected_paths = [files[i] for i in sorted(list(selected_indices))]
            
            print(f"\n4. 已选择 {len(selected_paths)} 个文件:")
            for p in selected_paths:
                print(f"  - {os.path.basename(p)}")
            
            return selected_paths

        except ValueError:
            print("输入错误：请确保输入的是数字编号、范围（如1-5）或 'all'。")
        except KeyboardInterrupt:
            print("\n用户取消操作。")
            sys.exit()

# ============= 确认函数 =============
def make_sure():
    response = input(f"\n5. 确认对上述选中的文件执行【清除标题序号】操作? (y/n): ").strip().lower()
    if response not in ['y', 'yes']:
        print("操作已取消。")
        return False
    return True

# ============= 核心清洗逻辑 =============
def process_html_file(filepath):
    """
    读取文件，清除 <h2> 和 <h3> 内部前置的手动序号，然后保存。
    """
    try:
        with open(filepath, 'r', encoding='utf-8') as f:
            content = f.read()
    except Exception as e:
        print(f"  [读取失败] {os.path.basename(filepath)}: {e}")
        return False

    def clean_heading(match):
        tag_start = match.group(1)  # 例如: <h2 id="xxx">
        inner_text = match.group(2) # 例如: 1. 引言
        tag_end = match.group(3)    # 例如: </h2>
        
        # 正则解析：
        # ^\s* : 忽略开头的空格
        # \d{1,2}     : 匹配 1-2 位数字（防止把年份 "2026" 误删）
        # (?:\.\d{1,2})* : 匹配后续的小数点及数字，如 ".1" 或 ".1.2"
        # \.?         : 匹配最后一个可能存在的小数点，如 "1."
        # \s+         : 必须以至少一个空格结尾（"1.标题" 会被放过，必须是 "1. 标题"）
        cleaned_text = re.sub(r'^\s*\d{1,2}(?:\.\d{1,2})*\.?\s+', '', inner_text)
        
        return f"{tag_start}{cleaned_text}{tag_end}"

    # 使用正则匹配 <h2>...</h2> 和 <h3>...</h3>
    new_content = re.sub(r'(<h[23][^>]*>)(.*?)(</h[23]>)', clean_heading, content, flags=re.IGNORECASE | re.DOTALL)
    
    if new_content != content:
        with open(filepath, 'w', encoding='utf-8') as f:
            f.write(new_content)
        return True # 表示有修改
    return False # 表示无需修改

# ============= 主函数 =============
def main():
    print("=========================================")
    print("   海钊知识港 - HTML 标题序号自动化清洗工具")
    print("=========================================\n")
    
    # 默认直接查找 .html 文件，省去手动输入的步骤
    target_ext = ".html"
    print(f"1. 正在搜索当前目录及子目录下的 {target_ext} 文件...")
    
    files = find_files(target_ext)
    chosen_list = choose_file(files, desc=f"{target_ext} 文章")
    
    if chosen_list and make_sure():
        print(f"\n--- 开始处理 {len(chosen_list)} 个文件 ---")
        
        modified_count = 0
        skipped_count = 0
        
        for filepath in chosen_list:
            filename = os.path.basename(filepath)
            is_modified = process_html_file(filepath)
            
            if is_modified:
                print(f"  [修改成功] {filename}")
                modified_count += 1
            else:
                print(f"  [无需修改] {filename}")
                skipped_count += 1
                
        print("\n--- 处理完成 ---")
        print(f"总计处理: {len(chosen_list)} 篇")
        print(f"成功清洗: {modified_count} 篇")
        print(f"无需改动: {skipped_count} 篇")

if __name__ == "__main__":
    main()