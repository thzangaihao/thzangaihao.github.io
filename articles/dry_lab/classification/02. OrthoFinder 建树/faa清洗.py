#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import os
import glob
import re
import time

def clean_fasta_strictly():
    print("=== FASTA 序列终极洁癖清洗工具 ===")
    
    # 获取当前目录下所有的 .faa, .pep, .fasta 文件
    extensions = ['*.faa', '*.pep', '*.fasta']
    files = []
    for ext in extensions:
        files.extend(glob.glob(ext))
        
    if not files:
        print("未找到序列文件，请确保脚本放在包含序列的文件夹中。")
        return

    out_dir = "Cleaned_faa"
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
        
    print(f"找到 {len(files)} 个文件，开始执行严格清洗（过滤所有含非法字符的日志行）...\n")
    
    start = time.time()
    
    for f in files:
        file_name = os.path.basename(f)
        out_path = os.path.join(out_dir, file_name)
        
        valid_lines_count = 0
        deleted_lines_count = 0
        
        with open(f, 'r', encoding='utf-8') as fin, open(out_path, 'w', encoding='utf-8') as fout:
            for line in fin:
                line = line.strip()
                if not line:
                    continue  # 跳过空行
                    
                if line.startswith(">"):
                    fout.write(line + "\n")
                    valid_lines_count += 1
                else:
                    # 【核心过滤】：纯字母校验。只要这一行包含数字、标点符号，re.match 就会返回 None
                    if re.match(r'^[A-Za-z]+$', line):
                        fout.write(line + "\n")
                        valid_lines_count += 1
                    else:
                        deleted_lines_count += 1
                        
        print(f"  √ {file_name} -> 保留了 {valid_lines_count} 行，删除了 {deleted_lines_count} 行垃圾日志")

    print("\n" + "="*40)
    print(f"🎉 全部清洗完毕！耗时: {time.time()-start:.2f} 秒。")
    print(f"干净的序列已保存在 '{out_dir}' 文件夹中。")

if __name__ == "__main__":
    clean_fasta_strictly()