#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import sys
import glob
import subprocess
import shutil

# ==============================================================================
# 0. 全局配置
# ==============================================================================
# 输出目录名称
OUTPUT_DIR_NAME = "02_Merged_Results"
# 输出文件名
OUTPUT_FILENAME = "Merged_Population.vcf.gz"
# 线程数
THREADS = 16

# BCFtools 可执行命令 (如果需要指定路径请修改这里)
BCFTOOLS_EXEC = "bcftools"

# ==============================================================================
# 1. 交互与搜索模块 (Cite Logic)
# ==============================================================================
def get_base_dir():
    """获取脚本所在目录"""
    return os.path.dirname(os.path.abspath(__file__))

def find_files(ext, path=None):
    """递归查找指定后缀文件"""
    if path is None: path = get_base_dir()
    if not ext.startswith('.'): ext = '.' + ext
    
    print(f"正在搜索 {ext} 文件...")
    # 排除掉输出目录，防止把上次合并的结果又加进来了
    all_files = glob.glob(os.path.join(path, '**', f'*{ext}'), recursive=True)
    filtered_files = [f for f in all_files if OUTPUT_DIR_NAME not in f]
    
    return sorted(filtered_files)

def choose_files(files, desc="文件"):
    """交互式多选逻辑"""
    if not files:
        print(f"[提示] 未找到任何 {desc}。")
        return []
    
    print(f"\n--- 找到 {len(files)} 个 {desc} ---")
    limit = 10
    if len(files) > limit * 2:
        for i, f in enumerate(files[:limit], 1):
            print(f"  [{i}] {os.path.basename(f)}")
        print(f"  ... (中间省略 {len(files) - limit*2} 个) ...")
        for i, f in enumerate(files[-limit:], len(files)-limit+1):
            print(f"  [{i}] {os.path.basename(f)}")
    else:
        for i, f in enumerate(files, 1):
            print(f"  [{i}] {os.path.basename(f)}")

    while True:
        try:
            prompt = f"\n请输入要合并的文件编号 (支持: all | 1-5 | 1,3,5): "
            user_input = input(prompt).strip().lower()
            if not user_input: continue
            
            if user_input in ['all', 'a']:
                return files
            
            selected_indices = set()
            parts = user_input.split(',')
            for part in parts:
                if '-' in part:
                    s, e = map(int, part.split('-'))
                    selected_indices.update(range(s-1, e))
                else:
                    selected_indices.add(int(part)-1)
            
            selected_files = [files[i] for i in sorted(selected_indices) if 0 <= i < len(files)]
            
            if selected_files:
                print(f"已选择 {len(selected_files)} 个文件。")
                return selected_files
            print("选择无效，请重试。")
        except ValueError:
            print("输入格式错误。")

# ==============================================================================
# 2. 核心处理模块
# ==============================================================================
def check_bcftools():
    if shutil.which(BCFTOOLS_EXEC) is None:
        print(f"[错误] 未找到 '{BCFTOOLS_EXEC}'。请先安装: conda install bcftools")
        sys.exit(1)

def ensure_index(vcf_files):
    """
    检查并构建索引 (.tbi 或 .csi)
    bcftools merge 要求所有输入文件必须有索引
    """
    print("\n--- 步骤 1/2: 检查文件索引 ---")
    files_to_index = []
    
    for f in vcf_files:
        # 检查是否存在 .tbi 或 .csi
        if not (os.path.exists(f + ".tbi") or os.path.exists(f + ".csi")):
            files_to_index.append(f)
    
    if not files_to_index:
        print("所有文件均已有索引，跳过构建步骤。")
        return

    print(f"发现 {len(files_to_index)} 个文件缺少索引，正在构建索引...")
    
    # 简单循环构建 (为了稳定性)
    for i, f in enumerate(files_to_index, 1):
        print(f"[{i}/{len(files_to_index)}] Indexing: {os.path.basename(f)}")
        try:
            # -t: 生成 .tbi (兼容性好), -f: 强制刷新
            subprocess.run([BCFTOOLS_EXEC, "index", "-t", "-f", f], check=True)
        except subprocess.CalledProcessError:
            print(f"  [警告] 索引构建失败: {os.path.basename(f)}")

def run_merge(file_list):
    """执行合并"""
    print("\n--- 步骤 2/2: 执行合并 ---")
    
    # 1. 准备输出目录
    base_dir = get_base_dir()
    out_dir = os.path.join(base_dir, OUTPUT_DIR_NAME)
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
        print(f"创建输出目录: {out_dir}")
    
    output_path = os.path.join(out_dir, OUTPUT_FILENAME)
    
    # 2. 生成文件列表 (避免命令行过长)
    list_file_path = os.path.join(out_dir, "merge_list.txt")
    with open(list_file_path, 'w') as f:
        for p in file_list:
            f.write(p + "\n")
            
    # 3. 构造命令
    # -m none: 解决 multiallelics 报错
    # -0: 缺失转参考 (SV 必须)
    # -O z: 输出压缩的 VCF (vcf.gz)
    cmd = [
        BCFTOOLS_EXEC, "merge",
        "-l", list_file_path,
        "-o", output_path,
        "-O", "z",
        "--threads", str(THREADS),
        "-m", "none",
        "-0"
    ]
    
    print(f"执行命令: {' '.join(cmd)}")
    
    try:
        subprocess.run(cmd, check=True)
        print("-" * 50)
        print(f"[成功] 合并完成！")
        print(f"结果文件: {output_path}")
        
        # 顺手给结果文件建个索引，方便后续查看
        print("正在为结果文件建立索引...")
        subprocess.run([BCFTOOLS_EXEC, "index", "-t", output_path], check=False)
        
        # 清理临时列表
        os.remove(list_file_path)
        
    except subprocess.CalledProcessError as e:
        print(f"\n[失败] 合并过程中出错: {e}")

# ==============================================================================
# 主函数
# ==============================================================================
def main():
    print("=== VCF.GZ 批量合并工具 (优化版) ===")
    check_bcftools()
    
    # 1. 搜索
    # 假设都在 01_Raw_VCFs 或当前目录下
    target_files = find_files("vcf.gz")
    
    # 2. 选择
    selected = choose_files(target_files, ".vcf.gz 文件")
    
    if selected:
        # 3. 索引
        ensure_index(selected)
        
        # 4. 合并
        run_merge(selected)

if __name__ == "__main__":
    main()