#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''
全自动 Kraken2 交互式鉴定脚本

功能特色：
1. 自动全盘扫描并配对双端 Fastq 文件 (支持 .fq.gz / .fastq.gz)
2. 自动定位 Kraken2 数据库 (通过寻找 hash.k2d)
3. 交互式样本选择 (参考 cite_v2 逻辑，支持 1, 2, 3-5, all)
4. 多线程全自动批量比对

使用方法：
在终端中直接运行: python run_kraken2.py

thz - 自动生成
'''

import os
import sys
import glob
import subprocess
import time
import re

# ============= 基础路径获取 =============
def get_base_dir():
    """获取脚本所在目录，作为搜索的起点"""
    if getattr(sys, 'frozen', False):
        return os.path.dirname(sys.executable)
    return os.path.dirname(os.path.abspath(__file__))

# ============= 工具检查 =============
def check_tools():
    """检查 kraken2 是否安装并在环境变量中"""
    if subprocess.call(['which', 'kraken2'], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL) != 0:
        print("错误: 未找到 kraken2 命令，请确保已激活包含 kraken2 的 conda 环境！")
        sys.exit(1)

# ============= 数据库自动查找 =============
def find_kraken_db():
    """在当前目录及其子目录寻找包含 hash.k2d 的 Kraken2 数据库目录"""
    print("\n🔍 正在扫描 Kraken2 数据库...")
    base_dir = get_base_dir()
    search_pattern = os.path.join(base_dir, '**', 'hash.k2d')
    db_files = glob.glob(search_pattern, recursive=True)
    
    if not db_files:
        print("错误: 未在当前目录及子目录中找到 Kraken2 数据库 (缺失 hash.k2d)。")
        sys.exit(1)
    
    # 获取唯一的数据库目录
    db_dirs = list(set([os.path.dirname(f) for f in db_files]))
    
    if len(db_dirs) == 1:
        print(f"✅ 自动锁定数据库: {os.path.basename(db_dirs[0])} ({db_dirs[0]})")
        return db_dirs[0]
    else:
        print("\n找到多个 Kraken2 数据库，请选择:")
        for i, d in enumerate(db_dirs, 1):
            print(f"  [{i}] {d}")
        while True:
            try:
                choice = int(input("请输入数据库编号: ").strip())
                if 1 <= choice <= len(db_dirs):
                    return db_dirs[choice-1]
                else:
                    print("编号无效，请重试。")
            except ValueError:
                print("请输入有效数字。")

# ============= 样本自动发现与配对 =============
def get_samples():
    """扫描所有的 .fq.gz 和 .fastq.gz 并按样本名配对"""
    print("\n🔍 正在扫描测序样本数据...")
    base_dir = get_base_dir()
    fq_files = glob.glob(os.path.join(base_dir, '**', '*.fq.gz'), recursive=True)
    fq_files.extend(glob.glob(os.path.join(base_dir, '**', '*.fastq.gz'), recursive=True))
    
    samples = {}
    for f in fq_files:
        basename = os.path.basename(f)
        # 使用正则匹配提取样本名 (兼容 _R1, _1 等命名习惯)
        match = re.search(r'(.+)(_R[12]|_[12])(\.fastq\.gz|\.fq\.gz)$', basename)
        if match:
            sample_name = match.group(1)
            read_type = match.group(2)
            
            if sample_name not in samples:
                samples[sample_name] = {'r1': None, 'r2': None}
                
            if '1' in read_type:
                samples[sample_name]['r1'] = f
            elif '2' in read_type:
                samples[sample_name]['r2'] = f
    
    # 过滤掉单端不完整的样本
    valid_samples = {k: v for k, v in samples.items() if v['r1'] and v['r2']}
    
    if not valid_samples:
        print("错误: 未找到完整的双端测序数据 (需要 R1 和 R2)。")
        sys.exit(1)
        
    return valid_samples

# ============= 交互式样本选择 (参考 cite_v2) =============
def select_samples(samples_dict):
    """支持多种模式选择样本"""
    sample_list = sorted(list(samples_dict.keys()))
    
    print(f"\n✅ 找到 {len(sample_list)} 个有效双端样本:")
    for i, name in enumerate(sample_list, 1):
        print(f"  [{i}] {name}")
        
    print("\n💡 选择模式说明:")
    print("  - 单选: 1")
    print("  - 多选: 1,3,5")
    print("  - 范围: 1-3")
    print("  - 全部: all")
    
    while True:
        choice_str = input("\n👉 请输入要处理的样本编号 (输入 q 退出): ").strip().lower()
        
        if choice_str == 'q':
            print("已退出程序。")
            sys.exit(0)
            
        if choice_str == 'all':
            return {name: samples_dict[name] for name in sample_list}
            
        selected_samples = {}
        try:
            parts = choice_str.replace(' ', '').split(',')
            indices = set()
            for part in parts:
                if '-' in part:
                    start, end = map(int, part.split('-'))
                    indices.update(range(start, end + 1))
                else:
                    indices.add(int(part))
            
            # 验证索引
            invalid_indices = [idx for idx in indices if not (1 <= idx <= len(sample_list))]
            if invalid_indices:
                print(f"⚠️ 编号 {invalid_indices} 超出范围，请重新输入！")
                continue
                
            for idx in sorted(list(indices)):
                name = sample_list[idx-1]
                selected_samples[name] = samples_dict[name]
                
            return selected_samples
            
        except ValueError:
            print("⚠️ 输入格式错误！请确保输入的是数字编号、范围或 'all'。")

# ============= 确认步骤 =============
def make_sure(selected_names):
    print("\n" + "="*40)
    print(f"您已选择以下 {len(selected_names)} 个样本进行鉴定:")
    for name in selected_names:
        print(f"  - {name}")
    print("="*40)
    
    response = input(f"👉 确认开始执行 Kraken2 鉴定吗? (y/n): ").strip().lower()
    if response not in ['y', 'yes']:
        print("操作已取消。")
        sys.exit(0)

# ============= 主运行逻辑 =============
def run_kraken2(db_dir, selected_samples, threads):
    """执行 Kraken2 命令"""
    timestamp = time.strftime("%Y%m%d_%H%M%S")
    out_dir = os.path.join(get_base_dir(), f"03_Kraken2_鉴定结果_{timestamp}")
    os.makedirs(out_dir, exist_ok=True)
    
    print(f"\n🚀 开始执行鉴定，结果将保存在: {os.path.basename(out_dir)}\n")
    
    for i, (sample_name, paths) in enumerate(selected_samples.items(), 1):
        print("-" * 50)
        print(f"⏳ 正在处理 [{i}/{len(selected_samples)}]: {sample_name} ...")
        
        report_file = os.path.join(out_dir, f"{sample_name}_report.txt")
        output_file = os.path.join(out_dir, f"{sample_name}_kraken.out")
        
        cmd = [
            "kraken2",
            "--db", db_dir,
            "--threads", str(threads),
            "--paired", paths['r1'], paths['r2'],
            "--use-names",  # 输出中包含具体的物种名称而不仅仅是 TaxID
            "--report", report_file
        ]
        
        try:
            # 运行命令并将标准输出重定向到 .out 文件
            with open(output_file, "w") as f_out:
                subprocess.check_call(cmd, stdout=f_out)
            print(f"✅ {sample_name} 鉴定完成！")
            print(f"   📄 报告文件: {os.path.basename(report_file)}")
        except subprocess.CalledProcessError as e:
            print(f"❌ {sample_name} 鉴定失败！错误信息: {e}")
            
    print("\n🎉 所有选中样本处理完毕！")
    print(f"📂 结果总目录: {out_dir}")

# ============= 主函数 =============
def main():
    print("="*50)
    print("      Kraken2 物种鉴定全自动交互流水线      ")
    print("="*50)
    
    # 1. 检查环境变量
    check_tools()
    
    # 2. 获取 Kraken2 数据库路径
    db_dir = find_kraken_db()
    
    # 3. 获取所有有效样本
    all_samples = get_samples()
    
    # 4. 交互式选择样本
    selected_samples = select_samples(all_samples)
    
    # 5. 获取线程数
    while True:
        try:
            threads = int(input("\n💻 请输入运行使用的 CPU 核心数 (推荐 8-32): ").strip())
            if threads > 0:
                break
            else:
                print("⚠️ 线程数必须大于 0")
        except ValueError:
            print("⚠️ 请输入有效的数字。")
            
    # 6. 最终确认
    make_sure(list(selected_samples.keys()))
    
    # 7. 运行分析
    run_kraken2(db_dir, selected_samples, threads)

if __name__ == "__main__":
    # 确保遇到 ctrl+c 能优雅退出
    try:
        main()
    except KeyboardInterrupt:
        print("\n\n程序被用户强制中断 (Ctrl+C)。退出运行。")
        sys.exit(0)