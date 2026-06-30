'''
fastq文件合并脚本

目录构型：
    自动识别脚本同级目录及子目录中的.fastq文件

用于DNA_annotation工作流
所有样本的R1端合并为一个总R1文件, R2端合并为一个总R2文件
支持.fastq, .fq及对应的.gz压缩格式

thz 2025/11/19
'''

import glob
import os
import subprocess

def find_fastq_files():
    """查找所有FASTQ文件"""
    patterns = [
        ".fastq", ".fq", 
        ".fastq.gz", ".fq.gz"
    ]
    
    fastq_files_R1 = []
    for pattern in patterns:
        # 搜索当前目录和所有子目录
        current_path = os.path.dirname(os.path.abspath(__file__))
        path = os.path.join(current_path, '**', f'*_R1{pattern}')
        files = glob.glob(path, recursive=True)
        fastq_files_R1.extend(files)

    fastq_files_R2 = []
    for pattern in patterns:
        # 搜索当前目录和所有子目录
        current_path = os.path.dirname(os.path.abspath(__file__))
        path = os.path.join(current_path, '**', f'*_R2{pattern}')
        files = glob.glob(path, recursive=True)
        fastq_files_R2.extend(files)
    
    find_sample(fastq_files_R1, fastq_files_R2)
    return fastq_files_R1, fastq_files_R2

def find_sample(R1, R2):
    """总结所有样品"""
    print(f'R1共发现{len(R1)}个样品：')
    for i in R1:
        print(f'  {i}')
    print(f'R2共发现{len(R2)}个样品：')
    for i in R2:
        print(f'  {i}')

def extract_sample_name(filename):
    """从文件名中提取样本名"""
    basename = os.path.basename(filename)
    # 去除_R1或_R2及后续部分
    if '_R1' in basename:
        return basename.split('_R1')[0]
    elif '_R2' in basename:
        return basename.split('_R2')[0]
    else:
        return os.path.splitext(basename)[0]
    
def group_files_by_sample(R1_files, R2_files):
    """按样本名分组文件"""
    samples = {}
    
    for r1_file in R1_files:
        sample_name = extract_sample_name(r1_file)
        if sample_name not in samples:
            samples[sample_name] = {'R1': [], 'R2': []}
        samples[sample_name]['R1'].append(r1_file)
    
    for r2_file in R2_files:
        sample_name = extract_sample_name(r2_file)
        if sample_name not in samples:
            samples[sample_name] = {'R1': [], 'R2': []}
        samples[sample_name]['R2'].append(r2_file)
    
    return samples

def select_samples(samples):
    """让用户选择要合并的样本"""
    print("\n" + "=" * 50)
    print("请选择要合并的样本:")
    print("-" * 50)
    
    # 显示所有样本列表
    sample_list = list(samples.keys())
    for i, sample_name in enumerate(sample_list, 1):
        r1_count = len(samples[sample_name]['R1'])
        r2_count = len(samples[sample_name]['R2'])
        print(f"{i}. {sample_name} (R1: {r1_count}个文件, R2: {r2_count}个文件)")
    
    print(f"{len(sample_list) + 1}. 合并所有样本")
    
    while True:
        try:
            choice = input("\n请选择样本编号(多个编号用逗号分隔，如: 1,3,5 或输入'all'合并所有): ").strip()
            
            if choice.lower() == 'all':
                return samples  # 返回所有样本
            
            # 处理用户输入
            selected_indices = [int(x.strip()) for x in choice.split(',')]
            
            # 验证输入
            valid_selection = True
            selected_samples = {}
            
            for idx in selected_indices:
                if 1 <= idx <= len(sample_list):
                    sample_name = sample_list[idx - 1]
                    selected_samples[sample_name] = samples[sample_name]
                else:
                    print(f"错误: 编号 {idx} 无效，请输入1到{len(sample_list) + 1}之间的数字")
                    valid_selection = False
                    break
            
            if valid_selection and selected_samples:
                return selected_samples
            elif not selected_samples:
                print("错误: 未选择任何样本，请重新选择")
                
        except ValueError:
            print("错误: 请输入有效的数字编号，如: 1,3,5 或输入'all'")

def merge_fastq_files(samples, output_dir="merged_fastq"):
    """合并FASTQ文件"""
    
    # 创建输出目录
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
        print(f"创建输出目录: {output_dir}")

    print("=" * 50)    
    print("开始合并FASTQ文件...")

    # 构建输出文件
    output_r1 = os.path.join(output_dir, "merged_R1.fastq.gz")
    output_r2 = os.path.join(output_dir, "merged_R2.fastq.gz")
    # 构建待合并列表
    merged_list_R1 = []
    merged_list_R2 = []
    
    for sample_name, files in samples.items():
        print(f"\n样本: {sample_name}")
        print("-" * 30)
        
        r1_files = files['R1']
        r2_files = files['R2']
        
        print(f"R1文件 ({len(r1_files)}个):")
        for f in r1_files:
            print(f"  {f}")
            merged_list_R1.append(f)
        
        print(f"R2文件 ({len(r2_files)}个):")
        for f in r2_files:
            print(f"  {f}")
            merged_list_R2.append(f)
        
        # 检查文件数量是否匹配
        if len(r1_files) != len(r2_files):
            print(f"⚠ 警告: R1和R2文件数量不匹配 ({len(r1_files)} vs {len(r2_files)})，但仍然执行合并")

    # 合并所有R1文件
    r1_cmd = f"cat {' '.join(merged_list_R1)} > {output_r1}"
    if run_command(r1_cmd, f"合并R1文件 -> {output_r1}"):
        # 合并所有R2文件
        r2_cmd = f"cat {' '.join(merged_list_R2)} > {output_r2}"
        run_command(r2_cmd, f"合并R2文件 -> {output_r2}")

    # 验证合并结果
    if os.path.exists(output_r1) and os.path.exists(output_r2):
        print("验证合并结果:")
        # 检查文件大小
        r1_size = os.path.getsize(output_r1) / (1024*1024)  # MB
        r2_size = os.path.getsize(output_r2) / (1024*1024)  # MB
        print(f"  合并后R1文件大小: {r1_size:.2f} MB")
        print(f"  合并后R2文件大小: {r2_size:.2f} MB")

def run_command(cmd, description):
    """运行命令行命令"""
    print(f"执行: {description}")
    print(f"命令: {cmd}")
    try:
        subprocess.run(cmd, shell=True, check=True, 
            stdout=subprocess.PIPE, stderr=subprocess.PIPE,
            text=True)
        print("✓ 成功完成\n")
        print("=" * 50)
        print('合并成功')
        return True
    except subprocess.CalledProcessError as e:
        print(f"✗ 执行失败: {e}")
        print(f"错误输出: {e.stderr}")
        print("=" * 50)
        print('合并失败, 请检查linux运行命令')        
        return False

def main():
    # 样本检测
    fastq_files_R1, fastq_files_R2 = find_fastq_files()
    samples = group_files_by_sample(fastq_files_R1, fastq_files_R2)
    print(f"\n共发现 {len(samples)} 个样本:")
    for sample_name in samples.keys():
        print(f"  {sample_name}")
    
    # 让用户选择要合并的样本
    selected_samples = select_samples(samples)
    
    print(f"\n已选择 {len(selected_samples)} 个样本进行合并:")
    for sample_name in selected_samples.keys():
        print(f"  {sample_name}")
    
    # 样本合并
    response = input("\n是否继续合并? (y/n): ").strip().lower()
    if response not in ['y', 'yes']:
        print("操作已取消")
        return
    merge_fastq_files(selected_samples)

main()