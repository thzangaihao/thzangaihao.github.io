'''
转录组分析全自动脚本 - HISAT2基因组比对版本（增强文件查找功能）

目录构型
    RNAseq
        ├── reference
        │   ├── fasta
        │   ├── gtf
        ├── samples
        │   ├── fastq
        ├── script
        │   ├── 此脚本

基于HISAT2流程，只进行基因组比对，输出BAM文件

关键HISAT2参数修改：
--no-spliced-alignment：完全禁用剪接比对，只进行基因组层面的连续比对
--no-softclip：禁用软剪辑，确保比对结果更干净

thz 2025/11/10
'''

import os
import glob
import subprocess
import sys
import re
import time

def get_base_dir():
    """获取基础目录, 兼容打包环境和普通Python环境"""
    if getattr(sys, 'frozen', False):
        # 打包后的情况 - 使用可执行文件所在目录
        return os.path.dirname(sys.executable)
    else:
        # 普通Python脚本运行的情况
        return os.path.dirname(os.path.abspath(__file__))

def find_files(file_type, search_dir=None, recursive=True):
    """灵活查找文件，支持当前目录、子目录和父级目录搜索"""
    if search_dir is None:
        search_dir = get_base_dir()
    
    # 在当前目录及其子目录下搜索
    pattern = os.path.join(search_dir, '**', f'*{file_type}') if recursive else os.path.join(search_dir, f'*{file_type}')
    file_list = glob.glob(pattern, recursive=recursive)
    
    if file_list:
        return file_list
    
    # 如果没有找到，询问是否在父级目录搜索
    print(f"在目录 {search_dir} 中未找到任何{file_type}文件")
    response = input("是否在父级目录搜索? (y/n): ").strip().lower()
    
    if response in ['y', 'yes']:
        parent_dir = os.path.dirname(search_dir)
        if parent_dir != search_dir:  # 避免无限递归
            return find_files(file_type, parent_dir, recursive)
    
    return []

def choose_file(file_list, file_type="文件"):
    """让用户从文件列表中选择文件"""
    if not file_list:
        print(f"没有找到任何{file_type}")
        return None
    
    if len(file_list) == 1:
        selected_file = file_list[0]
        print(f"自动选择: {os.path.basename(selected_file)}")
        return selected_file
    
    print(f"找到 {len(file_list)} 个{file_type}:")
    for i, file in enumerate(file_list, 1):
        print(f"  [{i}] {file}")
    
    while True:
        try:
            choice = input(f"请选择{file_type}编号 (输入q退出): ").strip()
            
            if choice.lower() == 'q':
                return None
            
            choice_index = int(choice) - 1
            
            if 0 <= choice_index < len(file_list):
                selected_file = file_list[choice_index]
                print(f"已选择: {os.path.basename(selected_file)}")
                return selected_file
            else:
                print(f"请输入有效的编号 (1-{len(file_list)})")
        
        except ValueError:
            print("请输入有效的数字")

def check_tools():
    """检查必要的工具是否安装"""
    tools = ['hisat2', 'samtools']
    missing_tools = []
    
    for tool in tools:
        if subprocess.call(['which', tool], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL) != 0:
            missing_tools.append(tool)
    
    if missing_tools:
        print(f"错误: 未找到以下工具: {', '.join(missing_tools)}，请确保已安装")
        sys.exit(1)

def get_reference_files():
    """获取参考基因组文件"""
    print("\n正在查找参考基因组文件...")
    
    # 查找FASTA文件
    fasta_files = find_files('.fasta')
    if not fasta_files:
        fasta_files = find_files('.fa')
        if not fasta_files:
            fasta_files = find_files('.fna')
    
    fasta_file = choose_file(fasta_files, "FASTA文件")
    if not fasta_file:
        print("未选择FASTA文件，程序退出")
        sys.exit(1)
    
    return fasta_file

def select_samples(fastq_type):
    """让用户选择要处理的样本"""
    print(f"\n正在查找{fastq_type}文件...")
    
    # 查找fastq文件
    fastq_files = find_files(fastq_type)
    
    if not fastq_files:
        print(f"错误: 未找到任何{fastq_type}文件")
        return []
    
    # 提取样本名和对应的文件
    sample_files = {}
    for f in fastq_files:
        basename = os.path.basename(f)
        # 匹配常见的双端命名模式
        if '_R1' in basename:
            sample_name = basename.split('_R1')[0]
        elif '_R2' in basename:
            sample_name = basename.split('_R2')[0]
        elif '_1' in basename:
            sample_name = basename.split('_1')[0]
        elif '_2' in basename:
            sample_name = basename.split('_2')[0]
        else:
            # 如果没有明确的R1/R2标记，使用文件名（不含扩展名）作为样本名
            sample_name = re.sub(fr'\{fastq_type}$', '', basename)
        
        if sample_name not in sample_files:
            sample_files[sample_name] = {'r1': None, 'r2': None}
        
        # 确定是R1还是R2文件
        if '_R1' in basename or '_1' in basename:
            sample_files[sample_name]['r1'] = f
        elif '_R2' in basename or '_2' in basename:
            sample_files[sample_name]['r2'] = f
    
    # 显示所有找到的样本
    sample_list = list(sample_files.keys())
    sample_list.sort()
    
    print(f"\n找到 {len(sample_list)} 个样本:")
    for i, sample in enumerate(sample_list, 1):
        r1_file = os.path.basename(sample_files[sample]['r1']) if sample_files[sample]['r1'] else "未找到"
        r2_file = os.path.basename(sample_files[sample]['r2']) if sample_files[sample]['r2'] else "未找到"
        print(f"{i}. {sample}")
        print(f"   R1: {r1_file}")
        print(f"   R2: {r2_file}")
    
    # 让用户选择要处理的样本
    print("\n请选择要处理的样本:")
    print("  [all] - 处理所有样本")
    print("  [1,3,5] - 处理指定编号的样本（用逗号分隔）")
    print("  [1-3] - 处理编号范围的样本（用连字符分隔）")
    
    while True:
        choice = input("\n请输入选择 (输入q退出): ").strip()
        
        if choice.lower() == 'q':
            print("程序已退出")
            exit(0)
        
        if choice.lower() == 'all':
            selected_samples = sample_list
            break
        
        # 处理范围选择 (如: 1-3)
        elif '-' in choice:
            try:
                start, end = map(int, choice.split('-'))
                if 1 <= start <= len(sample_list) and 1 <= end <= len(sample_list) and start <= end:
                    selected_samples = [sample_list[i-1] for i in range(start, end+1)]
                    break
                else:
                    print(f"请输入有效的范围 (1-{len(sample_list)})")
            except ValueError:
                print("请输入有效的范围格式 (如: 1-3)")
        
        # 处理多个选择 (如: 1,3,5)
        elif ',' in choice:
            try:
                indices = [int(x.strip()) for x in choice.split(',')]
                if all(1 <= i <= len(sample_list) for i in indices):
                    selected_samples = [sample_list[i-1] for i in indices]
                    break
                else:
                    print(f"请输入有效的编号 (1-{len(sample_list)})")
            except ValueError:
                print("请输入有效的编号格式 (如: 1,3,5)")
        
        # 处理单个选择
        else:
            try:
                index = int(choice)
                if 1 <= index <= len(sample_list):
                    selected_samples = [sample_list[index-1]]
                    break
                else:
                    print(f"请输入有效的编号 (1-{len(sample_list)})")
            except ValueError:
                print("请输入有效的选择")
    
    # 返回选中的样本及其文件信息
    result = []
    for sample in selected_samples:
        if sample_files[sample]['r1'] and sample_files[sample]['r2']:
            result.append({
                'name': sample,
                'r1': sample_files[sample]['r1'],
                'r2': sample_files[sample]['r2']
            })
        else:
            print(f"警告: 样本 {sample} 的fastq文件不完整，跳过")
    
    return result

def build_hisat2_index(fasta_file, out_dir, threads):
    """构建HISAT2参考基因组索引"""
    if not os.path.exists(fasta_file):
        print("错误: 参考基因组fasta文件不存在")
        sys.exit(1)
    
    # 创建索引目录
    index_dir = os.path.join(out_dir, 'hisat2_index')
    os.makedirs(index_dir, exist_ok=True)
    
    # 索引前缀
    index_prefix = os.path.join(index_dir, 'genome')
    
    # 构建HISAT2索引
    cmd = [
        'hisat2-build',
        '-p', str(threads),
        fasta_file,
        index_prefix
    ]
    
    print("正在构建HISAT2参考索引...")
    print(f"命令: {' '.join(cmd)}")
    
    try:
        subprocess.check_call(cmd)
        print("HISAT2索引构建完成!")
        return index_prefix
    except subprocess.CalledProcessError as e:
        print(f"构建索引时出错: {e}")
        sys.exit(1)

def process_selected_samples(selected_samples, out_dir, threads, index_prefix):
    """处理用户选择的样本 - 只进行基因组比对"""
    if not selected_samples:
        print("没有有效的样本需要处理")
        return
    
    print('='*30)
    print(f"开始处理 {len(selected_samples)} 个选中的样本")

    for sample_info in selected_samples:
        sample = sample_info['name']
        r1 = sample_info['r1']
        r2 = sample_info['r2']
        
        print(f"正在处理样本 {sample}...")
        print(f"  R1: {os.path.basename(r1)}")
        print(f"  R2: {os.path.basename(r2)}")
        
        # 创建样本输出目录
        sample_out_dir = os.path.join(out_dir, f'alignment_results{time.time()}', sample)
        os.makedirs(sample_out_dir, exist_ok=True)
        
        # 输出文件路径
        output_bam = os.path.join(sample_out_dir, f'{sample}.bam')
        sam_file = os.path.join(sample_out_dir, f'{sample}.sam')
        
        # 构建HISAT2比对命令 - 关键修改：禁用剪接比对，只进行基因组比对
        cmd = [
            'hisat2',
            '-x', index_prefix,
            '-1', r1,
            '-2', r2,
            '-p', str(threads),
            '--no-spliced-alignment',  # 禁用剪接比对
            '--no-softclip',           # 禁用软剪辑
            '--time',                  # 输出运行时间
            '--summary-file', os.path.join(sample_out_dir, f'{sample}_summary.txt'),
            '-S', sam_file             # 输出SAM文件
        ]
        
        print(f"  运行HISAT2基因组比对...")
        try:
            # 运行HISAT2比对
            subprocess.check_call(cmd)
            
            # 将SAM转换为BAM并排序
            print(f"  将SAM转换为排序的BAM...")
            
            # 第一步：SAM转BAM
            bam_unsorted = os.path.join(sample_out_dir, f'{sample}_unsorted.bam')
            view_cmd = [
                'samtools', 'view',
                '-@', str(threads),
                '-b',
                '-o', bam_unsorted,
                sam_file
            ]
            subprocess.check_call(view_cmd)
            
            # 第二步：排序BAM
            sort_cmd = [
                'samtools', 'sort',
                '-@', str(threads),
                '-o', output_bam,
                bam_unsorted
            ]
            subprocess.check_call(sort_cmd)
            
            # 第三步：建立索引
            index_cmd = [
                'samtools', 'index',
                output_bam
            ]
            subprocess.check_call(index_cmd)
            
            # 清理临时文件
            if os.path.exists(sam_file):
                os.remove(sam_file)
            if os.path.exists(bam_unsorted):
                os.remove(bam_unsorted)
            
            print(f"  样本 {sample} 处理完成!")
            print(f"  排序BAM文件: {output_bam}")
            print(f"  BAM索引文件: {output_bam}.bai")
            
        except subprocess.CalledProcessError as e:
            print(f"处理样本 {sample} 时出错: {e}")
            # 清理可能产生的临时文件
            if os.path.exists(sam_file):
                os.remove(sam_file)
            continue

def main():
    # 获取基础目录
    base_dir = get_base_dir()
    out_dir = os.path.join(base_dir, f'out{time.time()}')
    
    print(f"基础目录: {base_dir}")
    print(f"输出目录: {out_dir}")
    print("注意: 此版本只进行基因组比对，不进行剪接比对")
    
    # 创建输出目录
    os.makedirs(out_dir, exist_ok=True)
    
    # 获取用户输入
    print("\n请选择fastq文件格式:")
    print("  [1]: .fastq")
    print("  [2]: .fastq.gz")
    print("  [3]: .fq")
    print("  [4]: .fq.gz")
    
    fastq_type_in = input("请输入选择 (1-4): ").strip()
    if fastq_type_in == '1':
        fastq_type = '.fastq'
    elif fastq_type_in == '2':
        fastq_type = '.fastq.gz'
    elif fastq_type_in == '3':
        fastq_type = '.fq'
    elif fastq_type_in == '4':
        fastq_type = '.fq.gz'
    else:
        print('请输入正确序号 (1, 2, 3 或 4)')
        return
    
    # 获取线程数
    try:
        threads = int(input('请确定运算核心数: '))
        if threads <= 0:
            print("线程数必须大于0")
            return
    except ValueError:
        print("请输入有效的数字")
        return
    
    # 检查工具
    check_tools()
    
    # 获取参考基因组文件
    fasta_file = get_reference_files()
    
    # 让用户选择要处理的样本
    selected_samples = select_samples(fastq_type)
    
    if not selected_samples:
        print("没有选择任何样本，程序退出")
        return
    
    # 构建HISAT2索引
    print("\n正在构建基因组索引...")
    index_prefix = build_hisat2_index(fasta_file, out_dir, threads)
    
    # 处理用户选择的样本
    print("\n开始基因组比对...")
    process_selected_samples(selected_samples, out_dir, threads, index_prefix)
    
    print("\n" + "="*50)
    print("HISAT2基因组比对完成！")
    print("输出特点:")
    print("  - 只进行基因组层面比对")
    print("  - 禁用剪接比对")
    print("  - 输出排序的BAM文件和索引")
    print("  - 适合IGV可视化")
    print(f"索引文件位置: {index_prefix}.*")
    print(f"比对结果位置: {os.path.join(out_dir, 'alignment_results')}")
    print("="*50)

if __name__ == '__main__':
    main()