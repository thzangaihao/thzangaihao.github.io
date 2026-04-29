import os
import glob
import time
import subprocess
import sys
from datetime import datetime

# ==========================================
# 1. 基础配置
# ==========================================
THREADS = 128  # 默认线程数，可在此修改

def log_info(message):
    """打印带时间戳的日志信息"""
    current_time = time.strftime('%Y-%m-%d %H:%M:%S')
    print(f"[{current_time}] {message}")

def check_dependencies():
    """检查必要的生信软件是否在环境变量中"""
    for cmd in ["samtools", "hifiasm"]:
        if subprocess.call(f"type {cmd}", shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE) != 0:
            log_info(f"严重错误：系统未找到命令 '{cmd}'。请先激活你的 conda/uv 环境或加载对应模块。")
            sys.exit(1)

# ==========================================
# 2. 核心任务函数
# ==========================================
def bam_to_fastq(bam_file, fastq_file, threads):
    """将 BAM 转换为 FASTQ"""
    log_info(f"  -> 开始格式转换: BAM to FASTQ")
    # 使用 samtools 提取 fastq，-T 重定向一些可能需要的标签，通常直接提取即可
    cmd = ["samtools", "fastq", "-@", str(threads), bam_file]
    
    try:
        with open(fastq_file, "w") as out_fq:
            subprocess.run(cmd, stdout=out_fq, stderr=subprocess.DEVNULL, check=True)
        log_info(f"  -> 格式转换完成！已生成: {fastq_file}")
    except subprocess.CalledProcessError as e:
        log_info(f"  -> 错误：BAM 转换失败！{e}")
        raise

def run_hifiasm(fastq_file, prefix, threads):
    """运行 Hifiasm"""
    log_info(f"  -> 开始运行 hifiasm 进行序列拼接...")
    cmd = ["hifiasm", "-o", prefix, "-t", str(threads), fastq_file]
    
    try:
        subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        log_info(f"  -> hifiasm 组装核心步骤完成！")
    except subprocess.CalledProcessError as e:
        log_info(f"  -> 错误：hifiasm 运行失败！{e}")
        raise

def extract_fasta(gfa_file, fasta_file):
    """从 GFA 提取主组装 FASTA"""
    log_info(f"  -> 正在提取最终 FASTA 序列...")
    if not os.path.exists(gfa_file):
        log_info(f"  -> 警告：找不到图文件 {gfa_file}，跳过提取。")
        return

    try:
        contig_count = 0
        with open(gfa_file, 'r') as gfa, open(fasta_file, 'w') as fasta:
            for line in gfa:
                if line.startswith('S'):
                    parts = line.strip().split('\t')
                    if len(parts) >= 3:
                        fasta.write(f">{parts[1]}\n{parts[2]}\n")
                        contig_count += 1
        log_info(f"  -> 提取完成！共获得 {contig_count} 条 Contigs。")
    except Exception as e:
        log_info(f"  -> 错误：提取 FASTA 时发生异常：{e}")

# ==========================================
# 3. 主控制流
# ==========================================
if __name__ == "__main__":
    print("="*60)
    print("      自动化 HiFi 基因组组装队列系统 (BAM -> FASTQ -> FASTA)")
    print("="*60)

    check_dependencies()

    # 1. 扫描所有的 BAM 文件
    log_info("正在扫描当前目录及子目录下的 .bam 文件...")
    bam_files = glob.glob("**/*.bam", recursive=True)

    if not bam_files:
        log_info("未找到任何 .bam 文件，脚本退出。")
        sys.exit(0)

    # 2. 交互式选择菜单
    print("\n发现以下待处理文件：")
    for i, file_path in enumerate(bam_files, 1):
        print(f"  [{i}] {file_path}")
    
    print("\n请选择要组装的文件编号 (例如: 1,3,4)，或者输入 'all' 处理全部：")
    choice = input("你的选择 >>> ").strip().lower()

    selected_files = []
    if choice == 'all':
        selected_files = bam_files
    else:
        try:
            indices = [int(x.strip()) - 1 for x in choice.split(',')]
            selected_files = [bam_files[i] for i in indices if 0 <= i < len(bam_files)]
        except Exception:
            log_info("输入格式有误，请输入用逗号分隔的数字，脚本退出。")
            sys.exit(1)

    if not selected_files:
        log_info("未选择任何有效文件，脚本退出。")
        sys.exit(0)

    # 3. 创建带时间戳的顶级工作目录
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    work_dir = os.path.abspath(f"assemble_{timestamp}")
    os.makedirs(work_dir, exist_ok=True)
    log_info(f"创建总输出目录: {work_dir}\n")

    # 4. 执行队列任务
    for index, bam_path in enumerate(selected_files, 1):
        sample_name = os.path.splitext(os.path.basename(bam_path))[0]
        # 去掉可能的 .hifi_reads 后缀让名字更清爽
        sample_name = sample_name.replace(".hifi_reads", "") 
        
        print("-" * 60)
        log_info(f"任务 [{index}/{len(selected_files)}] 开始处理样本: {sample_name}")

        # 为该样本创建专属的三个子目录
        sample_dir = os.path.join(work_dir, sample_name)
        dir_fastq = os.path.join(sample_dir, "01_fastq")
        dir_hifiasm = os.path.join(sample_dir, "02_hifiasm")
        dir_fasta = os.path.join(sample_dir, "03_fasta")
        
        for d in [dir_fastq, dir_hifiasm, dir_fasta]:
            os.makedirs(d, exist_ok=True)

        # 配置文件路径
        out_fastq = os.path.join(dir_fastq, f"{sample_name}.fastq")
        prefix_hifiasm = os.path.join(dir_hifiasm, sample_name)
        primary_gfa = f"{prefix_hifiasm}.bp.p_ctg.gfa"
        final_fasta = os.path.join(dir_fasta, f"{sample_name}.p_ctg.fasta")

        # 依次执行管线
        try:
            bam_to_fastq(bam_path, out_fastq, THREADS)
            run_hifiasm(out_fastq, prefix_hifiasm, THREADS)
            extract_fasta(primary_gfa, final_fasta)
            log_info(f"样本 {sample_name} 全部处理完毕！结果位于: {dir_fasta}")
        except Exception as e:
            log_info(f"样本 {sample_name} 处理失败，跳过该样本进入下一个。")
            continue

    print("-" * 60)
    log_info(f"所有队列任务执行完毕！数据均保存在: {work_dir}")