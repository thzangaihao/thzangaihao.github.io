import os
import glob
import time
import subprocess
import sys
from datetime import datetime

# ==========================================
# 1. 基础配置
# ==========================================
THREADS = 128  # 充分利用你的服务器性能

def log_info(message):
    current_time = time.strftime('%Y-%m-%d %H:%M:%S')
    print(f"[{current_time}] {message}")

def check_samtools():
    """检查 samtools 是否可用"""
    if subprocess.call("type samtools", shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE) != 0:
        log_info("严重错误：未找到 'samtools'。请先加载环境。")
        sys.exit(1)

# ==========================================
# 2. 核心转换函数
# ==========================================
def convert_bam_to_fastq(bam_path, out_dir):
    """
    将 BAM 转换为 .fastq.gz
    """
    sample_name = os.path.splitext(os.path.basename(bam_path))[0]
    sample_name = sample_name.replace(".hifi_reads", "")
    
    # 设定输出路径
    output_fq = os.path.join(out_dir, f"{sample_name}.fastq.gz")
    
    log_info(f"正在处理: {sample_name}")
    
    # 构建命令：使用 samtools 提取并直接通过重定向/管道压缩
    # -c 6 代表压缩级别，-@ 代表线程
    cmd = f"samtools fastq -@ {THREADS} {bam_path} | gzip -c6 > {output_fq}"
    
    try:
        # 使用 shell=True 是为了支持管道操作
        subprocess.run(cmd, shell=True, check=True)
        log_info(f"转换成功 -> {output_fq}")
    except subprocess.CalledProcessError:
        log_info(f"错误：样本 {sample_name} 转换失败！")

# ==========================================
# 3. 交互流程
# ==========================================
if __name__ == "__main__":
    print("="*60)
    print("         PacBio BAM -> FASTQ.GZ 批量转换工具")
    print("="*60)

    check_samtools()

    # 1. 扫描文件
    bam_files = sorted(glob.glob("**/*.bam", recursive=True))
    if not bam_files:
        log_info("当前目录下没发现 .bam 文件。")
        sys.exit(0)

    print("\n发现以下待转换文件：")
    for i, file_path in enumerate(bam_files, 1):
        print(f"  [{i}] {file_path}")
    
    # 2. 交互选择
    choice = input("\n请选择编号 (例如 1,3) 或输入 'all' >>> ").strip().lower()
    selected = []
    if choice == 'all':
        selected = bam_files
    else:
        try:
            indices = [int(x.strip()) - 1 for x in choice.split(',')]
            selected = [bam_files[i] for i in indices]
        except:
            print("输入无效。")
            sys.exit(1)

    # 3. 创建输出目录
    timestamp = datetime.now().strftime("%Y%m%d_%H%M")
    work_dir = os.path.abspath(f"fastq_outputs_{timestamp}")
    os.makedirs(work_dir, exist_ok=True)
    log_info(f"总输出目录: {work_dir}\n")

    # 4. 循环处理
    start_time = time.time()
    for bam in selected:
        convert_bam_to_fastq(bam, work_dir)
    
    end_time = time.time()
    duration = (end_time - start_time) / 60
    print("-" * 60)
    log_info(f"任务全部完成！总耗时: {duration:.2f} 分钟。")