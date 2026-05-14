import os
import glob
import time
import subprocess
import sys
from datetime import datetime

# ==========================================
# 基础配置
# ==========================================
THREADS = 128         # 默认线程数
MEMORY = "256g"       # 允许 Canu 使用的最大内存，请根据服务器配置调整
GENOME_SIZE = "40m"   # 预估基因组大小。真菌通常在 30m-50m 左右，请根据具体菌株调整

def log_info(message):
    """打印带时间戳的日志信息"""
    current_time = time.strftime('%Y-%m-%d %H:%M:%S')
    print(f"[{current_time}] {message}")

def check_dependencies():
    """检查必要的生信软件"""
    for cmd in ["samtools", "canu"]:
        if subprocess.call(f"type {cmd}", shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE) != 0:
            log_info(f"严重错误：未找到命令 '{cmd}'。请确认已激活 conda 环境。")
            sys.exit(1)

# ==========================================
# 核心任务函数
# ==========================================
def bam_to_fastq(bam_file, fastq_file, threads):
    """将 BAM 转换为 FASTQ"""
    log_info(f"  -> 开始格式转换: BAM to FASTQ")
    cmd = ["samtools", "fastq", "-@", str(threads), bam_file]
    try:
        with open(fastq_file, "w") as out_fq:
            subprocess.run(cmd, stdout=out_fq, stderr=subprocess.DEVNULL, check=True)
        log_info(f"  -> 格式转换完成！已生成临时 FASTQ。")
    except subprocess.CalledProcessError as e:
        log_info(f"  -> 错误：BAM 转换失败！{e}")
        raise

def run_hicanu(fastq_file, output_dir, prefix, threads, memory, genome_size):
    """运行 HiCanu 进行组装"""
    log_info(f"  -> 开始运行 HiCanu 流程...")
    
    # Canu 的核心命令构建
    cmd = [
        "canu",
        f"-p", prefix,
        f"-d", output_dir,
        f"genomeSize={genome_size}",
        f"maxThreads={threads}",
        f"maxMemory={memory}",
        "useGrid=false",       # 强制在当前节点运行，避免自动提交到 Slurm/SGE 系统
        "-pacbio-hifi", fastq_file
    ]
    
    log_info(f"  -> 执行命令: {' '.join(cmd)}")
    try:
        # Canu 自身会输出非常详细的日志到控制台，这里保留输出以便监控进度
        subprocess.run(cmd, check=True)
        log_info(f"  -> HiCanu 组装流程顺利完成！")
    except subprocess.CalledProcessError as e:
        log_info(f"  -> 错误：HiCanu 运行失败！{e}")
        raise

# ==========================================
# 主控制流
# ==========================================
if __name__ == "__main__":
    print("="*60)
    print("   自动化 HiCanu 基因组组装队列系统 (支持 BAM/FASTQ 直投)")
    print("="*60)

    check_dependencies()

    log_info("正在扫描当前目录及子目录下的序列文件...")
    extensions = ["**/*.bam", "**/*.fastq", "**/*.fastq.gz", "**/*.fq", "**/*.fq.gz"]
    all_files = []
    for ext in extensions:
        all_files.extend(glob.glob(ext, recursive=True))
    
    all_files = sorted(list(set([f for f in all_files if not os.path.basename(f).startswith(".")])))

    if not all_files:
        log_info("未找到任何 .bam 或 FASTQ 相关文件，脚本退出。")
        sys.exit(0)

    print("\n发现以下待处理文件：")
    for i, file_path in enumerate(all_files, 1):
        print(f"  [{i}] {file_path}")
    
    choice = input("\n请选择要组装的文件编号 (例如: 1,3,4)，或者输入 'all' 处理全部：\n>>> ").strip().lower()

    selected_files = []
    if choice == 'all':
        selected_files = all_files
    else:
        try:
            indices = [int(x.strip()) - 1 for x in choice.split(',')]
            selected_files = [all_files[i] for i in indices if 0 <= i < len(all_files)]
        except Exception:
            log_info("输入格式有误，脚本退出。")
            sys.exit(1)

    if not selected_files:
        sys.exit(0)

    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    work_dir = os.path.abspath(f"hicanu_assemble_{timestamp}")
    os.makedirs(work_dir, exist_ok=True)
    log_info(f"创建总输出目录: {work_dir}\n")

    for index, file_path in enumerate(selected_files, 1):
        base_name = os.path.basename(file_path)
        for ext in [".bam", ".fastq.gz", ".fastq", ".fq.gz", ".fq", ".hifi_reads"]:
            base_name = base_name.replace(ext, "")
        sample_name = base_name
        
        print("-" * 60)
        log_info(f"任务 [{index}/{len(selected_files)}] 开始处理样本: {sample_name}")

        sample_dir = os.path.join(work_dir, sample_name)
        dir_fastq = os.path.join(sample_dir, "01_fastq")
        dir_canu_out = os.path.join(sample_dir, "02_hicanu_out")
        
        for d in [dir_fastq, dir_canu_out]:
            os.makedirs(d, exist_ok=True)

        try:
            input_for_canu = file_path
            
            if file_path.lower().endswith(".bam"):
                out_fastq = os.path.join(dir_fastq, f"{sample_name}.fastq")
                bam_to_fastq(file_path, out_fastq, THREADS)
                input_for_canu = out_fastq
            else:
                log_info("  -> 输入已是 FASTQ 格式。")
                symlink_path = os.path.join(dir_fastq, os.path.basename(file_path))
                try:
                    os.symlink(os.path.abspath(file_path), symlink_path)
                except FileExistsError:
                    pass

            # 运行 HiCanu
            run_hicanu(input_for_canu, dir_canu_out, sample_name, THREADS, MEMORY, GENOME_SIZE)
            
            # Canu 会自动生成 .contigs.fasta
            final_fasta = os.path.join(dir_canu_out, f"{sample_name}.contigs.fasta")
            if os.path.exists(final_fasta):
                log_info(f"样本 {sample_name} 成功组装！最终序列文件位于: {final_fasta}")
            else:
                log_info(f"警告：未在预期位置找到 {final_fasta}，请检查 Canu 日志。")
            
        except Exception as e:
            log_info(f"样本 {sample_name} 处理失败，跳过。错误信息：{e}")
            continue

    print("-" * 60)
    log_info(f"所有队列任务执行完毕！数据均保存在: {work_dir}")