import os
import glob
import time
import subprocess
import sys
from datetime import datetime

# ==========================================
# 1. 基础配置
# ==========================================
THREADS = 128  # 默认线程数，充分利用集群性能

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
    cmd = ["samtools", "fastq", "-@", str(threads), bam_file]
    try:
        with open(fastq_file, "w") as out_fq:
            subprocess.run(cmd, stdout=out_fq, stderr=subprocess.DEVNULL, check=True)
        log_info(f"  -> 格式转换完成！已生成临时 FASTQ。")
    except subprocess.CalledProcessError as e:
        log_info(f"  -> 错误：BAM 转换失败！{e}")
        raise

def run_hifiasm(fastq_file, prefix, threads, ploidy_mode):
    """根据倍性选择运行 Hifiasm"""
    log_info(f"  -> 开始运行 hifiasm 进行序列拼接...")
    
    # 根据用户选择动态生成参数
    if ploidy_mode in ['1', '2']:
        # 模式1和2使用 --primary 生成 p_ctg 和 a_ctg
        cmd = ["hifiasm", "-o", prefix, "-t", str(threads), "--primary", fastq_file]
    else:
        # 模式3默认行为，生成 hap1 和 hap2
        cmd = ["hifiasm", "-o", prefix, "-t", str(threads), fastq_file]
        
    try:
        log_info(f"  -> 执行命令: {' '.join(cmd)}")
        subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        log_info(f"  -> hifiasm 组装核心步骤完成！")
    except subprocess.CalledProcessError as e:
        log_info(f"  -> 错误：hifiasm 运行失败！{e}")
        raise

def extract_fasta(gfa_file, fasta_file):
    """从 GFA 提取主组装 FASTA"""
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
        log_info(f"  -> 成功从 {os.path.basename(gfa_file)} 提取 {contig_count} 条 Contigs。")
    except Exception as e:
        log_info(f"  -> 错误：提取 FASTA 时发生异常：{e}")

# ==========================================
# 3. 主控制流
# ==========================================
if __name__ == "__main__":
    print("="*60)
    print("   自动化 HiFi 基因组组装队列系统 (支持多倍型/单倍体)")
    print("="*60)

    check_dependencies()

    # 1. 扫描文件
    log_info("正在扫描当前目录及子目录下的序列文件...")
    extensions = ["**/*.bam", "**/*.fastq", "**/*.fastq.gz", "**/*.fq", "**/*.fq.gz"]
    all_files = []
    for ext in extensions:
        all_files.extend(glob.glob(ext, recursive=True))
    
    all_files = sorted(list(set([f for f in all_files if not os.path.basename(f).startswith(".")])))

    if not all_files:
        log_info("未找到任何 .bam 或 FASTQ 相关文件，脚本退出。")
        sys.exit(0)

    # 2. 交互式选择样本
    print("\n发现以下待处理文件：")
    for i, file_path in enumerate(all_files, 1):
        print(f"  [{i}] {file_path}")
    
    print("\n请选择要组装的文件编号 (例如: 1,3,4)，或输入 'all' 处理全部：")
    choice = input("你的选择 >>> ").strip().lower()

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
        log_info("未选择任何有效文件，脚本退出。")
        sys.exit(0)

    # 3. 交互式选择组装策略
    print("\n请选择物种倍性及组装策略：")
    print("  [1] 单倍体模式 (仅提取 Primary Contig，适合真菌等单倍体物种)")
    print("  [2] 主副组装模式 (生成 Primary/Alternate，适合高杂合二倍体构建单一参考基因组)")
    print("  [3] 默认双倍体分型模式 (生成 Hap1/Hap2，适合植物二倍体全相态组装)")
    ploidy_mode = input("策略选择 (1/2/3) >>> ").strip()
    
    if ploidy_mode not in ['1', '2', '3']:
        log_info("无效输入，默认回退至 [3] 双倍体分型模式。")
        ploidy_mode = '3'

    # 4. 创建工作目录
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    work_dir = os.path.abspath(f"assemble_{timestamp}")
    os.makedirs(work_dir, exist_ok=True)
    log_info(f"创建总输出目录: {work_dir}\n")

    # 5. 执行队列任务
    for index, file_path in enumerate(selected_files, 1):
        base_name = os.path.basename(file_path)
        for ext in [".bam", ".fastq.gz", ".fastq", ".fq.gz", ".fq", ".hifi_reads"]:
            base_name = base_name.replace(ext, "")
        sample_name = base_name
        
        print("-" * 60)
        log_info(f"任务 [{index}/{len(selected_files)}] 开始处理样本: {sample_name}")

        sample_dir = os.path.join(work_dir, sample_name)
        dir_fastq = os.path.join(sample_dir, "01_fastq")
        dir_hifiasm = os.path.join(sample_dir, "02_hifiasm")
        dir_fasta = os.path.join(sample_dir, "03_fasta")
        
        for d in [dir_fastq, dir_hifiasm, dir_fasta]:
            os.makedirs(d, exist_ok=True)

        prefix_hifiasm = os.path.join(dir_hifiasm, sample_name)

        try:
            input_for_hifiasm = file_path
            
            if file_path.lower().endswith(".bam"):
                out_fastq = os.path.join(dir_fastq, f"{sample_name}.fastq")
                bam_to_fastq(file_path, out_fastq, THREADS)
                input_for_hifiasm = out_fastq
            else:
                log_info("  -> 输入已是 FASTQ 格式，智能跳过格式转换。")
                symlink_path = os.path.join(dir_fastq, os.path.basename(file_path))
                try:
                    os.symlink(os.path.abspath(file_path), symlink_path)
                except FileExistsError:
                    pass

            # 运行 Hifiasm
            run_hifiasm(input_for_hifiasm, prefix_hifiasm, THREADS, ploidy_mode)
            
            # 动态提取 GFA
            log_info(f"  -> 正在根据选择的策略提取最终 FASTA 序列...")
            extract_targets = []
            
            if ploidy_mode == '1':
                extract_targets.append(("bp.p_ctg.gfa", "p_ctg.fasta"))
            elif ploidy_mode == '2':
                extract_targets.append(("bp.p_ctg.gfa", "p_ctg.fasta"))
                extract_targets.append(("bp.a_ctg.gfa", "a_ctg.fasta"))
            elif ploidy_mode == '3':
                extract_targets.append(("bp.hap1.p_ctg.gfa", "hap1.fasta"))
                extract_targets.append(("bp.hap2.p_ctg.gfa", "hap2.fasta"))

            for gfa_suffix, fasta_suffix in extract_targets:
                gfa_file = f"{prefix_hifiasm}.{gfa_suffix}"
                fasta_file = os.path.join(dir_fasta, f"{sample_name}.{fasta_suffix}")
                extract_fasta(gfa_file, fasta_file)

            log_info(f"样本 {sample_name} 全部处理完毕！结果位于: {dir_fasta}")
            
        except Exception as e:
            log_info(f"样本 {sample_name} 处理失败，跳过该样本进入下一个。错误信息：{e}")
            continue

    print("-" * 60)
    log_info(f"所有队列任务执行完毕！数据均保存在: {work_dir}")