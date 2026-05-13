import os
import subprocess
import sys
from datetime import datetime

# ==========================================
# 1. 基础配置
# ==========================================
THREADS = int(input("请输入线程数 (如: 128) >>> "))  # 充分利用你的服务器性能

def log_info(message):
    print(f"[{datetime.now().strftime('%H:%M:%S')}] {message}")

# ==========================================
# 2. 执行回比对
# ==========================================
if __name__ == "__main__":
    print("="*60)
    print("      HiFi Reads 回比对与覆盖度准备工具")
    print("="*60)

    # 输入：组装好的 FASTA 和 原始 FASTQ
    ref_fasta = input("请输入组装好的 FASTA 路径 (如: sample.p_ctg.fasta) >>> ").strip()
    reads_fq = input("请输入原始 HiFi Reads 路径 (如: merged.fastq.gz) >>> ").strip()
    
    if not os.path.exists(ref_fasta) or not os.path.exists(reads_fq):
        print("路径不存在，请检查！")
        sys.exit(1)

    sample_name = os.path.splitext(os.path.basename(ref_fasta))[0]
    output_bam = f"{sample_name}.sorted.bam"

    # 核心命令流：
    # 1. minimap2 比对 (-ax map-hifi 是针对 HiFi 的优化参数)
    # 2. samtools view 转为 BAM
    # 3. samtools sort 排序 (IGV 必须要求排序过的 BAM)
    log_info("开始进行 Minimap2 比对...")
    cmd = (
        f"minimap2 -ax map-hifi -t {THREADS} {ref_fasta} {reads_fq} | "
        f"samtools view -@ {THREADS} -bS - | "
        f"samtools sort -@ {THREADS} -o {output_bam}"
    )
    
    try:
        subprocess.run(cmd, shell=True, check=True)
        log_info(f"排序后的 BAM 已生成: {output_bam}")
        
        # 4. 创建索引 (IGV 必须要求 .bai 文件)
        log_info("正在创建 BAM 索引...")
        subprocess.run(f"samtools index -@ {THREADS} {output_bam}", shell=True, check=True)
        log_info("索引文件 .bai 已生成。")
        
    except subprocess.CalledProcessError as e:
        log_info(f"运行失败: {e}")
        sys.exit(1)

    print("-" * 60)
    log_info("恭喜！现在你可以将以下两个文件下载并拖入 IGV 了：")
    print(f"1. {ref_fasta} (作为 Genome 导入)")
    print(f"2. {output_bam} (作为 Data Track 导入)")