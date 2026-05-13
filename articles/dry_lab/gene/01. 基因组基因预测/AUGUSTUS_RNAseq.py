import os
import glob
import time
import subprocess
import sys
from datetime import datetime

# ==========================================
# 1. 基础配置与辅助函数
# ==========================================
def log_info(message):
    current_time = time.strftime('%Y-%m-%d %H:%M:%S')
    print(f"[{current_time}] {message}")

def check_dependencies():
    for cmd in ["augustus", "bam2hints"]:
        if subprocess.call(f"type {cmd}", shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE) != 0:
            log_info(f"严重错误：系统未找到命令 '{cmd}'。请确认已加载环境。")
            sys.exit(1)

# ==========================================
# 2. 交互式选择菜单
# ==========================================
def select_file(prompt_msg, extensions, allow_skip=False):
    """通用的交互式文件选择函数"""
    files = []
    for ext in extensions:
        files.extend(glob.glob(ext, recursive=True))
    files = sorted(list(set([f for f in files if not os.path.basename(f).startswith(".")])))

    if not files:
        print(f"\n未找到相关文件: {extensions}")
        if allow_skip:
            return None
        sys.exit(0)

    print(f"\n{prompt_msg}")
    if allow_skip:
        print("  [0] 不选择 (跳过此步骤 / 纯算法预测)")
        
    for i, f in enumerate(files, 1):
        # 如果是 BAM 文件，顺便打印下大小
        size_str = ""
        if f.endswith(".bam"):
            size_mb = os.path.getsize(f) / (1024 * 1024)
            size_str = f" ({size_mb:.1f} MB)"
        print(f"  [{i}] {f}{size_str}")

    while True:
        choice = input("\n请输入编号 >>> ").strip()
        if allow_skip and choice == '0':
            return None
        try:
            idx = int(choice) - 1
            if 0 <= idx < len(files):
                return files[idx]
            else:
                print("编号超出范围，请重试。")
        except ValueError:
            print("请输入有效的数字编号。")

# ==========================================
# 3. 主控制流
# ==========================================
if __name__ == "__main__":
    print("="*60)
    print("   AUGUSTUS 单样本精准交互预测 (支持 RNA-seq 共预测)")
    print("="*60)

    check_dependencies()

    # --- 步骤 1: 选择 FASTA 基因组 ---
    fasta_file = select_file(
        "发现以下基因组文件，请选择你要预测的 FASTA:",
        ["**/*.fasta", "**/*.fa", "**/*.fna", "**/*.p_ctg.fasta"]
    )
    
    # 提取干净的样本名
    sample_name = os.path.splitext(os.path.basename(fasta_file))[0]
    for ext in [".p_ctg", ".ctg", ".hifi"]:
        sample_name = sample_name.replace(ext, "")
    log_info(f"已选择目标基因组: {fasta_file}")

    # --- 步骤 2: 选择 BAM 转录组比对文件 (可选) ---
    bam_file = select_file(
        "发现以下 BAM 文件，请选择对应的转录组比对文件:",
        ["**/*.bam", "**/*.sorted.bam"],
        allow_skip=True
    )
    if bam_file:
        log_info(f"已选择转录组证据: {bam_file}")
    else:
        log_info("已跳过转录组证据，将执行纯算法从头预测。")

    # --- 步骤 3: 设定物种模型 ---
    print("\n" + "-"*40)
    print("请输入 AUGUSTUS 物种模型 (如: aspergillus_fumigatus, neurospora_crassa)")
    species = input("物种模型 [必填] >>> ").strip()
    if not species:
        log_info("物种模型不能为空，程序退出。")
        sys.exit(1)

    # --- 步骤 4: 准备输出目录 ---
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    out_dir = os.path.abspath(f"augustus_{sample_name}_{timestamp}")
    os.makedirs(out_dir, exist_ok=True)
    out_gff = os.path.join(out_dir, f"{sample_name}.augustus.gff")
    hints_gff = os.path.join(out_dir, f"{sample_name}.hints.gff")

    print("\n" + "="*60)
    log_info("任务配置确认:")
    log_info(f"  基因组: {fasta_file}")
    log_info(f"  转录组: {bam_file if bam_file else '无 (纯算法)'}")
    log_info(f"  模  型: {species}")
    log_info(f"  输出至: {out_dir}")
    print("="*60)
    
    confirm = input("按回车键开始运行 (或输入 q 退出) >>> ").strip().lower()
    if confirm == 'q':
        sys.exit(0)

    # --- 步骤 5: 执行核心逻辑 ---
    start_t = time.time()
    use_hints = False

    # 阶段 A: 提取 Hints
    if bam_file:
        log_info(f">>> [1/2] 正在利用 bam2hints 提取转录组内含子特征...")
        cmd_hints = ["bam2hints", f"--in={bam_file}", f"--out={hints_gff}"]
        try:
            subprocess.run(cmd_hints, check=True)
            use_hints = True
            log_info(f"    Hints 提取完成！")
        except subprocess.CalledProcessError:
            log_info(f"    严重警告: bam2hints 提取失败！可能是 BAM 未排序或损坏，将退回纯算法预测。")
            use_hints = False

    # 阶段 B: 执行 AUGUSTUS
    cmd_aug = ["augustus", f"--species={species}"]
    
    if use_hints:
        cmd_aug.append(f"--hintsfile={hints_gff}")
        cmd_aug.append("--allow_hinted_splicesites=atac") 
        log_info(f">>> [2/2] 开始基于 RNA-seq 的 AUGUSTUS 共预测 (请耐心等待)...")
    else:
        log_info(f">>> [1/1] 开始纯算法 AUGUSTUS 从头预测 (请耐心等待)...")
        
    cmd_aug.append(fasta_file)
    
    try:
        # 重定向标准输出到文件，将错误输出打印到屏幕以便监控
        with open(out_gff, "w") as f_out:
            # 去除 stdout=subprocess.DEVNULL，让错误和警告直接显示在你的屏幕上
            subprocess.run(cmd_aug, stdout=f_out, check=True, text=True)
            
        elapsed = (time.time() - start_t) / 60
        print("\n" + "="*60)
        log_info(f"🎉 预测圆满完成！耗时: {elapsed:.2f} 分钟。")
        log_info(f"📁 最终 GFF 注释文件路径: {out_gff}")
        if use_hints:
            log_info(f"📁 提取的提示文件路径: {hints_gff}")
            
    except subprocess.CalledProcessError as e:
        print("\n" + "="*60)
        log_info(f"❌ 运行过程中发生错误，AUGUSTUS 异常退出。")
        log_info("建议检查上方屏幕输出的错误日志。")