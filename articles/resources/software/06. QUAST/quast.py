import os
import glob
import time
import subprocess
import sys
from datetime import datetime

# ==========================================
# 1. 基础配置
# ==========================================
def log_info(message):
    current_time = time.strftime('%Y-%m-%d %H:%M:%S')
    print(f"[{current_time}] {message}")

def check_quast():
    """检查 quast.py 是否可用"""
    # 有些环境下是 quast，有些是 quast.py
    for cmd in ["quast.py", "quast"]:
        if subprocess.call(f"type {cmd}", shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE) == 0:
            return cmd
    log_info("严重错误：未找到 quast 或 quast.py，请确认已安装（conda install -c bioconda quast）。")
    sys.exit(1)

# ==========================================
# 2. 核心功能
# ==========================================
def run_quast_task(fasta_list, output_dir, label="evaluation"):
    """运行 QUAST 命令"""
    quast_cmd = check_quast()
    log_info(f"正在对 {len(fasta_list)} 个文件进行评估...")
    
    # 构建命令
    cmd = [quast_cmd] + fasta_list + ["-o", output_dir, "--threads", "16"]
    
    try:
        # QUAST 运行时输出较多，这里直接打印到屏幕以便监控进度
        subprocess.run(cmd, check=True)
        log_info(f"评估完成！报告路径: {os.path.abspath(output_dir)}/report.html")
    except subprocess.CalledProcessError as e:
        log_info(f"QUAST 运行出错: {e}")

# ==========================================
# 3. 主控制流
# ==========================================
if __name__ == "__main__":
    print("="*60)
    print("         QUAST 交互式基因组评估工具")
    print("="*60)

    # 1. 扫描 FASTA 文件
    log_info("正在扫描当前目录及子目录下的 .fasta 文件...")
    # 优先扫描你之前脚本生成的 p_ctg.fasta
    all_fastas = glob.glob("**/*.fasta", recursive=True) + glob.glob("**/*.fa", recursive=True)
    # 去重
    all_fastas = sorted(list(set(all_fastas)))

    if not all_fastas:
        log_info("未找到任何 .fasta 或 .fa 文件，请检查路径。")
        sys.exit(0)

    # 2. 交互式选择
    print("\n发现以下待评估文件：")
    for i, file_path in enumerate(all_fastas, 1):
        print(f"  [{i}] {file_path}")
    
    print("\n请输入编号选择文件 (例如: 1,2,5)，或输入 'all' 全选：")
    choice = input("你的选择 >>> ").strip().lower()

    selected_files = []
    if choice == 'all':
        selected_files = all_fastas
    else:
        try:
            indices = [int(x.strip()) - 1 for x in choice.split(',')]
            selected_files = [all_fastas[i] for i in indices if 0 <= i < len(all_fastas)]
        except:
            log_info("输入格式错误。")
            sys.exit(1)

    if not selected_files:
        log_info("未选择有效文件。")
        sys.exit(0)

    # 3. 选择运行模式
    print("\n请选择评估模式：")
    print("  [1] 合并评估 (所有选中的样本生成一份对比报告，推荐)")
    print("  [2] 独立评估 (每个样本生成各自的文件夹)")
    mode = input("模式选择 (1/2) >>> ").strip()

    # 4. 执行
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    base_output = f"quast_eval_{timestamp}"
    os.makedirs(base_output, exist_ok=True)

    if mode == "1":
        # 模式1：合并评估
        run_quast_task(selected_files, os.path.join(base_output, "combined_report"))
    else:
        # 模式2：独立评估
        for f in selected_files:
            sample_name = os.path.splitext(os.path.basename(f))[0]
            log_info(f"正在处理样本: {sample_name}")
            run_quast_task([f], os.path.join(base_output, sample_name))

    print("\n" + "="*60)
    log_info(f"所有任务已完成！结果存放在: {base_output}")