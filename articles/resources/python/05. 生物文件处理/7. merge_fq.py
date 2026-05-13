import os
import glob
import subprocess
import sys
from datetime import datetime
import shutil

# ==========================================
# 1. 配置
# ==========================================
THREADS = 128  # 充分利用你的服务器性能

def log_info(message):
    current_time = datetime.now().strftime('%Y-%m-%d %H:%M:%S')
    print(f"[{current_time}] {message}")

def get_tool(name):
    """检查是否有 pigz，没有则回退到 gzip"""
    if subprocess.call(f"type {name}", shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE) == 0:
        return name
    return None

# ==========================================
# 2. 交互逻辑
# ==========================================
if __name__ == "__main__":
    print("="*60)
    print("      极致速度：先解压-再合并-后压缩工具")
    print("="*60)

    # 1. 扫描文件
    extensions = ["**/*.fastq", "**/*.fastq.gz", "**/*.fq", "**/*.fq.gz"]
    all_files = []
    for ext in extensions:
        all_files.extend(glob.glob(ext, recursive=True))
    all_files = sorted(list(set([f for f in all_files if not os.path.basename(f).startswith(".")])))

    if not all_files:
        print("没找到文件，退出。")
        sys.exit(0)

    for i, f in enumerate(all_files, 1):
        print(f"  [{i}] {f} ({os.path.getsize(f)/(1024**3):.2f} GB)")

    try:
        choice = input("\n请选择编号 (例如 1,2,3) >>> ").strip()
        selected = [all_files[int(x)-1] for x in choice.split(',')]
    except:
        sys.exit(1)

    # 2. 选择最终状态
    print("\n合并后的文件处理方式：")
    print("  [1] 保持不压缩 (.fastq) - 速度最快，但占空间")
    print("  [2] 最终压缩 (.fastq.gz) - 使用多线程压缩")
    final_choice = input("你的选择 (1/2) >>> ").strip()

    # 3. 创建临时工作目录（解压用）
    work_dir = f"temp_merge_{datetime.now().strftime('%H%M%S')}"
    os.makedirs(work_dir, exist_ok=True)
    
    uncompressed_files = []
    
    # 第一阶段：多线程并行解压
    log_info(">>> 阶段 1: 正在将所有压缩文件解压到临时目录...")
    unzip_tool = get_tool("unpigz") or get_tool("gunzip")
    
    for f in selected:
        base_name = os.path.basename(f)
        target_path = os.path.join(work_dir, base_name.replace(".gz", ""))
        
        if f.endswith(".gz"):
            # 使用 -c 将解压流定向到新文件，保持原文件不动
            cmd = f"{unzip_tool} -c -p {THREADS} {f} > {target_path}" if "unpigz" in unzip_tool else f"gunzip -c {f} > {target_path}"
            log_info(f"正在解压: {base_name} ...")
            subprocess.run(cmd, shell=True, check=True)
            uncompressed_files.append(target_path)
        else:
            # 如果本身没压缩，直接创建软链接到临时目录
            os.symlink(os.path.abspath(f), target_path)
            uncompressed_files.append(target_path)

    # 第二阶段：合并
    out_name = "fully_merged.fastq"
    log_info(f">>> 阶段 2: 正在执行快速合并 (cat) ...")
    files_str = " ".join([f"'{f}'" for f in uncompressed_files])
    subprocess.run(f"cat {files_str} > {out_name}", shell=True, check=True)

    # 第三阶段：处理最终结果
    if final_choice == "2":
        log_info(">>> 阶段 3: 正在进行最终多线程压缩...")
        zip_tool = get_tool("pigz") or "gzip"
        zip_cmd = f"{zip_tool} -p {THREADS} {out_name}" if zip_tool == "pigz" else f"gzip {out_name}"
        subprocess.run(zip_cmd, shell=True, check=True)
        final_file = out_name + ".gz"
    else:
        final_file = out_name

    # 清理
    log_info(">>> 正在清理临时文件...")
    shutil.rmtree(work_dir)

    print("-" * 60)
    log_info(f"全部完成！最终文件: {final_file}")