import os
import glob
import time
import subprocess
import sys
from datetime import datetime
from concurrent.futures import ThreadPoolExecutor, as_completed

# ==========================================
# 1. 基础配置与辅助函数
# ==========================================
def log_info(message):
    """打印带时间戳的日志信息"""
    current_time = time.strftime('%Y-%m-%d %H:%M:%S')
    print(f"[{current_time}] {message}")

def check_dependencies():
    """检查 augustus 是否在环境变量中"""
    if subprocess.call("type augustus", shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE) != 0:
        log_info("严重错误：系统未找到命令 'augustus'。请先激活相应的 conda 环境。")
        sys.exit(1)

# ==========================================
# 2. 核心预测任务函数 (将被扔进线程池)
# ==========================================
def run_augustus(fasta_file, species, out_dir):
    """运行单样本的 AUGUSTUS 预测"""
    sample_name = os.path.splitext(os.path.basename(fasta_file))[0]
    
    # 清理常见的后缀名，让输出文件名更清爽
    for ext in [".p_ctg", ".ctg", ".hifi"]:
        sample_name = sample_name.replace(ext, "")
        
    out_gff = os.path.join(out_dir, f"{sample_name}.augustus.gff")
    
    # 构造 AUGUSTUS 命令 (标准算法预测，不借助任何外部提示)
    # augustus 默认将结果输出到标准输出 (STDOUT)，我们需要重定向到文件
    cmd = ["augustus", f"--species={species}", fasta_file]
    
    log_info(f"  -> [启动] 样本 {sample_name} 的预测任务...")
    start_t = time.time()
    
    try:
        with open(out_gff, "w") as f_out:
            # 捕获错误信息并隐藏标准输出流的滚动
            subprocess.run(cmd, stdout=f_out, stderr=subprocess.PIPE, check=True, text=True)
            
        elapsed = (time.time() - start_t) / 60
        log_info(f"  -> [完成] 样本 {sample_name} 预测完毕，耗时: {elapsed:.2f} 分钟。")
        return True, sample_name
        
    except subprocess.CalledProcessError as e:
        log_info(f"  -> [失败] 样本 {sample_name} 运行出错！")
        log_info(f"     错误详情: {e.stderr.strip()}")
        return False, sample_name
    except Exception as e:
        log_info(f"  -> [异常] 样本 {sample_name}: {e}")
        return False, sample_name

# ==========================================
# 3. 主控制流
# ==========================================
if __name__ == "__main__":
    print("="*60)
    print("      AUGUSTUS 纯算法基因预测 - 批量并发流水线")
    print("="*60)

    check_dependencies()

    # 1. 扫描所有的 FASTA 文件
    log_info("正在扫描当前目录及子目录下的序列文件...")
    extensions = ["**/*.fasta", "**/*.fa", "**/*.fna"]
    all_files = []
    for ext in extensions:
        all_files.extend(glob.glob(ext, recursive=True))
    
    # 过滤掉隐藏文件并去重排序
    all_files = sorted(list(set([f for f in all_files if not os.path.basename(f).startswith(".")])))

    if not all_files:
        log_info("未找到任何 FASTA 文件，脚本退出。")
        sys.exit(0)

    # 2. 交互式选择菜单
    print("\n发现以下待处理文件：")
    for i, file_path in enumerate(all_files, 1):
        print(f"  [{i}] {file_path}")
    
    print("\n请选择要预测的文件编号 (例如: 1,3,4)，或者输入 'all' 处理全部：")
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

    # 3. 参数设定：指定模型
    print("\n" + "-"*40)
    print("AUGUSTUS 需要指定物种模型 (Species Model)。")
    print("提示: 可用 'augustus --species=help' 查看所有内置模型。")
    print("常用的真菌模型例如: aspergillus_fumigatus, saccharomyces_cerevisiae_S288C, neurospora_crassa 等")
    
    species = input("请输入物种模型名称 [必填] >>> ").strip()
    if not species:
        log_info("物种模型不能为空，脚本退出。")
        sys.exit(1)

    # 4. 参数设定：指定并发数
    print("\n" + "-"*40)
    print(f"已选择 {len(selected_files)} 个样本。")
    max_workers_input = input("请输入最大并发样本数 (默认: 4) >>> ").strip()
    try:
        max_workers = int(max_workers_input) if max_workers_input else 4
    except ValueError:
        max_workers = 4

    # 5. 创建输出目录
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    work_dir = os.path.abspath(f"augustus_results_{timestamp}")
    os.makedirs(work_dir, exist_ok=True)
    
    print("\n" + "="*60)
    log_info(f"任务总数: {len(selected_files)}")
    log_info(f"并发线程: {max_workers}")
    log_info(f"物种模型: {species}")
    log_info(f"输出目录: {work_dir}")
    print("="*60 + "\n")

    # 6. 使用线程池执行并发任务
    start_total_time = time.time()
    success_count = 0
    fail_count = 0
    
    # ThreadPoolExecutor 负责任务的并发调度
    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        # 提交所有任务到线程池
        futures = {
            executor.submit(run_augustus, fasta_path, species, work_dir): fasta_path 
            for fasta_path in selected_files
        }
        
        # 实时监控完成情况
        for future in as_completed(futures):
            success, s_name = future.result()
            if success:
                success_count += 1
            else:
                fail_count += 1

    # 7. 汇总报告
    total_elapsed = (time.time() - start_total_time) / 60
    print("\n" + "="*60)
    log_info(f"全部预测队列执行完毕！总耗时: {total_elapsed:.2f} 分钟。")
    log_info(f"成功: {success_count} 个，失败: {fail_count} 个。")
    log_info(f"所有 GFF 结果文件保存在: {work_dir}")