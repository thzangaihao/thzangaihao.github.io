import os
import subprocess
import glob
import time

# ================= 配置区域 (Configuration) =================
# 1. 集群资源配置 (根据之前的成功测试设定)
PARTITION = "tmp_7d_gpu"             # 提交队列
GPU_RES   = "gpu:nvidia_h20_96gb:1"  # 显卡型号
CPU_CORES = 32                       # CPU资源配置
TIME_LIMIT = "24:00:00"              # 给足时间
MEM_LIMIT = "250G"                   # 内存也要给足

# 2. 文件夹路径配置
JSON_ROOT_DIR = "JSON"      # 输入: 存放 JSON 的根目录
RESULT_ROOT_DIR = "Results" # 输出: 存放结果的根目录
LOG_ROOT_DIR = "Logs"       # 日志: 存放 .out 和 .err 的根目录

# 3. AlphaFold 3 镜像与数据路径
AF3_DATA_DIR = "/data/AlphaFold3"
AF3_IMAGE    = f"{AF3_DATA_DIR}/images/alphafold3.sif"
AF3_CODE     = f"{AF3_DATA_DIR}/alphafold3"
AF3_WEIGHTS  = f"{AF3_DATA_DIR}/weights"
AF3_DB       = f"{AF3_DATA_DIR}/databases"

# ================= Slurm 脚本模板 (不要轻易修改花括号内容) =================
SLURM_TEMPLATE = """#!/bin/bash
#SBATCH --job-name="{job_name}"
#SBATCH --partition={partition}
#SBATCH --gres={gpu_res}
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task={cpu_cores}
#SBATCH --mem={mem_limit}
#SBATCH --time={time_limit}
#SBATCH --output={log_out}
#SBATCH --error={log_err}

# 加载环境
module load singularity/4.3.3

echo "Job Started at $(date)"
echo "Input JSON: {input_json}"
echo "Output Dir: {output_dir}"

# 创建输出目录 (双重保险)
mkdir -p {output_dir}

# 运行 AlphaFold 3
# 注意: 这里的路径映射是固定的，宿主机路径 -> 容器内 /root/...
singularity exec --nv \\
  --bind {input_json}:/root/af_input/input.json \\
  --bind {output_dir}:/root/af_output \\
  --bind {af3_weights}:/root/models \\
  --bind {af3_db}:/root/public_database \\
  --bind {af3_code}:/root/af3code \\
  {af3_image} \\
  python3 /root/af3code/run_alphafold.py \\
  --json_path=/root/af_input/input.json \\
  --model_dir=/root/models \\
  --db_dir=/root/public_database \\
  --output_dir=/root/af_output

echo "Job Finished at $(date)"
"""

# ================= 核心逻辑 =================

def submit_task(json_path, cas_id, task_name):
    """
    组装参数并提交单个任务
    """
    # 1. 构造输出和日志的具体路径
    # 结果路径: Results/CAS号/任务名/
    job_output_dir = os.path.join(RESULT_ROOT_DIR, cas_id, task_name)
    
    # 日志路径: Logs/CAS号/
    job_log_dir = os.path.join(LOG_ROOT_DIR, cas_id)
    
    # 确保目录存在 (Slurm 不会自动创建日志文件夹，必须先创建)
    os.makedirs(job_output_dir, exist_ok=True)
    os.makedirs(job_log_dir, exist_ok=True)
    
    log_out = os.path.join(job_log_dir, f"{task_name}.out")
    log_err = os.path.join(job_log_dir, f"{task_name}.err")

    # 2. 填充 Slurm 模板
    script_content = SLURM_TEMPLATE.format(
        job_name=task_name,
        partition=PARTITION,
        gpu_res=GPU_RES,
        cpu_cores=CPU_CORES,
        mem_limit=MEM_LIMIT,
        time_limit=TIME_LIMIT,
        log_out=log_out,
        log_err=log_err,
        input_json=os.path.abspath(json_path), # 使用绝对路径更安全
        output_dir=os.path.abspath(job_output_dir),
        af3_weights=AF3_WEIGHTS,
        af3_db=AF3_DB,
        af3_code=AF3_CODE,
        af3_image=AF3_IMAGE
    )

    # 3. 调用 sbatch 提交
    try:
        # 将脚本内容通过 stdin 传给 sbatch，不产生中间文件
        process = subprocess.Popen(
            ['sbatch'], 
            stdin=subprocess.PIPE, 
            stdout=subprocess.PIPE, 
            stderr=subprocess.PIPE,
            text=True
        )
        stdout, stderr = process.communicate(input=script_content)
        
        if process.returncode == 0:
            # 成功输出: Submitted batch job 12345
            job_id = stdout.strip().split()[-1]
            return True, job_id
        else:
            return False, stderr.strip()
            
    except Exception as e:
        return False, str(e)

def main():
    # 查找所有的 JSON 文件
    # 结构预期: JSON/CAS_ID/filename.json
    search_pattern = os.path.join(JSON_ROOT_DIR, "*", "*.json")
    files = glob.glob(search_pattern)
    files.sort()
    
    total_files = len(files)
    if total_files == 0:
        print(f"错误: 在 {JSON_ROOT_DIR} 下未找到任何 .json 文件。")
        return

    print(f"准备提交 {total_files} 个任务...")
    print(f"结果目录: {RESULT_ROOT_DIR}/<CAS_ID>/<Task_Name>")
    print(f"日志目录: {LOG_ROOT_DIR}/<CAS_ID>/")
    
    # 简单的确认
    confirm = input("确认开始提交? (y/n): ").strip().lower()
    if confirm not in ['y', 'yes']:
        print("已取消。")
        return

    success_count = 0
    
    for i, json_file in enumerate(files, 1):
        # 解析路径信息
        # json_file 类似于: JSON/328968-36-1/328968-36-1_Fs000286.t01.json
        parent_dir = os.path.dirname(json_file)     # JSON/328968-36-1
        cas_id = os.path.basename(parent_dir)       # 328968-36-1
        filename = os.path.basename(json_file)      # 328968-36-1_Fs000286.t01.json
        task_name = os.path.splitext(filename)[0]   # 328968-36-1_Fs000286.t01
        
        print(f"[{i}/{total_files}] 提交: {task_name} ... ", end="", flush=True)
        
        is_success, msg = submit_task(json_file, cas_id, task_name)
        
        if is_success:
            print(f"成功 (JobID: {msg})")
            success_count += 1
        else:
            print(f"失败! 原因: {msg}")
        
        # 稍微暂停一下，避免把调度器冲垮
        time.sleep(0.1)

    print("-" * 50)
    print(f"完成! 成功提交 {success_count}/{total_files} 个任务。")
    print("使用 'squeue -u guoxj' 查看队列状态。")

if __name__ == "__main__":
    main()