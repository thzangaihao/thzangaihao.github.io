import os
import glob
import subprocess
import sys
import shlex
from datetime import datetime

# ==========================================
# 1. 配置
# ==========================================
THREADS = 128  # 充分利用你的服务器性能
OUTPUT_FASTQ = "fully_merged.fastq"
OUTPUT_FASTQ_GZ = OUTPUT_FASTQ + ".gz"


def log_info(message):
    current_time = datetime.now().strftime('%Y-%m-%d %H:%M:%S')
    print(f"[{current_time}] {message}")


def get_tool(name):
    """检查是否有 pigz，没有则回退到 gzip"""
    if subprocess.call(f"type {name}", shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE) == 0:
        return name
    return None


def parse_selection(choice, total):
    """解析 1,2,3、1-5、all 或 1-5,7,19 格式的选择。"""
    choice = choice.strip().lower()
    if not choice:
        raise ValueError("选择不能为空")
    if choice == "all":
        return list(range(total))

    selected = []
    seen = set()
    for part in choice.split(","):
        part = part.strip()
        if not part:
            raise ValueError("逗号之间缺少编号")

        if "-" in part:
            bounds = [value.strip() for value in part.split("-")]
            if len(bounds) != 2 or not all(bounds):
                raise ValueError(f"无效区间: {part}")
            try:
                start, end = map(int, bounds)
            except ValueError:
                raise ValueError(f"无效区间: {part}") from None
            if start > end:
                raise ValueError(f"区间起始编号不能大于结束编号: {part}")
            numbers = range(start, end + 1)
        else:
            try:
                numbers = [int(part)]
            except ValueError:
                raise ValueError(f"无效编号: {part}") from None

        for number in numbers:
            if not 1 <= number <= total:
                raise ValueError(f"编号 {number} 超出范围（1-{total}）")
            index = number - 1
            if index not in seen:
                seen.add(index)
                selected.append(index)

    return selected


def build_fastq_stream(files, unzip_tool):
    """构建按选择顺序输出未压缩 FASTQ 内容的 shell 数据流。"""
    commands = []
    for file_path in files:
        quoted_path = shlex.quote(file_path)
        if file_path.lower().endswith(".gz"):
            if unzip_tool == "unpigz":
                commands.append(f"unpigz -c -p {THREADS} -- {quoted_path}")
            else:
                commands.append(f"gzip -cd -- {quoted_path}")
        else:
            commands.append(f"cat -- {quoted_path}")
    return "{ " + " && ".join(commands) + "; }"

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
    output_paths = {
        os.path.abspath(OUTPUT_FASTQ),
        os.path.abspath(OUTPUT_FASTQ_GZ),
    }
    all_files = sorted(set(
        f for f in all_files
        if not os.path.basename(f).startswith(".")
        and os.path.abspath(f) not in output_paths
    ))

    if not all_files:
        print("没找到文件，退出。")
        sys.exit(0)

    for i, f in enumerate(all_files, 1):
        print(f"  [{i}] {f} ({os.path.getsize(f)/(1024**3):.2f} GB)")

    while True:
        choice = input(
            "\n请选择编号（支持 1,2,3 / 1-5 / all / 1-5,7,19）>>> "
        )
        try:
            selected_indices = parse_selection(choice, len(all_files))
            selected = [all_files[index] for index in selected_indices]
            break
        except ValueError as error:
            print(f"输入错误: {error}，请重新选择。")

    # 2. 选择最终状态
    print("\n合并后的文件处理方式：")
    print("  [1] 保持不压缩 (.fastq) - 速度最快，但占空间")
    print("  [2] 最终压缩 (.fastq.gz) - 使用多线程压缩")
    final_choice = input("你的选择 (1/2) >>> ").strip()

    # 3. 流式解压、合并并按需重新压缩，不生成临时 FASTQ 文件
    unzip_tool = "unpigz" if get_tool("unpigz") else "gzip"
    input_stream = build_fastq_stream(selected, unzip_tool)

    if final_choice == "2":
        log_info(">>> 正在流式解压、合并并压缩，不生成临时文件...")
        zip_tool = get_tool("pigz") or "gzip"
        final_file = OUTPUT_FASTQ_GZ
        compressor = (
            f"pigz -c -p {THREADS}"
            if zip_tool == "pigz"
            else "gzip -c"
        )
        command = f"{input_stream} | {compressor} > {shlex.quote(final_file)}"
    else:
        final_file = OUTPUT_FASTQ
        log_info(">>> 正在流式解压并合并，不生成临时文件...")
        command = f"{input_stream} > {shlex.quote(final_file)}"

    # pipefail 可确保解压、读取或压缩中的任一步失败时脚本立即报错。
    subprocess.run(["bash", "-o", "pipefail", "-c", command], check=True)

    print("-" * 60)
    log_info(f"全部完成！最终文件: {final_file}")
