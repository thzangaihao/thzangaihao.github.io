#!/usr/bin/env python3
"""交互式调用 Flye 组装 Oxford Nanopore 测序数据。"""

import glob
import os
import shlex
import shutil
import subprocess
import sys
import time
from datetime import datetime


DEFAULT_THREADS = 128
READ_EXTENSIONS = (".fastq.gz", ".fq.gz", ".fasta.gz", ".fa.gz", ".fastq", ".fq", ".fasta", ".fa")


def log_info(message):
    print(f"[{time.strftime('%Y-%m-%d %H:%M:%S')}] {message}", flush=True)


def ask_text(prompt, default=None, allow_empty=False):
    suffix = f"（直接回车使用 {default}）" if default is not None else ""
    try:
        value = input(f"{prompt}{suffix} >>> ").strip()
    except EOFError:
        if default is not None:
            return str(default)
        if allow_empty:
            return ""
        raise RuntimeError("当前环境不支持交互式输入。")
    if value:
        return value
    if default is not None:
        return str(default)
    if allow_empty:
        return ""
    raise ValueError("输入不能为空。")


def ask_yes_no(prompt, default=False):
    hint = "Y/n" if default else "y/N"
    while True:
        try:
            value = input(f"{prompt} [{hint}] >>> ").strip().lower()
        except EOFError:
            return default
        if not value:
            return default
        if value in {"y", "yes", "是"}:
            return True
        if value in {"n", "no", "否"}:
            return False
        print("请输入 y 或 n。")


def ask_positive_int(prompt, default):
    while True:
        value = ask_text(prompt, default)
        try:
            number = int(value)
            if number > 0:
                return number
        except ValueError:
            pass
        print("请输入大于 0 的整数。")


def select_threads(default=DEFAULT_THREADS):
    detected = os.cpu_count()
    print(f"\n检测到当前系统有 {detected if detected is not None else '未知'} 个逻辑核心。")
    threads = ask_positive_int("本次运行使用的核心数", default)
    if detected is not None and threads > detected:
        print(f"提示：选择的核心数 ({threads}) 超过检测到的逻辑核心数 ({detected})。")
    return threads


def check_dependencies():
    if shutil.which("flye") is None:
        log_info("严重错误：找不到 flye 命令，请先安装 Flye 并将其加入 PATH。")
        sys.exit(1)


def scan_read_files():
    files = []
    for extension in READ_EXTENSIONS:
        files.extend(glob.glob(f"**/*{extension}", recursive=True))
    return sorted({path for path in files if not os.path.basename(path).startswith(".")})


def select_files(files):
    print("\n请选择用于同一次组装的 Nanopore reads：")
    for number, path in enumerate(files, 1):
        print(f"  [{number}] {path}")
    print("输入编号（多个编号用逗号分隔），或输入 all 选择全部。")

    choice = ask_text("你的选择").lower()
    if choice == "all":
        return files
    try:
        indices = [int(value.strip()) - 1 for value in choice.split(",")]
    except ValueError as error:
        raise ValueError("请输入编号、逗号分隔的编号或 all。") from error
    if not indices or any(index < 0 or index >= len(files) for index in indices):
        raise ValueError("选择中包含无效编号。")
    return [files[index] for index in dict.fromkeys(indices)]


def ask_optional_number(prompt, number_type=float):
    while True:
        value = ask_text(prompt, allow_empty=True)
        if not value:
            return None
        try:
            number = number_type(value)
            if number > 0:
                return value
        except ValueError:
            pass
        print("请输入大于 0 的数值，或直接回车跳过。")


def collect_options():
    threads = select_threads()
    nano_hq = ask_yes_no("是否使用 --nano-hq 高质量 Nanopore reads 模式（通常适合 Q20+ 数据）")
    genome_size = ask_text(
        "预估基因组大小，例如 500m、2.8g；未知可留空",
        allow_empty=True,
    )
    iterations = ask_positive_int("Flye polishing 迭代次数", 1)
    meta = ask_yes_no("是否启用 --meta 宏基因组/覆盖度不均一模式", default=False)
    min_overlap = ask_optional_number("最小 overlap 长度（bp；留空使用 Flye 默认值）", int)
    asm_coverage = ask_optional_number("初始组装使用的最高覆盖度（留空不设置 --asm-coverage）")
    extra_text = ask_text(
        "其它 Flye 参数（可留空，例如：--keep-haplotypes）",
        allow_empty=True,
    )
    return {
        "threads": threads,
        "nano_hq": nano_hq,
        "genome_size": genome_size,
        "iterations": iterations,
        "meta": meta,
        "min_overlap": min_overlap,
        "asm_coverage": asm_coverage,
        "extra_args": shlex.split(extra_text),
    }


def build_command(read_files, output_dir, options):
    read_mode = "--nano-hq" if options["nano_hq"] else "--nano-raw"
    cmd = ["flye", read_mode, *[os.path.abspath(path) for path in read_files]]
    cmd.extend([
        "--out-dir", output_dir,
        "--threads", str(options["threads"]),
        "--iterations", str(options["iterations"]),
    ])
    if options["genome_size"]:
        cmd.extend(["--genome-size", options["genome_size"]])
    if options["meta"]:
        cmd.append("--meta")
    if options["min_overlap"]:
        cmd.extend(["--min-overlap", options["min_overlap"]])
    if options["asm_coverage"]:
        cmd.extend(["--asm-coverage", options["asm_coverage"]])
    cmd.extend(options["extra_args"])
    return cmd


def main():
    print("=" * 68)
    print("Nanopore 基因组组装（Flye）")
    print("=" * 68)
    check_dependencies()

    log_info("正在扫描当前目录及子目录中的 FASTQ/FASTA 文件……")
    files = scan_read_files()
    if not files:
        log_info("未找到 FASTQ/FQ/FASTA/FA 测序文件，脚本退出。")
        return

    try:
        read_files = select_files(files)
        options = collect_options()
        default_output = f"flye_assemble_{datetime.now().strftime('%Y%m%d_%H%M%S')}"
        output_name = ask_text("输出目录", default_output)
    except (ValueError, RuntimeError) as error:
        log_info(f"参数设置失败：{error}")
        sys.exit(1)

    output_dir = os.path.abspath(output_name)
    resume_requested = any(
        argument == "--resume" or argument.startswith("--resume-from")
        for argument in options["extra_args"]
    )
    if os.path.exists(output_dir) and os.listdir(output_dir) and not resume_requested:
        log_info(f"错误：输出目录已存在且非空：{output_dir}")
        log_info("请更换输出目录；如需续跑，可在其它参数中加入 --resume 并使用原目录。")
        sys.exit(1)

    cmd = build_command(read_files, output_dir, options)
    print("\n运行配置：")
    print(f"  Reads 模式：{'--nano-hq' if options['nano_hq'] else '--nano-raw'}")
    print(f"  输入文件数：{len(read_files)}")
    print(f"  核心数：{options['threads']}")
    print(f"  输出目录：{output_dir}")
    print(f"  命令：{shlex.join(cmd)}")
    if not ask_yes_no("确认开始组装", default=True):
        log_info("用户取消运行。")
        return

    try:
        log_info("开始运行 Flye……")
        subprocess.run(cmd, check=True)
    except subprocess.CalledProcessError as error:
        log_info(f"Flye 运行失败，退出码：{error.returncode}")
        sys.exit(error.returncode)
    except KeyboardInterrupt:
        log_info("运行被用户中断。")
        sys.exit(130)

    assembly = os.path.join(output_dir, "assembly.fasta")
    log_info(f"组装完成：{assembly if os.path.exists(assembly) else output_dir}")


if __name__ == "__main__":
    main()
