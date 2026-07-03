#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import sys
import glob
import subprocess
import time
import multiprocessing
from collections import OrderedDict

'''
🧬 RSEM 转录组定量批量队列自动化流程 (Batch Queue Version)
修复：完美适配 .gz 文件，利用 bash 进程替换动态传入数据流，防崩溃。
升级：采用严格的尾部后缀匹配算法，彻底解决样本名自带 _1 或 .R2 被误认的 Bug。
'''

def get_base_dir():
    return os.path.dirname(os.path.abspath(__file__))

def check_dependency():
    deps = ['rsem-prepare-reference', 'rsem-calculate-expression', 'bowtie2', 'gunzip']
    missing = [dep for dep in deps if subprocess.call(['which', dep], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL) != 0]
    if missing:
        print(f"❌ 错误：未找到以下依赖软件：{', '.join(missing)}")
        sys.exit(1)

def interactive_select(all_files, desc, valid_exts):
    valid_files = [f for f in all_files if f.lower().endswith(valid_exts)]
    if not valid_files:
        path = input(f"👉 未扫描到{desc}，请手动输入绝对路径 (q 退出): ").strip()
        if path.lower() == 'q': sys.exit(0)
        if os.path.isfile(path): return path
        print("❌ 文件不存在。")
        sys.exit(1)

    print(f"\n📂 扫描到以下 {len(valid_files)} 个候选 {desc}:")
    for i, f in enumerate(valid_files, 1):
        print(f"  [{i}] {os.path.relpath(f, get_base_dir())}")
        
    while True:
        choice = input(f"👉 请选择 {desc} (输入编号或绝对路径，q 退出): ").strip()
        if choice.lower() == 'q': sys.exit(0)
        if choice.isdigit() and 1 <= int(choice) <= len(valid_files):
            return valid_files[int(choice)-1]
        if os.path.isfile(choice): return choice
        print("⚠️ 输入无效，请重新输入。")

def select_retention(file_type, first_label, second_label, default='4'):
    """交互式选择需要保留的两类输出文件。"""
    print(f"\n📦 请选择需要保留的 {file_type}：")
    print(f"  [1] 仅保留{first_label}")
    print(f"  [2] 仅保留{second_label}")
    print("  [3] 都不保存")
    print("  [4] 都保存")

    while True:
        choice = input(f"👉 请选择 (1/2/3/4，默认 {default}): ").strip() or default
        if choice in {'1', '2', '3', '4'}:
            return choice in {'1', '4'}, choice in {'2', '4'}
        print("⚠️ 输入无效，请输入 1、2、3 或 4。")

def remove_if_exists(path):
    """删除用户选择不保留的输出文件。"""
    if os.path.exists(path):
        os.remove(path)
        print(f"🧹 已按设置清理: {os.path.basename(path)}")

def ask_yes_no(prompt, default=True):
    """读取 y/n 交互选项。"""
    default_hint = "Y/n" if default else "y/N"
    while True:
        choice = input(f"👉 {prompt} ({default_hint}): ").strip().lower()
        if not choice:
            return default
        if choice in {'y', 'yes'}:
            return True
        if choice in {'n', 'no'}:
            return False
        print("⚠️ 输入无效，请输入 y 或 n。")

def sort_and_index_bam(input_bam, threads):
    """使用 samtools 对 BAM 按坐标排序并创建 BAI 索引。"""
    sorted_bam = input_bam[:-4] + ".sorted.bam"
    print(f"\n🔃 正在排序 BAM: {os.path.basename(input_bam)}")
    subprocess.run(
        ['samtools', 'sort', '-@', str(threads), '-o', sorted_bam, input_bam],
        check=True,
    )
    print(f"🗂️ 正在建立索引: {os.path.basename(sorted_bam)}.bai")
    subprocess.run(
        ['samtools', 'index', '-@', str(threads), sorted_bam],
        check=True,
    )
    remove_if_exists(input_bam)
    return sorted_bam

def build_directory_tree(root_dir, excluded=None):
    """生成输出目录的文本树。"""
    excluded = set(excluded or [])
    lines = [os.path.basename(root_dir) + "/"]
    for current_dir, dirnames, filenames in os.walk(root_dir):
        dirnames.sort()
        filenames = sorted(name for name in filenames if name not in excluded)
        relative = os.path.relpath(current_dir, root_dir)
        depth = 0 if relative == '.' else relative.count(os.sep) + 1
        if relative != '.':
            lines.append(f"{'    ' * (depth - 1)}└── {os.path.basename(current_dir)}/")
        for filename in filenames:
            lines.append(f"{'    ' * depth}├── {filename}")
    return "\n".join(lines)

def write_readme(out_dir, settings):
    """记录本次任务配置及实际输出文件层次结构。"""
    readme_path = os.path.join(out_dir, 'readme.txt')
    tree = build_directory_tree(out_dir, excluded={'readme.txt'})
    with open(readme_path, 'w', encoding='utf-8') as handle:
        handle.write("RSEM 定量分析任务说明\n")
        handle.write("=" * 60 + "\n")
        handle.write(f"生成时间：{time.strftime('%Y-%m-%d %H:%M:%S')}\n")
        for label, value in settings:
            handle.write(f"{label}：{value}\n")
        handle.write("\n输出文件层次结构\n")
        handle.write("-" * 60 + "\n")
        handle.write(tree + "\n")
    print(f"📝 已生成任务说明: {readme_path}")

def build_sample_queue(fastq_dir, is_paired):
    """核心队列逻辑：扫描文件夹并严格配对样本"""
    all_fastqs = glob.glob(os.path.join(fastq_dir, "**", "*.fq"), recursive=True) + \
                 glob.glob(os.path.join(fastq_dir, "**", "*.fastq"), recursive=True) + \
                 glob.glob(os.path.join(fastq_dir, "**", "*.fq.gz"), recursive=True) + \
                 glob.glob(os.path.join(fastq_dir, "**", "*.fastq.gz"), recursive=True)
                 
    sample_queue = OrderedDict()
    
    # 扩展的后缀匹配列表，优先匹配最长的后缀
    valid_extensions = ['.raw.fastq.gz', '.raw.fq.gz', '.fastq.gz', '.fq.gz', '.raw.fastq', '.raw.fq', '.fastq', '.fq']
    
    if not is_paired:
        for fq in sorted(all_fastqs):
            basename = os.path.basename(fq)
            for ext in valid_extensions:
                if basename.endswith(ext):
                    sample_name = basename[:-len(ext)]
                    sample_queue[sample_name] = [fq]
                    break
        return sample_queue

    # 双端测序严格配对逻辑 (修复名称中含 _1 的 bug)
    for fq in sorted(all_fastqs):
        basename = os.path.basename(fq)
        
        # 1. 匹配并去除文件扩展名
        matched_ext = ""
        for ext in valid_extensions:
            if basename.endswith(ext):
                matched_ext = ext
                break
        
        if not matched_ext:
            continue
            
        name_without_ext = basename[:-len(matched_ext)]
        
        # 2. 严格检查 Read 1 标识符（必须在去除后缀后的最末尾）
        r1_suffix, r2_suffix = "", ""
        if name_without_ext.endswith('.R1'):
            r1_suffix, r2_suffix = '.R1', '.R2'
        elif name_without_ext.endswith('_R1'):
            r1_suffix, r2_suffix = '_R1', '_R2'
        elif name_without_ext.endswith('_1'):
            r1_suffix, r2_suffix = '_1', '_2'
        else:
            # 如果不是 R1 文件（例如是 R2 文件，或者未识别格式），直接跳过
            # 我们只需要找到 R1，然后去反推 R2 的位置即可，这样可以避免 R2 文件报错
            continue
            
        # 3. 提取干净的样本名并寻找 R2
        sample_name = name_without_ext[:-len(r1_suffix)]
        expected_r2_basename = sample_name + r2_suffix + matched_ext
        expected_r2_path = os.path.join(os.path.dirname(fq), expected_r2_basename)
        
        if os.path.exists(expected_r2_path):
            sample_queue[sample_name] = [fq, expected_r2_path]
        else:
            print(f"⚠️ 警告: 找不到 {basename} 对应的 Read 2 文件 ({expected_r2_basename})，已跳过该样本。")
            
    return sample_queue

def run_cmd(cmd, desc):
    """执行系统命令，强制使用 bash 以支持高级语法"""
    print(f"\n▶️ [{time.strftime('%H:%M:%S')}] 正在执行: {desc}")
    print(f"   💻 命令: {cmd}") 
    try:
        subprocess.run(cmd, shell=True, check=True, executable='/bin/bash')
        print(f"✅ {desc} 完成！")
    except subprocess.CalledProcessError:
        print(f"\n❌ 错误：{desc} 失败。请检查上方的报错信息。")
        sys.exit(1)

def main():
    print("=" * 60)
    print(" 🧬 RSEM 批量队列自动化流程 (严格配对版)")
    print("=" * 60)

    check_dependency()
    base_dir = get_base_dir()
    
    print("🔍 正在扫描当前目录文件...")
    all_files = glob.glob(os.path.join(base_dir, "**", "*.*"), recursive=True)

    # 1. 配置参考基因组
    print("\n" + "-"*40)
    print(" 🛠️ 步骤 1：配置参考基因组与注释")
    genome_fasta = interactive_select(all_files, "参考基因组文件", ('.fa', '.fasta', '.fna'))
    annotation_gtf = interactive_select(all_files, "基因组注释文件", ('.gtf', '.gff', '.gff3'))

    ref_basename = os.path.basename(genome_fasta).rsplit('.', 1)[0]
    ref_dir = os.path.join(base_dir, f"{ref_basename}_rsem_idx")
    ref_prefix = os.path.join(ref_dir, "ref")

    # 2. 配置测序数据目录与队列
    print("\n" + "-"*40)
    print(" 🛠️ 步骤 2：配置测序数据与队列")
    fastq_dir = input(f"👉 请输入存放测序数据的文件夹路径 (直接回车扫描当前目录): ").strip()
    if not fastq_dir: fastq_dir = base_dir
    
    is_paired = input("👉 数据是双端测序 (Paired-end) 吗？(y/n, 默认 y): ").strip().lower() != 'n'
    
    # 构建任务队列
    sample_queue = build_sample_queue(fastq_dir, is_paired)
    
    if not sample_queue:
        print("❌ 在指定目录中未找到任何匹配的测序文件。")
        sys.exit(1)

    print(f"\n📋 成功构建任务队列，共扫描到 \033[92m{len(sample_queue)}\033[0m 个样本：")
    for i, (sample, files) in enumerate(sample_queue.items(), 1):
        print(f"  [{i}] {sample} -> {len(files)} 个文件")

    # 3. 性能设置与输出目录
    print("\n" + "-"*40)
    threads = max(1, int(multiprocessing.cpu_count() * 0.8))
    user_threads = input(f"👉 请输入 RSEM 调用的线程数 (直接回车使用 {threads}): ").strip()
    if user_threads.isdigit(): threads = int(user_threads)

    keep_genome_bam, keep_transcriptome_bam = select_retention(
        "BAM 文件", "比对到基因组上的 BAM 文件", "比对到转录本上的 BAM 文件"
    )
    keep_gene_results, keep_isoform_results = select_retention(
        "results 文件", "基因层面的 genes.results", "转录本层面的 isoforms.results"
    )
    if not keep_gene_results and not keep_isoform_results:
        print("⚠️ 提示：两类 results 文件都不保留，将不会留下 RSEM 定量结果。")

    sort_index_bam = False
    if keep_genome_bam or keep_transcriptome_bam:
        sort_index_bam = ask_yes_no("是否对保留的 BAM 文件按坐标排序并自动建立 BAI 索引？", default=True)
        if sort_index_bam and subprocess.call(
            ['which', 'samtools'], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL
        ) != 0:
            print("❌ 错误：BAM 排序与索引需要 samtools，但当前环境未找到该命令。")
            sys.exit(1)

    out_dir = os.path.join(base_dir, "RSEM_Results")
    if not os.path.exists(out_dir): os.makedirs(out_dir)

    print("\n" + "="*60)
    print("🚀 队列执行即将开始")
    print("="*60)

    # 4. 构建参考库
    if not os.path.exists(ref_dir):
        os.makedirs(ref_dir)
        gtf_flag = "--gff3" if annotation_gtf.lower().endswith(('.gff', '.gff3')) else "--gtf"
        build_cmd = f"rsem-prepare-reference --num-threads {threads} {gtf_flag} '{annotation_gtf}' --bowtie2 '{genome_fasta}' '{ref_prefix}'"
        run_cmd(build_cmd, "构建 RSEM 参考基因组库")
    else:
        print(f"\n✅ 检测到参考基因组索引已存在，跳过建库步骤。")

    # 5. 执行任务队列
    for idx, (sample_name, files) in enumerate(sample_queue.items(), 1):
        sample_out_dir = os.path.join(out_dir, sample_name)
        os.makedirs(sample_out_dir, exist_ok=True)
        output_prefix = os.path.join(sample_out_dir, sample_name)
        output_files = {
            'genome_bam': f"{output_prefix}.genome.bam",
            'transcriptome_bam': f"{output_prefix}.transcript.bam",
            'gene_results': f"{output_prefix}.genes.results",
            'isoform_results': f"{output_prefix}.isoforms.results",
        }
        expected_outputs = []
        for keep, path in (
            (keep_genome_bam, output_files['genome_bam']),
            (keep_transcriptome_bam, output_files['transcriptome_bam']),
        ):
            if keep:
                bam_path = path[:-4] + '.sorted.bam' if sort_index_bam else path
                expected_outputs.append(bam_path)
                if sort_index_bam:
                    expected_outputs.append(bam_path + '.bai')
        if keep_gene_results:
            expected_outputs.append(output_files['gene_results'])
        if keep_isoform_results:
            expected_outputs.append(output_files['isoform_results'])
        
        print("\n" + "-"*40)
        print(f"⏳ 正在处理队列任务 [{idx}/{len(sample_queue)}]: \033[93m{sample_name}\033[0m")
        
        if expected_outputs and all(os.path.exists(path) for path in expected_outputs):
            print(f"⏭️ 检测到 {sample_name} 所有所选输出均已存在，跳过该样本。")
            continue

        rsem_params = f"--num-threads {threads} --bowtie2"
        if is_paired: rsem_params += " --paired-end"
        # RSEM 默认输出 transcript.bam；genome.bam 需要显式开启。
        if not keep_genome_bam and not keep_transcriptome_bam:
            rsem_params += " --no-bam-output"
        elif keep_genome_bam:
            rsem_params += " --output-genome-bam"
        
        if files[0].endswith('.gz'):
            fq_inputs = " ".join([f"<(gunzip -c '{f}')" for f in files])
        else:
            fq_inputs = " ".join([f"'{f}'" for f in files])
            
        quant_cmd = f"rsem-calculate-expression {rsem_params} {fq_inputs} '{ref_prefix}' '{output_prefix}'"
        
        run_cmd(quant_cmd, f"定量计算: {sample_name}")

        if sort_index_bam:
            if keep_genome_bam:
                sort_and_index_bam(output_files['genome_bam'], threads)
            if keep_transcriptome_bam:
                sort_and_index_bam(output_files['transcriptome_bam'], threads)

        # 使用 --output-genome-bam 时 RSEM 仍会默认生成 transcript.bam。
        # 若用户选择仅保留基因组 BAM，则在任务成功后清理转录本 BAM。
        if not keep_transcriptome_bam:
            remove_if_exists(output_files['transcriptome_bam'])

        # results 文件由 RSEM 默认产生；按用户选择在成功完成后清理不需要的文件。
        if not keep_gene_results:
            remove_if_exists(output_files['gene_results'])
        if not keep_isoform_results:
            remove_if_exists(output_files['isoform_results'])

    write_readme(out_dir, [
        ('参考基因组', os.path.abspath(genome_fasta)),
        ('注释文件', os.path.abspath(annotation_gtf)),
        ('测序数据目录', os.path.abspath(fastq_dir)),
        ('样本数量', len(sample_queue)),
        ('测序类型', '双端' if is_paired else '单端'),
        ('线程数', threads),
        ('保留基因组 BAM', '是' if keep_genome_bam else '否'),
        ('保留转录本 BAM', '是' if keep_transcriptome_bam else '否'),
        ('BAM 排序并建立索引', '是' if sort_index_bam else '否'),
        ('保留 genes.results', '是' if keep_gene_results else '否'),
        ('保留 isoforms.results', '是' if keep_isoform_results else '否'),
    ])

    print("\n" + "="*60)
    print("🎉 所有队列任务处理完毕！")
    print(f"📂 结果文件保存在: {out_dir}/")
    print("="*60 + "\n")

if __name__ == "__main__":
    main()
