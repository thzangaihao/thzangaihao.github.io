#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import sys
import glob
import subprocess
import time
import multiprocessing
from datetime import datetime

'''
👑 全基因组共线性分析终极自动化管线 (Master Synteny Pipeline)
整合模块：
1. 蛋白提取与最长转录本过滤 (gffread_pro + filter_longest)
2. GFF 转 BED 并强制注入前缀 (gff_to_bed 终极版)
3. Diamond 全对全比对 (智能分配算力)
4. JCVI 共线性区块识别 (软链接防伪装术)
'''

def get_base_dir():
    return os.path.dirname(os.path.abspath(__file__))

def check_dependencies():
    """检查核心依赖库"""
    deps = ['gffread', 'diamond']
    for dep in deps:
        if subprocess.call(['which', dep], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL) != 0:
            print(f"❌ 致命错误：未找到系统组件 '{dep}'！请先安装。")
            sys.exit(1)
    try:
        subprocess.run(["python", "-m", "jcvi"], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    except Exception:
        print("❌ 致命错误：未检测到 jcvi！请确认已激活对应的 conda 环境。")
        sys.exit(1)

def interactive_select(files, desc):
    """交互式文件选择器"""
    if not files:
        print(f"❌ 严重错误：未在当前目录及其子目录下找到任何 {desc}！")
        sys.exit(1)
    
    print(f"\n📂 扫描到以下 {len(files)} 个 {desc}:")
    for i, f in enumerate(files, 1):
        print(f"  [{i}] {os.path.relpath(f, get_base_dir())}")
        
    while True:
        choice = input(f"👉 请选择 {desc} 编号 (q退出): ").strip().lower()
        if choice == 'q': sys.exit(0)
        try:
            return os.path.abspath(files[int(choice)-1])
        except:
            print("⚠️ 输入无效，请重新输入。")

def process_proteins(ref_fasta, gff_file, prefix, out_dir):
    """提取蛋白 -> 清洗 -> 筛选最长转录本 -> 注入前缀"""
    print(f"\n⚙️ [{prefix}] 正在提取全基因组蛋白序列...")
    temp_faa = os.path.join(out_dir, f"temp_{prefix}.faa")
    final_faa = os.path.join(out_dir, f"{prefix}.faa")
    
    # 1. 运行 gffread
    cmd = f"gffread '{gff_file}' -g '{ref_fasta}' -y '{temp_faa}'"
    subprocess.run(cmd, shell=True, check=True)
    
    # 2. 读取并筛选最长转录本
    print(f"⚙️ [{prefix}] 正在筛选最长转录本并注入物种前缀 '{prefix}_'...")
    sequences, current_header, current_seq = [], None, []
    with open(temp_faa, 'r') as f:
        for line in f:
            line = line.strip()
            if not line: continue
            if line.startswith('>'):
                if current_header: sequences.append((current_header, ''.join(current_seq)))
                current_header = line
                current_seq = []
            else:
                current_seq.append(line)
        if current_header: sequences.append((current_header, ''.join(current_seq)))
        
    longest_transcripts = {}
    for header, seq in sequences:
        # 提取转录本纯 ID
        raw_id = header.split()[0][1:] 
        gene_id = raw_id.split('.')[0] if '.' in raw_id else raw_id # 简单粗暴找 Gene 归属
        
        # 处理 gffread 的 gene= 标签
        for part in header.split():
            if part.startswith('gene='): gene_id = part.split('=')[1]
            
        clean_seq = "".join(filter(str.isalpha, seq))
        if gene_id not in longest_transcripts or len(clean_seq) > len(longest_transcripts[gene_id][1]):
            # 🌟 核心：在此处直接将前缀注入 FASTA Header
            new_header = f">{prefix}_{raw_id}"
            longest_transcripts[gene_id] = (new_header, clean_seq)
            
    # 3. 写入最终蛋白文件
    with open(final_faa, 'w') as out:
        for _, (header, seq) in longest_transcripts.items():
            out.write(f"{header}\n")
            for i in range(0, len(seq), 80): out.write(f"{seq[i:i+80]}\n")
            
    os.remove(temp_faa)
    print(f"✅ [{prefix}] 蛋白提取完成！代表转录本数量: {len(longest_transcripts)}")
    return final_faa

def process_bed(gff_file, prefix, out_dir):
    """GFF 转 BED 并强制注入前缀"""
    print(f"⚙️ [{prefix}] 正在生成简化基因位置文件 (BED)...")
    bed_path = os.path.join(out_dir, f"{prefix}.bed")
    
    cmd = f"python -m jcvi.formats.gff bed '{gff_file}' -o '{bed_path}' --type=mRNA"
    try:
        subprocess.run(cmd, shell=True, check=True)
    except subprocess.CalledProcessError:
        print(f"⚠️ [{prefix}] mRNA 提取失败，尝试使用 transcript 类型...")
        subprocess.run(f"python -m jcvi.formats.gff bed '{gff_file}' -o '{bed_path}' --type=transcript", shell=True, check=True)
        
    # 🌟 核心：强制注入前缀
    lines = []
    with open(bed_path, 'r') as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) >= 4 and not parts[3].startswith(f"{prefix}_"):
                parts[3] = f"{prefix}_{parts[3]}"
            lines.append('\t'.join(parts))
            
    with open(bed_path, 'w') as f:
        f.write('\n'.join(lines) + '\n')
    print(f"✅ [{prefix}] BED 文件生成完毕，前缀注入成功！")
    return bed_path

def run_diamond_and_synteny(ref_prefix, ref_faa, qry_prefix, qry_faa, out_dir, threads):
    """运行 Diamond 双向比对并直接触发 JCVI"""
    print("\n" + "="*50)
    print(f" 🚀 开始分析核心对: [{ref_prefix}] ⚔️ [{qry_prefix}]")
    print("="*50)
    
    # 1. 建库
    print(f"▶️ 正在构建 Diamond 数据库...")
    ref_db = os.path.join(out_dir, ref_prefix)
    qry_db = os.path.join(out_dir, qry_prefix)
    subprocess.run(f"diamond makedb --threads {threads} --in '{ref_faa}' -d '{ref_db}'", shell=True, check=True, stdout=subprocess.DEVNULL)
    subprocess.run(f"diamond makedb --threads {threads} --in '{qry_faa}' -d '{qry_db}'", shell=True, check=True, stdout=subprocess.DEVNULL)
    
    # 2. 双向比对
    print(f"▶️ 正在运行全对全双向比对 (Diamond BLASTP)...")
    blast_r_vs_q = os.path.join(out_dir, f"temp_{ref_prefix}_{qry_prefix}.blast")
    blast_q_vs_r = os.path.join(out_dir, f"temp_{qry_prefix}_{ref_prefix}.blast")
    params = f"--threads {threads} --evalue 1e-5 --outfmt 6 --max-target-seqs 5"
    
    subprocess.run(f"diamond blastp -q '{ref_faa}' -d '{qry_db}' -o '{blast_r_vs_q}' {params}", shell=True, check=True, stderr=subprocess.DEVNULL)
    subprocess.run(f"diamond blastp -q '{qry_faa}' -d '{ref_db}' -o '{blast_q_vs_r}' {params}", shell=True, check=True, stderr=subprocess.DEVNULL)
    
    # 3. 合并为 .last (完美欺骗 JCVI)
    print(f"▶️ 正在合并比对结果并封装...")
    final_last = os.path.join(out_dir, f"{ref_prefix}.{qry_prefix}.last")
    with open(final_last, 'w') as outfile:
        for fname in [blast_r_vs_q, blast_q_vs_r]:
            with open(fname) as infile: outfile.write(infile.read())
            os.remove(fname)
            
    # 4. 召唤 JCVI
    print(f"\n🧠 正在调用 JCVI MCscan 核心算法寻找 Anchors...")
    cmd_jcvi = f"python -m jcvi.compara.catalog ortholog {ref_prefix} {qry_prefix} --no_strip_names"
    try:
        subprocess.run(cmd_jcvi, shell=True, check=True, cwd=out_dir)
        print(f"🎉 成功！[{ref_prefix}] 与 [{qry_prefix}] 的共线性区块识别完毕！")
    except subprocess.CalledProcessError:
        print(f"❌ [{ref_prefix} vs {qry_prefix}] JCVI 运行失败，请检查报错日志。")

def main():
    print("=" * 60)
    print(" 👑 全基因组共线性分析终极自动化管线")
    print("   (One-Click Synteny Master Pipeline)")
    print("=" * 60)

    check_dependencies()
    base_dir = get_base_dir()

    # 1. 扫描可用文件
    all_fasta = glob.glob(os.path.join(base_dir, "**", "*.fasta"), recursive=True) + \
                glob.glob(os.path.join(base_dir, "**", "*.fna"), recursive=True)
    all_gff = glob.glob(os.path.join(base_dir, "**", "*.gff*"), recursive=True)

    # 2. 注册参考物种 (Reference)
    print("\n" + "🌟 步骤一：配置【参考物种】(Reference)")
    ref_fasta = interactive_select(all_fasta, "参考物种基因组序列 (.fasta/.fna)")
    ref_gff = interactive_select(all_gff, "参考物种注释文件 (.gff/.gff3)")
    ref_prefix = input(f"👉 请输入参考物种前缀代号 (例如 {os.path.basename(ref_fasta).split('.')[0]}): ").strip()

    # 3. 注册查询物种 (Queries)
    print("\n" + "🌟 步骤二：配置【待比较物种】(Query)")
    queries = []
    while True:
        print(f"\n--- 正在添加第 {len(queries)+1} 个查询物种 ---")
        q_fasta = interactive_select([f for f in all_fasta if f != ref_fasta], "待比较物种基因组序列")
        q_gff = interactive_select([f for f in all_gff if f != ref_gff], "待比较物种注释文件")
        q_prefix = input(f"👉 请输入该物种前缀代号: ").strip()
        queries.append((q_fasta, q_gff, q_prefix))
        
        cont = input("\n[?] 是否继续添加下一个待比较物种？(y/n): ").strip().lower()
        if cont != 'y': break

    # 4. 配置算力
    max_cores = multiprocessing.cpu_count()
    threads = max(1, int(max_cores * 0.8))
    user_threads = input(f"\n⚡ 请设置并发线程数 (当前系统 {max_cores} 核，默认使用 {threads}): ").strip()
    if user_threads.isdigit(): threads = int(user_threads)

    # 5. 创建主控工作区
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    work_dir = os.path.join(base_dir, f"Master_Synteny_Workspace_{timestamp}")
    os.makedirs(work_dir, exist_ok=True)
    print(f"\n📁 成功创建独立主控工作区: {os.path.relpath(work_dir, base_dir)}")

    # 6. 数据预处理池
    start_time = time.time()
    print("\n" + "="*50)
    print(" 🛠️ 阶段 A：数据标准化预处理 (蛋白提取与 BED 格式化)")
    print("="*50)
    
    # 处理参考种
    ref_faa_path = process_proteins(ref_fasta, ref_gff, ref_prefix, work_dir)
    process_bed(ref_gff, ref_prefix, work_dir)
    
    # 处理查询种
    qry_faa_paths = {}
    for q_fasta, q_gff, q_prefix in queries:
        qry_faa_paths[q_prefix] = process_proteins(q_fasta, q_gff, q_prefix, work_dir)
        process_bed(q_gff, q_prefix, work_dir)

    # 7. 比对与共线性核心池
    print("\n" + "="*50)
    print(" ⚔️ 阶段 B：全基因组比对与共线性搜寻")
    print("="*50)
    
    for _, _, q_prefix in queries:
        run_diamond_and_synteny(ref_prefix, ref_faa_path, q_prefix, qry_faa_paths[q_prefix], work_dir, threads)

    # 8. 总结
    total_time = time.time() - start_time
    print("\n" + "🏆"*20)
    print(" 全自动化管线运行圆满结束！")
    print(f" ⏱️ 总耗时: {total_time/60:.2f} 分钟")
    print(f" 📂 所有结果都已安全存放在:\n   \033[92m{os.path.relpath(work_dir, base_dir)}/\033[0m")
    print("\n 💡 接下来，你只需进入该文件夹，直接执行绘图命令：")
    for _, _, q_prefix in queries:
        print(f"    python -m jcvi.graphics.dotplot {ref_prefix}.{q_prefix}.anchors")
    print("🏆"*20 + "\n")

if __name__ == "__main__":
    try:
        main()
    except KeyboardInterrupt:
        print("\n🛑 用户强制终止。")