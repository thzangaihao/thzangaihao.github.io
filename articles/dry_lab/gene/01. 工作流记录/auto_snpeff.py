#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import sys
import glob
import shutil
import subprocess

'''
SnpEff 变异注释大师版 (建库与批处理分离)
功能 1：一键编译自建基因组数据库 (只需运行一次)
功能 2：调用已有数据库，批量对多个 VCF 进行注释与靶标过滤
'''

def get_base_dir():
    return os.path.dirname(os.path.abspath(__file__))

def get_snpeff_dir():
    """自动寻找 Conda 环境中的 SnpEff 安装目录"""
    if 'CONDA_PREFIX' not in os.environ:
        print("❌ 错误: 未检测到 Conda 环境，请先执行: conda activate snpeff_env")
        sys.exit(1)
    
    prefix = os.environ['CONDA_PREFIX']
    snpeff_paths = glob.glob(os.path.join(prefix, 'share', 'snpeff*'))
    
    if not snpeff_paths:
        print("❌ 错误: 在当前环境中未找到 SnpEff，请确认是否安装成功。")
        sys.exit(1)
        
    return snpeff_paths[0]

def interactive_select(files, desc, multiple=False):
    """交互式选择 (支持单选或多选)"""
    if not files:
        print(f"⚠️ 未找到任何 {desc}！")
        return [] if multiple else None
    
    print(f"\n📂 扫描到以下 {len(files)} 个 {desc}:")
    for i, f in enumerate(files, 1):
        # 如果是提取出来的数据库名（字符串），直接打印
        if isinstance(f, str) and not os.path.exists(f) and not f.endswith(('.fasta', '.vcf', '.gff3')):
            print(f"  [{i}] {f}")
        else:
            print(f"  [{i}] {os.path.relpath(f, get_base_dir())}")
        
    while True:
        prompt = f"\n👉 请选择 (单选输入数字" + (", 多选如1,3 或all" if multiple else "") + ", q退出): "
        choice = input(prompt).strip().lower()
        
        if choice == 'q': sys.exit(0)
        if multiple and choice == 'all': return files
        
        try:
            if multiple:
                selected = []
                parts = choice.replace(' ', '').split(',')
                for part in parts:
                    if '-' in part:
                        start, end = map(int, part.split('-'))
                        selected.extend(files[start-1:end])
                    else:
                        selected.append(files[int(part)-1])
                return list(set(selected))
            else:
                idx = int(choice)
                if 1 <= idx <= len(files):
                    return files[idx-1]
        except:
            pass
        print("⚠️ 输入无效，请重新输入。")

def get_custom_genomes(snpeff_dir):
    """从配置文件中读取所有自建的数据库名"""
    config_file = os.path.join(snpeff_dir, "snpEff.config")
    genomes = []
    try:
        with open(config_file, 'r') as f:
            for line in f:
                if ".genome : " in line and "Custom Genome" not in line: # 排除掉自带注释
                    genome_name = line.split(".genome :")[0].strip()
                    genomes.append(genome_name)
    except:
        pass
    return genomes

def build_database(snpeff_dir):
    """建库模式：一次性工作"""
    print("\n" + "="*45)
    print(" 🛠️ 模式一：构建自定义参考基因组数据库")
    print("="*45)
    
    base_dir = get_base_dir()
    
    fastas = glob.glob(os.path.join(base_dir, "**", "*.fasta"), recursive=True)
    gffs = glob.glob(os.path.join(base_dir, "**", "*.gff3"), recursive=True)
    
    ref_fasta = interactive_select(fastas, "【建库】参考基因组序列 (.fasta)", multiple=False)
    if not ref_fasta: return
    
    ref_gff = interactive_select(gffs, "【建库】对应的注释文件 (.gff3) [Liftoff生成]", multiple=False)
    if not ref_gff: return
    
    genome_name = input("\n👉 请为这个数据库起个简短的名字 (例如 CK-3-5): ").strip()
    if not genome_name: genome_name = "Custom_Ref"
        
    data_dir = os.path.join(snpeff_dir, "data", genome_name)
    os.makedirs(data_dir, exist_ok=True)
    
    print(f"\n⚙️ 正在装载文件...")
    shutil.copy(ref_fasta, os.path.join(data_dir, "sequences.fa"))
    shutil.copy(ref_gff, os.path.join(data_dir, "genes.gff"))
    
    config_file = os.path.join(snpeff_dir, "snpEff.config")
    with open(config_file, 'r') as f:
        config_content = f.read()
        
    if f"{genome_name}.genome" not in config_content:
        with open(config_file, 'a') as f:
            f.write(f"\n# {genome_name} Custom Genome\n")
            f.write(f"{genome_name}.genome : {genome_name}\n")

    # 包含跳过检查的强力编译命令
    print(f"\n🚀 开始强力编译 {genome_name} 数据库...")
    cmd_build = f"snpEff build -gff3 -v -noCheckCds -noCheckProtein {genome_name}"
    
    try:
        subprocess.run(cmd_build, shell=True, check=True)
        print(f"\n✅ 成功！数据库 [{genome_name}] 已永久加入 SnpEff 核心系统！")
    except subprocess.CalledProcessError:
        print("\n❌ 编译失败，请检查上面终端的日志。")

def annotate_vcf_batch(snpeff_dir):
    """批量注释模式"""
    print("\n" + "="*45)
    print(" 🔬 模式二：批量变异注释与靶标锁定")
    print("="*45)
    
    # 1. 选择已建好的库
    genomes = get_custom_genomes(snpeff_dir)
    if not genomes:
        print("⚠️ 尚未建立任何自定义数据库！请先运行模式一建库。")
        return
        
    target_genome = interactive_select(genomes, "你要使用哪个数据库进行注释？", multiple=False)
    if not target_genome: return
    
    # 2. 批量选择 VCF
    base_dir = get_base_dir()
    vcfs = glob.glob(os.path.join(base_dir, "**", "*_SNP_filtered.vcf"), recursive=True) + \
           glob.glob(os.path.join(base_dir, "**", "*_SV_delly.vcf"), recursive=True)
           
    target_vcfs = interactive_select(vcfs, f"待分析的变异结果 (.vcf) [请选择对应 {target_genome} 的样本]", multiple=True)
    if not target_vcfs: return
    
    # 3. 循环批处理
    print(f"\n🔥 开始使用 [{target_genome}] 数据库批量处理 {len(target_vcfs)} 个样本...")
    
    for i, target_vcf in enumerate(target_vcfs, 1):
        out_dir = os.path.dirname(target_vcf)
        base_name = os.path.basename(target_vcf).replace('.vcf', '')
        
        annotated_vcf = os.path.join(out_dir, f"{base_name}_annotated.vcf")
        final_vcf = os.path.join(out_dir, f"{base_name}_HIGH_MODERATE.vcf")
        
        print(f"\n▶️ [任务 {i}/{len(target_vcfs)}] 正在处理: {base_name}")
        
        cmd_eff = f"snpEff -v {target_genome} {target_vcf} > {annotated_vcf}"
        subprocess.run(cmd_eff, shell=True, stderr=subprocess.DEVNULL)
        
        cmd_sift = f"cat {annotated_vcf} | SnpSift filter \"( ANN[*].IMPACT = 'HIGH' ) | ( ANN[*].IMPACT = 'MODERATE' )\" > {final_vcf}"
        subprocess.run(cmd_sift, shell=True, stderr=subprocess.DEVNULL)
        
        print(f"   ✅ 提纯完成 -> {final_vcf}")
        
    print(f"\n🎉 恭喜！所有样本均已完成变异注释与提纯，请去各个文件夹查收结果！")

if __name__ == "__main__":
    try:
        snpeff_dir = get_snpeff_dir()
        
        print("\n🌟 欢迎使用 SnpEff 自动变异注释大师")
        print("  [1] 构建新数据库 (建库)")
        print("  [2] 批量注释 VCF (分析)")
        print("  [q] 退出程序")
        
        choice = input("\n👉 请选择工作模式 (1或2): ").strip().lower()
        if choice == '1':
            build_database(snpeff_dir)
        elif choice == '2':
            annotate_vcf_batch(snpeff_dir)
        elif choice == 'q':
            sys.exit(0)
        else:
            print("⚠️ 输入有误，退出程序。")
            
    except KeyboardInterrupt:
        print("\n🛑 用户强制退出。")