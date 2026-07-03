#!/usr/bin/env python3
import os
import sys

def parse_fasta(fasta_path):
    """解析FASTA文件，返回一个列表，每个元素为 (header, sequence)"""
    sequences = []
    current_header = None
    current_seq = []
    
    with open(fasta_path, 'r') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith('>'):
                if current_header:
                    sequences.append((current_header, ''.join(current_seq)))
                current_header = line
                current_seq = []
            else:
                current_seq.append(line)
        if current_header:
            sequences.append((current_header, ''.join(current_seq)))
            
    return sequences

def extract_gene_id(header):
    """从FASTA表头中提取基因ID"""
    header_content = header[1:]  # 去掉 '>'
    parts = header_content.split()
    transcript_id = parts[0]
    
    # 1. 优先寻找 gffread 默认生成的 'gene=' 标签
    for part in parts:
        if part.startswith('gene='):
            return part.split('=')[1]
            
    # 2. 如果没有 'gene=' 标签，尝试根据点号 '.' 分隔符猜测（如 Gene001.1 -> Gene001）
    if '.' in transcript_id:
        return '.'.join(transcript_id.split('.')[:-1])
        
    # 3. 如果都不符合，则将转录本ID直接视为基因ID
    return transcript_id

def main():
    print("=" * 60)
    print(" 🧬  基因组蛋白序列最长转录本提取工具 (自动命名版)  🧬")
    print("=" * 60)

    # 1. 获取并验证输入文件
    while True:
        input_fasta = input("\n[?] 请输入原始 FASTA 文件路径 (或输入 q 退出): ").strip()
        
        if input_fasta.lower() == 'q':
            print("[INFO] 程序已退出。")
            sys.exit(0)
            
        if not input_fasta:
            print("❌ 路径不能为空，请重新输入。")
            continue
            
        if os.path.exists(input_fasta):
            break
        else:
            print(f"❌ 错误: 找不到文件 '{input_fasta}'，请检查路径是否正确！")

    # 2. 自动生成输出文件名 (形如: speciesA_all.pep.fa -> speciesA_all.pep_longest.fa)
    file_root, file_ext = os.path.splitext(input_fasta)
    output_fasta = f"{file_root}_longest{file_ext}"
    print(f"[💡] 已自动生成输出路径: {output_fasta}")

    # 3. 安全检查：处理潜在的覆盖风险
    if os.path.exists(output_fasta):
        overwrite = input(f"⚠️  警告: 自动生成的文件 '{output_fasta}' 已存在！是否覆盖？(y/n): ").strip().lower()
        if overwrite != 'y':
            print("[INFO] 已取消操作，程序退出。")
            sys.exit(0)

    # 4. 开始执行过滤逻辑
    print("\n[运行中] 正在读取原始序列...")
    try:
        records = parse_fasta(input_fasta)
    except Exception as e:
        print(f"❌ 读取文件失败，错误原因: {e}")
        sys.exit(1)
    
    print(f"[运行中] 正在分析同源转录本 (共 {len(records)} 条序列)...")
    longest_transcripts = {}
    
    for header, seq in records:
        gene_id = extract_gene_id(header)
        seq_len = len(seq)
        
        if gene_id not in longest_transcripts:
            longest_transcripts[gene_id] = (header, seq)
        else:
            existing_header, existing_seq = longest_transcripts[gene_id]
            if seq_len > len(existing_seq):
                longest_transcripts[gene_id] = (header, seq)
                
    # 5. 写入输出文件
    print("[运行中] 正在写入最长转录本序列...")
    try:
        with open(output_fasta, 'w') as out:
            for gene_id, (header, seq) in longest_transcripts.items():
                out.write(f"{header}\n")
                # 规整格式，每80个氨基酸换一行
                for i in range(0, len(seq), 80):
                    out.write(f"{seq[i:i+80]}\n")
    except Exception as e:
        print(f"❌ 写入文件失败，错误原因: {e}")
        sys.exit(1)
                
    print("\n" + "=" * 60)
    print(f"🎉 恭喜！处理成功！")
    print(f"   🔹 原始总转录本数: {len(records)}")
    print(f"   🔹 过滤后代表基因数: {len(longest_transcripts)}")
    print(f"   💾 结果已保存至: {os.path.abspath(output_fasta)}")
    print("=" * 60 + "\n")

if __name__ == "__main__":
    main()