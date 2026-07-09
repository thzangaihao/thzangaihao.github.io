#!/usr/bin/env python3
import os
import sys
import glob
import csv
import argparse
from collections import Counter

def reverse_complement(seq):
    """获取序列的反向互补序列"""
    mapping = str.maketrans('ATCGN', 'TAGCN')
    return seq.upper().translate(mapping)[::-1]

def get_circular_shifts(seq):
    """获取一个序列的所有循环位移"""
    return set(seq[i:] + seq[:i] for i in range(len(seq)))

def parse_fasta(filepath):
    """简单的FASTA解析器 (按条目生成，节约内存)"""
    with open(filepath, 'r') as f:
        header = ""
        seq = []
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if header:
                    yield header, "".join(seq)
                header = line[1:].split()[0]
                seq = []
            else:
                seq.append(line.upper())
        if header:
            yield header, "".join(seq)

def count_motif(sequence, motif):
    """统计序列中该motif及其所有循环位移出现的总次数（非重叠）"""
    shifts = get_circular_shifts(motif)
    best_count = 0
    for shift in shifts:
        best_count = max(best_count, sequence.count(shift))
    return best_count

def find_global_motif(sequences, k):
    """【新逻辑】在所有端部序列汇总中寻找最常见的串联k-mer"""
    counts = Counter()
    for seq in sequences:
        kmers = [seq[i:i+k] for i in range(len(seq) - k + 1)]
        counts.update(kmers)
        
    for kmer, _ in counts.most_common():
        if len(set(kmer)) > 1: # 排除单碱基重复，如 AAAAAA
            return kmer
    return None

def process_fasta(filepath, mode, window_size, min_count, motif_5=None, k_len=None):
    """处理单个FASTA文件：先寻找全局Motif，再统一扫描"""
    results = []
    
    # 用于存储提取出来的端部信息，避免重复读取大文件
    # 格式: (contig_id, seq_len, left_window, right_window)
    contigs_data = [] 
    
    left_windows_for_discovery = []
    right_windows_for_discovery = []
    
    # 步骤 1: 遍历 FASTA，提取所有端部序列
    for contig_id, seq in parse_fasta(filepath):
        seq_len = len(seq)
        
        if seq_len < window_size * 2:
            print(f"  [提示] Contig {contig_id} 长度({seq_len})小于两倍搜索窗口，仍会扫描，但不参与全局 Motif 寻找。")
            # 序列太短，直接把全长当做窗口存起来
            contigs_data.append((contig_id, seq_len, seq, seq))
            continue
            
        left_win = seq[:window_size]
        right_win = seq[-window_size:]
        
        contigs_data.append((contig_id, seq_len, left_win, right_win))
        left_windows_for_discovery.append(left_win)
        right_windows_for_discovery.append(right_win)

    if not contigs_data:
        return []

    # 步骤 2: 确定该 FASTA 文件的全局 Motif
    global_motif_5 = motif_5
    global_motif_3 = None
    
    if mode == 'denovo':
        # 把所有长 contig 的左端汇总，寻找最强信号
        global_motif_5 = find_global_motif(left_windows_for_discovery, k_len)
        if global_motif_5:
            global_motif_3 = reverse_complement(global_motif_5)
        else:
            # 如果左端什么都没找到，尝试去所有右端找
            global_motif_3 = find_global_motif(right_windows_for_discovery, k_len)
            if global_motif_3:
                global_motif_5 = reverse_complement(global_motif_3)
            else:
                print(f"  [警告] 在该文件中未能自动找到有效的长度为 {k_len} 的重复序列！跳过。")
                return []
    else:
        # fixed 模式
        if global_motif_5:
            global_motif_3 = reverse_complement(global_motif_5)

    print(f"  [成功] 锁定该文件全局 Motif -> 5'端: {global_motif_5} | 3'端: {global_motif_3}")

    # 步骤 3: 使用全局唯一的 Motif 统一扫描之前存下来的所有 Contig 端部
    for contig_id, seq_len, left_win, right_win in contigs_data:
        left_count = count_motif(left_win, global_motif_5)
        right_count = count_motif(right_win, global_motif_3)
        
        has_left = "Yes" if left_count >= min_count else "No"
        has_right = "Yes" if right_count >= min_count else "No"
        has_both = "Yes" if (has_left == "Yes" and has_right == "Yes") else "No"
        
        results.append({
            'contig_id': contig_id,
            'contig_len': seq_len,
            'left_telo_count': left_count,
            'right_telo_count': right_count,
            'has_left_telo': has_left,
            'has_right_telo': has_right,
            'has_both': has_both,
            'detected_5_motif': global_motif_5,
            'detected_3_motif': global_motif_3
        })
        
    return results

def main():
    parser = argparse.ArgumentParser(description="端粒序列批量扫描工具 (全局唯一Motif优化版)")
    parser.add_argument('-i', '--input', nargs='+', help="指定的FASTA文件列表 (支持多个，空格分隔)")
    parser.add_argument('-m', '--mode', choices=['fixed', 'denovo'], help="模式：fixed 或 denovo")
    parser.add_argument('-s', '--seq5', help="5端固定序列 (例如 AATCCC)")
    parser.add_argument('-k', '--kmer', type=int, help="De novo模式的k-mer长度 (例如 6)")
    parser.add_argument('-w', '--window', type=int, default=1000, help="两端搜索范围(bp)，默认1000")
    parser.add_argument('-c', '--min_count', type=int, default=3, help="判定为存在的最小重复次数，默认3")
    
    args = parser.parse_args()
    
    is_interactive = len(sys.argv) == 1
    target_files = []
    script_dir = os.path.dirname(os.path.abspath(__file__))
    
    local_fastas = (
        glob.glob(os.path.join(script_dir, "*.fa")) + 
        glob.glob(os.path.join(script_dir, "*.fasta")) + 
        glob.glob(os.path.join(script_dir, "*.fna"))
    )

    if is_interactive:
        print("====== 端粒批量扫描程序 (交互模式) ======")
        print(f"当前搜索目录: {script_dir}")
        
        if not local_fastas:
            print("当前目录下未找到任何 FASTA 文件！")
            return
            
        print("\n找到以下 FASTA 文件：")
        for i, f in enumerate(local_fastas):
            filename = os.path.basename(f)
            print(f"[{i+1}] {filename}")
        
        file_choice = input("\n请选择要处理的文件编号 (输入 'all' 处理全部，输入如 '1,3' 处理特定多个): ").strip().lower()
        
        if file_choice == 'all':
            target_files = local_fastas
        else:
            try:
                indices = [int(idx.strip()) - 1 for idx in file_choice.split(',')]
                target_files = [local_fastas[i] for i in indices if 0 <= i < len(local_fastas)]
            except ValueError:
                print("输入格式错误！请确保输入的是数字和逗号。")
                return
            
        mode_choice = input("请选择模式 (1. 指定5'端序列  2. 全局自动寻找最短重复): ").strip()
        if mode_choice == '1':
            args.mode = 'fixed'
            args.seq5 = input("请输入5'端端粒重复序列 (如 CCCTAA): ").strip().upper()
        else:
            args.mode = 'denovo'
            args.kmer = int(input("请输入最短重复序列长度 (如 5, 6, 7): ").strip())
            
        user_win = input("请输入两端搜索范围 (默认 1000，直接回车使用默认值): ").strip()
        args.window = int(user_win) if user_win else 1000

    else:
        if args.input:
            for f_pattern in args.input:
                full_pattern = os.path.join(script_dir, f_pattern) if not os.path.isabs(f_pattern) else f_pattern
                target_files.extend(glob.glob(full_pattern))
            target_files = list(set(target_files)) 
        else:
            target_files = local_fastas
            
        if not args.mode:
            print("错误：非交互模式下必须通过 -m 指定模式！")
            return

    if not target_files:
        print("未找到需要处理的文件！")
        return

    print(f"\n即将处理 {len(target_files)} 个文件...")
    
    for fasta_file in target_files:
        filename_only = os.path.basename(fasta_file)
        print(f"\n>>> 正在分析: {filename_only}")
        
        out_csv = f"{os.path.splitext(fasta_file)[0]}_telomere.csv"
        
        results = process_fasta(
            filepath=fasta_file,
            mode=args.mode,
            window_size=args.window,
            min_count=args.min_count,
            motif_5=args.seq5,
            k_len=args.kmer
        )
        
        if not results:
            print(f"  [跳过] 未提取到有效 contig 或未找到 Motif，不生成 csv。")
            continue
            
        fieldnames = ['contig_id', 'contig_len', 'left_telo_count', 'right_telo_count', 
                      'has_left_telo', 'has_right_telo', 'has_both', 
                      'detected_5_motif', 'detected_3_motif']
                      
        with open(out_csv, 'w', newline='', encoding='utf-8') as csvfile:
            writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
            writer.writeheader()
            writer.writerows(results)
            
        print(f"  [完成] 结果已保存至: {os.path.basename(out_csv)}")

if __name__ == "__main__":
    main()