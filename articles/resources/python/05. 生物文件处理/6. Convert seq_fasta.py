'''
seq合并fasta脚本

目录构型
    自动检测脚本同级目录及子目录中的.seq文件

以seq文件名作为fasta头序列合并为一个fasta文件

thz 2025/10/18
'''

import os
import glob
import sys

def seq_to_fasta():
    # 获取exe文件所在目录（适用于打包环境）
    if getattr(sys, 'frozen', False):
        script_dir = os.path.dirname(sys.executable)
    else:
        script_dir = os.path.dirname(os.path.abspath(__file__))
    
    # 查找所有.seq文件
    seq_files = glob.glob(os.path.join(script_dir, "*.seq"))
    
    if not seq_files:
        print(f"未在目录 {script_dir} 中找到任何.seq文件")
        print(f"当前工作目录: {os.getcwd()}")
        return
    
    # 创建输出FASTA文件名
    output_file = os.path.join(script_dir, "combined_sequences.fasta")
    
    # 处理每个.seq文件
    processed_count = 0
    with open(output_file, 'w') as fasta_out:
        for seq_file in seq_files:
            # 获取文件名（不含扩展名）作为序列头
            file_name = os.path.splitext(os.path.basename(seq_file))[0]
            
            try:
                # 读取序列内容
                with open(seq_file, 'r') as f:
                    sequence = f.read().strip()
                
                # 写入FASTA格式
                fasta_out.write(f">{file_name}\n")
                
                # 按每行80个字符格式化序列（FASTA标准格式）
                for i in range(0, len(sequence), 80):
                    fasta_out.write(sequence[i:i+80] + "\n")
                
                processed_count += 1
                    
            except Exception as e:
                print(f"处理文件 {seq_file} 时出错: {str(e)}")
    
    print(f"FASTA文件已生成: {output_file}")
    print(f"共处理了 {processed_count} 个序列文件")

if __name__ == "__main__":
    seq_to_fasta()