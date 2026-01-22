import os
import glob
import sys
import re

def current_path_function():
    # 获取基础路径，兼容打包模式
    if getattr(sys, 'frozen', False):
        # 打包后的情况 - 使用可执行文件所在目录
        current_path = os.path.dirname(sys.executable)
    else:
        # 普通Python脚本运行的情况
        current_path = os.path.dirname(os.path.abspath(__file__))
    return current_path

def find_file(file_type):
    current_path = current_path_function()

    # 在当前及其子目录下搜索
    current_file_path = os.path.join(current_path, '**', f'*{file_type}')
    current_file_list = glob.glob(current_file_path, recursive=True)

    print('='*50)
    if len(current_file_list) >= 1:
        print(f'在当前目录及子目录下共发现{len(current_file_list)}个{file_type}文件: ')
        for file in current_file_list:
            print(file)
        return current_file_list
    else:
        response = input(f'在当前目录及子目录下没有发现{file_type}文件, 是否在父级目录搜索(y/n): \n').strip().lower()
        if response not in ['y', 'yes']:
            print("操作已取消, 程序已退出")
            exit()
        else:
            parent_path = os.path.dirname(current_path)
            # 在父级目录下搜索
            parent_file_path = os.path.join(parent_path, '**', f'*{file_type}')
            parent_file_list = glob.glob(parent_file_path, recursive=True)

            if len(parent_file_list) >= 1:
                print(f'在父目录下共发现{len(parent_file_list)}个{file_type}文件: ')
                for file in parent_file_list:
                    print(file)
                return parent_file_list
            else:
                print('未发现任何文件, 程序已退出')
                exit()

def choose_file(file_list):
    number = 1
    print('='*50)
    print(f'共发现{len(file_list)}个文件: ')
    for file in file_list:
        print(f'  [{number}]{file}')
        number += 1
    print('='*50)
    response = int(input('请选择目标文件: '))
    choosing_file = file_list[response - 1]
    return choosing_file

def extract_protein_sequences(gff_file):
    """
    从GFF文件中提取gene ID和蛋白序列
    """
    gene_proteins = {}
    current_gene = None
    current_sequence = []
    in_protein_section = False
    
    with open(gff_file, 'r', encoding='utf-8') as file:
        for line in file:
            line = line.strip()
            
            # 检测gene开始
            if line.startswith('# start gene'):
                match = re.search(r'g\d+', line)
                if match:
                    current_gene = match.group()
                    in_protein_section = False
                    current_sequence = []
            
            # 检测protein sequence开始
            elif line.startswith('# protein sequence = ['):
                in_protein_section = True
                # 提取第一行序列，去掉开头的"["
                seq_start = line.find('[') + 1
                if seq_start > 0:
                    seq_part = line[seq_start:].strip()
                    if seq_part:
                        current_sequence.append(seq_part)
            
            # 在protein section中继续读取序列
            elif in_protein_section and line.startswith('#'):
                # 检查是否到达protein section的结束
                if line.startswith('# Evidence for and against this transcript:'):
                    # 保存当前gene的蛋白序列
                    if current_gene and current_sequence:
                        full_sequence = ''.join(current_sequence)
                        # 移除可能的换行符和空格，以及末尾的"]"
                        full_sequence = full_sequence.replace('\n', '').replace(' ', '')
                        # 移除末尾的"]"如果存在
                        if full_sequence.endswith(']'):
                            full_sequence = full_sequence[:-1]
                        gene_proteins[current_gene] = full_sequence
                    
                    in_protein_section = False
                    current_sequence = []
                else:
                    # 提取序列行（去除#和空格）
                    seq_line = line[1:].strip()
                    if seq_line and not seq_line.startswith('protein sequence'):
                        # 移除末尾的"]"如果存在
                        if seq_line.endswith(']'):
                            seq_line = seq_line[:-1]
                        current_sequence.append(seq_line)
    
    return gene_proteins

def save_fasta(gene_proteins, output_file):
    """
    将gene ID和蛋白序列保存为FASTA格式
    """
    with open(output_file, 'w', encoding='utf-8') as file:
        for gene_id, protein_seq in gene_proteins.items():
            # 写入FASTA格式
            file.write(f'>{gene_id}\n')
            # 每行80个字符，符合FASTA格式标准
            for i in range(0, len(protein_seq), 80):
                file.write(protein_seq[i:i+80] + '\n')
    
    print(f'FASTA文件已保存: {output_file}')
    print(f'共提取了 {len(gene_proteins)} 个基因的蛋白序列')

def make_sure():
    print('='*50)
    response = input("是否继续? (y/n): \n").strip().lower()
    if response not in ['y', 'yes']:
        print("操作已取消, 程序已退出")
        return False
    else:
        return True

def main():
    # 查找并选择GFF文件
    gff_file = choose_file(find_file('.gff'))
    
    print(f'正在处理文件: {gff_file}')
    
    # 提取蛋白序列
    gene_proteins = extract_protein_sequences(gff_file)
    
    if not gene_proteins:
        print("未在文件中找到任何蛋白序列")
        return
    
    print(f'成功提取了 {len(gene_proteins)} 个基因的蛋白序列')
    
    # 确认是否继续
    if not make_sure():
        return
    
    # 获取保存文件名
    save_name = input('请输入保存文件名(不含扩展名): \n').strip()
    if not save_name:
        save_name = 'protein_sequences'
    
    # 构建输出文件路径
    current_path = current_path_function()
    output_file = os.path.join(current_path, f'{save_name}.fasta')
    
    # 保存FASTA文件
    save_fasta(gene_proteins, output_file)
    
    print('='*50)
    print("程序执行完成")

if __name__ == '__main__':
    main()