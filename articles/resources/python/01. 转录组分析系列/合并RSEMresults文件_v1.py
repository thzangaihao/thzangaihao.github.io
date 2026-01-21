'''
.results expect_count 矩阵构建器（修改版）
支持根据 .results 文件名作为样本名称进行矩阵合并
thz 修改：2025/10/20 → 你要求的更新
'''

import os
import glob
import sys
import pandas as pd

def find_file(file_type):
    '''文件获得函数'''
    if getattr(sys, 'frozen', False):
        current_path = os.path.dirname(sys.executable)
    else:
        current_path = os.path.dirname(os.path.abspath(__file__))

    current_file_path = os.path.join(current_path, '**', f'*{file_type}')
    current_file_list = glob.glob(current_file_path, recursive=True)

    print('='*50)
    if len(current_file_list) >= 1:
        print(f'在当前目录及子目录下共发现{len(current_file_list)}个{file_type}文件:')
        for file in current_file_list:
            print(file)
        return current_file_list
    else:
        response = input(f'未找到{file_type}文件，是否在父级目录搜索？(y/n): ').strip().lower()
        if response not in ['y','yes']:
            print("操作已取消")
            exit()

        parent_path = os.path.dirname(current_path)
        parent_file_path = os.path.join(parent_path, '**', f'*{file_type}')
        parent_file_list = glob.glob(parent_file_path, recursive=True)

        if len(parent_file_list) >= 1:
            print(f'在父级目录共发现{len(parent_file_list)}个{file_type}文件:')
            for file in parent_file_list:
                print(file)
            return parent_file_list
        else:
            print("仍未找到文件，程序退出")
            exit()

def merge_results(results_files, datatype):
    '''合并多个.results 文件'''
    expression_data = {}

    for file_path in results_files:
        # 样本名 = 文件名（去扩展）
        filename = os.path.basename(file_path)
        sample_name = filename.replace('.genes.results', '').replace('.isoforms.results', '').replace('.results', '')

        try:
            df = pd.read_csv(file_path, sep='\t')

            # 清理列名（防止前后空格）
            df.columns = df.columns.str.strip()

            if 'gene_id' not in df.columns or datatype not in df.columns:
                print(f"警告: 文件 {file_path} 缺少 gene_id 或 {datatype} 列，跳过")
                continue

            df.set_index('gene_id', inplace=True)
            expression_data[sample_name] = df[datatype]

        except Exception as e:
            print(f"错误: 无法读取文件 {file_path}: {e}")
            continue

    if not expression_data:
        print("错误: 所有文件均未成功读取")
        return

    expression_matrix = pd.DataFrame(expression_data)

    # 列按字母排序
    expression_matrix = expression_matrix.reindex(sorted(expression_matrix.columns), axis=1)

    # ★★ 输出文件到脚本目录 ★★
    if getattr(sys, 'frozen', False):
        script_dir = os.path.dirname(sys.executable)
    else:
        script_dir = os.path.dirname(os.path.abspath(__file__))

    output_filename = os.path.join(script_dir, f"{datatype}_matrix.tsv")
    expression_matrix.to_csv(output_filename, sep='\t', index=True)

    print(f"\n{datatype} 矩阵已保存至：{output_filename}")
    return expression_matrix

def make_sure():
    print('='*50)
    response = input("是否继续？(y/n): ").strip().lower()
    if response in ['y', 'yes']:
        return True
    else:
        print("操作已取消")
        exit()

def main():
    print('='*50)
    response = int(input("请选择合并类型:\n [1] genes.results\n [2] isoforms.results\n [3] 自定义\n"))
    if response == 1:
        file_list = find_file('.genes.results')
    elif response == 2:
        file_list = find_file('.isoforms.results')
    elif response == 3:
        file_list = find_file(input('请输入自定义后缀（例如 .txt）:\n'))
    else:
        print("无效输入"); exit()

    print('='*50)
    response = int(input("请选择合并数据类型:\n [1] expected_count\n [2] TPM\n [3] FPKM\n [4] 自定义\n"))
    if response == 1:
        datatype = 'expected_count'
    elif response == 2:
        datatype = 'TPM'
    elif response == 3:
        datatype = 'FPKM'
    elif response == 4:
        datatype = input("请输入列名:\n")
    else:
        print("无效输入"); exit()

    if make_sure():
        merge_results(file_list, datatype)

main()
