'''
.results expect_count 矩阵构建器

目录构型
    自动识别脚本父级目录下所有目标文件

支持.genes.results, .isoforms.results, 自定义文件格式
支持expected_count, TPM, FPKN, 自定义数据类型合并

列索引为.results所处文件夹名称

thz 2025/10/20
'''

import os
import glob
import sys
import pandas as pd

def find_file(file_type):
    '''文件获得函数'''
    # 获取基础路径，兼容打包模式
    if getattr(sys, 'frozen', False):
        # 打包后的情况 - 使用可执行文件所在目录
        current_path = os.path.dirname(sys.executable)
    else:
        # 普通Python脚本运行的情况
        current_path = os.path.dirname(os.path.abspath(__file__))

    # 在当前及其子目录下搜索
    current_file_path = os.path.join(current_path, '**', f'*{file_type}')
    current_file_list = glob.glob(current_file_path, recursive=True)

    print('='*50)
    if len(current_file_list) >= 1:
        print(f'在当前目录及子目录下共发现{len(current_file_list)}个{file_type}文件: ')
        for file in current_file_list:
            print(file)
        return(current_file_list)

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
                for file in parent_file_path:
                    print(file)
                    return parent_file_list
            else:
                print('未发现任何文件, 程序已退出')
                exit()

def merge_expected_counts(results_files, datatype):
    '''合并多个.results文件中的expected_count数据'''
    # 创建一个空字典来存储所有数据
    expression_data = {}
    
    # 遍历所有.results文件
    for file_path in results_files:
        # 获取文件所在文件夹的名称作为列名
        folder_name = os.path.basename(os.path.dirname(file_path))
        
        try:
            # 读取TSV文件
            df = pd.read_csv(file_path, sep='\t')
            
            # 检查必要的列是否存在
            if 'gene_id' not in df.columns or datatype not in df.columns:
                print(f"警告: 文件 {file_path} 缺少必要的列 'gene_id' 或 {datatype}")
                continue
            
            # 设置gene_id为索引
            df.set_index('gene_id', inplace=True)
            
            # 将expected_count列添加到字典中
            expression_data[folder_name] = df[datatype]
            
        except Exception as e:
            print(f"错误: 无法读取文件 {file_path}: {e}")
            continue
    
    if not expression_data:
        print("错误: 没有成功读取任何文件")
        return
    
    # 创建DataFrame对象并转置，使基因在行中，样本在列中
    expression_matrix = pd.DataFrame(expression_data)
    
    # 对列索引按英文字母升序排列
    expression_matrix = expression_matrix.reindex(sorted(expression_matrix.columns), axis=1)
    
    # 保存到TSV文件
    output_filename=f"{datatype}_matrix.tsv"
    expression_matrix.to_csv(output_filename, sep='\t', index=True)
    print(f"{datatype}矩阵已保存到脚本同级目录下: {output_filename}")
    
    return expression_matrix

def make_sure():
    '''确认函数'''
    print('='*50)
    response = input("是否继续? (y/n): ").strip().lower()
    if response in ['y', 'yes']:
        return True
    else:
        print("操作已取消, 程序已退出")
        exit()

def main():
    '''主函数'''
    print('='*50)
    response = int(input(f'请选择合并results类型: \n  [1].genes.results \n  [2].isoforms.results \n  [3]自定义 \n'))
    if response == 1:
        file_list = find_file('.genes.results')
    elif response == 2:
        file_list = find_file('.isoforms.results')
    elif response == 3:
        file_list = find_file(input('请输入文件格式: \n'))
    else:
        print('请选择正确的序号')

    print('='*50)
    response = int(input(f'请选择合并数据类型: \n  [1]expected_count \n  [2]TPM \n  [3]FPKM \n  [4]自定义 \n'))
    if response == 1:
        datatype = 'expected_count'
    elif response == 2:
        datatype = 'TPM'
    elif response == 3:
        datatype = 'FPKM'
    elif response == 4:
        datatype = input('请输入数据类型: ')
    else:
        print('请选择正确的序号')

    response = make_sure()
    if response == True:
        merge_expected_counts(file_list, datatype=datatype)

main()