import os
import glob
import sys
import pandas as pd

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
        return current_file_list

    else:
        response = input(f'在当前目录及子目录下没有发现{file_type}文件, 是否在父级目录搜索(y/n): ').strip().lower()
        if response not in ['y', 'yes']:
            print("操作已取消, 程序已退出")
            exit()
        else:
            parent_path = os.path.dirname(current_path)

            # 在父级目录下搜索
            parent_file_path = os.path.join(parent_path, '**', f'*{file_type}')
            parent_file_list = glob.glob(parent_file_path, recursive=True)

            if len(parent_file_list) >= 1:
                return parent_file_list
            else:
                print('未发现任何文件, 程序已退出')
                exit()

def make_sure():
    print('='*50)
    response = input("是否保存? (y/n): \n").strip().lower()
    if response not in ['y', 'yes']:
        print("操作已取消, 程序已退出")
        exit()
    else:
        return True

def choose_file(file_list):
    number = 1
    print(f'共发现{len(file_list)}个文件: ')
    for file in file_list:
        print(f'  [{number}]{file}')
        number += 1
    response = int(input('请选择目标文件: '))
    choosing_file = file_list[response - 1]

    return choosing_file

def save_dataframe(select_df, save_name='default'):
    # 获取当前路径
    current_path = current_path_function()
    
    # 构建完整的文件路径
    output_filename = os.path.join(current_path, save_name)
    
    # 保存到TSV文件
    select_df.to_csv(output_filename, sep='\t', index=True)
    print('='*50)
    print(f"矩阵已保存到脚本同级目录下: {output_filename}")

import pandas as pd

def round_dataframe(df):
    # 创建DataFrame的副本，避免修改原始数据
    rounded_df = df.copy()
    
    # 选择数值类型的列
    numeric_columns = rounded_df.select_dtypes(include=['number']).columns
    
    # 对数值列进行四舍五入并转换为整数
    rounded_df[numeric_columns] = rounded_df[numeric_columns].round().astype(int)
    
    return rounded_df

def main():
    print('='*50)
    print('矩阵数据量子化\n四舍五入法')

    file_path = choose_file(find_file('.tsv'))
    df = round_dataframe(pd.DataFrame(pd.read_csv(file_path, sep='\t')))
    print('='*50)
    print(f'已成功量子化: \n  {df}')
    df.set_index('gene_id', inplace=True)
    
    response = make_sure()
    if response == True:
        # 保存到TSV文件
        save_dataframe(select_df=df, save_name='count_martix_quantization.tsv')


main()