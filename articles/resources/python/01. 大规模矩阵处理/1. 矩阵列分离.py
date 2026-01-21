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
        print(f'在当前目录及子目录下共发现{len(current_file_list)}个{file_type}文件: ')
        i = 0
        for file in current_file_list:
            print(f'  {[i]}{file}')
            i += 1
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
                print(f'在父目录下共发现{len(parent_file_list)}个{file_type}文件: ')
                i = 0
                for file in parent_file_list:
                    print(f'  {[i]}{file}')
                    i += 1
                return parent_file_list
            else:
                print('未发现任何文件, 程序已退出')
                exit()

def select_columns_interactively(df):
    # 显示所有可用的列
    print('='*50)
    print("可用的列:")
    for i, col in enumerate(df.columns):
        print(f"[{i}]: {col}")
    
    # 获取用户输入
    while True:
        try:
            user_input = input("\n请选择列索引(多个索引用逗号分隔，例如: 0,1,2,4): ")
            
            # 处理用户输入
            selected_indices = [int(idx.strip()) for idx in user_input.split(',')]
            
            # 验证索引是否有效
            for idx in selected_indices:
                if idx < 0 or idx >= len(df.columns):
                    print(f"错误: 索引 {idx} 超出范围。请使用 0 到 {len(df.columns)-1} 之间的索引。")
                    break
            else:
                # 所有索引都有效，创建新的DataFrame
                selected_columns = [df.columns[idx] for idx in selected_indices]
                new_df = df[selected_columns].copy()

                print('='*50)
                print(f"已选择列: {selected_columns}")
                return new_df
                
        except ValueError:
            print("错误: 请输入有效的数字索引，用逗号分隔。")
        except KeyboardInterrupt:
            print("\n用户中断操作。")
            return None
        except Exception as e:
            print(f"发生错误: {e}")

def choose_file(file_list):
    response = int(input('\n请选择目标文件: '))
    choosing_file = file_list[response]

    return choosing_file

def make_sure():
    print('='*50)
    response = input("是否保存? (y/n): \n").strip().lower()
    if response not in ['y', 'yes']:
        print("操作已取消, 程序已退出")
        exit()
    else:
        return True

def save_dataframe(select_df, save_name='default'):
    # 获取当前路径
    current_path = current_path_function()
    
    # 构建完整的文件路径
    output_filename = os.path.join(current_path, save_name)
    
    # 保存到TSV文件
    select_df.to_csv(output_filename, sep='\t', index=False)
    print(f"矩阵已保存到脚本同级目录下: {output_filename}")

def set_index_interactively(df):
    print('='*50)
    print("是否要将某一列设置为索引?")
    print("可用的列:")
    for i, col in enumerate(df.columns):
        print(f"[{i}]: {col}")
    
    while True:
        try:
            user_input = input("\n请选择要设置为索引的列索引(输入-1表示不设置索引): ").strip()
            
            if user_input == '-1':
                print("不设置任何列为索引")
                return df
                
            idx = int(user_input)
            if idx < 0 or idx >= len(df.columns):
                print(f"错误: 索引 {idx} 超出范围。请使用 0 到 {len(df.columns)-1} 之间的索引，或输入-1表示不设置索引。")
                continue
                
            # 设置索引
            index_column = df.columns[idx]
            df.set_index(index_column, inplace=True)
            print(f"已将列 '{index_column}' 设置为索引")
            return df
            
        except ValueError:
            print("错误: 请输入有效的数字索引。")
        except KeyboardInterrupt:
            print("\n用户中断操作。")
            return df
        except Exception as e:
            print(f"发生错误: {e}")

def main():
    file_list = find_file('.tsv')
    file = choose_file(file_list)
    with open(file, 'r') as file:
        df = pd.read_csv(file, sep='\t')
        martix_obj = pd.DataFrame(df)
        select_df = pd.DataFrame(select_columns_interactively(df=martix_obj))
        print(f'已选择列: \n  {select_df}')
        
        # 询问用户是否要设置索引
        select_df = set_index_interactively(select_df)
        
        response = make_sure()
        if response == True:
            # 保存到TSV文件
            save_name = input(f'请输入保存文件名: \n')
            save_dataframe(select_df, save_name=f'{save_name} segmentation_martix.tsv')

main()