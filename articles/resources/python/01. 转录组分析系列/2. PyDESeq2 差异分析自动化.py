import os
import glob
import sys
import pandas as pd

from pydeseq2.dds import DeseqDataSet
from pydeseq2.default_inference import DefaultInference
from pydeseq2.ds import DeseqStats

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

    if len(current_file_list) >= 1:
        number = 1
        print('='*50)
        print(f'共发现{len(current_file_list)}个文件: ')
        for file in current_file_list:
            print(f'  [{number}]{file}')
            number += 1
        return current_file_list

    else:
        print('='*50)
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
                number = 1
                print('='*50)
                print(f'共发现{len(parent_file_list)}个文件: ')
                for file in parent_file_list:
                    print(f'  [{number}]{file}')
                    number += 1
                return parent_file_list
            else:
                print('未发现任何文件, 程序已退出')
                exit()

def choose_file(file_list, text):
    print('='*50)
    response = int(input(f'请选择{text}目标文件: \n'))
    if response >0 and response <= len(file_list):
        choosing_file = file_list[response - 1]
        print(f'已选择: {choosing_file}')
        tsv_file = pd.read_csv(choosing_file, sep='\t', index_col=0)
        print(tsv_file)
        return choosing_file
    else:
        print('请输入正确序号, 程序已退出')
        exit()
    
def save_dataframe(select_df, save_name='default'):
    # 获取当前路径
    current_path = current_path_function()
    
    # 构建完整的文件路径
    output_filename = os.path.join(current_path, save_name)
    
    # 保存到TSV文件
    select_df.to_csv(output_filename, sep='\t', index=True)
    print(f"矩阵已保存到脚本同级目录下: {output_filename}")

def make_sure():
    print('='*50)
    response = input("是否保存? (y/n): \n").strip().lower()
    if response not in ['y', 'yes']:
        print("操作已取消, 程序已退出")
        exit()
    else:
        return True

def main():
    count_martix_file_path_list = find_file('.tsv')
    count_martix_file_path = choose_file(file_list=count_martix_file_path_list, text='计数矩阵count_martix')
    metadata_martix_file_path = choose_file(file_list=count_martix_file_path_list, text='元数据metadata_martix')

    # 读取TSV格式数据
    counts_df = pd.read_csv(count_martix_file_path, sep='\t', index_col=0)
    metadata = pd.read_csv(metadata_martix_file_path, sep='\t', index_col=0)

    # 过滤掉表达量过低的基因(总counts < 10)
    print('='*50)
    print('开始执行数据过滤')
    genes_to_keep = counts_df.columns[counts_df.sum(axis=0) >= 10]
    counts_df = counts_df[genes_to_keep]
    print(f"过滤后基因数量: {len(counts_df.columns)}")

    # 创建推断对象(设置CPU数量)
    inference = DefaultInference(n_cpus=8)
    # 创建DeseqDataSet对象
    dds = DeseqDataSet(
        counts=counts_df,
        metadata=metadata,
        design_factors=['condition'],  # 使用condition作为设计因子
        refit_cooks=True,              # 推荐：重新拟合Cook's异常值
        inference=inference
    )

    # 运行负二项分析(核心步骤)
    print('='*50)
    print('开始执行分析')
    dds.deseq2()

    print('='*50)
    print('分析完毕, DeseqDataSet已拓展为AnnData类: ')
    print(dds)
    print('='*50)
    print('下面输出线性拟合数据: ')
    # 查看离散度估计
    print("\n基因离散度:")
    print(dds.var["dispersions"])
    # 查看LFCs (log fold changes)
    print("\nLog Fold Changes (LFC):")
    print(dds.varm["LFC"])

    # 运行Wald检验
    # 创建统计检验DeseqStats对象
    ds = DeseqStats(dds, contrast=['condition', 'A', 'B'], inference=inference)
    print('='*50)
    print('Wald检验')
    ds.summary()
    ds_result = ds.results_df

    # 对LFC进行收缩(减少噪音)
    print('='*50)
    print('LFC收缩')
    ds.lfc_shrink(coeff="condition[T.B]")
    ds_result_shrink = ds.results_df

    # 保存输出结果
    response = make_sure()
    if response == True:
        # 保存到TSV文件
        print('='*50)
        save_dataframe(select_df=ds_result, save_name='ds_result.tsv')
        save_dataframe(select_df=ds_result_shrink, save_name='ds_result_shrink.tsv')

main()