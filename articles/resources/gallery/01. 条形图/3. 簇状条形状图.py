import os
import glob
import sys
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

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

def make_sure():
    print('='*50)
    response = input("是否继续? (y/n): \n").strip().lower()
    if response not in ['y', 'yes']:
        print("操作已取消, 程序已退出")
        return False
    else:
        return True

def choose_file(file_list, text):
    number = 1
    print('='*50)
    print(f'共发现{len(file_list)}个文件: ')
    for file in file_list:
        print(f'  [{number}]{file}')
        number += 1
    print('='*50)
    response = int(input(f'请选择{text}文件: '))
    choosing_file = file_list[response - 1]

    return choosing_file

def plot_grouped_barchart(df, t_test_df_ck1, t_test_df_ck2):
    print(df)

    # 设置中文字体
    plt.rcParams['font.sans-serif'] = ['SimHei', 'Arial Unicode MS', 'DejaVu Sans']
    plt.rcParams['axes.unicode_minus'] = False

    # 1. 数据准备
    x_labels = df.iloc[:, 0].tolist()

    # 1.1 柱状图数据
    y_data_aver_dict = {'CK1':[], 'T':[], 'CK2':[]}
    for i in range(0, len(x_labels)):
        y_axes = df.iloc[i, :].tolist()[1:]
        # 计算均值
        ck1_average = float(np.mean(y_axes[0:3]))
        t_average = float(np.mean(y_axes[3:6]))
        ck2_average = float(np.mean(y_axes[6:]))

        y_data_aver_dict['CK1'].append(ck1_average)
        y_data_aver_dict['T'].append(t_average)
        y_data_aver_dict['CK2'].append(ck2_average)

    # 1.2 散点图数据
    y_scatter_dict = {'CK1':[], 'T':[], 'CK2':[]}
    for i in range(0, len(x_labels)):
        y_axes = df.iloc[i, :].tolist()[1:]

        y_scatter_dict['CK1'].append(y_axes[0:3])
        y_scatter_dict['T'].append(y_axes[3:6])
        y_scatter_dict['CK2'].append(y_axes[6:])

    # 1.3 标准差数据
    y_err_dict = {'CK1':[], 'T':[], 'CK2':[]}
    for i in range(0, len(x_labels)):
        y_axes = df.iloc[i, :].tolist()[1:]
        # 计算标准差
        ck1_std = float(np.std(y_axes[0:3]))
        t_std = float(np.std(y_axes[3:6]))
        ck2_std = float(np.std(y_axes[6:]))

        y_err_dict['CK1'].append(ck1_std)
        y_err_dict['T'].append(t_std)
        y_err_dict['CK2'].append(ck2_std)

    # 1.4 显著性数据
    p_values_ck1_vs_t = t_test_df_ck1.iloc[: ,2].tolist() # 对t检验tsv的DF第三列切片为列表
    p_values_ck2_vs_t = t_test_df_ck2.iloc[: ,2].tolist()

    # 2. 创建画布并分离对象
    fig, ax = plt.subplots(figsize=(16, 10), dpi=100)
    x = np.arange(len(x_labels)) # 设置x轴初始刻度，初始x为[0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10...]
    
    # 3. 设置图像
    # 3.1 柱状图绘制
    width = 0.25  # 柱子宽度
    multiplier = 0
    colors = ['#f5785f', '#f7a964', '#4ac0d5']
    for group, measurement in y_data_aver_dict.items():
        offset = width * multiplier # 偏移量
        ax.bar(x=x + offset, height=measurement, width=width, label=group, color=colors[multiplier], # 柱状图设置
               yerr=y_err_dict[group], capsize=5, error_kw={'elinewidth': 1, 'markeredgewidth': 1} # 误差线设置
               )
        multiplier += 1
    
    # 3.2 绘制散点图
    multiplier = 0
    x = np.arange(len(x_labels))   # 重新定义x，因为之前被修改过
    for group, measurement in y_scatter_dict.items():
        offset = width * multiplier # 偏移量
        for i, me_list in enumerate(measurement):
            # 为每个柱子的三个重复值生成x坐标
            x_positions = [x[i] + offset - width/4, x[i] + offset, x[i] + offset + width/4]
            # 绘制三个散点
            ax.scatter(x_positions, me_list, color="#000000", s=5, zorder=5)
        multiplier += 1

    # 3.3 绘制显著性标记    
    # 为每个时间点添加显著性标记
    for i in range(min(len(x_labels), len(p_values_ck1_vs_t))): # 当数据不足时，选用p_values_ck1_vs_t的长度只画有数据的地方
        # CK1 vs T 的显著性标记
        x1 = x[i]                  # CK1的水平位置，作为基准点
        x2 = x[i] + width - 0.05   # T的水平位置，偏移1个柱宽,减去0.05是为了和另一个标注线分开
        # y_max 两个柱子间最大的高度（柱高+误差线高）再偏移两个单位（+2）
        y_max = max(y_data_aver_dict['CK1'][i]+y_err_dict['CK1'][i], y_data_aver_dict['T'][i]+y_err_dict['T'][i]) + 2
        # y_1 左边柱子的高度（柱高+误差线高）再偏移两个单位（+2）
        y_1 = y_data_aver_dict['CK1'][i] + y_err_dict['CK1'][i] + 2
        # y_2 左边柱子的高度（柱高+误差线高）再偏移两个单位（+2）
        y_2 = y_data_aver_dict['T'][i] + y_err_dict['T'][i] + 2
        add_significance_bar(ax, x1, x2, y_max, y_1, y_2, p_values_ck1_vs_t[i])
        # CK2 vs T 的显著性标记
        x1 = x[i] + width + 0.05    # T的水平位置，偏移1个柱宽
        x2 = x[i] + 2*width         # CK2的水平位置，偏移2个柱宽
        y_max = max(y_data_aver_dict['CK2'][i]+y_err_dict['CK2'][i], y_data_aver_dict['T'][i]+y_err_dict['T'][i]) + 2
        y_1 = y_data_aver_dict['T'][i] + y_err_dict['T'][i] + 2
        y_2 = y_data_aver_dict['CK2'][i] + y_err_dict['CK2'][i] + 2
        add_significance_bar(ax, x1, x2, y_max, y_1, y_2, p_values_ck2_vs_t[i])

    # 3.4 其它设置
    ax.set_ylabel('生长直径 (mm)')
    ax.set_xlabel('接种日期')

    ax.set_title(df.columns.tolist()[0], loc='left') # 以数据tsv的列索引第一个作为图标题

    ax.set_xticks(x + width, x_labels, rotation=45) # 因为有三个，把刻度放中间加一个width就行了
    
    ax.legend(loc='upper right')

    ax.set_ylim(0, 100) # 设置y轴范围
    ax.set_yticks(range(0, 101, 10)) # 在0到100（左闭右开）手动设置步进为10的刻度

    # 4. 显示图像
    plt.show()

def add_significance_bar(ax, x1, x2, y_max, y_1, y_2, p_value):
    # 根据p值确定星号数量
    if p_value < 0.001:
        stars = '***'
    elif p_value < 0.01:
        stars = '**'
    elif p_value < 0.05:
        stars = '*'
    else:
        stars = 'ns'
        
    # 计算线条和文本位置
    bar_height = 2
    text_y = y_max + bar_height
    
    # 绘制横线
    ax.plot([x1, x1, x1, x2, x2], [y_1, y_max, y_max+bar_height, y_max+bar_height, y_2], 
            lw=1.5, c='black')
    
    # 添加星号
    ax.text((x1+x2)*0.5, text_y, stars, ha='center', va='bottom', color='black')

def main():
    # 选择并读取数据文件
    file_path = choose_file(find_file('.tsv'), '数据')
    # 读取文件，将第一列作为字符串处理
    df = pd.read_csv(file_path, sep='\t', dtype={0: str})
    
    # 选择t检验输出文件
    file_path = choose_file(find_file('.tsv'), 'CK1-T t检验')
    t_test_df_ck1 = pd.read_csv(file_path, sep='\t')
    file_path = choose_file(find_file('.tsv'), 'CK2-T t检验')
    t_test_df_ck2 = pd.read_csv(file_path, sep='\t')

    # 绘制柱状图
    plot_grouped_barchart(df, t_test_df_ck1, t_test_df_ck2)

if __name__ == "__main__":
    main()