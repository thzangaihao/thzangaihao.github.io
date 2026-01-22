import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import os

# ================= 参数设置区 =================
# 1. 画布设置
FIG_WIDTH = 6
FIG_HEIGHT = 200
FIG_DPI = 100        

# 2. X轴设置 (横坐标)
X_LABEL_SIZE = 12    
X_LABEL_ROT = 90     # 建议长文字保持 45 度或 90 度

# 3. Y轴设置 (纵坐标)
Y_LABEL_SIZE = 10    
Y_LABEL_ROT = 0      

# 4. 【新增】边距设置 (关键参数)
# 这是一个比例值 (0~1)，0.1 代表底部留 10% 空白，0.3 代表留 30%。
# 如果字还是被切掉，请把这个数字改大，比如改成 0.25 或 0.3
BOTTOM_MARGIN = 0.2
# ============================================

# 数据准备
current_dir = os.path.dirname(os.path.abspath(__file__))
excel_dir = os.path.join(current_dir, 'figure.xlsx')

# 容错：如果没有文件生成测试数据
if not os.path.exists(excel_dir):
    print("生成长名字测试数据...")
    # 模拟很长的横坐标名字
    long_names = [f'Very_Long_Sample_Name_{i}_Condition_XY' for i in range(1, 11)]
    data = pd.DataFrame(np.random.rand(10, 10), 
                        columns=long_names,
                        index=[f'Gene_{i}' for i in range(1, 11)])
else:
    data = pd.read_excel(excel_dir, index_col=0)

data_tuple = data.values
data_index = data.index 
data_columns = data.columns 

# 主图及子图创建
figure_screen = plt.figure(figsize=(FIG_WIDTH, FIG_HEIGHT), dpi=FIG_DPI)
ax = figure_screen.subplots()

# 热图及颜色棒绘制
imshow = ax.imshow(data_tuple, cmap='coolwarm', aspect='auto')
figure_screen.colorbar(mappable=imshow)

# 坐标轴设置
ax.set_xlabel('Experimental Group', fontsize=X_LABEL_SIZE + 2)
ax.set_ylabel('Physiological Parameters', fontsize=Y_LABEL_SIZE + 2)

ax.set_xticks(np.arange(len(data_columns)))
ax.set_xticklabels(data_columns, fontsize=X_LABEL_SIZE, rotation=X_LABEL_ROT, ha='right') 
# ha='right' 让旋转后的文字右对齐，看起来更整齐，且不容易挡住下面的字

ax.set_yticks(np.arange(len(data_index)))
ax.set_yticklabels(data_index, fontsize=Y_LABEL_SIZE, rotation=Y_LABEL_ROT)

# ================= 核心修改 =================
# 以前用的是 tight_layout()，现在改用 subplots_adjust 手动强行留白
# bottom=BOTTOM_MARGIN 意思就是把图的下边缘往上提，给字留出位置
figure_screen.subplots_adjust(bottom=BOTTOM_MARGIN)
# ============================================

print("当前底部留白比例:", BOTTOM_MARGIN)
plt.show()