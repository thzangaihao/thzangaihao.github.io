import pandas as pd
import os
import sys

def main():
    # 1. 将脚本运行环境自动切换为脚本所在目录下
    try:
        # 兼容打包后的执行环境
        if getattr(sys, 'frozen', False):
            script_dir = os.path.dirname(sys.executable)
        else:
            script_dir = os.path.dirname(os.path.abspath(__file__))
        os.chdir(script_dir)
        print(f"当前工作目录已切换至: {script_dir}")
    except Exception as e:
        print(f"切换目录失败: {e}")
        return

    # 2. 在终端让用户输入文件路径
    file_path = input("请输入qPCR下机数据文件(CSV格式)的路径: ").strip()
    
    # 移除路径两边可能带有的引号（直接拖拽文件进终端常带有引号）
    file_path = file_path.strip('"').strip("'")

    if not os.path.exists(file_path):
        print(f"错误：找不到文件 '{file_path}'，请检查路径是否正确！")
        return

    try:
        # 3. 读取CSV数据
        print("正在读取数据...")
        df = pd.read_csv(file_path)
        
        # 检查是否包含必需的列
        if 'Well' not in df.columns or 'Cq' not in df.columns:
            print("错误：数据中未找到 'Well' 或 'Cq' 列。请检查该CSV文件是否为标准的qPCR下机格式。")
            return

        # 4. 初始化 384孔板 DataFrame (16行 x 24列)
        rows = [chr(i) for i in range(ord('A'), ord('P') + 1)] # 行标签 A 到 P
        cols = list(range(1, 25)) # 列标签 1 到 24
        plate_df = pd.DataFrame(index=rows, columns=cols)

        # 5. 遍历数据并填入 384孔板 中
        for index, row in df.iterrows():
            well = str(row['Well']).strip()
            cq_value = row['Cq']
            
            # 跳过空孔或格式不正确的孔
            if pd.notna(well) and len(well) >= 2:
                row_letter = well[0].upper() # 提取字母 (A-P)
                try:
                    col_num = int(well[1:])  # 提取数字 (1-24)
                    
                    # 确保提取的行列在 384 孔板范围内
                    if row_letter in rows and col_num in cols:
                        plate_df.at[row_letter, col_num] = cq_value
                except ValueError:
                    continue # 忽略无法解析列数字的异常数据

        # 6. 输出结果到脚本同级目录
        output_filename = "Cq_384_plate_output.csv"
        output_path = os.path.join(script_dir, output_filename)
        
        plate_df.to_csv(output_path)
        print(f"处理成功！Cq值矩阵已保存至: {output_path}")

    except Exception as e:
        print(f"处理文件时发生错误: {e}")

if __name__ == "__main__":
    main()