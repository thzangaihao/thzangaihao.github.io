import os
import sys
import glob
import re
import pandas as pd

# ============= 基础路径与交互逻辑 (保持一致) =============
def get_base_dir():
    if getattr(sys, 'frozen', False):
        return os.path.dirname(sys.executable)
    return os.path.dirname(os.path.abspath(__file__))

def find_files(exts, path=None):
    """支持搜索多种后缀，如 ['.xlsx', '.xls']"""
    if path is None:
        path = get_base_dir()
    all_files = []
    for ext in exts:
        search_pattern = os.path.join(path, '**', f'*{ext}')
        all_files.extend(glob.glob(search_pattern, recursive=True))
    return sorted(all_files)

def choose_one_file(files, desc="文件"):
    if not files:
        print(f"提示：在当前目录及子目录下未找到 {desc}")
        return None
    print(f"\n找到以下 {desc}:")
    for i, f in enumerate(files, 1):
        print(f"  [{i}] {os.path.basename(f)}")
    
    while True:
        try:
            val = input(f"\n请选择 {desc} 的编号: ").strip()
            idx = int(val) - 1
            if 0 <= idx < len(files):
                return files[idx]
            print("编号超出范围，请重新输入。")
        except ValueError:
            print("请输入数字编号。")

# ============= 核心逻辑：从 Excel 映射并替换 =============
def run_excel_replacement():
    print("--- Newick 标签替换工具 (Excel 映射版) ---")
    
    # 1. 选择 Excel 映射表
    excel_files = find_files([".xlsx", ".xls"])
    excel_path = choose_one_file(excel_files, "Excel 映射表")
    if not excel_path: return

    # 2. 读取 Excel
    try:
        # 默认读取第一个 Sheet
        df = pd.read_excel(excel_path)
        print(f"\n✅ 成功读取 Excel。表格列名如下：")
        cols = list(df.columns)
        for i, col in enumerate(cols, 1):
            print(f"  {i}. {col}")
        
        gca_idx = int(input("\n请输入包含 GCA 编号的列序号: ")) - 1
        name_idx = int(input("请输入要替换成的名称列序号 (如 NAME 列): ")) - 1
        
        gca_col = cols[gca_idx]
        name_col = cols[name_idx]
        
        # 建立映射字典，去掉 GCA 编号可能的版本号后缀差异，统一转为字符串
        mapping = {}
        for _, row in df.iterrows():
            key = str(row[gca_col]).strip()
            val = str(row[name_col]).strip()
            if key != 'nan' and val != 'nan':
                mapping[key] = val
        
        print(f"✅ 已载入 {len(mapping)} 组映射关系。")
    except Exception as e:
        print(f"❌ 读取 Excel 出错 (请确保已安装 openpyxl): {e}")
        return

    # 3. 选择 Newick 树文件
    nwk_files = find_files([".txt", ".nwk"])
    nwk_path = choose_one_file(nwk_files, "Newick 树文件")
    if not nwk_path: return

    with open(nwk_path, 'r', encoding='utf-8') as f:
        content = f.read()

    # 4. 识别并替换
    # 精准匹配：前面是 ( 或 , 后面是 : 的内容
    leaf_pattern = r'(?<=[,\(])([^,:\(\)]+)(?=:)'
    found_labels = re.findall(leaf_pattern, content)
    unique_labels = set(found_labels)

    new_content = content
    # 长度降序替换，防止子字符串误伤
    for label in sorted(unique_labels, key=len, reverse=True):
        clean_label = label.strip("'\"")
        
        # 提取标签中的 GCA 核心编号 (例如 GCA_001.1)
        gca_match = re.search(r'(GCA_\d+\.\d+)', clean_label)
        
        target_name = None
        if gca_match:
            gca_id = gca_match.group(1)
            # 尝试完全匹配或核心编号匹配
            target_name = mapping.get(gca_id) or mapping.get(clean_label)
        
        if target_name:
            # 如果拉丁名有空格，加单引号
            if ' ' in target_name:
                target_name = f"'{target_name}'"
            
            # 执行精准替换
            new_content = new_content.replace(f"({label}:", f"({target_name}:")
            new_content = new_content.replace(f",{label}:", f",{target_name}:")
        else:
            # 处理非 GCA 样品
            if not clean_label.startswith("GCA_") and not re.match(r'^\d+\.?\d*$', clean_label):
                print(f"\n发现未知样品: [ {clean_label} ]")
                u_input = input("请输入新名称 (回车跳过): ").strip()
                if u_input:
                    if ' ' in u_input: u_input = f"'{u_input}'"
                    new_content = new_content.replace(f"({label}:", f"({u_input}:")
                    new_content = new_content.replace(f",{label}:", f",{u_input}:")

    # 5. 保存
    output_path = nwk_path.rsplit('.', 1)[0] + "_LatinName.nwk"
    with open(output_path, 'w', encoding='utf-8') as f:
        f.write(new_content)
    
    print(f"\n✨ 完成！处理后的进化树已保存至：\n{output_path}")

if __name__ == "__main__":
    run_excel_replacement()