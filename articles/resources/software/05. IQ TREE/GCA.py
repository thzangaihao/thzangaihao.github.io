import os
import sys
import glob
import re

# ============= 基础路径与搜索 (保持一致) =============
def get_base_dir():
    if getattr(sys, 'frozen', False):
        return os.path.dirname(sys.executable)
    return os.path.dirname(os.path.abspath(__file__))

def find_files(ext, path=None):
    if not ext.startswith('.'):
        ext = '.' + ext
    if path is None:
        path = get_base_dir()
    search_pattern = os.path.join(path, '**', f'*{ext}')
    return sorted(glob.glob(search_pattern, recursive=True))

def choose_file(files, desc="文件"):
    if not files:
        print(f"提示：未找到 {desc}")
        return []
    print(f"\n找到 {len(files)} 个 {desc}:")
    for i, f in enumerate(files, 1):
        print(f"  [{i}] {os.path.basename(f)}")
    
    user_input = input(f"\n请输入欲操作的编号 (或输入 'all'): ").strip().lower()
    selected_indices = []
    if user_input in ['all', 'a']:
        selected_indices = range(len(files))
    else:
        try:
            selected_indices = [int(i.strip())-1 for i in user_input.split(',')]
        except:
            return []
    return [files[i] for i in selected_indices if 0 <= i < len(files)]

# ============= 核心逻辑改进：精准定位叶节点 =============
def process_newick(file_path):
    with open(file_path, 'r', encoding='utf-8') as f:
        content = f.read()

    # 改进的正则：匹配前面是 ( 或 , 的字符，直到碰到 :
    # 这样就能完美避开 )0.953: 这种支持率格式
    leaf_pattern = r'(?<=[,\(])([^,:\(\)]+)(?=:)'
    labels = re.findall(leaf_pattern, content)
    unique_labels = sorted(list(set(labels)))

    mapping = {}
    print(f"\n--- 正在处理文件: {os.path.basename(file_path)} ---")
    
    for label in unique_labels:
        # 清洗标签（去掉两端引号，如果有的话）
        clean_label = label.strip("'\"")
        
        # 模式1：NCBI GCA 编号
        if "GCA_" in clean_label:
            match = re.search(r'(GCA_\d+\.\d+)', clean_label)
            if match:
                mapping[label] = match.group(1)
            else:
                mapping[label] = clean_label
        
        # 模式2：自定义样品 (如 CK-5-4)
        else:
            print(f"\n发现自定义样品: [ {clean_label} ]")
            choice = input(f"请输入新名称 (直接回车保持原样, 输入 's' 跳过所有非GCA): ").strip()
            if choice == 's':
                mapping[label] = clean_label
            elif choice:
                mapping[label] = choice
            else:
                mapping[label] = clean_label

    # 执行精准替换
    new_content = content
    # 按长度降序排列，防止短字符串替换了长字符串的一部分
    for old_name in sorted(mapping.keys(), key=len, reverse=True):
        new_name = mapping[old_name]
        # 确保只替换作为标签出现的字符串
        # 我们寻找 (old_name: 或 ,old_name:
        new_content = new_content.replace(f"({old_name}:", f"({new_name}:")
        new_content = new_content.replace(f",{old_name}:", f",{new_name}:")

    # 保存文件
    new_file_path = file_path.rsplit('.', 1)[0] + "_clean.nwk"
    with open(new_file_path, 'w', encoding='utf-8') as f:
        f.write(new_content)
    print(f"\n✅ 成功！节点支持率已保留，叶节点已更名。")
    print(f"结果路径: {new_file_path}")

def main():
    print("--- Newick 标签清洗工具 v2.0 (已修正支持率识别错误) ---")
    files = find_files("txt") + find_files("nwk") # 自动搜寻 txt 和 nwk
    chosen_files = choose_file(files, "进化树文件")

    if chosen_files:
        for f in chosen_files:
            process_newick(f)

if __name__ == "__main__":
    main()