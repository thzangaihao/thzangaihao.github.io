import sys
from pathlib import Path

def select_from_list(items, prompt, is_file=False, base_dir=None):
    """交互式命令行菜单"""
    print(f"\n{prompt}")
    for i, item in enumerate(items):
        if is_file and base_dir:
            # 如果是文件，显示相对于脚本所在目录的路径，方便区分不同子文件夹中的同名文件
            try:
                display_name = item.relative_to(base_dir)
            except ValueError:
                display_name = item.name
        else:
            display_name = item
        print(f"[{i + 1}] {display_name}")
    
    while True:
        try:
            choice = int(input(f"\n请输入对应数字进行选择 (1-{len(items)}，输入 0 退出): "))
            if choice == 0:
                print("已退出程序。")
                sys.exit(0)
            if 1 <= choice <= len(items):
                return items[choice - 1]
            else:
                print("⚠️ 数字超出范围，请重新输入。")
        except ValueError:
            print("⚠️ 无效输入，请输入数字。")

def get_chromosomes_from_fasta(filepath):
    """从 FASTA 文件中快速提取所有染色体 ID"""
    chroms = []
    with open(filepath, 'r', encoding='utf-8') as f:
        for line in f:
            if line.startswith('>'):
                chrom_id = line[1:].strip().split()[0]
                chroms.append(chrom_id)
    return chroms

def get_chromosomes_from_gff(filepath):
    """从 GFF3 文件中提取所有染色体/Scaffold ID"""
    chroms = []
    seen = set()
    with open(filepath, 'r', encoding='utf-8') as f:
        for line in f:
            if line.startswith('#'):
                continue
            parts = line.split('\t')
            if len(parts) > 1:
                chrom_id = parts[0]
                if chrom_id not in seen:
                    seen.add(chrom_id)
                    chroms.append(chrom_id)
    return chroms

def extract_fasta(input_file, target_chrom, output_dir):
    """逐行提取目标染色体的 FASTA 序列"""
    # 严格将输出路径指定为 output_dir (即脚本同级目录)
    output_file = output_dir / f"{target_chrom}.fasta"
    keep = False
    found = False
    
    with open(input_file, 'r', encoding='utf-8') as fin, \
         open(output_file, 'w', encoding='utf-8') as fout:
        for line in fin:
            if line.startswith('>'):
                current_chrom = line[1:].strip().split()[0]
                if current_chrom == target_chrom:
                    keep = True
                    found = True
                    fout.write(line)
                else:
                    keep = False
            elif keep:
                fout.write(line)
                
    return output_file if found else None

def extract_gff(input_file, target_chrom, output_dir):
    """逐行提取目标染色体的 GFF3 注释信息"""
    # 严格将输出路径指定为 output_dir (即脚本同级目录)
    output_file = output_dir / f"{target_chrom}.gff3"
    found = False
    
    with open(input_file, 'r', encoding='utf-8') as fin, \
         open(output_file, 'w', encoding='utf-8') as fout:
        for line in fin:
            if line.startswith('##gff-version'):
                fout.write(line)
            elif line.startswith('#'):
                continue
            else:
                parts = line.split('\t')
                if len(parts) > 1 and parts[0] == target_chrom:
                    found = True
                    fout.write(line)
                    
    return output_file if found else None

def main():
    print("="*55)
    print(" 🧬 染色体序列/注释提取工具 (FASTA/GFF3)")
    print("="*55)

    # 获取脚本自身的绝对路径所在目录
    script_dir = Path(__file__).resolve().parent
    print(f"📂 脚本运行及输出目录固定为: {script_dir}\n")

    # 1. 搜寻脚本同级及子级文件夹下的目标文件
    files = []
    for ext in ['*.fasta', '*.fa', '*.fna', '*.gff', '*.gff3']:
        files.extend(script_dir.rglob(ext))

    if not files:
        print("❌ 在脚本所在的目录及子文件夹下未找到任何 FASTA 或 GFF3 文件。")
        return

    # 2. 交互式选择文件
    selected_file = select_from_list(files, "🔍 找到以下文件，请选择你要处理的文件：", is_file=True, base_dir=script_dir)
    file_ext = selected_file.suffix.lower()
    is_fasta = file_ext in ['.fasta', '.fa', '.fna']
    is_gff = file_ext in ['.gff', '.gff3']

    # 3. 解析文件中的染色体列表
    print(f"\n⏳ 正在扫描文件中的染色体信息，大文件请稍候...")
    if is_fasta:
        chroms = get_chromosomes_from_fasta(selected_file)
    elif is_gff:
        chroms = get_chromosomes_from_gff(selected_file)
    else:
        print("❌ 不支持的文件格式。")
        return

    if not chroms:
        print("❌ 未在文件中解析到有效的染色体或序列 ID。")
        return

    # 4. 交互式选择染色体
    selected_chrom = select_from_list(chroms, "🎯 扫描到以下染色体/Scaffold，请选择你要提取的目标：")

    # 5. 提取并输出 (强制输出到 script_dir)
    print(f"\n⏳ 正在提取染色体 '{selected_chrom}' 的数据...")
    if is_fasta:
        out_file = extract_fasta(selected_file, selected_chrom, script_dir)
    elif is_gff:
        out_file = extract_gff(selected_file, selected_chrom, script_dir)

    if out_file:
        print(f"✅ 提取完成！文件已精确输出到脚本同级目录: \n📁 {out_file}")
    else:
        print("❌ 提取失败，未找到匹配的数据内容。")

if __name__ == '__main__':
    main()