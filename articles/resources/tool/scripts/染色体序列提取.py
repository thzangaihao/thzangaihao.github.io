import re
import sys
from pathlib import Path


def select_from_list(items, prompt, is_file=False, base_dir=None):
    """交互式命令行菜单，返回单个条目。"""
    print(f"\n{prompt}")
    for i, item in enumerate(items):
        if is_file and base_dir:
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
            print("⚠️ 数字超出范围，请重新输入。")
        except ValueError:
            print("⚠️ 无效输入，请输入数字。")


def select_multiple_from_list(items, prompt):
    """支持 all、逗号分隔、范围连接的多选菜单。"""
    print(f"\n{prompt}")
    for i, item in enumerate(items):
        print(f"[{i + 1}] {item}")

    print("\n选择方式示例：")
    print("  1,3,5        提取第 1、3、5 条")
    print("  2-6          提取第 2 到 6 条")
    print("  1,4-7,chr10  混合使用数字、范围和染色体 ID")
    print("  all          提取全部")

    while True:
        raw = input("请输入要提取的染色体/Scaffold (输入 0 退出): ").strip()
        if raw == "0":
            print("已退出程序。")
            sys.exit(0)

        try:
            selected = parse_chromosome_selection(raw, items)
        except ValueError as exc:
            print(f"⚠️ {exc}")
            continue

        if selected:
            print(f"\n✅ 已选择 {len(selected)} 条：{', '.join(selected)}")
            return selected
        print("⚠️ 未选择任何条目，请重新输入。")


def parse_chromosome_selection(raw, items):
    """解析 all、逗号列表、数字范围和染色体 ID。"""
    if not raw:
        raise ValueError("输入为空，请重新输入。")

    if raw.strip().lower() == "all":
        return list(items)

    selected = []
    seen = set()
    tokens = [part.strip() for part in raw.split(",") if part.strip()]

    for token in tokens:
        matches = resolve_selection_token(token, items)
        for item in matches:
            if item not in seen:
                selected.append(item)
                seen.add(item)

    return selected


def resolve_selection_token(token, items):
    """解析单个选择片段：单项或范围。"""
    exact = find_item_by_name(token, items)
    if exact is not None:
        return [exact]

    if "-" in token:
        left, right = [part.strip() for part in token.split("-", 1)]
        if not left or not right:
            raise ValueError(f"范围 '{token}' 不完整，请使用类似 2-6 的格式。")

        start = resolve_item_index(left, items)
        end = resolve_item_index(right, items)
        if start is None or end is None:
            raise ValueError(f"无法识别范围 '{token}'，请确认两端是序号或染色体 ID。")
        if start > end:
            start, end = end, start
        return items[start:end + 1]

    index = resolve_item_index(token, items)
    if index is None:
        raise ValueError(f"无法识别 '{token}'，请使用列表序号、染色体 ID、范围或 all。")
    return [items[index]]


def resolve_item_index(value, items):
    """将序号或染色体 ID 解析为 0-based index。"""
    if value.isdigit():
        index = int(value) - 1
        if 0 <= index < len(items):
            return index
        raise ValueError(f"序号 {value} 超出范围 1-{len(items)}。")

    item = find_item_by_name(value, items)
    if item is not None:
        return items.index(item)
    return None


def find_item_by_name(value, items):
    """优先精确匹配染色体 ID，同时兼容大小写不敏感匹配。"""
    if value in items:
        return value

    lowered = value.lower()
    for item in items:
        if item.lower() == lowered:
            return item
    return None


def select_output_mode(selected_count):
    """选择合并输出或分别输出。"""
    if selected_count <= 1:
        return "separate"

    print("\n📦 多条染色体输出方式：")
    print("[1] 合并输出到一个文件")
    print("[2] 分别输出为多个文件")

    while True:
        choice = input("请选择输出方式 (1-2，默认 1): ").strip() or "1"
        if choice == "1":
            return "merge"
        if choice == "2":
            return "separate"
        print("⚠️ 无效输入，请输入 1 或 2。")


def get_chromosomes_from_fasta(filepath):
    """从 FASTA 文件中快速提取所有染色体 ID。"""
    chroms = []
    with open(filepath, "r", encoding="utf-8") as f:
        for line in f:
            if line.startswith(">"):
                chrom_id = line[1:].strip().split()[0]
                chroms.append(chrom_id)
    return chroms


def get_chromosomes_from_gff(filepath):
    """从 GFF3 文件中提取所有染色体/Scaffold ID。"""
    chroms = []
    seen = set()
    with open(filepath, "r", encoding="utf-8") as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) > 1:
                chrom_id = parts[0]
                if chrom_id not in seen:
                    seen.add(chrom_id)
                    chroms.append(chrom_id)
    return chroms


def safe_filename(text):
    """生成适合 Windows/Linux 文件名的安全字符串。"""
    return re.sub(r'[\\/:*?"<>|\s]+', "_", text).strip("_") or "selected"


def build_merged_output_path(input_file, selected_chroms, output_dir, extension):
    """构建合并输出文件名。"""
    if len(selected_chroms) > 5:
        label = "all"
    else:
        label = "_".join(safe_filename(chrom) for chrom in selected_chroms)
    return output_dir / f"{input_file.stem}_{label}{extension}"


def extract_fasta(input_file, target_chroms, output_dir, output_mode):
    """逐行提取一个或多个染色体的 FASTA 序列。"""
    target_set = set(target_chroms)
    found = set()

    if output_mode == "merge":
        output_file = build_merged_output_path(input_file, target_chroms, output_dir, ".fasta")
        keep = False
        with open(input_file, "r", encoding="utf-8") as fin, open(output_file, "w", encoding="utf-8") as fout:
            for line in fin:
                if line.startswith(">"):
                    current_chrom = line[1:].strip().split()[0]
                    keep = current_chrom in target_set
                    if keep:
                        found.add(current_chrom)
                        fout.write(line)
                elif keep:
                    fout.write(line)
        return ([output_file] if found else []), found

    handles = {}
    output_files = {}
    keep_handle = None
    try:
        with open(input_file, "r", encoding="utf-8") as fin:
            for line in fin:
                if line.startswith(">"):
                    current_chrom = line[1:].strip().split()[0]
                    if current_chrom in target_set:
                        found.add(current_chrom)
                        if current_chrom not in handles:
                            output_file = output_dir / f"{safe_filename(current_chrom)}.fasta"
                            output_files[current_chrom] = output_file
                            handles[current_chrom] = open(output_file, "w", encoding="utf-8")
                        keep_handle = handles[current_chrom]
                        keep_handle.write(line)
                    else:
                        keep_handle = None
                elif keep_handle:
                    keep_handle.write(line)
    finally:
        for handle in handles.values():
            handle.close()

    return [output_files[chrom] for chrom in target_chroms if chrom in found], found


def extract_gff(input_file, target_chroms, output_dir, output_mode):
    """逐行提取一个或多个染色体的 GFF3 注释信息。"""
    target_set = set(target_chroms)
    found = set()

    if output_mode == "merge":
        output_file = build_merged_output_path(input_file, target_chroms, output_dir, ".gff3")
        with open(input_file, "r", encoding="utf-8") as fin, open(output_file, "w", encoding="utf-8") as fout:
            wrote_version = False
            for line in fin:
                if line.startswith("##gff-version"):
                    if not wrote_version:
                        fout.write(line)
                        wrote_version = True
                elif line.startswith("#"):
                    continue
                else:
                    parts = line.split("\t")
                    if len(parts) > 1 and parts[0] in target_set:
                        found.add(parts[0])
                        fout.write(line)
        return ([output_file] if found else []), found

    handles = {}
    output_files = {}
    try:
        with open(input_file, "r", encoding="utf-8") as fin:
            for line in fin:
                if line.startswith("#"):
                    continue

                parts = line.split("\t")
                if len(parts) <= 1 or parts[0] not in target_set:
                    continue

                chrom = parts[0]
                found.add(chrom)
                if chrom not in handles:
                    output_file = output_dir / f"{safe_filename(chrom)}.gff3"
                    output_files[chrom] = output_file
                    handles[chrom] = open(output_file, "w", encoding="utf-8")
                    handles[chrom].write("##gff-version 3\n")
                handles[chrom].write(line)
    finally:
        for handle in handles.values():
            handle.close()

    return [output_files[chrom] for chrom in target_chroms if chrom in found], found


def main():
    print("=" * 55)
    print(" 🧬 染色体序列/注释提取工具 (FASTA/GFF3)")
    print("=" * 55)

    script_dir = Path(__file__).resolve().parent
    print(f"📂 脚本运行及输出目录固定为: {script_dir}\n")

    files = []
    for ext in ["*.fasta", "*.fa", "*.fna", "*.gff", "*.gff3"]:
        files.extend(script_dir.rglob(ext))

    if not files:
        print("❌ 在脚本所在的目录及子文件夹下未找到任何 FASTA 或 GFF3 文件。")
        return

    selected_file = select_from_list(files, "🔍 找到以下文件，请选择你要处理的文件：", is_file=True, base_dir=script_dir)
    file_ext = selected_file.suffix.lower()
    is_fasta = file_ext in [".fasta", ".fa", ".fna"]
    is_gff = file_ext in [".gff", ".gff3"]

    print("\n⏳ 正在扫描文件中的染色体信息，大文件请稍候...")
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

    selected_chroms = select_multiple_from_list(chroms, "🎯 扫描到以下染色体/Scaffold，请选择你要提取的目标：")
    output_mode = select_output_mode(len(selected_chroms))
    mode_text = "合并输出" if output_mode == "merge" else "分别输出"

    print(f"\n⏳ 正在提取 {len(selected_chroms)} 条染色体/Scaffold 的数据（{mode_text}）...")
    if is_fasta:
        out_files, found = extract_fasta(selected_file, selected_chroms, script_dir, output_mode)
    elif is_gff:
        out_files, found = extract_gff(selected_file, selected_chroms, script_dir, output_mode)

    if out_files:
        print("✅ 提取完成！文件已输出到脚本同级目录：")
        for out_file in out_files:
            print(f"📁 {out_file}")

        missing = [chrom for chrom in selected_chroms if chrom not in found]
        if missing:
            print(f"⚠️ 以下条目未找到匹配内容：{', '.join(missing)}")
    else:
        print("❌ 提取失败，未找到匹配的数据内容。")


if __name__ == "__main__":
    main()
