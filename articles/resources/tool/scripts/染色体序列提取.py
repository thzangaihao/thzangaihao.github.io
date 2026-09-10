#!/usr/bin/env python3
"""从 FASTA/GFF3 提取指定序列；支持每行一个 ID 的批量清单。"""

import argparse
import re
import sys
from pathlib import Path

FASTA_SUFFIXES = {".fasta", ".fa", ".fna"}
GFF_SUFFIXES = {".gff", ".gff3"}
LIST_SUFFIXES = {".txt", ".list", ".tsv", ".csv"}


def unique(items):
    """按首次出现顺序去重。"""
    return list(dict.fromkeys(items))


def select_from_list(items, prompt, base_dir=None):
    print(f"\n{prompt}")
    for number, item in enumerate(items, 1):
        display = item
        if base_dir and isinstance(item, Path):
            try:
                display = item.relative_to(base_dir)
            except ValueError:
                pass
        print(f"[{number}] {display}")
    while True:
        value = input(f"请选择 (1-{len(items)}，0 退出): ").strip()
        if value == "0":
            raise SystemExit(0)
        if value.isdigit() and 1 <= int(value) <= len(items):
            return items[int(value) - 1]
        print("输入无效，请重新输入。")


def read_id_list(filepath):
    """读取清单首列；空行和 # 注释行会忽略，重复 ID 自动去重。"""
    ids = []
    with open(filepath, "r", encoding="utf-8-sig") as handle:
        for line in handle:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            # 兼容每行一个 ID，以及以空白、制表符或逗号分隔的表格首列。
            sequence_id = re.split(r"[\s,]+", line, maxsplit=1)[0]
            if sequence_id:
                ids.append(sequence_id)
    return unique(ids)


def get_ids(filepath, is_fasta):
    ids = []
    with open(filepath, "r", encoding="utf-8-sig") as handle:
        for line in handle:
            if is_fasta and line.startswith(">"):
                ids.append(line[1:].strip().split()[0])
            elif not is_fasta and not line.startswith("#"):
                fields = line.rstrip("\n").split("\t")
                if len(fields) >= 2:
                    ids.append(fields[0])
    return unique(ids)


def resolve_manual_token(token, available):
    """解析单个 ID、菜单序号或菜单范围。"""
    by_lower = {item.lower(): item for item in available}
    if token.lower() in by_lower:
        return [by_lower[token.lower()]]
    if "-" in token:
        left, right = token.split("-", 1)
        if left.isdigit() and right.isdigit():
            start, end = int(left), int(right)
            if 1 <= start <= len(available) and 1 <= end <= len(available):
                start, end = sorted((start, end))
                return available[start - 1:end]
    if token.isdigit() and 1 <= int(token) <= len(available):
        return [available[int(token) - 1]]
    raise ValueError(f"无法识别“{token}”")


def choose_targets(available, script_dir):
    print("\n目标 ID 输入方式：\n[1] 屏幕列表手工选择\n[2] 批量清单（每行一个 ID）")
    while True:
        choice = input("请选择 (1-2，默认 2): ").strip() or "2"
        if choice in {"1", "2"}:
            break
        print("请输入 1 或 2。")
    if choice == "2":
        files = sorted(p for p in script_dir.rglob("*") if p.is_file() and p.suffix.lower() in LIST_SUFFIXES)
        if files:
            filepath = select_from_list(files, "请选择 ID 清单：", script_dir)
        else:
            filepath = Path(input("清单路径: ").strip().strip('"')).expanduser()
        ids = read_id_list(filepath)
        if not ids:
            raise ValueError(f"清单为空：{filepath}")
        return ids

    print("\n可用 ID：")
    for number, item in enumerate(available, 1):
        print(f"[{number}] {item}")
    print("示例：1,3,5 / 2-6 / query001,query002 / all")
    while True:
        raw = input("请输入目标: ").strip()
        if raw.lower() == "all":
            return list(available)
        try:
            selected = []
            for token in (part.strip() for part in raw.split(",")):
                if token:
                    selected.extend(resolve_manual_token(token, available))
            if selected:
                return unique(selected)
            raise ValueError("输入不能为空")
        except ValueError as exc:
            print(f"输入有误：{exc}")


def safe_filename(text):
    return re.sub(r'[\\/:*?"<>|\s]+', "_", text).strip("_") or "selected"


def extract_fasta(input_file, targets, output_dir, mode):
    target_set, found = set(targets), set()
    handles, paths = {}, {}
    merged = output_dir / f"{input_file.stem}_batch.fasta"
    merged_handle = open(merged, "w", encoding="utf-8") if mode == "merge" else None
    try:
        with open(input_file, "r", encoding="utf-8-sig") as source:
            current_handle = None
            for line in source:
                if line.startswith(">"):
                    current = line[1:].strip().split()[0]
                    current_handle = None
                    if current in target_set:
                        found.add(current)
                        if merged_handle:
                            current_handle = merged_handle
                        else:
                            path = output_dir / f"{safe_filename(current)}.fasta"
                            paths[current] = path
                            if current not in handles:
                                handles[current] = open(path, "w", encoding="utf-8")
                            current_handle = handles[current]
                if current_handle:
                    current_handle.write(line)
    finally:
        if merged_handle:
            merged_handle.close()
        for handle in handles.values():
            handle.close()
    if mode == "merge":
        if not found:
            merged.unlink(missing_ok=True)
            return [], found
        return [merged], found
    return [paths[item] for item in targets if item in found], found


def extract_gff(input_file, targets, output_dir, mode):
    target_set, found = set(targets), set()
    handles, paths = {}, {}
    merged = output_dir / f"{input_file.stem}_batch.gff3"
    merged_handle = open(merged, "w", encoding="utf-8") if mode == "merge" else None
    if merged_handle:
        merged_handle.write("##gff-version 3\n")
    try:
        with open(input_file, "r", encoding="utf-8-sig") as source:
            for line in source:
                if line.startswith("#"):
                    continue
                fields = line.rstrip("\n").split("\t")
                if len(fields) < 2 or fields[0] not in target_set:
                    continue
                current = fields[0]
                found.add(current)
                if merged_handle:
                    merged_handle.write(line)
                else:
                    if current not in handles:
                        path = output_dir / f"{safe_filename(current)}.gff3"
                        paths[current] = path
                        handles[current] = open(path, "w", encoding="utf-8")
                        handles[current].write("##gff-version 3\n")
                    handles[current].write(line)
    finally:
        if merged_handle:
            merged_handle.close()
        for handle in handles.values():
            handle.close()
    if mode == "merge":
        if not found:
            merged.unlink(missing_ok=True)
            return [], found
        return [merged], found
    return [paths[item] for item in targets if item in found], found


def parse_args(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("-i", "--input", type=Path, help="输入 FASTA/GFF3")
    parser.add_argument("-l", "--list", dest="id_list", type=Path, help="批量 ID 清单")
    parser.add_argument("-o", "--output-dir", type=Path, help="输出目录")
    parser.add_argument("--mode", choices=("merge", "separate"), default="merge")
    return parser.parse_args(argv)


def main(argv=None):
    args = parse_args(argv)
    script_dir = Path(__file__).resolve().parent
    input_file = args.input
    interactive = input_file is None
    if interactive:
        files = sorted(p for p in script_dir.rglob("*") if p.is_file() and p.suffix.lower() in FASTA_SUFFIXES | GFF_SUFFIXES)
        if not files:
            print("未找到 FASTA/GFF3 文件。", file=sys.stderr)
            return 1
        input_file = select_from_list(files, "请选择要处理的文件：", script_dir)
    input_file = input_file.resolve()
    if not input_file.is_file() or input_file.suffix.lower() not in FASTA_SUFFIXES | GFF_SUFFIXES:
        print(f"输入文件不存在或格式不支持：{input_file}", file=sys.stderr)
        return 1

    is_fasta = input_file.suffix.lower() in FASTA_SUFFIXES
    available = get_ids(input_file, is_fasta)
    try:
        targets = read_id_list(args.id_list) if args.id_list else choose_targets(available, script_dir)
    except (OSError, ValueError) as exc:
        print(f"读取 ID 清单失败：{exc}", file=sys.stderr)
        return 1
    output_dir = (args.output_dir or script_dir).resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    mode = args.mode
    if interactive and len(targets) > 1:
        mode = "separate" if input("合并输出到一个文件？(Y/n): ").strip().lower() in {"n", "no"} else "merge"

    extractor = extract_fasta if is_fasta else extract_gff
    outputs, found = extractor(input_file, targets, output_dir, mode)
    missing = [item for item in targets if item not in found]
    print(f"\n完成：清单 {len(targets)} 个唯一 ID，找到 {len(found)} 个，缺失 {len(missing)} 个。")
    for path in outputs:
        print(f"输出：{path}")
    if missing:
        path = output_dir / f"{input_file.stem}_missing_ids.txt"
        path.write_text("\n".join(missing) + "\n", encoding="utf-8")
        print(f"缺失 ID 清单：{path}")
    return 0 if found else 2


if __name__ == "__main__":
    raise SystemExit(main())
