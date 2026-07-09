#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import re
import sys
from datetime import datetime
from pathlib import Path


DEFAULT_FEATURES = ["gene", "mRNA", "transcript", "exon", "CDS"]


def script_dir():
    return Path(__file__).resolve().parent


def find_annotation_files(base_dir):
    files = []
    for pattern in ("*.gff3", "*.gff", "*.gtf"):
        files.extend(base_dir.rglob(pattern))
    return sorted(set(files), key=lambda p: str(p).lower())


def select_file(files, base_dir):
    print("\n--- 选择注释文件 ---")
    for index, path in enumerate(files, 1):
        try:
            display = path.relative_to(base_dir)
        except ValueError:
            display = path
        print(f"  [{index}] {display}")
    print("  [p] 手动输入文件路径")
    print("  [q] 退出")

    while True:
        choice = input(f"\n请选择文件 (1-{len(files)}/p/q): ").strip()
        if choice.lower() == "q":
            sys.exit(0)
        if choice.lower() == "p":
            manual = input("请输入 gff3/gff/gtf 文件完整路径: ").strip().strip('"').strip("'")
            path = Path(manual).expanduser()
            if path.exists() and path.is_file():
                return path
            print("文件不存在，请重新输入。")
            continue
        if choice.isdigit():
            index = int(choice)
            if 1 <= index <= len(files):
                return files[index - 1]
        print("输入无效，请重新选择。")


def parse_feature_selection():
    print("\n--- 选择保留的 feature 类型 ---")
    print("默认保留: gene, mRNA, transcript, exon, CDS")
    print("说明: transcript 用于兼容 GTF；如需严格四类，可输入 gene,mRNA,exon,CDS")

    raw = input("请输入要保留的 feature，逗号分隔；直接回车使用默认: ").strip()
    if not raw:
        features = DEFAULT_FEATURES
    else:
        features = [item.strip() for item in re.split(r"[,，;\s]+", raw) if item.strip()]

    if not features:
        print("未指定任何 feature，程序退出。")
        sys.exit(0)

    # GFF/GTF 第 3 列 feature 大小写不完全统一，这里做大小写不敏感匹配。
    normalized = {feature.lower() for feature in features}
    print("将保留: " + ", ".join(features))
    return normalized


def ask_keep_comments():
    print("\n--- 注释行处理 ---")
    print("  [1] 保留所有 # 开头的注释/头部行")
    print("  [2] 仅保留 ## 开头的 GFF 指令行")
    print("  [3] 不保留注释行")
    choice = input("请选择 (1/2/3) [默认 1]: ").strip()
    if choice not in {"1", "2", "3"}:
        choice = "1"
    return choice


def build_output_path(input_file):
    timestamp = datetime.now().strftime("%Y%m%d_%H%M")
    return input_file.with_name(f"{input_file.stem}_filtered_{timestamp}{input_file.suffix}")


def should_keep_comment(line, comment_mode):
    if comment_mode == "1":
        return True
    if comment_mode == "2":
        return line.startswith("##")
    return False


def filter_annotation(input_file, output_file, keep_features, comment_mode):
    total_records = 0
    kept_records = 0
    skipped_malformed = 0
    feature_counts = {}

    with input_file.open("r", encoding="utf-8", errors="replace") as fin, output_file.open("w", encoding="utf-8") as fout:
        for line in fin:
            if not line.strip():
                continue

            if line.startswith("#"):
                if should_keep_comment(line, comment_mode):
                    fout.write(line)
                continue

            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9:
                skipped_malformed += 1
                continue

            total_records += 1
            feature = parts[2].strip()
            if feature.lower() in keep_features:
                kept_records += 1
                feature_counts[feature] = feature_counts.get(feature, 0) + 1
                fout.write(line)

    return total_records, kept_records, skipped_malformed, feature_counts


def main():
    print("=" * 60)
    print(" GFF3/GFF/GTF 注释文件过滤工具")
    print("=" * 60)
    print("功能: 按 feature 类型过滤注释文件，默认只保留 gene/mRNA/exon/CDS 相关记录。")

    base_dir = script_dir()
    files = find_annotation_files(base_dir)

    if files:
        input_file = select_file(files, base_dir)
    else:
        print("\n脚本目录及子目录下未找到 .gff3/.gff/.gtf 文件。")
        manual = input("请输入注释文件完整路径，或直接回车退出: ").strip().strip('"').strip("'")
        if not manual:
            return
        input_file = Path(manual).expanduser()
        if not input_file.exists() or not input_file.is_file():
            print("文件不存在，程序退出。")
            return

    keep_features = parse_feature_selection()
    comment_mode = ask_keep_comments()
    output_file = build_output_path(input_file)

    print("\n正在过滤注释文件，请稍候...")
    total, kept, skipped, feature_counts = filter_annotation(input_file, output_file, keep_features, comment_mode)

    print("\n过滤完成。")
    print(f"输入文件: {input_file}")
    print(f"输出文件: {output_file}")
    print(f"注释记录总数: {total}")
    print(f"保留记录数: {kept}")
    if skipped:
        print(f"跳过格式异常行: {skipped}")
    if feature_counts:
        print("保留类型统计:")
        for feature, count in sorted(feature_counts.items(), key=lambda item: item[0].lower()):
            print(f"  {feature}: {count}")
    else:
        print("未保留任何记录，请检查 feature 名称是否与文件第 3 列一致。")


if __name__ == "__main__":
    try:
        main()
    except KeyboardInterrupt:
        print("\n用户中断，程序退出。")
