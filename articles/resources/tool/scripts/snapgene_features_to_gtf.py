#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""将 SnapGene 导出的特征 CSV 转换为 IGV 可读取的 GTF 文件。

SnapGene CSV 的常见无表头格式：
    特征名称,起止位置,长度,方向,特征类型
    EGFP,2200..2919,720,=>,CDS

坐标按 SnapGene 与 GTF 共有的 1-based 闭区间原样转换；=> 为正链，<= 为负链。
"""

import argparse
import csv
import re
import sys
from collections import defaultdict
from pathlib import Path


HEADER_WORDS = {
    "name", "feature", "feature name", "名称", "特征", "特征名称",
    "range", "location", "位置", "坐标", "起止位置",
}


def parse_args():
    parser = argparse.ArgumentParser(
        description="将 SnapGene 导出的特征 CSV 转换为 GTF。"
    )
    parser.add_argument("input_csv", type=Path, help="SnapGene 导出的 CSV 文件")
    parser.add_argument(
        "output_gtf", nargs="?", type=Path,
        help="输出 GTF；默认与输入文件同名、后缀改为 .gtf",
    )
    parser.add_argument(
        "--seqname", default="sequence",
        help="GTF 第一列的参考序列/染色体名称，必须与 IGV 中参考序列名一致（默认: sequence）",
    )
    parser.add_argument(
        "--source", default="SnapGene", help="GTF source 字段（默认: SnapGene）"
    )
    return parser.parse_args()


def open_csv(path):
    """尝试常见中文 CSV 编码，并返回文本句柄和实际编码。"""
    for encoding in ("utf-8-sig", "gb18030"):
        try:
            handle = path.open("r", encoding=encoding, newline="")
            handle.read(4096)
            handle.seek(0)
            return handle, encoding
        except UnicodeDecodeError:
            handle.close()
    raise UnicodeError("文件既不是 UTF-8，也无法按 GB18030 解码")


def is_header(row):
    if not row:
        return False
    normalized = {cell.strip().lower() for cell in row}
    return bool(normalized & HEADER_WORDS)


def parse_location(text):
    """解析 33..49、33-49 或单点 33，返回 1-based 闭区间。"""
    cleaned = text.strip().replace(",", "")
    cleaned = cleaned.replace("<", "").replace(">", "")
    match = re.fullmatch(r"(\d+)\s*(?:\.\.|-)\s*(\d+)", cleaned)
    if match:
        start, end = map(int, match.groups())
    elif cleaned.isdigit():
        start = end = int(cleaned)
    else:
        raise ValueError(f"无法识别坐标 {text!r}")
    if start < 1 or end < 1:
        raise ValueError("GTF 坐标必须大于或等于 1")
    if start > end:
        raise ValueError(
            f"起点 {start} 大于终点 {end}；环状序列跨原点特征需先拆分为两行"
        )
    return start, end


def normalize_strand(text):
    value = text.strip()
    if value in {"=>", "+", "forward", "Forward", "正向"}:
        return "+"
    if value in {"<=", "-", "reverse", "Reverse", "反向"}:
        return "-"
    return "."


def normalize_feature_type(text):
    value = text.strip() or "feature"
    # GTF 第三列不能含空白；保留 SnapGene 原始类型的语义。
    return re.sub(r"\s+", "_", value)


def escape_attribute(value):
    return value.replace("\\", "\\\\").replace('"', '\\"').replace("\n", " ")


def safe_id(name):
    value = re.sub(r"[^A-Za-z0-9_.:-]+", "_", name.strip()).strip("_")
    return value or "feature"


def convert(input_csv, output_gtf, seqname, source):
    if not input_csv.is_file():
        raise FileNotFoundError(f"找不到输入文件: {input_csv}")
    if not seqname.strip() or any(char.isspace() for char in seqname):
        raise ValueError("--seqname 不能为空或包含空白字符")

    records = []
    warnings = []
    id_counts = defaultdict(int)

    handle, encoding = open_csv(input_csv)
    with handle:
        reader = csv.reader(handle, skipinitialspace=True)
        for line_no, row in enumerate(reader, 1):
            if not row or not any(cell.strip() for cell in row):
                continue
            if line_no == 1 and is_header(row):
                continue
            if len(row) < 5:
                warnings.append(f"第 {line_no} 行不足 5 列，已跳过")
                continue

            name, location, exported_length, direction, feature_type = (
                cell.strip() for cell in row[:5]
            )
            try:
                start, end = parse_location(location)
            except ValueError as error:
                warnings.append(f"第 {line_no} 行已跳过: {error}")
                continue

            actual_length = end - start + 1
            if exported_length.isdigit() and int(exported_length) != actual_length:
                warnings.append(
                    f"第 {line_no} 行长度不一致: CSV={exported_length}, 坐标计算={actual_length}"
                )

            feature_id_base = safe_id(name)
            id_counts[feature_id_base] += 1
            feature_id = feature_id_base
            if id_counts[feature_id_base] > 1:
                feature_id = f"{feature_id_base}_{id_counts[feature_id_base]}"

            attributes = (
                f'gene_id "{escape_attribute(feature_id)}"; '
                f'transcript_id "{escape_attribute(feature_id)}"; '
                f'feature_name "{escape_attribute(name or feature_id)}"; '
                f'snapgene_type "{escape_attribute(feature_type or "feature")}";'
            )
            fields = [
                seqname, source, normalize_feature_type(feature_type),
                str(start), str(end), ".", normalize_strand(direction), ".", attributes,
            ]
            records.append("\t".join(fields))

    if not records:
        raise ValueError("没有可转换的有效特征记录")

    output_gtf.parent.mkdir(parents=True, exist_ok=True)
    with output_gtf.open("w", encoding="utf-8", newline="\n") as handle:
        handle.write("# SnapGene features converted to GTF\n")
        handle.write(f"# source_csv: {input_csv.name}\n")
        for record in records:
            handle.write(record + "\n")

    return len(records), warnings, encoding


def main():
    args = parse_args()
    output_gtf = args.output_gtf or args.input_csv.with_suffix(".gtf")
    try:
        count, warnings, encoding = convert(
            args.input_csv, output_gtf, args.seqname, args.source
        )
    except (OSError, UnicodeError, ValueError) as error:
        print(f"[错误] 转换失败: {error}", file=sys.stderr)
        return 1

    print(f"[完成] 共写入 {count} 个特征")
    print(f"[信息] 输入编码: {encoding}")
    print(f"[输出] {output_gtf.resolve()}")
    for warning in warnings:
        print(f"[警告] {warning}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
