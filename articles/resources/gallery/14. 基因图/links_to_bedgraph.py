#!/usr/bin/env python3
"""Sum links scores over both anchor intervals to make a 1D bedGraph."""

from __future__ import annotations

import argparse
import math
import sys
from collections import defaultdict
from pathlib import Path


def read_events(path: Path):
    events = defaultdict(lambda: defaultdict(list))
    count = 0
    with path.open("r", encoding="utf-8-sig") as handle:
        for number, line in enumerate(handle, 1):
            if not line.strip() or line.startswith(("#", "track", "browser")):
                continue
            fields = line.split()
            if len(fields) < 7:
                raise ValueError(f"第 {number} 行必须包含两个区间和第7列互作分值")
            try:
                score = float(fields[6])
                anchors = (
                    (fields[0], int(fields[1]), int(fields[2])),
                    (fields[3], int(fields[4]), int(fields[5])),
                )
            except ValueError as exc:
                raise ValueError(f"第 {number} 行的坐标或分值无效") from exc
            if not math.isfinite(score) or score < 0:
                raise ValueError(f"第 {number} 行的互作分值必须是有限非负数")
            for chrom, start, end in anchors:
                if start < 0 or end <= start:
                    raise ValueError(f"第 {number} 行必须满足 0 <= start < end")
                events[chrom][start].append(score)
                events[chrom][end].append(-score)
            count += 1
    if not count:
        raise ValueError("输入文件没有有效互作记录")
    return events, count


def bedgraph_rows(events):
    """Sweep interval boundaries; include zero-valued gaps within each chromosome."""
    for chrom in sorted(events):
        positions = sorted(events[chrom])
        value = 0.0
        pending = None
        for index, start in enumerate(positions[:-1]):
            value = math.fsum((value, *events[chrom][start]))
            if abs(value) < 1e-12:
                value = 0.0
            end = positions[index + 1]
            formatted = f"{value:.10g}"
            if pending and pending[2] == start and pending[3] == formatted:
                pending = (chrom, pending[1], end, formatted)
            else:
                if pending:
                    yield pending
                pending = (chrom, start, end, formatted)
        if pending:
            yield pending


def main() -> int:
    parser = argparse.ArgumentParser(
        description="将 links 第7列分值累加到两个端点，生成一维互作强度 bedGraph。"
    )
    parser.add_argument("links", type=Path, help="输入 dp.links 或映射后的 .links 文件")
    parser.add_argument("-o", "--output", type=Path, help="默认输出为输入文件同名 .bedgraph")
    parser.add_argument("--force", action="store_true", help="覆盖已有输出文件")
    args = parser.parse_args()
    output = args.output or args.links.with_suffix(".bedgraph")
    try:
        if args.links.resolve() == output.resolve():
            raise ValueError("输出路径不能与输入文件相同")
        if output.exists() and not args.force:
            raise ValueError(f"输出文件已存在：{output}；使用 --force 覆盖")
        events, count = read_events(args.links)
        output.parent.mkdir(parents=True, exist_ok=True)
        rows = 0
        with output.open("w", encoding="utf-8", newline="\n") as handle:
            for row in bedgraph_rows(events):
                handle.write("\t".join(map(str, row)) + "\n")
                rows += 1
        print(f"已将 {count} 条互作转换为 {rows} 行 bedGraph：{output.resolve()}")
        print("坐标原样保留；强度为所有相连互作的分值之和。")
    except (OSError, ValueError) as exc:
        print(f"错误：{exc}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
