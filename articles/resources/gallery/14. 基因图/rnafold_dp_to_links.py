#!/usr/bin/env python3
"""Convert an RNAfold dot-plot EPS to transcript-relative links."""

from __future__ import annotations

import argparse
import math
import re
import sys
from pathlib import Path


UBOX_PATTERN = re.compile(
    r"^\s*(\d+)\s+(\d+)\s+"
    r"([+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[Ee][+-]?\d+)?)\s+ubox\b"
)


def parse_dot_plot(
    path: Path, minimum_probability: float = 0.0
) -> list[tuple[int, int, float]]:
    """Read 1-based RNA positions and restore p(i,j) from sqrt(p)."""
    pairs: dict[tuple[int, int], float] = {}
    with path.open("r", encoding="latin-1") as handle:
        for line in handle:
            match = UBOX_PATTERN.match(line)
            if not match:
                continue

            left, right = int(match.group(1)), int(match.group(2))
            sqrt_probability = float(match.group(3))
            probability = sqrt_probability * sqrt_probability
            if left < 1 or right < 1 or left == right or not math.isfinite(probability):
                continue

            probability = min(max(probability, 0.0), 1.0)
            if probability < minimum_probability:
                continue
            pair = tuple(sorted((left, right)))
            pairs[pair] = max(probability, pairs.get(pair, 0.0))

    if not pairs:
        raise ValueError(
            "没有找到达到最低概率要求的 `i j sqrt(p) ubox` 记录；"
            "请检查文件或降低 --min-probability"
        )

    return [(left, right, pairs[(left, right)]) for left, right in sorted(pairs)]


def make_links(pairs: list[tuple[int, int, float]], chrom: str) -> list[str]:
    """Create links whose interval starts retain RNAfold's 1-based positions."""
    lines: list[str] = []
    for left, right, probability in pairs:
        lines.append(
            "\t".join(
                (
                    chrom,
                    str(left),
                    str(left + 1),
                    chrom,
                    str(right),
                    str(right + 1),
                    f"{probability:.8g}",
                )
            )
        )
    return lines


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "将 RNAfold/ViennaRNA dp.eps 转换为从位置 1 开始的 pyGenomeTracks .links 中间文件。"
        )
    )
    parser.add_argument("dot_plot", type=Path, help="RNAfold dot-plot EPS 文件")
    parser.add_argument("chrom", help="写入 links 文件的染色体名称，例如 c05")
    parser.add_argument(
        "--min-probability",
        type=float,
        default=0.0,
        metavar="P",
        help="保留的最低真实配对概率，范围 0～1（默认：0，即不过滤）",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    if not args.dot_plot.is_file():
        print(f"错误：找不到 dot plot：{args.dot_plot}", file=sys.stderr)
        return 1
    if not args.chrom.strip() or any(character.isspace() for character in args.chrom):
        print("错误：染色体名称不能为空或包含空白字符", file=sys.stderr)
        return 1
    if not 0.0 <= args.min_probability <= 1.0:
        print("错误：--min-probability 必须在 0～1 之间", file=sys.stderr)
        return 1

    output = args.dot_plot.with_suffix(".links")
    try:
        pairs = parse_dot_plot(args.dot_plot, args.min_probability)
        lines = make_links(pairs, args.chrom)
        output.write_text("\n".join(lines) + "\n", encoding="utf-8", newline="\n")
    except (OSError, ValueError) as exc:
        print(f"错误：{exc}", file=sys.stderr)
        return 1

    print(f"已输出 {len(lines)} 对碱基配对：{output.resolve()}")
    print("坐标为转录本相对位置（从 1 开始），尚未映射到真实基因组位置。")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
