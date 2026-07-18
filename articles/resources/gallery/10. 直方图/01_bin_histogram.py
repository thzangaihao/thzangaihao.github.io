#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
交互式等宽分箱脚本。

运行方式：
  python 01_bin_histogram.py

输出 CSV 可直接交给 02_fit_histogram.py 使用。
"""

from __future__ import annotations

import csv
import math
import re
from pathlib import Path


SCRIPT_DIR = Path.cwd()


config = {
    # 输入文件。无表头单列数值文件，支持空格、逗号、制表符和换行分隔。
    "input": "input_cis.tsv",

    # 输出分箱 CSV。
    "output_csv": "histogram_bins.csv",

    # 等宽分箱宽度，必须大于 0。
    "bin_width": 50.0,

    # 分箱统计范围。留空/None 表示按数据自动计算。
    "x_min": None,
    "x_max": None,

    # 是否把分箱边界对齐到 bin_width 的整数倍。
    # True: 例如 bin_width=50 时从 0、50、100... 这类边界开始。
    # False: 直接从 x_min 开始。
    "align_to_bin_width": True,

    # 区间闭合方式。推荐 "left": [start, end)，最后一箱包含右边界。
    # 可选 "left" 或 "right"；"right" 表示 (start, end]，第一箱包含左边界。
    "closed": "left",

    # 是否在交互式运行时逐项询问参数。False 时直接使用上面的 config。
    "interactive": True,
}


def as_script_path(value: str | Path) -> Path:
    path = Path(value).expanduser()
    if path.is_absolute():
        return path
    return path


def prompt_text(label: str, default: str) -> str:
    raw = input(f"{label} [{default}]: ").strip()
    return raw if raw else default


def prompt_float(label: str, default: float | None, allow_none: bool = True) -> float | None:
    shown = "" if default is None else str(default)
    while True:
        raw = input(f"{label} [{shown}]: ").strip()
        if not raw:
            return default
        if allow_none and raw.lower() in {"na", "none", "null"}:
            return None
        try:
            return float(raw)
        except ValueError:
            print("请输入数值；留空表示使用默认值。")


def prompt_bool(label: str, default: bool) -> bool:
    shown = "y" if default else "n"
    while True:
        raw = input(f"{label} y/n [{shown}]: ").strip().lower()
        if not raw:
            return default
        if raw in {"y", "yes", "true", "1"}:
            return True
        if raw in {"n", "no", "false", "0"}:
            return False
        print("请输入 y 或 n。")


def prompt_choice(label: str, default: str, choices: set[str]) -> str:
    while True:
        raw = input(f"{label} {sorted(choices)} [{default}]: ").strip().lower()
        if not raw:
            return default
        if raw in choices:
            return raw
        print(f"可选值：{', '.join(sorted(choices))}")


def collect_config() -> dict:
    cfg = dict(config)
    if not cfg["interactive"]:
        return cfg

    print("等宽分箱参数设置；直接回车使用默认值。")
    cfg["input"] = prompt_text("输入文件", str(cfg["input"]))
    cfg["output_csv"] = prompt_text("输出分箱 CSV", str(cfg["output_csv"]))
    cfg["bin_width"] = prompt_float("分箱宽度 bin_width", cfg["bin_width"], allow_none=False)
    cfg["x_min"] = prompt_float("统计范围 x_min，留空自动", cfg["x_min"])
    cfg["x_max"] = prompt_float("统计范围 x_max，留空自动", cfg["x_max"])
    cfg["align_to_bin_width"] = prompt_bool("边界是否对齐到 bin_width 整数倍", cfg["align_to_bin_width"])
    cfg["closed"] = prompt_choice("区间闭合方式", cfg["closed"], {"left", "right"})
    return cfg


def read_values(path: Path) -> list[float]:
    text = path.read_text(encoding="utf-8-sig", errors="ignore")
    tokens = [x for x in re.split(r"[\s,;]+", text.strip()) if x]
    values: list[float] = []
    bad = 0
    for token in tokens:
        try:
            value = float(token)
        except ValueError:
            bad += 1
            continue
        if math.isfinite(value):
            values.append(value)
        else:
            bad += 1
    if bad:
        print(f"已忽略 {bad} 个非数值或非有限值。")
    if not values:
        raise SystemExit("输入文件中没有可用数值。")
    return values


def make_breaks(values: list[float], cfg: dict) -> list[float]:
    width = float(cfg["bin_width"])
    if width <= 0:
        raise SystemExit("bin_width 必须大于 0。")

    x_min = min(values) if cfg["x_min"] is None else float(cfg["x_min"])
    x_max = max(values) if cfg["x_max"] is None else float(cfg["x_max"])
    if x_min >= x_max:
        raise SystemExit("x_min 必须小于 x_max。")

    if cfg["align_to_bin_width"]:
        start = math.floor(x_min / width) * width
        end = math.ceil(x_max / width) * width
    else:
        start = x_min
        end = x_max

    breaks = []
    current = start
    # 加半个步长容忍浮点误差。
    while current <= end + width * 1e-10:
        breaks.append(current)
        current += width
    if len(breaks) < 2 or breaks[-1] < end:
        breaks.append(breaks[-1] + width)
    return breaks


def count_bins(values: list[float], breaks: list[float], closed: str) -> list[int]:
    counts = [0] * (len(breaks) - 1)
    first, last = breaks[0], breaks[-1]
    width = breaks[1] - breaks[0]
    for value in values:
        if value < first or value > last:
            continue
        if closed == "left":
            if value == last:
                idx = len(counts) - 1
            else:
                idx = int(math.floor((value - first) / width))
        else:
            if value == first:
                idx = 0
            else:
                idx = int(math.ceil((value - first) / width)) - 1
        if 0 <= idx < len(counts):
            counts[idx] += 1
    return counts


def interval_label(left: float, right: float, index: int, total: int, closed: str) -> str:
    if closed == "left":
        right_bracket = "]" if index == total - 1 else ")"
        return f"[{left:g}, {right:g}{right_bracket}"
    left_bracket = "[" if index == 0 else "("
    return f"{left_bracket}{left:g}, {right:g}]"


def main() -> None:
    cfg = collect_config()
    input_path = as_script_path(cfg["input"])
    output_path = as_script_path(cfg["output_csv"])
    if not input_path.exists():
        raise SystemExit(f"输入文件不存在：{input_path}")

    values = read_values(input_path)
    breaks = make_breaks(values, cfg)
    counts = count_bins(values, breaks, cfg["closed"])

    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", newline="", encoding="utf-8-sig") as handle:
        writer = csv.writer(handle)
        writer.writerow(["interval", "bin_midpoint", "count"])
        total = len(counts)
        for i, count in enumerate(counts):
            left, right = breaks[i], breaks[i + 1]
            writer.writerow([
                interval_label(left, right, i, total, cfg["closed"]),
                (left + right) / 2,
                count,
            ])

    print(f"分箱 CSV 已输出：{output_path}")
    print(f"有效原始数值：{len(values)}；分箱数：{len(counts)}；总频数：{sum(counts)}")
    print(f"分箱范围：{breaks[0]:g} 至 {breaks[-1]:g}；bin_width={cfg['bin_width']:g}")


if __name__ == "__main__":
    main()
