#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
交互式线性拟合脚本。

读取 01_bin_histogram.py 输出的 CSV，手动选择自变量和因变量字段，只做
普通最小二乘线性拟合。输出一个曲线 CSV，供 03_plot_histogram_fit.R 使用。
"""

from __future__ import annotations

import csv
import math
from pathlib import Path
from statistics import NormalDist


config = {
    # 输入 CSV，通常是 01_bin_histogram.py 生成的 histogram_bins.csv。
    "input_csv": "histogram_bins.csv",

    # 输出拟合曲线 CSV。
    "output_csv": "histogram_fit_curve.csv",

    # 默认自变量和因变量字段；交互式运行时可修改。
    "x_field": "bin_midpoint",
    "y_field": "count",

    # 曲线点数量。
    "curve_points": 500,

    # 是否交互式询问参数。False 时直接使用 config。
    "interactive": True,
}


def as_path(value: str | Path) -> Path:
    return Path(value).expanduser()


def prompt_text(label: str, default: str) -> str:
    raw = input(f"{label} [{default}]: ").strip()
    return raw if raw else default


def prompt_int(label: str, default: int) -> int:
    while True:
        raw = input(f"{label} [{default}]: ").strip()
        if not raw:
            return default
        try:
            value = int(raw)
        except ValueError:
            print("请输入整数。")
            continue
        if value > 1:
            return value
        print("数值必须大于 1。")


def collect_config() -> dict:
    cfg = dict(config)
    if not cfg["interactive"]:
        return cfg
    print("线性拟合参数设置；直接回车使用默认值。")
    cfg["input_csv"] = prompt_text("输入 CSV", cfg["input_csv"])
    cfg["output_csv"] = prompt_text("输出拟合 CSV", cfg["output_csv"])
    cfg["x_field"] = prompt_text("自变量字段 x", cfg["x_field"])
    cfg["y_field"] = prompt_text("因变量字段 y", cfg["y_field"])
    cfg["curve_points"] = prompt_int("曲线点数量", cfg["curve_points"])
    return cfg


def read_xy(path: Path, x_field: str, y_field: str) -> tuple[list[float], list[float]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle)
        if reader.fieldnames is None:
            raise SystemExit("输入 CSV 没有表头。")
        missing = [field for field in (x_field, y_field) if field not in reader.fieldnames]
        if missing:
            raise SystemExit(f"输入 CSV 缺少字段：{', '.join(missing)}")
        x_values: list[float] = []
        y_values: list[float] = []
        skipped = 0
        for row in reader:
            try:
                x = float(row[x_field])
                y = float(row[y_field])
            except (TypeError, ValueError):
                skipped += 1
                continue
            if math.isfinite(x) and math.isfinite(y):
                x_values.append(x)
                y_values.append(y)
            else:
                skipped += 1
    if skipped:
        print(f"已忽略 {skipped} 行非数值或非有限值。")
    if len(x_values) < 3:
        raise SystemExit("线性拟合至少需要 3 个有效点。")
    if len(set(x_values)) < 2:
        raise SystemExit("自变量至少需要 2 个不同取值。")
    return x_values, y_values


def linear_fit(x_values: list[float], y_values: list[float]) -> dict:
    n = len(x_values)
    x_mean = sum(x_values) / n
    y_mean = sum(y_values) / n
    sxx = sum((x - x_mean) ** 2 for x in x_values)
    sxy = sum((x - x_mean) * (y - y_mean) for x, y in zip(x_values, y_values))
    slope = sxy / sxx
    intercept = y_mean - slope * x_mean
    fitted = [intercept + slope * x for x in x_values]
    residuals = [y - y_hat for y, y_hat in zip(y_values, fitted)]
    sse = sum(r * r for r in residuals)
    sst = sum((y - y_mean) ** 2 for y in y_values)
    r_squared = 1 - sse / sst if sst > 0 else float("nan")
    df = n - 2
    mse = sse / df
    se_slope = math.sqrt(mse / sxx)
    t_value = slope / se_slope if se_slope > 0 else float("inf")
    p_value = two_sided_t_p_value(t_value, df)
    return {
        "intercept": intercept,
        "slope": slope,
        "r_squared": r_squared,
        "p_value": p_value,
        "n": n,
    }


def two_sided_t_p_value(t_value: float, df: int) -> float:
    if not math.isfinite(t_value):
        return 0.0
    t_abs = abs(t_value)
    if df <= 0:
        return float("nan")
    # 优先使用 scipy；没有 scipy 时使用正态近似。对常见分箱数量而言近似足够稳定。
    try:
        from scipy import stats  # type: ignore

        return float(2 * stats.t.sf(t_abs, df))
    except Exception:
        return 2 * (1 - NormalDist().cdf(t_abs))


def fmt(value: float, digits: int = 8) -> str:
    return f"{value:.{digits}g}"


def make_equation(intercept: float, slope: float, x_name: str, y_name: str) -> str:
    sign = "+" if slope >= 0 else "-"
    return f"{y_name} = {fmt(intercept)} {sign} {fmt(abs(slope))} {x_name}"


def write_curve(path: Path, fit: dict, x_min: float, x_max: float, points: int, x_name: str, y_name: str) -> None:
    equation = make_equation(fit["intercept"], fit["slope"], x_name, y_name)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["x", "predicted_count", "equation", "R^2", "p_value"])
        if points == 1:
            xs = [x_min]
        else:
            step = (x_max - x_min) / (points - 1)
            xs = [x_min + i * step for i in range(points)]
        for x in xs:
            y_hat = fit["intercept"] + fit["slope"] * x
            writer.writerow([x, y_hat, equation, fit["r_squared"], fit["p_value"]])


def main() -> None:
    cfg = collect_config()
    input_path = as_path(cfg["input_csv"])
    output_path = as_path(cfg["output_csv"])
    if not input_path.exists():
        raise SystemExit(f"输入 CSV 不存在：{input_path}")
    x_values, y_values = read_xy(input_path, cfg["x_field"], cfg["y_field"])
    fit = linear_fit(x_values, y_values)
    write_curve(output_path, fit, min(x_values), max(x_values), int(cfg["curve_points"]), cfg["x_field"], cfg["y_field"])
    print(f"拟合 CSV 已输出：{output_path}")
    print(make_equation(fit["intercept"], fit["slope"], cfg["x_field"], cfg["y_field"]))
    print(f"R^2={fit['r_squared']:.8g}; p_value={fit['p_value']:.8g}; n={fit['n']}")


if __name__ == "__main__":
    main()
