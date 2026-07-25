#!/usr/bin/env python3
"""根据 Circos 颜色配置和图例表生成独立 SVG 图例。

图例表使用制表符分隔，每行格式：

    label<TAB>color<TAB>symbol

例如：

    GC content              circle_color_03    line
    Sequencing coverage     circle_color_01    box
    Gene density            circle_color_04    box

支持的 symbol：line、box、circle、hollow_circle、ribbon。
脚本只使用 Python 标准库，不修改任何 Circos 文件。
"""

from __future__ import annotations

import argparse
import html
import math
import re
import sys
import unicodedata
from dataclasses import dataclass
from pathlib import Path


SUPPORTED_SYMBOLS = {"line", "box", "circle", "hollow_circle", "ribbon"}


@dataclass(frozen=True)
class LegendItem:
    """一条图例记录。"""

    label: str
    color: str
    symbol: str
    minimum: str = "0"
    maximum: str = "1"


def positive_float(text: str) -> float:
    value = float(text)
    if value <= 0:
        raise argparse.ArgumentTypeError("必须是大于 0 的数值")
    return value


def positive_integer(text: str) -> int:
    value = int(text)
    if value <= 0:
        raise argparse.ArgumentTypeError("必须是大于 0 的整数")
    return value


def normalize_rgb(red: int, green: int, blue: int) -> str:
    """校验 RGB，并返回 SVG 使用的十六进制颜色。"""
    if any(channel < 0 or channel > 255 for channel in (red, green, blue)):
        raise ValueError(f"RGB 分量必须在 0–255：{red},{green},{blue}")
    return f"#{red:02X}{green:02X}{blue:02X}"


def parse_color_value(value: str) -> str | None:
    """解析 R,G,B 或 #RRGGBB；不是直接颜色时返回 None。"""
    value = value.strip()
    if re.fullmatch(r"#[0-9A-Fa-f]{6}", value):
        return value.upper()
    rgb_match = re.fullmatch(
        r"\s*(\d{1,3})\s*,\s*(\d{1,3})\s*,\s*(\d{1,3})\s*", value
    )
    if rgb_match:
        return normalize_rgb(*(int(channel) for channel in rgb_match.groups()))
    return None


def read_circos_colors(path: Path) -> dict[str, str]:
    """读取 Circos <colors> 中的 name = R,G,B 定义。

    也支持一个颜色名称引用另一个已定义的颜色名称。
    """
    raw_definitions: dict[str, str] = {}
    inside_colors = False

    with path.open("r", encoding="utf-8-sig") as handle:
        for line_number, raw_line in enumerate(handle, start=1):
            line = raw_line.strip()
            if not line or line.startswith("#"):
                continue
            if line.lower() == "<colors>":
                inside_colors = True
                continue
            if line.lower() == "</colors>":
                inside_colors = False
                continue
            if not inside_colors:
                continue
            if "=" not in line:
                raise ValueError(
                    f"颜色配置第 {line_number} 行缺少 '='：{line}"
                )
            name, value = (part.strip() for part in line.split("=", 1))
            if not re.fullmatch(r"[A-Za-z0-9_.-]+", name):
                raise ValueError(f"颜色名称无效：{name}")
            raw_definitions[name] = value

    if not raw_definitions:
        raise ValueError(f"没有在 <colors> 中找到颜色定义：{path}")

    resolved: dict[str, str] = {}

    def resolve(name: str, chain: set[str]) -> str:
        if name in resolved:
            return resolved[name]
        if name in chain:
            raise ValueError(f"颜色定义形成循环引用：{name}")
        if name not in raw_definitions:
            raise ValueError(f"引用了未定义颜色：{name}")
        value = raw_definitions[name]
        direct = parse_color_value(value)
        if direct is not None:
            resolved[name] = direct
        else:
            resolved[name] = resolve(value, chain | {name})
        return resolved[name]

    for color_name in raw_definitions:
        resolve(color_name, set())
    return resolved


def read_legend_items(path: Path) -> list[tuple[str, str, str, str, str]]:
    """读取手动 TSV。

    支持以下格式：
    label, color
    label, color, symbol
    label, color, min, max
    label, color, min, max, symbol
    """
    items: list[tuple[str, str, str, str, str]] = []
    with path.open("r", encoding="utf-8-sig") as handle:
        for line_number, raw_line in enumerate(handle, start=1):
            line = raw_line.rstrip("\r\n")
            if not line.strip() or line.lstrip().startswith("#"):
                continue
            fields = line.split("\t")
            if len(fields) not in (2, 3, 4, 5):
                raise ValueError(
                    f"图例表第 {line_number} 行必须是 2–5 个制表符分隔字段"
                )
            label, color = fields[0].strip(), fields[1].strip()
            minimum, maximum, symbol = "0", "1", "box"
            if len(fields) == 3:
                symbol = fields[2].strip().lower()
            elif len(fields) >= 4:
                minimum, maximum = fields[2].strip(), fields[3].strip()
                if len(fields) == 5:
                    symbol = fields[4].strip().lower()
            if not label or not color:
                raise ValueError(f"图例表第 {line_number} 行标签或颜色为空")
            if symbol not in SUPPORTED_SYMBOLS:
                choices = ", ".join(sorted(SUPPORTED_SYMBOLS))
                raise ValueError(
                    f"图例表第 {line_number} 行 symbol 无效：{symbol}；可选 {choices}"
                )
            if not minimum or not maximum:
                raise ValueError(f"图例表第 {line_number} 行 min 或 max 为空")
            items.append((label, color, symbol, minimum, maximum))
    if not items:
        raise ValueError("图例表中没有有效记录")
    return items


def resolve_item_colors(
    raw_items: list[tuple[str, str, str, str, str]],
    color_definitions: dict[str, str],
) -> list[LegendItem]:
    """将颜色名称、#RRGGBB 或 R,G,B 转换为 SVG 十六进制颜色。"""
    resolved: list[LegendItem] = []
    for label, color_text, symbol, minimum, maximum in raw_items:
        direct = parse_color_value(color_text)
        if direct is not None:
            color = direct
        elif color_text in color_definitions:
            color = color_definitions[color_text]
        else:
            raise ValueError(f"图例 '{label}' 使用了未定义颜色：{color_text}")
        resolved.append(
            LegendItem(
                label=label,
                color=color,
                symbol=symbol,
                minimum=minimum,
                maximum=maximum,
            )
        )
    return resolved


def expand_track_config(
    path: Path,
    project_dir: Path,
    active_files: set[Path] | None = None,
) -> str:
    """读取 tracks.conf，并递归展开其中的 <<include ...>>。

    Circos 通常从项目根目录解析 ``config/example.conf``。同时也兼容相对于
    当前配置文件目录的 include 路径。
    """
    active_files = set() if active_files is None else active_files
    resolved_path = path.resolve()
    if resolved_path in active_files:
        raise ValueError(f"配置文件形成循环 include：{resolved_path}")
    if not resolved_path.is_file():
        raise ValueError(f"找不到轨道配置：{resolved_path}")

    active_files.add(resolved_path)
    expanded_lines: list[str] = []
    include_pattern = re.compile(r"^\s*<<include\s+(.+?)>>\s*$")
    with resolved_path.open("r", encoding="utf-8-sig") as handle:
        for raw_line in handle:
            match = include_pattern.match(raw_line.rstrip("\r\n"))
            if not match:
                expanded_lines.append(raw_line)
                continue
            include_text = match.group(1).strip().strip("\"'")
            candidate_from_project = project_dir / include_text
            candidate_from_parent = resolved_path.parent / include_text
            include_path = (
                candidate_from_project
                if candidate_from_project.is_file()
                else candidate_from_parent
            )
            expanded_lines.append(
                expand_track_config(
                    include_path, project_dir, active_files=active_files
                )
            )
    active_files.remove(resolved_path)
    return "".join(expanded_lines)


def infer_legend_label(data_file: str, plot_type: str) -> str:
    """在没有 legend_label 注释时，根据数据文件名推断标签。"""
    stem = Path(data_file).name
    for suffix in (".txt", ".tsv", ".bedgraph", ".bed", ".conf"):
        if stem.lower().endswith(suffix):
            stem = stem[: -len(suffix)]
            break
    lower = stem.lower()
    if lower.endswith("_gc") or lower == "gc_content" or "_gc_" in lower:
        return "GC content"
    if "gene_density" in lower:
        return "Gene density"
    if "coverage" in lower:
        return "Sequencing coverage"
    label = re.sub(r"[_-]+", " ", stem).strip()
    return label or f"{plot_type} track"


def infer_symbol(plot_type: str) -> str:
    """将 Circos plot type 映射为 SVG 图例符号。"""
    return {
        "line": "line",
        "scatter": "circle",
        "text": "circle",
        "histogram": "box",
        "heatmap": "box",
        "tile": "box",
    }.get(plot_type, "box")


def read_project_legend_items(
    project_dir: Path,
    tracks_path: Path,
) -> list[tuple[str, str, str, str, str]]:
    """从项目唯一 tracks.conf 及其 include 中提取所有 <plot>。

    可在每个 plot 前或 plot 内使用 Circos 注释元数据：

    ``# legend_label = GC content``
    ``# legend_symbol = line``
    ``# legend_order = 1``
    ``# legend_show = no``
    """
    text = expand_track_config(tracks_path, project_dir)
    lines = text.splitlines()
    metadata_pattern = re.compile(
        r"^\s*#\s*legend_(label|symbol|order|show)\s*=\s*(.*?)\s*$",
        re.IGNORECASE,
    )
    parameter_pattern = re.compile(
        r"^\s*([A-Za-z_][A-Za-z0-9_]*)\s*=\s*(.*?)\s*$"
    )

    pending_metadata: dict[str, str] = {}
    inside_plot = False
    plot_lines: list[str] = []
    plot_metadata: dict[str, str] = {}
    extracted: list[tuple[float, int, str, str, str, str, str]] = []

    for line in lines:
        metadata_match = metadata_pattern.match(line)
        if metadata_match:
            key, value = metadata_match.groups()
            if inside_plot:
                plot_metadata[key.lower()] = value.strip()
            else:
                pending_metadata[key.lower()] = value.strip()
            continue

        stripped = line.strip().lower()
        if stripped == "<plot>":
            if inside_plot:
                raise ValueError("检测到嵌套的 <plot>，请检查 tracks.conf")
            inside_plot = True
            plot_lines = []
            plot_metadata = dict(pending_metadata)
            pending_metadata.clear()
            continue
        if stripped == "</plot>":
            if not inside_plot:
                raise ValueError("检测到没有起始 <plot> 的 </plot>")
            parameters: dict[str, str] = {}
            nested_depth = 0
            for plot_line in plot_lines:
                stripped_plot_line = plot_line.strip()
                # 只提取 <plot> 的直接参数。axes、axis、rules、rule 等
                # 嵌套块中的 color 不能覆盖轨道自身的 color/fill_color。
                if re.fullmatch(r"</[A-Za-z_][A-Za-z0-9_]*>", stripped_plot_line):
                    nested_depth = max(0, nested_depth - 1)
                    continue
                if re.fullmatch(r"<[A-Za-z_][A-Za-z0-9_]*>", stripped_plot_line):
                    nested_depth += 1
                    continue
                if nested_depth != 0:
                    continue
                parameter_match = parameter_pattern.match(plot_line)
                if parameter_match:
                    key, value = parameter_match.groups()
                    parameters[key.lower()] = value.strip()

            if plot_metadata.get("show", "yes").lower() not in {"no", "false", "0"}:
                plot_type = parameters.get("type", "histogram").lower()
                data_file = parameters.get("file", "")
                label = plot_metadata.get("label") or infer_legend_label(
                    data_file, plot_type
                )
                symbol = (
                    plot_metadata.get("symbol", infer_symbol(plot_type)).lower()
                )
                if symbol not in SUPPORTED_SYMBOLS:
                    raise ValueError(
                        f"轨道 '{label}' 的 legend_symbol 无效：{symbol}"
                    )
                color = parameters.get("fill_color") or parameters.get("color")
                if not color:
                    raise ValueError(
                        f"轨道 '{label}' 没有 color 或 fill_color，无法生成图例"
                    )
                # 对逗号分隔的 Circos 颜色取第一种作为图例代表色。
                color = color.split(",", 1)[0].strip()
                minimum = parameters.get("min", "0")
                maximum = parameters.get("max", "1")
                try:
                    order = float(plot_metadata.get("order", len(extracted) + 1))
                except ValueError:
                    raise ValueError(
                        f"轨道 '{label}' 的 legend_order 不是数字"
                    )
                extracted.append(
                    (
                        order,
                        len(extracted),
                        label,
                        color,
                        symbol,
                        minimum,
                        maximum,
                    )
                )
            inside_plot = False
            plot_lines = []
            plot_metadata = {}
            continue
        if inside_plot:
            plot_lines.append(line)

    if inside_plot:
        raise ValueError("存在未闭合的 <plot> 块")
    if not extracted:
        raise ValueError(f"没有在轨道配置中找到可用于图例的 <plot>：{tracks_path}")
    extracted.sort(key=lambda item: (item[0], item[1]))
    return [
        (label, color, symbol, minimum, maximum)
        for _, _, label, color, symbol, minimum, maximum in extracted
    ]


def display_width(text: str) -> float:
    """估算文字宽度：中日韩宽字符按两个拉丁字符计算。"""
    width = 0.0
    for character in text:
        width += (
            2.0
            if unicodedata.east_asian_width(character) in {"W", "F"}
            else 1.0
        )
    return width


def draw_symbol(
    item: LegendItem,
    x: float,
    center_y: float,
    symbol_width: float,
    symbol_height: float,
    line_width: float,
) -> str:
    """返回单个图例符号的 SVG 元素。"""
    color = item.color
    if item.symbol == "line":
        return (
            f'<line x1="{x:g}" y1="{center_y:g}" '
            f'x2="{x + symbol_width:g}" y2="{center_y:g}" '
            f'stroke="{color}" stroke-width="{line_width:g}" '
            'stroke-linecap="round"/>'
        )
    if item.symbol == "circle":
        radius = min(symbol_width, symbol_height) / 2
        return (
            f'<circle cx="{x + symbol_width / 2:g}" cy="{center_y:g}" '
            f'r="{radius:g}" fill="{color}"/>'
        )
    if item.symbol == "hollow_circle":
        radius = min(symbol_width, symbol_height) / 2 - line_width / 2
        return (
            f'<circle cx="{x + symbol_width / 2:g}" cy="{center_y:g}" '
            f'r="{max(radius, 1):g}" fill="none" stroke="{color}" '
            f'stroke-width="{line_width:g}"/>'
        )
    top = center_y - symbol_height / 2
    if item.symbol == "ribbon":
        # 中间收窄的多边形，用于表示 Circos 内部 ribbon/link。
        points = [
            (x, top),
            (x + symbol_width, top),
            (x + symbol_width * 0.62, center_y),
            (x + symbol_width, top + symbol_height),
            (x, top + symbol_height),
            (x + symbol_width * 0.38, center_y),
        ]
        points_text = " ".join(f"{px:g},{py:g}" for px, py in points)
        return f'<polygon points="{points_text}" fill="{color}" fill-opacity="0.75"/>'
    return (
        f'<rect x="{x:g}" y="{top:g}" width="{symbol_width:g}" '
        f'height="{symbol_height:g}" rx="1" fill="{color}"/>'
    )


def make_svg(
    items: list[LegendItem],
    title: str | None,
    columns: int,
    font_size: float,
    title_size: float,
    font_family: str,
    symbol_width: float,
    symbol_height: float,
    line_width: float,
    item_gap: float,
    column_gap: float,
    padding: float,
    background: str,
) -> str:
    """根据内容自动计算画布尺寸并构建 SVG。"""
    columns = min(columns, len(items))
    rows = math.ceil(len(items) / columns)
    row_height = max(font_size * 1.35, symbol_height) + item_gap
    symbol_text_gap = font_size * 0.65

    # 使用按列排布：先填满第一列，再进入下一列。
    column_items: list[list[LegendItem]] = []
    for column in range(columns):
        start = column * rows
        column_items.append(items[start : min(start + rows, len(items))])

    column_widths: list[float] = []
    for group in column_items:
        label_width = max(
            (display_width(item.label) * font_size * 0.55 for item in group),
            default=0,
        )
        column_widths.append(symbol_width + symbol_text_gap + label_width)

    content_width = sum(column_widths) + column_gap * (columns - 1)
    title_height = title_size * 1.5 if title else 0.0
    width = padding * 2 + content_width
    height = padding * 2 + title_height + rows * row_height - item_gap

    escaped_font = html.escape(font_family, quote=True)
    elements: list[str] = []
    if background != "none":
        elements.append(
            f'<rect width="100%" height="100%" fill="{html.escape(background)}"/>'
        )
    if title:
        elements.append(
            f'<text x="{padding:g}" y="{padding + title_size:g}" '
            f'font-family="{escaped_font}" font-size="{title_size:g}" '
            f'font-weight="bold" fill="#222222">{html.escape(title)}</text>'
        )

    x = padding
    first_center_y = padding + title_height + row_height / 2 - item_gap / 2
    for column_index, group in enumerate(column_items):
        for row_index, item in enumerate(group):
            center_y = first_center_y + row_index * row_height
            elements.append(
                draw_symbol(
                    item,
                    x,
                    center_y,
                    symbol_width,
                    symbol_height,
                    line_width,
                )
            )
            text_x = x + symbol_width + symbol_text_gap
            # dominant-baseline=middle 使文字垂直中心与图例符号对齐。
            elements.append(
                f'<text x="{text_x:g}" y="{center_y:g}" '
                f'font-family="{escaped_font}" font-size="{font_size:g}" '
                'dominant-baseline="middle" fill="#222222">'
                f"{html.escape(item.label)}</text>"
            )
        x += column_widths[column_index] + column_gap

    body = "\n  ".join(elements)
    return f"""<?xml version="1.0" encoding="UTF-8"?>
<svg xmlns="http://www.w3.org/2000/svg"
     width="{width:g}" height="{height:g}"
     viewBox="0 0 {width:g} {height:g}">
  {body}
</svg>
"""


def make_scale_svg(
    items: list[LegendItem],
    title: str | None,
    columns: int,
    font_size: float,
    title_size: float,
    font_family: str,
    line_width: float,
    item_gap: float,
    column_gap: float,
    padding: float,
    background: str,
) -> str:
    """生成“彩色括号 + 最大/最小值 + 标签”的比例尺图例。"""
    columns = min(columns, len(items))
    rows = math.ceil(len(items) / columns)
    scale_height = font_size * 2.25
    tick_width = font_size * 0.45
    value_gap = font_size * 0.22
    label_gap = font_size * 0.9
    row_height = scale_height + item_gap

    column_items: list[list[LegendItem]] = []
    for column in range(columns):
        start = column * rows
        column_items.append(items[start : min(start + rows, len(items))])

    column_widths: list[float] = []
    for group in column_items:
        maximum_value_width = max(
            (
                max(display_width(item.minimum), display_width(item.maximum))
                * font_size
                * 0.55
                for item in group
            ),
            default=0,
        )
        label_width = max(
            (display_width(item.label) * font_size * 0.55 for item in group),
            default=0,
        )
        column_widths.append(
            tick_width + value_gap + maximum_value_width + label_gap + label_width
        )

    content_width = sum(column_widths) + column_gap * (columns - 1)
    title_height = title_size * 1.5 if title else 0.0
    width = padding * 2 + content_width
    height = padding * 2 + title_height + rows * row_height - item_gap
    escaped_font = html.escape(font_family, quote=True)
    elements: list[str] = []

    if background != "none":
        elements.append(
            f'<rect width="100%" height="100%" fill="{html.escape(background)}"/>'
        )
    if title:
        elements.append(
            f'<text x="{padding:g}" y="{padding + title_size:g}" '
            f'font-family="{escaped_font}" font-size="{title_size:g}" '
            f'font-weight="bold" fill="#222222">{html.escape(title)}</text>'
        )

    x = padding
    first_top = padding + title_height
    for column_index, group in enumerate(column_items):
        maximum_value_width = max(
            (
                max(display_width(item.minimum), display_width(item.maximum))
                * font_size
                * 0.55
                for item in group
            ),
            default=0,
        )
        for row_index, item in enumerate(group):
            top = first_top + row_index * row_height
            bottom = top + scale_height
            middle = (top + bottom) / 2
            # 三段折线构成向右开口的彩色比例尺括号。
            path = (
                f"M {x + tick_width:g},{top:g} H {x:g} "
                f"V {bottom:g} H {x + tick_width:g}"
            )
            elements.append(
                f'<path d="{path}" fill="none" stroke="{item.color}" '
                f'stroke-width="{line_width:g}" stroke-linecap="square" '
                'stroke-linejoin="miter"/>'
            )
            value_x = x + tick_width + value_gap
            elements.append(
                f'<text x="{value_x:g}" y="{top + font_size * 0.78:g}" '
                f'font-family="{escaped_font}" font-size="{font_size:g}" '
                f'fill="#222222">{html.escape(item.maximum)}</text>'
            )
            elements.append(
                f'<text x="{value_x:g}" y="{bottom:g}" '
                f'font-family="{escaped_font}" font-size="{font_size:g}" '
                f'fill="#222222">{html.escape(item.minimum)}</text>'
            )
            label_x = value_x + maximum_value_width + label_gap
            elements.append(
                f'<text x="{label_x:g}" y="{middle:g}" '
                f'font-family="{escaped_font}" font-size="{font_size:g}" '
                'dominant-baseline="middle" fill="#222222">'
                f"{html.escape(item.label)}</text>"
            )
        x += column_widths[column_index] + column_gap

    body = "\n  ".join(elements)
    return f"""<?xml version="1.0" encoding="UTF-8"?>
<svg xmlns="http://www.w3.org/2000/svg"
     width="{width:g}" height="{height:g}"
     viewBox="0 0 {width:g} {height:g}">
  {body}
</svg>
"""


def make_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="读取 Circos colors.conf 和 TSV 图例表，生成独立 SVG 图例。",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "legend_table", nargs="?", type=Path,
        help="可选的手动 TSV；省略时自动扫描 Circos 项目",
    )
    parser.add_argument(
        "--project-dir", type=Path, default=Path("."),
        help="Circos 项目目录；自动模式默认使用当前目录",
    )
    parser.add_argument(
        "-c", "--colors", type=Path,
        help="colors.conf；默认 <project-dir>/config/colors.conf",
    )
    parser.add_argument(
        "--tracks", type=Path,
        help="tracks.conf；默认 <project-dir>/config/tracks.conf",
    )
    parser.add_argument(
        "-o", "--output", type=Path, default=Path("legend.svg"),
        help="输出 SVG 文件",
    )
    parser.add_argument("--title", help="可选的图例标题")
    parser.add_argument(
        "--legend-style", choices=("scale", "symbol"), default="scale",
        help="scale 绘制彩色比例尺；symbol 使用传统符号图例",
    )
    parser.add_argument(
        "--orientation", choices=("vertical", "horizontal"), default="vertical",
        help="默认排版方向；指定 --columns 后以 columns 为准",
    )
    parser.add_argument(
        "--columns", type=positive_integer,
        help="图例列数；默认纵向 1 列、横向每项 1 列",
    )
    parser.add_argument("--font-size", type=positive_float, default=18)
    parser.add_argument("--title-size", type=positive_float, default=22)
    parser.add_argument(
        "--font-family", default="Arial, Noto Sans CJK SC, sans-serif",
        help="SVG 字体族；中文图例建议保留 Noto Sans CJK SC",
    )
    parser.add_argument("--symbol-width", type=positive_float, default=30)
    parser.add_argument("--symbol-height", type=positive_float, default=14)
    parser.add_argument("--line-width", type=positive_float, default=3)
    parser.add_argument("--item-gap", type=positive_float, default=10)
    parser.add_argument("--column-gap", type=positive_float, default=28)
    parser.add_argument("--padding", type=positive_float, default=16)
    parser.add_argument(
        "--background", default="none",
        help="背景色，如 white、#FFFFFF；none 表示透明",
    )
    parser.add_argument("--overwrite", action="store_true", help="覆盖已有 SVG")
    return parser


def main() -> int:
    args = make_parser().parse_args()
    project_dir = args.project_dir.resolve()
    colors_path = (
        args.colors.resolve()
        if args.colors
        else project_dir / "config" / "colors.conf"
    )
    tracks_path = (
        args.tracks.resolve()
        if args.tracks
        else project_dir / "config" / "tracks.conf"
    )
    output = args.output.resolve()

    if not colors_path.is_file():
        print(f"错误：找不到颜色配置：{colors_path}", file=sys.stderr)
        return 2
    if args.legend_table and not args.legend_table.is_file():
        print(f"错误：找不到图例表：{args.legend_table}", file=sys.stderr)
        return 2
    if not args.legend_table and not tracks_path.is_file():
        print(f"错误：找不到轨道配置：{tracks_path}", file=sys.stderr)
        return 2
    if output.exists() and not args.overwrite:
        print(f"错误：输出已存在：{output}；请添加 --overwrite。", file=sys.stderr)
        return 2

    try:
        color_definitions = read_circos_colors(colors_path)
        raw_items = (
            read_legend_items(args.legend_table.resolve())
            if args.legend_table
            else read_project_legend_items(project_dir, tracks_path)
        )
        items = resolve_item_colors(raw_items, color_definitions)
        columns = args.columns or (
            1 if args.orientation == "vertical" else len(items)
        )
        if args.legend_style == "scale":
            svg = make_scale_svg(
                items=items,
                title=args.title,
                columns=columns,
                font_size=args.font_size,
                title_size=args.title_size,
                font_family=args.font_family,
                line_width=args.line_width,
                item_gap=args.item_gap,
                column_gap=args.column_gap,
                padding=args.padding,
                background=args.background,
            )
        else:
            svg = make_svg(
                items=items,
                title=args.title,
                columns=columns,
                font_size=args.font_size,
                title_size=args.title_size,
                font_family=args.font_family,
                symbol_width=args.symbol_width,
                symbol_height=args.symbol_height,
                line_width=args.line_width,
                item_gap=args.item_gap,
                column_gap=args.column_gap,
                padding=args.padding,
                background=args.background,
            )
        output.parent.mkdir(parents=True, exist_ok=True)
        output.write_text(svg, encoding="utf-8", newline="\n")
    except (OSError, UnicodeError, ValueError) as error:
        output.unlink(missing_ok=True)
        print(f"错误：生成图例失败：{error}", file=sys.stderr)
        return 1

    print(f"完成：生成 {len(items)} 项图例")
    print(f"SVG：{output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
