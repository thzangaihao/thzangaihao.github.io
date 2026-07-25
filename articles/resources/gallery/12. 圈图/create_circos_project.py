#!/usr/bin/env python3
"""创建一个带详细中文注释的基础 Circos 绘图目录。

本脚本完全独立，不读取 FASTA，也不计算 GC。创建目录后，请将数据生成
脚本产生的 gc_content.txt 和 karyotype.txt 放入项目的 data/ 目录。
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path


# 主配置将各功能拆分到 config/，便于以后增加覆盖度、基因密度等轨道。
CIRCOS_CONF = """# Circos 主配置文件
# 优先载入本项目的自定义颜色名称，供 karyotype 和其他轨道引用。
<<include config/colors.conf>>

# karyotype 定义每条染色体/contig 的名称、长度和颜色。
karyotype = data/karyotype.txt

# u 是 Circos 的相对单位；这里规定 1u = 1,000,000 bp。
chromosomes_units = 1000000
# yes 表示展示 karyotype 中的所有序列。
chromosomes_display_default = yes

# 分别载入染色体外观、刻度、统一轨道容器和输出图像设置。
<<include config/ideogram.conf>>
<<include config/ticks.conf>>
<<include config/tracks.conf>>
<<include config/image.conf>>

# 载入 Circos 安装目录自带的颜色、字体、图案和基础运行参数。
<<include etc/colors_fonts_patterns.conf>>
<<include etc/housekeeping.conf>>
"""

COLORS_CONF = """# 圈图统一配色表
# Circos 在 karyotype.txt 中使用颜色名称，不直接使用带 # 的十六进制值。
# 每个十六进制颜色已转换为等价的十进制 RGB：R,G,B。
# karyotype.txt 最后一列可填写 circle_color_01 至 circle_color_12。
# 染色体多于 12 条时，prepare_circos_gc.py 会从第一种颜色开始循环。
<colors>
circle_color_01 = 59,110,168
circle_color_02 = 211,78,78
circle_color_03 = 99,163,117
circle_color_04 = 224,164,88
circle_color_05 = 108,91,123
circle_color_06 = 42,157,143
circle_color_07 = 231,111,81
circle_color_08 = 242,132,130
circle_color_09 = 132,165,157
circle_color_10 = 38,70,83
circle_color_11 = 168,218,220
circle_color_12 = 244,201,93
</colors>
"""

IDEOGRAM_CONF = """# 染色体圆环的外观参数
<ideogram>
<spacing>
# 相邻染色体间隔，占圆半径的 0.8%。
default = 0.008r
</spacing>

# 染色体圆环位于半径 90% 的位置。
radius = 0.90r
# 圆环宽度为 24 像素。
thickness = 24p
# 使用 karyotype.txt 指定的颜色填充染色体。
fill = yes
# 染色体边框颜色和线宽。
stroke_color = dgrey
stroke_thickness = 1p

# 显示染色体名称。
show_label = yes
label_font = default
# 标签位于染色体外侧，再向外偏移半径的 6%。
label_radius = dims(ideogram,radius_outer) + 0.06r
label_size = 24p
# 标签方向与染色体弧线平行。
label_parallel = yes
</ideogram>
"""

TICKS_CONF = """# 染色体坐标刻度参数
show_ticks = yes
show_tick_labels = yes

<ticks>
# 刻度从染色体圆环外边缘开始。
radius = dims(ideogram,radius_outer)
color = black
thickness = 1p
# 原始坐标乘以 1e-6，因此标签单位显示为 Mb。
multiplier = 1e-6
format = %d

<tick>
# 因为 1u = 1 Mb，所以每 1 Mb 绘制一个小刻度。
spacing = 1u
size = 4p
</tick>

<tick>
# 每 5 Mb 绘制一个带数字标签的大刻度。
spacing = 5u
size = 9p
show_label = yes
label_size = 14p
label_offset = 6p
suffix = Mb
</tick>
</ticks>
"""

TRACKS_CONF = """# 所有数据轨道的唯一容器
# 将 GC、bigWig、基因密度脚本生成的 <plot>...</plot> 片段粘贴到这里，
# 或者在本 <plots> 块内部使用 <<include config/某个_plot.conf>>。
# 整个 Circos 项目只应存在一个 <plots>...</plots> 容器。
<plots>

</plots>
"""

IMAGE_CONF = """# 输出图像参数
<image>
# 图片保存到项目中的 output/ 目录。
dir = output
file = circos_gc.png
# 同时输出 PNG 位图和 SVG 矢量图。
png = yes
svg = yes
# 画布半径；增大该值可提高 PNG 分辨率。
radius = 1500p
background = white
# 将圆形布局旋转 -90 度，让起点位于顶部。
angle_offset = -90
# 自动生成带透明度的颜色，供 *_aN 颜色名称使用。
auto_alpha_colors = yes
auto_alpha_steps = 5
</image>
"""

def make_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="创建独立的 Circos 基础绘图目录和中文注释配置。",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "project_dir", nargs="?", type=Path, default=Path("circos_project"),
        help="要创建的项目目录",
    )
    return parser


def main() -> int:
    args = make_parser().parse_args()
    project_dir = args.project_dir.resolve()

    # 为避免覆盖用户配置，项目目录必须不存在；也不提供强制覆盖选项。
    if project_dir.exists():
        print(f"错误：目录已存在：{project_dir}", file=sys.stderr)
        print("请指定新的目录名称。", file=sys.stderr)
        return 2

    config_dir = project_dir / "config"
    data_dir = project_dir / "data"
    output_dir = project_dir / "output"

    try:
        config_dir.mkdir(parents=True)
        data_dir.mkdir()
        output_dir.mkdir()
        (project_dir / "circos.conf").write_text(CIRCOS_CONF, encoding="utf-8")
        (config_dir / "colors.conf").write_text(COLORS_CONF, encoding="utf-8")
        (config_dir / "ideogram.conf").write_text(IDEOGRAM_CONF, encoding="utf-8")
        (config_dir / "ticks.conf").write_text(TICKS_CONF, encoding="utf-8")
        (config_dir / "tracks.conf").write_text(TRACKS_CONF, encoding="utf-8")
        (config_dir / "image.conf").write_text(IMAGE_CONF, encoding="utf-8")
    except OSError as error:
        print(f"错误：无法创建项目：{error}", file=sys.stderr)
        return 1

    print(f"Circos 项目已创建：{project_dir}")
    print(f"数据目录：{data_dir}")
    print(f"轨道容器：{config_dir / 'tracks.conf'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
