# Circos 圈图数据与配置工具

## 设计原则

- 整个项目只允许存在一个 `<plots>...</plots>`，位于 `config/tracks.conf`。
- 各数据脚本生成的配置都只有一个 `<plot>...</plot>`。
- 数据脚本不修改主配置，也不自动修改 `tracks.conf`。
- 将需要的 plot 片段复制到 `tracks.conf` 的 `<plots>` 内即可。

## 1. 创建工作目录

```bash
python create_circos_project.py circos_project
```

生成：

```text
circos_project/
├── circos.conf
├── config/
│   ├── colors.conf
│   ├── image.conf
│   ├── ticks.conf
│   ├── ideogram.conf
│   └── tracks.conf
├── data/
└── output/
```

其中 `tracks.conf` 默认为：

```text
<plots>

</plots>
```

## 2. 生成核型数据

```bash
python prepare_karyotype.py genome.fasta \
  -o circos_project/data/karyotype.txt
```

该脚本只输出 `karyotype.txt`。

## 3. 生成 GC 数据和 plot

```bash
python prepare_gc_track.py genome.fasta \
  -o gc_result \
  --prefix genome \
  --window-size 100000
```

输出：

```text
genome_gc.txt
genome_gc_plot.conf
```

## 4. 生成 bigWig 覆盖度数据和 plot

依赖：

```bash
conda install -c conda-forge -c bioconda pybigwig
```

运行：

```bash
python prepare_bigwig_track.py coverage.bw \
  -o coverage_result \
  --prefix sequencing \
  --window-size 100000 \
  --exact
```

输出：

```text
sequencing_coverage.txt
sequencing_coverage_plot.conf
```

## 5. 生成基因密度数据和 plot

```bash
python prepare_gene_density_track.py genome.gff3 \
  -o gene_result \
  --prefix genome \
  --window-size 100000 \
  --chrom-sizes circos_project/data/karyotype.txt
```

输出：

```text
genome_gene_density.txt
genome_gene_density_plot.conf
```

## 组合轨道

把数据文件复制到项目的 `data/`，然后把各 `_plot.conf` 的内容放入：

```text
config/tracks.conf
```

最终结构：

```text
<plots>

<plot>
# GC
</plot>

<plot>
# bigWig coverage
</plot>

<plot>
# gene density
</plot>

</plots>
```

最后运行：

```bash
cd circos_project
circos -conf circos.conf
```

## 6. 生成独立 SVG 图例

进入已经配置好轨道的 Circos 项目目录，直接运行：

```bash
cd circos_project
python /path/to/create_circos_legend.py
```

脚本默认自动读取：

```text
config/colors.conf
config/tracks.conf
```

并递归读取 `tracks.conf` 中的 `_plot.conf` include，输出当前目录下的
`legend.svg`。默认图例样式为彩色比例尺：颜色来自 `fill_color`（若没有则
使用 `color`），上下限来自 plot 的 `max` 和 `min`。

每个 plot 可以使用不会影响 Circos 的注释控制图例：

```text
# legend_label = GC content
# legend_symbol = line
# legend_order = 1
# legend_show = yes
```

如果不写这些注释，脚本会根据 `type`、`file`、`color` 和 `fill_color`
自动推断。

生成类似论文中的多列比例尺图例：

```bash
python /path/to/create_circos_legend.py \
  --columns 3 \
  --font-size 18 \
  --line-width 3
```

如需恢复传统线条/色块图例：

```bash
python /path/to/create_circos_legend.py --legend-style symbol
```

仍然可以准备制表符分隔的 `legend.tsv` 进行手动控制：

```text
GC content	circle_color_03	0	0.58
Sequencing coverage	circle_color_01	0	150
Gene density	circle_color_04	0	45
```

运行：

```bash
python create_circos_legend.py legend.tsv \
  --colors circos_project/config/colors.conf \
  --output legend.svg \
  --title "Tracks"
```

支持的图例符号为 `line`、`box`、`circle`、`hollow_circle` 和 `ribbon`。
SVG 默认使用透明背景。
