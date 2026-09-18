#!/usr/bin/env python3
"""Create a minimal, runnable pyGenomeTracks workspace."""

from __future__ import annotations

import argparse
import shutil
from pathlib import Path


FILES = {
    "tracks.ini": """# Plot settings read by plot.py
# plot_region = chr1:0-100000
# plot_output = output/example.png
# plot_dpi = 300

# pyGenomeTracks track configuration
[genes]
file = data/genes.bed
title = Example genes
height = 0.5
file_type = bed
color = #790000
border_color = black
labels = true
display = stacked
style = UCSC

[x-axis]
where = bottom
""",
    "plot.py": '''#!/usr/bin/env python3
"""Plot this workspace using settings stored in tracks.ini comments."""

from __future__ import annotations

import re
import shlex
import shutil
import subprocess
import sys
from pathlib import Path


SETTING_PATTERN = re.compile(
    r"^\\s*#\\s*plot_(region|output|dpi|width)\\s*=\\s*(.+?)\\s*$"
)


def read_plot_settings(config_path: Path) -> dict[str, str]:
    settings: dict[str, str] = {}
    for line in config_path.read_text(encoding="utf-8-sig").splitlines():
        match = SETTING_PATTERN.match(line)
        if match:
            settings[match.group(1)] = match.group(2)

    missing = [key for key in ("region", "output") if not settings.get(key)]
    if missing:
        names = ", ".join(f"# plot_{key} = ..." for key in missing)
        raise ValueError(f"tracks.ini 缺少绘图设置：{names}")
    return settings


def main() -> int:
    workspace = Path(__file__).resolve().parent
    config_path = workspace / "tracks.ini"
    if not config_path.is_file():
        print(f"错误：找不到配置文件：{config_path}", file=sys.stderr)
        return 1

    executable = shutil.which("pyGenomeTracks")
    if executable is None:
        print("错误：当前环境未找到 pyGenomeTracks，请先安装并激活相应环境。", file=sys.stderr)
        return 1

    try:
        settings = read_plot_settings(config_path)
        dpi = int(settings.get("dpi", "300"))
        if dpi <= 0:
            raise ValueError("plot_dpi 必须是正整数")
        width = settings.get("width")
        if width is not None and float(width) <= 0:
            raise ValueError("plot_width 必须是正数")
    except (OSError, ValueError) as exc:
        print(f"错误：{exc}", file=sys.stderr)
        return 1

    output_path = Path(settings["output"])
    if not output_path.is_absolute():
        output_path = workspace / output_path
    output_path.parent.mkdir(parents=True, exist_ok=True)

    command = [
        executable,
        "--tracks",
        str(config_path),
        "--region",
        settings["region"],
        "--outFileName",
        str(output_path),
        "--dpi",
        str(dpi),
    ]
    if width is not None:
        command.extend(["--width", width])

    print("执行：", shlex.join(command))
    try:
        completed = subprocess.run(command, cwd=workspace, check=False)
    except OSError as exc:
        print(f"错误：无法启动 pyGenomeTracks：{exc}", file=sys.stderr)
        return 1
    if completed.returncode != 0:
        print(f"绘图失败，pyGenomeTracks 退出码：{completed.returncode}", file=sys.stderr)
        return completed.returncode

    print(f"绘图完成：{output_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
''',
    "data/genes.bed": """chr1\t10000\t42000\tGeneA|Transcript1\t0\t+\t13000\t38000\t0,102,204\t3\t5000,7000,4000,\t0,10000,28000,
chr1\t56000\t90000\tGeneB|Transcript1\t0\t-\t60000\t85000\t0,102,204\t3\t6000,8000,5000,\t0,12000,29000,
""",
    "region.bed": "chr1\t0\t100000\n",
    "requirements.txt": "pyGenomeTracks\n",
    "README.md": """# pyGenomeTracks 最小工作区

## 安装

建议使用独立的 Conda 环境：

```bash
conda create -n pygenometracks -c conda-forge -c bioconda pygenometracks
conda activate pygenometracks
```

也可以在已配置好编译依赖的 Python 环境中运行：

```bash
python -m pip install -r requirements.txt
```

## 绘图

先在 `tracks.ini` 顶部设置绘图区间、输出路径等参数，然后运行：

```bash
python3 plot.py
```

生成的图片位于 `output/example.png`。将 `data/genes.bed` 替换为自己的
BED12 文件，并相应修改 `tracks.ini` 即可。
""",
    "output/.gitkeep": "",
}


def create_workspace(destination: Path, force: bool = False) -> None:
    destination = destination.expanduser().resolve()

    if destination.exists() and any(destination.iterdir()):
        if not force:
            raise FileExistsError(
                f"目标目录不是空目录：{destination}\n"
                "如需覆盖同名文件，请添加 --force（不会删除其他文件）。"
            )

    destination.mkdir(parents=True, exist_ok=True)

    for relative_path, content in FILES.items():
        path = destination / relative_path
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(content, encoding="utf-8", newline="\n")

    executable = shutil.which("pyGenomeTracks")
    print(f"已创建 pyGenomeTracks 工作区：{destination}")
    print("绘图命令：")
    print("  python3 plot.py")
    if executable is None:
        print("提示：当前环境未找到 pyGenomeTracks，请先按照 README.md 安装。")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="生成一个包含配置和示例数据的 pyGenomeTracks 工作目录。"
    )
    parser.add_argument(
        "directory",
        nargs="?",
        default="pygenometracks_workspace",
        type=Path,
        help="工作目录路径（默认：./pygenometracks_workspace）",
    )
    parser.add_argument(
        "--force",
        action="store_true",
        help="允许覆盖目标目录中的同名文件，但不会删除其他文件",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    try:
        create_workspace(args.directory, args.force)
    except (FileExistsError, OSError) as exc:
        print(f"错误：{exc}")
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
