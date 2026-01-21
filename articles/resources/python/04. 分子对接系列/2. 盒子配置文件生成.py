#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# 盒子计算，生成配置文件 #
import os
import sys
import glob

# ============= 0. 参数设置 =============
BUFFER_ANGSTROM = 10.0  # 缓冲距离 (埃米): 确保盒子比蛋白大一圈，防止边缘效应

# ============= 1. 基础工具 =============
def get_base_dir():
    if getattr(sys, 'frozen', False):
        return os.path.dirname(sys.executable)
    return os.path.dirname(os.path.abspath(__file__))

def find_files_recursive(directory, ext):
    if not ext.startswith('.'): ext = '.' + ext
    search_pattern = os.path.join(directory, '**', f'*{ext}')
    return sorted(glob.glob(search_pattern, recursive=True))

# ============= 2. 核心：计算盒子坐标 =============
def calculate_bounding_box(pdbqt_file):
    """
    读取 PDBQT，计算所有原子的最小/最大坐标
    返回: center_x, center_y, center_z, size_x, size_y, size_z
    """
    x_coords = []
    y_coords = []
    z_coords = []

    with open(pdbqt_file, 'r') as f:
        for line in f:
            # 只读取 ATOM 和 HETATM 行
            if line.startswith("ATOM") or line.startswith("HETATM"):
                try:
                    # PDBQT 格式固定列宽解析 (比 split 更稳健)
                    # X: 30-38, Y: 38-46, Z: 46-54
                    x = float(line[30:38])
                    y = float(line[38:46])
                    z = float(line[46:54])
                    
                    x_coords.append(x)
                    y_coords.append(y)
                    z_coords.append(z)
                except ValueError:
                    continue

    if not x_coords:
        return None

    # 计算极值
    min_x, max_x = min(x_coords), max(x_coords)
    min_y, max_y = min(y_coords), max(y_coords)
    min_z, max_z = min(z_coords), max(z_coords)

    # 计算中心 (Center)
    center_x = (min_x + max_x) / 2.0
    center_y = (min_y + max_y) / 2.0
    center_z = (min_z + max_z) / 2.0

    # 计算尺寸 (Size) + 缓冲
    # Vina 的 size 是指盒子边缘长度
    size_x = (max_x - min_x) + BUFFER_ANGSTROM
    size_y = (max_y - min_y) + BUFFER_ANGSTROM
    size_z = (max_z - min_z) + BUFFER_ANGSTROM

    return (center_x, center_y, center_z, size_x, size_y, size_z)

# ============= 3. 生成 Vina 配置文件 =============
def write_config_file(receptor_path, output_dir, box_params):
    """
    生成 conf.txt 文件
    """
    cx, cy, cz, sx, sy, sz = box_params
    
    # 获取相对文件名
    rec_name = os.path.basename(receptor_path)
    base_name = os.path.splitext(rec_name)[0]
    
    config_path = os.path.join(output_dir, f"conf_{base_name}.txt")
    
    # Vina 配置内容模板
    # 注意：receptor 和 ligand 的路径我们在后续运行脚本里动态指定，
    # 这里只存盒子参数，或者写入相对路径。
    # 为了通用性，这里只写入盒子参数。
    
    content = f"""
        receptor = {receptor_path}
        # ligand = (将在运行命令中指定)

        center_x = {cx:.3f}
        center_y = {cy:.3f}
        center_z = {cz:.3f}

        size_x = {sx:.3f}
        size_y = {sy:.3f}
        size_z = {sz:.3f}

        cpu = 64
        exhaustiveness = 64
        """
    with open(config_path, 'w') as f:
        f.write(content)
    
    return config_path

# ============= 4. 主程序 =============
def main():
    base_path = get_base_dir()
    
    # 输入：处理好的受体 PDBQT
    rec_in_dir = os.path.join(base_path, "1. Receptor_protein_pdbqt")
    
    # 输出：配置文件目录
    config_out_dir = os.path.join(base_path, "2. Docking_Config")
    
    if not os.path.exists(rec_in_dir):
        print(f"错误: 未找到目录 {rec_in_dir}，请先完成 Step 1。")
        return
        
    if not os.path.exists(config_out_dir):
        os.makedirs(config_out_dir)

    print("=== Step 2: 自动计算 Grid Box 并生成配置文件 ===")
    
    # 查找所有受体
    receptors = find_files_recursive(rec_in_dir, '.pdbqt')
    
    if not receptors:
        print("未找到受体 PDBQT 文件。")
        return

    # 简单交互
    print(f"找到 {len(receptors)} 个受体文件。")
    print(f"将为每个受体计算 Blind Docking 盒子 (Buffer = {BUFFER_ANGSTROM} Å)")
    
    if input("确认开始生成? (y/n): ").strip().lower() not in ['y', 'yes']:
        return

    # 批量处理
    print("\n[开始计算]...")
    count = 0
    for i, rec in enumerate(receptors, 1):
        try:
            box = calculate_bounding_box(rec)
            if box:
                write_config_file(rec, config_out_dir, box)
                count += 1
                sys.stdout.write(f"\r  进度: {i}/{len(receptors)} | 生成配置: conf_{os.path.basename(rec)[:-6]}.txt")
                sys.stdout.flush()
            else:
                print(f"\n  [跳过] 无法读取坐标: {os.path.basename(rec)}")
        except Exception as e:
            print(f"\n  [错误] {os.path.basename(rec)}: {e}")

    print(f"\n\n任务完成！")
    print(f"成功生成 {count} 个配置文件。")
    print(f"配置文件保存在: {config_out_dir}")
    print("你可以打开其中一个 conf_*.txt 检查 center 和 size 是否合理。")
    input("按回车键退出...")

if __name__ == "__main__":
    main()