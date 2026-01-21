#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import os
import sys
import glob
import shutil  # 用于文件复制

# ============= 基础路径获取（复用 cite_v1 逻辑） =============
def get_base_dir():
    """获取脚本所在目录"""
    if getattr(sys, 'frozen', False):
        return os.path.dirname(sys.executable)
    return os.path.dirname(os.path.abspath(__file__))

# ============= 高级文件查找（复用 cite_v1 逻辑） =============
def find_files(filename_pattern, path=None):
    """递归查找文件"""
    if path is None:
        path = get_base_dir()
    return glob.glob(os.path.join(path, '**', filename_pattern), recursive=True)

# ============= 核心处理逻辑 =============
def copy_cif_files(file_list, dest_dir):
    """
    复制文件到目标目录，并重命名以防冲突
    """
    print(f"准备处理 {len(file_list)} 个文件...")
    count = 0
    
    for i, src_path in enumerate(file_list, 1):
        try:
            # 获取原始文件名 (例如 t01_model_0.cif)
            original_name = os.path.basename(src_path)
            
            # 组合新文件名: 父文件夹名_原始文件名
            # 结果示例: fold_fs000286_t01_t01_model_0.cif
            new_filename = f"{original_name}"
            
            # 目标完整路径
            dest_path = os.path.join(dest_dir, new_filename)
            
            # 复制文件 (copy2 保留文件元数据如时间戳)
            shutil.copy2(src_path, dest_path)
            count += 1
            
            # 打印进度
            if i % 50 == 0:
                print(f"已复制 {i}/{len(file_list)}...")
                
        except Exception as e:
            print(f"[Error] 复制 {src_path} 失败: {e}")

    return count

# ============= 主函数 =============
def main():
    # 1. 设置参数
    target_pattern = '*_model_0.cif'  # 搜索模式
    output_folder_name = 'collected_cif_files' # 存放结果的文件夹名
    
    base_dir = get_base_dir()
    dest_dir = os.path.join(base_dir, output_folder_name)
    
    print(f"正在目录 {base_dir} 下搜索 {target_pattern} ...")
    
    # 2. 搜索文件
    found_files = find_files(target_pattern, base_dir)
    
    if not found_files:
        print("未找到任何符合条件的文件。")
        return

    print(f"共找到 {len(found_files)} 个 CIF 文件。")
    
    # 3. 创建存放目录
    if not os.path.exists(dest_dir):
        os.makedirs(dest_dir)
        print(f"已创建存放目录: {dest_dir}")
    else:
        print(f"存放目录已存在: {dest_dir}")

    # 4. 执行复制
    success_count = copy_cif_files(found_files, dest_dir)
    
    print("=" * 50)
    print(f"处理完成！")
    print(f"成功复制: {success_count} 个文件")
    print(f"文件存放位置: {dest_dir}")
    print("=" * 50)

if __name__ == "__main__":
    main()