'''
AUGUSTUS GFF 转 IGV 兼容 GFF3 转换器

目录构型：
    自动检测脚本同级目录及子目录中的.gff文件

用于DNA_annotation工作流
去除AUGUSYUS生成gff文件中所有的注释行, 只保留基因结构

thz 2025/10/18
'''

from pathlib import Path

def get_script_directory():
    """获取脚本所在目录的绝对路径"""
    return Path(__file__).parent.absolute()

def find_augustus_gff_files():
    """
    自动查找目录下的AUGUSTUS生成的GFF文件
    使用模糊匹配，查找所有.gff文件
    """
    script_dir = get_script_directory()
    print(script_dir)
    
    # 查找所有.gff文件
    gff_files = list(script_dir.glob("*.gff"))
    
    # 排除已经转换过的文件（以_igv.gff3结尾的）
    gff_files = [f for f in gff_files if not f.name.endswith('_igv.gff3')]
    
    return gff_files

def detect_sequence_id(gff_file):
    """
    检测GFF文件中的序列ID（如ptg000001l）
    """
    try:
        with open(gff_file, 'r') as f:
            for line in f:
                if line.startswith('ptg'):
                    # 提取序列ID（第一个字段）
                    parts = line.strip().split('\t')
                    if len(parts) > 0:
                        return parts[0]
                elif line.startswith('#'):
                    # 在注释中查找序列信息
                    if 'sequence number 1' in line or 'name =' in line:
                        # 尝试提取序列名称
                        for part in line.split():
                            if part.startswith('ptg'):
                                return part.strip(',;')
    except Exception as e:
        print(f"警告: 检测序列ID时出错: {e}")

def convert_augustus_to_igv_gff3(input_file, output_file):
    """
    将AUGUSTUS GFF文件转换为IGV兼容的GFF3格式
    """
    print(f"正在转换: {input_file.name}")
    
    try:
        # 检测序列ID
        seq_id = detect_sequence_id(input_file)
        print(f"检测到序列ID: {seq_id}")
        
        with open(input_file, 'r') as f_in, open(output_file, 'w') as f_out:
            # 写入GFF3头部
            f_out.write("##gff-version 3\n")
            
            line_count = 0
            feature_count = 0
            gene_count = 0
            current_gene = None
            
            for line in f_in:
                line_count += 1
                stripped_line = line.strip()
                
                if not stripped_line:
                    continue  # 跳过空行
                
                # 处理注释行
                if stripped_line.startswith('#'):
                    # 保留基因边界注释
                    if stripped_line.startswith('# start gene') or stripped_line.startswith('# end gene'):
                        f_out.write(line)
                        if 'start gene' in stripped_line:
                            gene_count += 1
                    continue
                
                # 跳过错误和日志信息
                if (stripped_line.startswith('Error') or 
                    stripped_line.startswith('Delete') or
                    stripped_line.startswith('Forced') or
                    'HintGroup' in stripped_line):
                    continue
                
                # 处理特征行
                parts = stripped_line.split('\t')
                if len(parts) >= 9:
                    # 确保使用检测到的序列ID
                    if not parts[0].startswith('ptg'):
                        parts[0] = seq_id
                    
                    # 修正特征类型以符合GFF3标准
                    feature_type = parts[2]
                    if feature_type == "transcript":
                        parts[2] = "mRNA"
                    
                    # 处理属性字段（第9列）
                    attributes = parts[8]
                    if 'Parent=' in attributes or 'ID=' in attributes:
                        # 已经是类似GFF3的属性格式，直接使用
                        pass
                    else:
                        # AUGUSTUS旧格式，需要转换
                        # 这里可以添加更复杂的属性转换逻辑
                        pass
                    
                    # 写入清理后的行
                    cleaned_line = '\t'.join(parts) + '\n'
                    f_out.write(cleaned_line)
                    feature_count += 1
                else:
                    print(f"跳过无效行 ({len(parts)}字段): {stripped_line[:50]}...")
            
            print(f"转换统计: {line_count}行读取, {feature_count}个特征, {gene_count}个基因")
            return True
            
    except Exception as e:
        print(f"转换文件时出错: {e}")
        return False

def main():
    """主函数"""
    print("=" * 60)
    print("AUGUSTUS GFF 转 GFF3")
    print("=" * 60)
    
    # 自动查找GFF文件
    gff_files = find_augustus_gff_files()
    
    if not gff_files:
        print("在脚本同级目录及子目录中未找到GFF文件")
        print("请确保AUGUSTUS已运行并生成了.gff文件")
        return
    
    print(f"找到 {len(gff_files)} 个GFF文件:")
    for i, file in enumerate(gff_files, 1):
        print(f"  {i}. {file.name}")
    
    # 处理所有找到的文件
    success_count = 0
    for input_file in gff_files:
        # 生成输出文件名
        output_name = input_file.stem + "_igv.gff3"
        output_file = input_file.parent / output_name
        
        print(f"\n处理文件 {success_count + 1}/{len(gff_files)}:")
        
        # 检查输出文件是否已存在
        if output_file.exists():
            overwrite = input(f"输出文件已存在，是否覆盖? (y/n): ").strip().lower()
            if overwrite not in ['y', 'yes']:
                print(f"跳过文件: {input_file.name}")
                continue
        
        # 执行转换
        if convert_augustus_to_igv_gff3(input_file, output_file):
            success_count += 1
            print(f"✓ 成功转换: {output_file.name}")
        else:
            print(f"✗ 转换失败: {input_file.name}")
    
    # 总结
    print("\n" + "=" * 60)
    print(f"转换完成: {success_count}/{len(gff_files)} 个文件成功转换")

if __name__ == "__main__":
    main()