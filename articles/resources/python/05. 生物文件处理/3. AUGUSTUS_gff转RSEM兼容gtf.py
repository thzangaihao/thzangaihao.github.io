import os
import glob
import sys
import pandas as pd

# This script converts an AUGUSTUS GFF3 file into a simplified gene-level GTF
# suitable for RSEM. Gene-level only (no transcript isoforms).

############################################################
# Utility functions copied from the user's cite.py pattern
############################################################
def current_path_function():
    if getattr(sys, 'frozen', False):
        return os.path.dirname(sys.executable)
    else:
        return os.path.dirname(os.path.abspath(__file__))


def find_file(file_type):
    current_path = current_path_function()
    search_path = os.path.join(current_path, '**', f'*{file_type}')
    file_list = glob.glob(search_path, recursive=True)

    print('=' * 50)
    if len(file_list) >= 1:
        print(f'在当前目录及子目录下共发现 {len(file_list)} 个 {file_type} 文件:')
        for f in file_list:
            print(f)
        return file_list
    else:
        response = input(f'未找到 {file_type} 文件，是否搜索父级目录? (y/n): ').strip().lower()
        if response not in ['y', 'yes']:
            print('操作已取消')
            sys.exit()
        parent_path = os.path.dirname(current_path)
        parent_search = os.path.join(parent_path, '**', f'*{file_type}')
        parent_list = glob.glob(parent_search, recursive=True)

        if len(parent_list) >= 1:
            print(f'在父级目录找到 {len(parent_list)} 个 {file_type} 文件:')
            for f in parent_list:
                print(f)
            return parent_list
        else:
            print('仍未找到文件，程序退出')
            sys.exit()


def choose_file(file_list):
    print('=' * 50)
    print(f'共发现 {len(file_list)} 个文件:')
    for i, f in enumerate(file_list, start=1):
        print(f' [{i}] {f}')
    idx = int(input('请选择目标文件编号: ')) - 1
    return file_list[idx]

############################################################
# Core conversion
############################################################
def convert_augustus_gff_to_gtf(gff_file, out_gtf):
    print(f"正在处理文件: {gff_file}")
    output = []

    with open(gff_file, 'r', encoding='utf-8', errors='replace') as f:
        for line in f:
            if line.startswith('#') or line.strip() == '':
                continue

            parts = line.strip().split('\t')
            if len(parts) < 9:
                continue

            seqid, source, feature, start, end, score, strand, frame, attrs = parts

            # We only keep gene-level entries
            if feature != 'gene':
                continue

            # Parse ID from attributes
            id_dict = {}
            for field in attrs.split(';'):
                if '=' in field:
                    k, v = field.split('=', 1)
                    id_dict[k] = v

            gene_id = id_dict.get('ID', None)
            if not gene_id:
                continue

            # Create GTF-like line: RSEM needs gene_id and transcript_id.
            # For gene-level quantification, we can set transcript_id = gene_id.
            gtf_attrs = f'gene_id "{gene_id}"; transcript_id "{gene_id}";'

            gtf_line = '\t'.join([
                seqid,
                'AUGUSTUS',
                'exon',  # RSEM requires exon features
                start,
                end,
                score if score != '.' else '0',
                strand,
                frame if frame != '.' else '0',
                gtf_attrs
            ])

            output.append(gtf_line)

    with open(out_gtf, 'w') as out:
        for l in output:
            out.write(l + '\n')

    print(f"转换完成！GTF 已保存到: {out_gtf}")

############################################################
# Main
############################################################
def main():
    gff_files = find_file('.gff')
    selected = choose_file(gff_files)

    out_name = input('请输入输出 GTF 文件名（无需扩展名）: ').strip()
    current_path = current_path_function()
    out_gtf = os.path.join(current_path, f'{out_name}.gtf')

    convert_augustus_gff_to_gtf(selected, out_gtf)


if __name__ == '__main__':
    main()
