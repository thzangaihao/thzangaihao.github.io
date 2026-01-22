'''
ab1转txt脚本

目录构型：
    自动识别脚本同级及子目录下的所有.ab1文件

本脚本依赖于第三方库

thz 2025/20/18
'''

from Bio import SeqIO
import os
import glob

def find_ab1():
    current_path = os.path.dirname(os.path.abspath(__file__))
    path = os.path.join(current_path, '**', f'*.ab1')
    ab1_list = glob.glob(path, recursive=True)
    print('='*50)
    print(f'在当前及子目录下共发现{len(ab1_list)}个ab1文件: ')
    for ab1 in ab1_list:
        print(ab1)

    return ab1_list

def convert_ab1(ab1_list):
    print('='*50)
    print('转换完整: ')
    for ab1 in ab1_list:
        seqence = SeqIO.read(ab1, 'abi')
        with open(f"{ab1}_output.txt", "w", encoding="utf-8") as file:
            file.write(str(seqence))
            print(f"{ab1}_output.txt")


def main():
    ab1_list = find_ab1()

    response = input("\n是否继续转化? (y/n): ").strip().lower()
    if response not in ['y', 'yes']:
        print("操作已取消")
        return
    convert_ab1(ab1_list)
    

main()