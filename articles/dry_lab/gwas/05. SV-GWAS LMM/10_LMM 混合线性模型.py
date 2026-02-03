#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import sys
import glob
import subprocess
import shutil

# ==============================================================================
# 0. 核心配置参数
# ==============================================================================
FINAL_OUTPUT_DIR = "10_LMM_GWAS_Results"
GEMMA_EXEC = "gemma"

# ==============================================================================
# 1. Cite2 交互逻辑模块 (复用版)
# ==============================================================================
def get_base_dir():
    if getattr(sys, 'frozen', False):
        return os.path.dirname(sys.executable)
    return os.path.dirname(os.path.abspath(__file__))

def find_files(exts, path=None):
    if path is None: path = get_base_dir()
    if isinstance(exts, str): exts = [exts]
    all_files = []
    for ext in exts:
        if not ext.startswith('.'): ext = '.' + ext
        search_pattern = os.path.join(path, '**', f'*{ext}')
        all_files.extend(glob.glob(search_pattern, recursive=True))
    # 排除输出目录
    return sorted(list(set([f for f in all_files if FINAL_OUTPUT_DIR not in f and "output" not in f])))

def choose_file(files, desc="文件"):
    if not files:
        print(f"[提示] 未找到任何 {desc}")
        return []

    print(f"\n--- 找到 {len(files)} 个 {desc}: ---")
    limit = 15
    for i, f in enumerate(files[:limit], 1):
        print(f"  [{i}] {os.path.basename(f)}  ({os.path.relpath(f, get_base_dir())})")
    if len(files) > limit:
        print(f"  ... (共 {len(files)} 个)")

    while True:
        try:
            prompt = f"\n请输入 {desc} 编号 (支持 1,2 | all): "
            user_input = input(prompt).strip().lower()
            if not user_input: continue

            selected_indices = set()
            if user_input in ['all', 'a']:
                selected_indices = set(range(len(files)))
            else:
                parts = user_input.split(',')
                for part in parts:
                    if '-' in part:
                        s, e = map(int, part.split('-'))
                        selected_indices.update(range(s-1, e))
                    else:
                        selected_indices.add(int(part)-1)
            
            selected_files = [files[i] for i in sorted(selected_indices) if 0 <= i < len(files)]
            if selected_files: return selected_files
        except: pass

def make_sure(action_name):
    response = input(f"\n[选项] 是否{action_name}? (y/n): ").strip().lower()
    return response in ['y', 'yes']

# ==============================================================================
# 2. 核心处理模块 (动态模型构建)
# ==============================================================================
def run_gwas_batch(bfile_path, trait_files, kinship_path=None, pca_path=None):
    base_dir = get_base_dir()
    final_dir = os.path.join(base_dir, FINAL_OUTPUT_DIR)
    gemma_default_out_dir = os.path.join(base_dir, "output")

    if not os.path.exists(gemma_default_out_dir): os.makedirs(gemma_default_out_dir)
    if not os.path.exists(final_dir): os.makedirs(final_dir)

    bfile_prefix = os.path.splitext(bfile_path)[0]

    # --- 确定统计模型 ---
    # 如果有 Kinship -> LMM (-lmm 1)
    # 如果无 Kinship -> LM  (-lm 1)
    if kinship_path:
        model_type = "LMM (混合线性模型)"
        model_flag = ["-lmm", "1"]
        print(f"\n[模型确认] 检测到 Kinship，将使用 {model_type}")
    else:
        model_type = "LM (普通线性模型)"
        model_flag = ["-lm", "1"]
        print(f"\n[模型确认] 未检测到 Kinship，自动切换为 {model_type}")

    if pca_path:
        print(f"[协变量] 已启用 PCA 校正")

    print(f"\n--- 开始批量 GWAS 分析 (共 {len(trait_files)} 个性状) ---")
    
    for idx, trait_file in enumerate(trait_files, 1):
        trait_name = os.path.splitext(os.path.basename(trait_file))[0]
        
        # 构造输出文件名，带上模型标记，防止混淆
        # 例如: PlantHeight_LMM_PCA
        suffix = "LMM" if kinship_path else "LM"
        if pca_path: suffix += "_PCA"
        output_name = f"{trait_name}_{suffix}"
        
        print(f"\n>>> [{idx}/{len(trait_files)}] 分析性状: {trait_name} (模式: {suffix})")
        
        # --- 动态构建命令 ---
        # 基础命令
        cmd = [GEMMA_EXEC, "-bfile", bfile_prefix, "-p", trait_file, "-o", output_name]
        
        # 加上模型参数 (-lmm 1 或 -lm 1)
        cmd.extend(model_flag)
        
        # 加上 Kinship (-k)
        if kinship_path:
            cmd.extend(["-k", kinship_path])
            
        # 加上 PCA/协变量 (-c)
        if pca_path:
            cmd.extend(["-c", pca_path])
        
        print(f"执行命令: {' '.join(cmd)}")
        
        try:
            subprocess.run(cmd, check=True)
            
            # 移动结果
            src_assoc = os.path.join(gemma_default_out_dir, f"{output_name}.assoc.txt")
            src_log   = os.path.join(gemma_default_out_dir, f"{output_name}.log.txt")
            
            dst_assoc = os.path.join(final_dir, f"{output_name}.assoc.txt")
            dst_log   = os.path.join(final_dir, f"{output_name}.log.txt")
            
            if os.path.exists(src_assoc):
                shutil.move(src_assoc, dst_assoc)
                print(f"  [成功] 结果已保存: {os.path.basename(dst_assoc)}")
            else:
                print(f"  [警告] 结果文件缺失。")
            
            if os.path.exists(src_log): shutil.move(src_log, dst_log)

        except subprocess.CalledProcessError as e:
            print(f"  [失败] 运行出错: {e}")

    print("\n" + "="*50)
    print(f"所有分析完成！")
    print(f"结果存放于: {final_dir}")
    print("="*50)

# ==============================================================================
# 主函数
# ==============================================================================
def main():
    print("==============================================")
    print("   Step 10: 全能 GWAS 分析工具 (LMM/LM + PCA)")
    print("==============================================")

    # 1. 必选: 基因型文件
    print("\n>>> 第一步: 选择基因型数据 (.fam)")
    fam_files = find_files('.fam')
    sel_fam = choose_file(fam_files, "基因型文件")
    if not sel_fam: return
    bfile_path = sel_fam[0]

    # 2. 可选: Kinship 矩阵
    kinship_path = None
    if make_sure("使用 Kinship 矩阵 (推荐: 是 -> LMM模型)"):
        kin_files = find_files('.cXX.txt')
        sel_kin = choose_file(kin_files, "Kinship 矩阵")
        if sel_kin: kinship_path = sel_kin[0]
    else:
        print("  -> 已跳过 Kinship。模式将切换为普通线性模型 (LM)。")

    # 3. 可选: PCA 协变量
    pca_path = None
    if make_sure("使用 PCA 协变量文件 (推荐: 是 -> 校正群体分层)"):
        # 自动搜索包含 covar 的 txt 文件
        cov_files = find_files(['.txt'])
        cov_candidates = [f for f in cov_files if "covar" in os.path.basename(f) or "pca" in os.path.basename(f)]
        sel_pca = choose_file(cov_candidates, "PCA/协变量文件")
        if sel_pca: pca_path = sel_pca[0]
    else:
        print("  -> 已跳过 PCA 校正。")

    # 4. 必选: 性状文件
    print("\n>>> 第四步: 选择性状文件 (支持多选)")
    txt_files = find_files('.txt')
    # 排除掉 kinship, covar, log 等非性状文件
    trait_candidates = [
        f for f in txt_files 
        if "covar" not in f and "kinship" not in f and "cXX" not in f and "log" not in f
    ]
    
    sel_traits = choose_file(trait_candidates, "性状文件")
    if not sel_traits: return

    # 5. 执行
    print("\n" + "-"*30)
    print(f"配置确认:")
    print(f"  基因型 : {os.path.basename(bfile_path)}")
    print(f"  Kinship: {'[启用] ' + os.path.basename(kinship_path) if kinship_path else '[禁用]'}")
    print(f"  PCA    : {'[启用] ' + os.path.basename(pca_path) if pca_path else '[禁用]'}")
    print(f"  性状数 : {len(sel_traits)} 个")
    print("-" * 30)
    
    if make_sure("开始批量分析"):
        run_gwas_batch(bfile_path, sel_traits, kinship_path, pca_path)

if __name__ == "__main__":
    main()