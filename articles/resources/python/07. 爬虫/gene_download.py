import pandas as pd
import os
from Bio import Entrez, SeqIO
import time

def fetch_sequences():
    # 1. 路径逻辑：强制切换到脚本所在目录
    script_dir = os.path.dirname(os.path.abspath(__file__))
    os.chdir(script_dir)
    
    input_csv = 'accession.csv'
    output_fasta = 'all_sequences.fasta'
    
    # 2. 配置 NCBI Entrez (请务必填写你的邮箱)
    Entrez.email = "your_email@example.com" # NCBI 要求提供邮箱以识别请求者
    
    if not os.path.exists(input_csv):
        print(f"错误: 找不到文件 {input_csv}")
        return

    # 3. 读取 Accession IDs
    try:
        df = pd.read_csv(input_csv)
        # 清洗列名，防止空格干扰
        df.columns = df.columns.str.strip()
        
        if 'Accession' not in df.columns:
            print(f"错误: CSV文件中未找到 'Accession' 列。请检查表头。")
            return
            
        accessions = df['Accession'].dropna().unique().tolist()
        print(f"从表格中读取到 {len(accessions)} 个唯一的 Accession ID。")
    except Exception as e:
        print(f"读取 CSV 时出错: {e}")
        return

    # 4. 批量下载并汇总
    print("开始从 NCBI 下载序列，请稍候...")
    
    # 为了避免请求过于频繁被封IP，我们分批下载（每批20个）
    batch_size = 20
    count = 0
    
    with open(output_fasta, "w") as out_handle:
        for i in range(0, len(accessions), batch_size):
            batch = accessions[i:i+batch_size]
            try:
                # 使用 efetch 获取 FASTA 格式
                # db="protein" 如果你下载的是蛋白序列；如果是核酸请改为 db="nucleotide"
                handle = Entrez.efetch(
                    db="protein", 
                    id=",".join(batch), 
                    rettype="fasta", 
                    retmode="text"
                )
                fasta_data = handle.read()
                out_handle.write(fasta_data)
                handle.close()
                
                count += len(batch)
                print(f"已进度: {count}/{len(accessions)}...")
                
                # 遵守 NCBI 频率限制：每秒不超过3个请求
                time.sleep(0.5) 
                
            except Exception as e:
                print(f"批次 {i//batch_size + 1} 下载失败: {e}")
                continue

    print(f"\n任务完成！所有序列已保存至: {os.path.join(script_dir, output_fasta)}")

if __name__ == "__main__":
    fetch_sequences()