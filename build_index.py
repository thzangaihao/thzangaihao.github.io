import glob
import os
import json
import re
import time

# 配置扫描根目录
base_dir = 'articles'  # 注意：这里最好用相对路径，不要加 ./ 前缀以免 glob 匹配出问题
# 输出文件路径
output_file = 'js/article_data.js'

articles = []

# 正则表达式 (保持不变)
# 提取 <title> 中的内容，并去掉后缀
title_pattern = re.compile(r'<title>(.*?) - 海钊知识港</title>', re.IGNORECASE)
# 提取 meta collection
collection_pattern = re.compile(r'<meta\s+name=["\']collection["\']\s+content=["\'](.*?)["\']', re.IGNORECASE)
# 提取 meta date
date_pattern = re.compile(r'<meta\s+name=["\']date["\']\s+content=["\'](.*?)["\']', re.IGNORECASE)
# 提取 meta description
desc_pattern = re.compile(r'<meta\s+name=["\']description["\']\s+content=["\'](.*?)["\']', re.IGNORECASE)

print(f"正在使用 glob 递归扫描 '{base_dir}' 目录下的所有 .html 文章...")

# [关键修改] 使用 glob 的递归模式
# ** 代表中间可以有任意多层目录
# recursive=True 是必须的
search_pattern = os.path.join(base_dir, '**', '*.html')
files = glob.glob(search_pattern, recursive=True)

count = 0

for file_path in files:
    file_name = os.path.basename(file_path)
    
    # [过滤] 排除索引页(index.html) 和 模板文件(模板.html)
    if 'index.html' in file_name or '模板' in file_name:
        continue

    try:
        with open(file_path, 'r', encoding='utf-8') as f:
            content = f.read()
            
            # 提取元数据
            coll_match = collection_pattern.search(content)
            
            # 只有包含 <meta name="collection" ...> 的文章才会被收录
            if coll_match:
                title_match = title_pattern.search(content)
                date_match = date_pattern.search(content)
                desc_match = desc_pattern.search(content)

                # 优先使用 title 标签，没找到就用文件名
                title = title_match.group(1) if title_match else file_name.replace('.html', '')
                
                # [关键修改] 路径标准化
                # 1. 计算相对于网站根目录的路径
                rel_path = os.path.relpath(file_path, '.')
                # 2. Windows 下路径是 \, 网页必须用 /，强制替换
                rel_path = rel_path.replace('\\', '/')

                article_info = {
                    'title': title,
                    'collection': coll_match.group(1), # 对应 HTML 中的 id="pcr"
                    'date': date_match.group(1) if date_match else 'Unknown',
                    'summary': desc_match.group(1) if desc_match else '暂无简介',
                    'path': rel_path  # 例如: "articles/dry_lab/data/GWAS.html"
                }
                
                articles.append(article_info)
                print(f"  [收录] {title} ({article_info['collection']})")
                count += 1
            else:
                # 这是一个 html 文件但不是文章（没有collection标记），跳过
                pass

    except Exception as e:
        print(f"  [读取错误] {file_path}: {e}")

# [排序] 按照文章标题排序 (如需按日期排序，把 x['title'] 改为 x['date'])
articles.sort(key=lambda x: x['title'])

# 生成 JS 文件
js_content = f"""// 自动生成的文件，请勿手动修改
// 生成时间: {time.strftime("%Y/%m/%d_%H:%M%S", time.localtime(time.time()))}
// 文章总数: {len(articles)}

window.ARTICLE_DATABASE = {json.dumps(articles, ensure_ascii=False, indent=4)};
"""

# 确保 js 目录存在
os.makedirs(os.path.dirname(output_file), exist_ok=True)

with open(output_file, 'w', encoding='utf-8') as f:
    f.write(js_content)

print(f"\n成功！共收录 {count} 篇文章。")
print(f"数据库已更新至: {output_file}")