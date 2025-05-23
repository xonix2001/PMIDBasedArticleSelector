import os
import re
import requests

def get_pubmed_pmids(query, count):
    base_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"  
    params = {
        "db": "pubmed",        # 指定数据库为PubMed
        "term": query,         # 查询关键词
        "retmax": count,       # 返回的PMID最大数量
        "sort": "relevance",   # 按相关性排序
        "retmode": "json"      # 返回结果格式为JSON
    }
    try:
        response = requests.get(base_url, params=params)  
        response.raise_for_status()  # 如果状态码非200，抛出异常
        data = response.json()  # 将返回的JSON解析为字典
        
        # 获取PMID列表和检索结果总数
        pmid_list = data.get("esearchresult", {}).get("idlist", [])
        total_count = int(data.get("esearchresult", {}).get("count", 0))
        return pmid_list, total_count
    except requests.exceptions.RequestException as e:
        print(f"请求失败：{e}")
        return [], 0

def sanitize_filename(query):
    sanitized = re.sub(r'[\\/:*?"<>|]', ' ', query)
    sanitized = re.sub(r'\s+', ' ', sanitized).strip()
    return sanitized[:50] if len(sanitized) > 50 else sanitized

def save_pmids_to_file(pmids, query, output_folder):
    os.makedirs(output_folder, exist_ok=True)
    filename = sanitize_filename(query)
    file_path = os.path.join(output_folder, f"{filename}.txt")
    try:
        with open(file_path, 'w', encoding='utf-8') as f:
            for pmid in pmids:
                f.write(f"{pmid}\n")
        print(f"PMID列表已保存到：{file_path}")
    except Exception as e:
        print(f"保存文件失败：{e}")

if __name__ == "__main__":
#——————————————————————Enter here————————————————————————————————

    query = '''
    #查询关键词（多行字符串也可以）
    '''
    count = 100  # 希望返回的PMID数量

    output_folder = "#path"# 指定保存PMID列表的文件夹路径

#————————————————————————————————————————————————————————————————
    pmids, total_count = get_pubmed_pmids(query, count)
    if pmids:
        print(f"前{min(count, total_count)}/{total_count}篇相关文章的PMID：")
        for pmid in pmids:
            print(pmid)
        save_pmids_to_file(pmids, query, output_folder)
    else:
        print("未找到相关PMID。")
