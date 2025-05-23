from base import PubMedDataProcessorBase
from Bio import Entrez
import pandas as pd
import openpyxl
from datetime import datetime
import os
from concurrent.futures import ThreadPoolExecutor, as_completed
from translator import translator

class PubMedDataProcessorExtended(PubMedDataProcessorBase):
    def fetch_abstract_info(self, pmid):
        """获取摘要信息"""
        try:
            handle = Entrez.efetch(db="pubmed", id=str(pmid), retmode="xml")
            records = Entrez.read(handle)
            handle.close()
            article = records['PubmedArticle'][0]
            abstract = article['MedlineCitation']['Article']\
                              .get('Abstract', {})\
                              .get('AbstractText', [])
            return " ".join(abstract) if abstract else "No abstract available"
        except Exception as e:
            print(f"[ERROR] PMID {pmid} 抓取摘要失败: {e}")
            return "No abstract available"

    def process_data(self, input_file_path):
        with open(input_file_path, 'r', encoding='utf-8') as f:
            pmid_list = [line.strip() for line in f if line.strip()]
        total = len(pmid_list)
        print(f"开始处理 {total} 条 PMID …\n")

        data = {
            "PMID": [], "Title": [], "Title_CN": [],
            "Year": [], "Journal": [], "Crossref-Cites": [],
            "Abstract": [], "Abstract_CN": [], "URL": []
        }

        completed = 0
        with ThreadPoolExecutor(max_workers=20) as executor:
            future_to_pmid = {
                executor.submit(self.fetch_and_process_pmid, pmid): pmid
                for pmid in pmid_list
            }
            for future in as_completed(future_to_pmid):
                pmid = future_to_pmid[future]
                result = None
                try:
                    result = future.result()
                except Exception as e:
                    print(f"[ERROR] PMID {pmid} 处理线程异常: {e}")
                if result:
                    for key, val in result.items():
                        data[key].append(val)
                completed += 1
                print(f"[{completed}/{total}] 已处理 PMID {pmid}", end="\r", flush=True)

        # 构建 DataFrame 并保存
        df = pd.DataFrame(data)
        df.sort_values(
            by=['Year', 'Crossref-Cites'],
            ascending=[False, False],
            inplace=True
        )

        fname = f"ExtendedInfo_{datetime.now():%Y%m%d_%H%M%S}.xlsx"
        out_path = os.path.join(os.path.dirname(input_file_path), fname)
        df.to_excel(out_path, index=False, engine='openpyxl')

        wb = openpyxl.load_workbook(out_path)
        sheet = wb.active
        self.format_excel(sheet)
        # 假设 URL 列是第 9 列（从 1 开始计数）
        self.add_hyperlink_to_URL(sheet, url_column_index=9)
        wb.save(out_path)

        print(f"\n所有 PMID 处理完毕，文件已保存至：\n{out_path}")

    def fetch_and_process_pmid(self, pmid):
        try:
            title, pub_date, journal, doi = super().fetch_pubmed_info(pmid)
            title_cn = translator(title)
            abstract = self.fetch_abstract_info(pmid)
            abstract_cn = (
                translator(abstract)
                if abstract != "No abstract available"
                else "无摘要"
            )
            cites = super().fetch_citation_count(doi) if doi else "No Citation Data"
            url = f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/"
            return {
                "PMID": pmid,
                "Title": title,
                "Title_CN": title_cn,
                "Year": pub_date,
                "Journal": journal,
                "Crossref-Cites": cites,
                "Abstract": abstract,
                "Abstract_CN": abstract_cn,
                "URL": url,
            }
        except Exception as e:
            print(f"[ERROR] PMID {pmid} 处理失败: {e}")
            return None

    def run(self, input_file_path):
        if not input_file_path:
            print("未提供输入文件路径。")
            return
        self.process_data(input_file_path)


if __name__ == "__main__":
    email = ""
    api_key = ""
    processor_extended = PubMedDataProcessorExtended(email, api_key)
    processor_extended.run('')  # 需要提供正确的输入文件路径