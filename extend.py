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
        """
        专门获取摘要信息的方法。这是一个扩展功能，不在基类中。
        """
        handle = Entrez.efetch(db="pubmed", id=str(pmid), retmode="xml")
        records = Entrez.read(handle)
        handle.close()
        article = records['PubmedArticle'][0]
        abstract = article['MedlineCitation']['Article'].get('Abstract', {}).get('AbstractText', [])
        if abstract:
            abstract = " ".join(abstract)
        else:
            abstract = "No abstract available"
        return abstract
        
    def process_data(self, input_file_path):
        with open(input_file_path, 'r', encoding='utf-8') as file:
            pmid_list = file.read().splitlines()

        data = {
            "PMID": [],
            "Title": [],
            "Title_CN": [],
            "Year": [],
            "Journal": [],
            "Crossref-Cites": [],
            "Abstract": [],
            "Abstract_CN": [],
            "URL": []
        }
        
        with ThreadPoolExecutor(max_workers=10) as executor:
            future_to_pmid = {executor.submit(self.fetch_and_process_pmid, pmid): pmid for pmid in pmid_list}
            
            for future in as_completed(future_to_pmid):
                result = future.result()
                if result is not None:
                    for key, value in result.items():
                        data[key].append(value)
        
        df = pd.DataFrame(data)
        df.sort_values(by=['Year', 'Crossref-Cites'], ascending=[False, False], inplace=True)

        excel_filename = f"ExtendedInfo_" + datetime.now().strftime('%Y%m%d_%H%M%S') + ".xlsx"
        excel_path = os.path.join(os.path.dirname(input_file_path), excel_filename)

        df.to_excel(excel_path, index=False, engine='openpyxl')
        workbook = openpyxl.load_workbook(excel_path)
        sheet = workbook.active
        self.format_excel(sheet)
        self.add_hyperlink_to_URL(sheet, url_column_index=9)  # Assuming URLs are in the 9th column
        workbook.save(excel_path)
        print(f"信息已经全部保存在 {excel_path}")

    def fetch_and_process_pmid(self, pmid):
        try:
            title, pub_date, journal, doi = super().fetch_pubmed_info(pmid)
            title_cn = translator(title)
            abstract = self.fetch_abstract_info(pmid)
            abstract_cn = translator(abstract) if abstract != "No abstract available" else "无摘要"
            citation_count = super().fetch_citation_count(doi) if doi else "No Citation Data"
            url = f'https://pubmed.ncbi.nlm.nih.gov/{pmid}/'
            return {
                "PMID": pmid,
                "Title": title,
                "Title_CN": title_cn,
                "Year": pub_date,
                "Journal": journal,
                "Crossref-Cites": citation_count,
                "Abstract": abstract,
                "Abstract_CN": abstract_cn,
                "URL": url
            }
        except Exception as e:
            print(f"Error processing PMID: {pmid}. {e}")
            return None

    def run(self, input_file_path):
        if not input_file_path:
            print("未提供输入文件路径。")
            return
        
        self.process_data(input_file_path)

    
if __name__ == "__main__":
    email = "263026560@qq.com"
    api_key = "c440fd575a5e0b26840b67cf3280c9b9bf08"
    processor_extended = PubMedDataProcessorExtended(email, api_key)
    processor_extended.run('D:\python code\塔吉科研\8\pmids_list\epilepsy [Title Abstract] AND neuroscience [Title .txt')
