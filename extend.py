from base import PubMedDataProcessorBase
from Bio import Entrez
import pandas as pd
import openpyxl
from datetime import datetime
import os
from concurrent.futures import ThreadPoolExecutor, as_completed
from translator import translator
from tqdm import tqdm  

class PubMedDataProcessorExtended(PubMedDataProcessorBase):
    def fetch_abstract_info(self, pmid):
        """
        获取摘要信息
        """
        try:
            handle = Entrez.efetch(db="pubmed", id=str(pmid), retmode="xml")
            records = Entrez.read(handle)
            handle.close()
            article = records['PubmedArticle'][0]
            abstract = article['MedlineCitation']['Article'].get('Abstract', {}).get('AbstractText', [])
            return " ".join(abstract) if abstract else "No abstract available"
        except Exception as e:
            print(f"Error fetching abstract for PMID {pmid}: {e}")
            return "No abstract available"

    def process_data(self, input_file_path):
        with open(input_file_path, 'r', encoding='utf-8') as file:
            pmid_list = file.read().splitlines()

        data = {"PMID": [], "Title": [], "Title_CN": [], "Year": [], "Journal": [],
                "Crossref-Cites": [], "Abstract": [], "Abstract_CN": [], "URL": []}

        with ThreadPoolExecutor(max_workers=10) as executor, tqdm(total=len(pmid_list), desc="Processing PMIDs") as pbar:
            future_to_pmid = {executor.submit(self.fetch_and_process_pmid, pmid): pmid for pmid in pmid_list}
            for future in as_completed(future_to_pmid):
                result = future.result()
                if result:
                    for key, value in result.items():
                        data[key].append(value)
                pbar.update(1)

        df = pd.DataFrame(data)
        df.sort_values(by=['Year', 'Crossref-Cites'], ascending=[False, False], inplace=True)

        excel_filename = f"ExtendedInfo_{datetime.now().strftime('%Y%m%d_%H%M%S')}.xlsx"
        excel_path = os.path.join(os.path.dirname(input_file_path), excel_filename)

        df.to_excel(excel_path, index=False, engine='openpyxl')
        workbook = openpyxl.load_workbook(excel_path)
        sheet = workbook.active
        self.format_excel(sheet)
        self.add_hyperlink_to_URL(sheet, url_column_index=9)
        workbook.save(excel_path)
        print(f"信息已保存至 {excel_path}")

    def fetch_and_process_pmid(self, pmid):
        try:
            title, pub_date, journal, doi = super().fetch_pubmed_info(pmid)
            title_cn = translator(title)
            abstract = self.fetch_abstract_info(pmid)
            abstract_cn = translator(abstract) if abstract != "No abstract available" else "无摘要"
            citation_count = super().fetch_citation_count(doi) if doi else "No Citation Data"
            url = f'https://pubmed.ncbi.nlm.nih.gov/{pmid}/'
            return {"PMID": pmid, "Title": title, "Title_CN": title_cn, "Year": pub_date,
                    "Journal": journal, "Crossref-Cites": citation_count,
                    "Abstract": abstract, "Abstract_CN": abstract_cn, "URL": url}
        except Exception as e:
            print(f"Error processing PMID {pmid}: {e}")
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