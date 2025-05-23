import pandas as pd
import os
from datetime import datetime

class PmidFilter:
    def __init__(self, file_path1, bad_journal_list):
        self.file_path1 = file_path1
        # 载入数据并将'Journal'列转换为小写
        self.df = pd.read_excel(file_path1, engine='openpyxl')
        self.df['Journal'] = self.df['Journal'].str.lower()  # 转换为小写
        self.exclude_bad_journals(bad_journal_list)  # 在初始化时执行，排除国际预警期刊

    def exclude_bad_journals(self, bad_journal_list):
        """排除国际预警期刊，这个操作是必须的"""
        # 将bad_journal_list中的元素转换为小写
        bad_journal_list_lower = [journal.lower() for journal in bad_journal_list]
        self.df_filtered = self.df[~self.df['Journal'].isin(bad_journal_list_lower)]

    def filter_by_journal(self, interested_journals):
        """根据期刊筛选文献，基于已排除国际预警期刊的数据集"""
        # 将感兴趣的期刊列表转换为小写
        interested_journals_lower = [journal.lower() for journal in interested_journals]
        self.df_interested_journals = self.df_filtered[self.df_filtered['Journal'].isin(interested_journals_lower)]
        return self.df_interested_journals

    def filter_by_citation(self, custom_filter_percentage, base_df=None):
        """根据被引次数筛选文献，基于已排除国际预警期刊的数据集"""
        if base_df is None:
            base_df = self.df_filtered
        custom_filter_multiplier = custom_filter_percentage / 100.0
        sum_F = base_df['Crossref-Cites'].sum()
        self.df_final = base_df[base_df['Crossref-Cites'] >= sum_F * custom_filter_multiplier]
        return self.df_final
    
    def output_result(self, df_to_output):
        """输出用户指定的DataFrame中包含PMID号的列到文件"""
        # 获取输入文件所在的目录
        directory = os.path.dirname(self.file_path1)
        # 创建新的文件名，包括路径
        sorted_pmid_filename = f"SortedPmid_{datetime.now().strftime('%Y%m%d_%H%M%S')}.txt"
        sorted_pmid_path = os.path.join(directory, sorted_pmid_filename)

        try:
            # 假设PMID号在名为"PMID"的列中
            df_to_output['PMID'].to_csv(sorted_pmid_path, index=False, header=False, encoding='utf-8')
            print(f"结果已保存到 {sorted_pmid_path}")
        except KeyError:
            # 如果找不到"PMID"列，则抛出错误
            print("DataFrame中没有找到'PMID'列，请检查列名。")

if __name__ == "__main__":
    file_path1 = r'D:\1\PubMedDataOutput.xlsx'  # Example file path
    bad_journal_list = ['CANCERS', 
                    'DIAGNOSTICS', 
                    'ENVIRONMENTAL SCIENCE AND POLLUTION RESEARCH', 
                    'FUEL', 
                    'JOURNAL OF CLINICAL MEDICINE', 
                    'JOURNAL OF PERSONALIZED MEDICINE', 
                    'RADIOLOGIA MEDICA',
                    'BIOENGINEERED', 
                    'CONNECTION SCIENCE', 
                    'MULTIMEDIA TOOLS AND APPLICATIONS', 
                    'PSYCHIATRIA DANUBINA', 
                    'JOURNAL OF BIOBASED MATERIALS AND BIOENERGY', 
                    'JOURNAL OF BIOMATERIALS AND TISSUE ENGINEERING', 
                    'JOURNAL OF BIOMEDICAL NANOTECHNOLOGY', 
                    'JOURNAL OF NANOELECTRONICS AND OPTOELECTRONICS', 
                    'JOURNAL OF SENSORS', 
                    'MATERIALS EXPRESS', 
                    'SCIENCE OF ADVANCED MATERIALS', 
                    'ALTERNATIVE THERAPIES IN HEALTH AND MEDICINE', 
                    'CMES-COMPUTER MODELING IN ENGINEERING & SCIENCES', 
                    'EXPERIMENTAL AND THERAPEUTIC MEDICINE', 
                    'FRONTIERS IN ENERGY RESEARCH', 
                    'MATHEMATICAL BIOSCIENCES AND ENGINEERING', 
                    'TROPICAL JOURNAL OF PHARMACEUTICAL RESEARCH']
    pmid_filter = PmidFilter(file_path1, bad_journal_list)

    # Example usage of the class for testing purposes
    print("Loaded data frame:")
    print(pmid_filter.df.head())  