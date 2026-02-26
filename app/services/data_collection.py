import sys
import os
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from utils.utils import ChemBLTargetCollector, BindingDB_df, collect_data
import pandas as pd

class DataCollectionService:
    def __init__(self):
        self.collector = ChemBLTargetCollector()
    
    def collect_target_data(self, target_list: list, standard_type: str = "IC50", 
                           binding_db_folder: str = None):
        """
        타겟 유전자에 대한 ChEMBL과 BindingDB 데이터 수집
        """
        try:
            # ChEMBL 데이터 수집
            print(f"Collecting ChEMBL data for targets: {target_list}")
            chembl_df = self.collector.get_target_data(target_list, standard_type=standard_type)
            chembl_df = chembl_df.dropna()
            
            # BindingDB 데이터 수집 (폴더가 제공된 경우)
            bindingdb_df = pd.DataFrame()
            if binding_db_folder:
                print(f"Collecting BindingDB data from folder: {binding_db_folder}")
                bindingdb_df = BindingDB_df(binding_db_folder)
                bindingdb_df = bindingdb_df.dropna()
            
            # 결과 저장
            output_path = f'saved_data/IC50/{target_list[0]}_summary.xlsx'
            os.makedirs(os.path.dirname(output_path), exist_ok=True)
            
            with pd.ExcelWriter(output_path, mode='w', engine='openpyxl') as writer:
                chembl_df.to_excel(writer, sheet_name='ChemBL', index=False)
                if not bindingdb_df.empty:
                    bindingdb_df.to_excel(writer, sheet_name='BindingDB', index=False)
            
            return {
                'chembl_df': chembl_df,
                'bindingdb_df': bindingdb_df,
                'chembl_count': len(chembl_df),
                'bindingdb_count': len(bindingdb_df),
                'total_compounds': len(chembl_df) + len(bindingdb_df),
                'output_path': output_path
            }
        
        except Exception as e:
            print(f"Error in data collection: {str(e)}")
            raise e
