import sys
import os
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from utils.molecular_descriptor import (
    calc_descriptors, 
    calc_specific_descriptors,
    smiles_to_fingerprint,
    remove_invalid_descriptors,
    remove_low_variance_descriptors,
    remove_correlated_descriptors
)
from utils.utils import convert_protein_data
import pandas as pd
import numpy as np
from rdkit import Chem
from rdkit.Chem import MACCSkeys, DataStructs

class FeatureTransformService:
    
    def __init__(self):
        # descriptor_selection.csv 로드
        self.descriptor_selection = self._load_descriptor_selection()
    
    def _load_descriptor_selection(self):
        """descriptor_selection.csv 파일 로드"""
        csv_path = "descriptor_selection.csv"
        if os.path.exists(csv_path):
            df = pd.read_csv(csv_path)
            file_md_list = {}
            for column in df.columns:
                cols = df[column].dropna().tolist()
                if cols:
                    file_md_list[column] = cols
            return file_md_list
        return {}
    
    def transform_to_fingerprint(self, smiles_list: list, fp_type: str = "ECFP4", 
                                 fp_size: int = 1024, radius: int = 2):
        """
        SMILES를 Fingerprint로 변환
        """
        fingerprints = []
        
        if fp_type == "ECFP4" or fp_type == "MORGAN":
            for smiles in smiles_list:
                fp = smiles_to_fingerprint(smiles, fp_size=fp_size, radius=radius)
                fingerprints.append(fp)
        
        elif fp_type == "MACCS":
            for smiles in smiles_list:
                mol = Chem.MolFromSmiles(smiles)
                if mol is None:
                    fingerprints.append(np.zeros(167))
                else:
                    fp = MACCSkeys.GenMACCSKeys(mol)
                    arr = np.zeros(167)
                    DataStructs.ConvertToNumpyArray(fp, arr)
                    fingerprints.append(arr)
        
        return np.array(fingerprints)
    
    def get_selected_descriptors(self, dataset_key: str):
        """
        특정 데이터셋에 대한 선택된 Descriptor 리스트 반환
        
        dataset_key 예시:
        - "descriptors_filtered_FTO_training_5x_ignore3D_False.csv"
        - "descriptors_filtered_FTO_training_10x_ignore3D_True.csv"
        """
        return self.descriptor_selection.get(dataset_key, [])
    
    def transform_to_fingerprint_with_descriptors(self, smiles_list: list, 
                                                   fp_type: str = "ECFP4",
                                                   dataset_ratio: str = "5x",
                                                   ignore3D: bool = True):
        """
        SMILES를 Fingerprint + 선택된 Descriptor로 변환
        (원본 노트북 코드 방식)
        
        Args:
            smiles_list: SMILES 문자열 리스트
            fp_type: "ECFP4", "MACCS", "MORGAN"
            dataset_ratio: "5x", "10x", "20x"
            ignore3D: True (2D descriptor), False (3D descriptor)
        
        Returns:
            DataFrame with fingerprint + selected descriptors
        """
        # 1. Fingerprint 생성
        fp_array = self.transform_to_fingerprint(smiles_list, fp_type=fp_type)
        
        # Fingerprint 크기 결정
        if fp_type == "MACCS":
            fp_size = 167
            fp_cols = [f'X{i}' for i in range(167)]
        elif fp_type == "ECFP4" or fp_type == "MORGAN":
            fp_size = 1024
            fp_cols = [f'X{i}' for i in range(1024)]
        
        fp_df = pd.DataFrame(fp_array, columns=fp_cols)
        
        # 2. 선택된 Descriptor 가져오기
        ignore3D_str = "True" if ignore3D else "False"
        dataset_key = f"descriptors_filtered_FTO_training_{dataset_ratio}_ignore3D_{ignore3D_str}.csv"
        selected_descriptors = self.get_selected_descriptors(dataset_key)
        
        if not selected_descriptors:
            print(f"Warning: No selected descriptors found for {dataset_key}")
            return fp_df
        
        # 3. 선택된 Descriptor만 계산
        try:
            desc_df = calc_specific_descriptors(
                smiles_list, 
                descriptor_list=selected_descriptors,
                ignore3D=ignore3D
            )
        except Exception as e:
            print(f"Error in calc_specific_descriptors: {e}")
            print(f"Selected descriptors: {selected_descriptors}")
            # Descriptor 계산 실패 시 Fingerprint만 반환
            return fp_df
        
        # SMILES 컬럼 제거
        if 'canonical_SMILES' in desc_df.columns:
            desc_df = desc_df.drop('canonical_SMILES', axis=1)
        
        # 4. Fingerprint + Descriptor 결합
        result_df = pd.concat([fp_df, desc_df], axis=1)
        
        print(f"Total features: {len(result_df.columns)} ({fp_size} fingerprint + {len(desc_df.columns)} descriptors)")
        
        return result_df
    
    def transform_to_descriptors(self, smiles_list: list, descriptor_type: str = "MORDRED_2D",
                                 descriptor_list: list = None):
        """
        SMILES를 Molecular Descriptors로 변환
        """
        ignore3D = True if descriptor_type == "MORDRED_2D" else False
        
        if descriptor_list:
            desc_df = calc_specific_descriptors(smiles_list, descriptor_list=descriptor_list, 
                                                ignore3D=ignore3D)
        else:
            desc_df = calc_descriptors(smiles_list, ignore3D=ignore3D)
        
        # 데이터 정제
        desc_df, _ = remove_invalid_descriptors(desc_df)
        desc_df_numeric = desc_df.select_dtypes(include=[np.number])
        desc_df_cleaned, _ = remove_low_variance_descriptors(desc_df_numeric)
        desc_df_final, _ = remove_correlated_descriptors(desc_df_cleaned)
        
        return desc_df_final
    
    def prepare_training_data(self, chembl_df: pd.DataFrame, bindingdb_df: pd.DataFrame,
                             protein_name: str, dataset_type: str = "fewshot",
                             pos_threshold: float = 10000, neg_threshold: float = 20000):
        """
        학습용 데이터 준비
        """
        result_df = convert_protein_data(
            chembl_df=chembl_df,
            BDB_df=bindingdb_df,
            protein_name=protein_name,
            dataset_type=dataset_type,
            pos_threshold=pos_threshold,
            neg_threshold=neg_threshold
        )
        
        return result_df
