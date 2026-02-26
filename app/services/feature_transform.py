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
