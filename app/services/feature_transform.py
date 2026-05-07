import sys
import os
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from utils.molecular_descriptor import (
    calc_descriptors,
    calc_specific_descriptors,
    remove_invalid_descriptors,
    remove_low_variance_descriptors,
    remove_correlated_descriptors
)
from utils.utils import convert_protein_data
import pandas as pd
import numpy as np
from rdkit import Chem
from rdkit.Chem import MACCSkeys, RDKFingerprint, rdMolDescriptors, rdFingerprintGenerator


def make_get_fingerprint(fp_type, n_bits=1024):
    """노트북과 동일한 fingerprint 생성 함수 (Pycaret/TabPFN 공통)"""
    def get_fingerprint(smiles):
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return [None] * (167 if fp_type == 'maccs' else n_bits)
        if fp_type == 'maccs':
            return list(MACCSkeys.GenMACCSKeys(mol))
        elif fp_type == 'ecfp4':
            gen = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=n_bits)
            return list(gen.GetFingerprintAsNumPy(mol))
        elif fp_type == 'ecfp6':
            gen = rdFingerprintGenerator.GetMorganGenerator(radius=3, fpSize=n_bits)
            return list(gen.GetFingerprintAsNumPy(mol))
        elif fp_type == 'fcfp4':
            invgen = rdFingerprintGenerator.GetMorganFeatureAtomInvGen()
            gen = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=n_bits, atomInvariantsGenerator=invgen)
            return list(gen.GetFingerprintAsNumPy(mol))
        elif fp_type == 'fcfp6':
            invgen = rdFingerprintGenerator.GetMorganFeatureAtomInvGen()
            gen = rdFingerprintGenerator.GetMorganGenerator(radius=3, fpSize=n_bits, atomInvariantsGenerator=invgen)
            return list(gen.GetFingerprintAsNumPy(mol))
        elif fp_type == 'rdkit':
            return list(RDKFingerprint(mol, fpSize=n_bits))
        elif fp_type == 'atompair':
            return list(rdMolDescriptors.GetHashedAtomPairFingerprintAsBitVect(mol, nBits=n_bits))
        elif fp_type == 'torsion':
            return list(rdMolDescriptors.GetHashedTopologicalTorsionFingerprintAsBitVect(mol, nBits=n_bits))
        else:
            raise ValueError(f'지원하지 않는 FP_TYPE: {fp_type}')
    return get_fingerprint


def get_fp_cols(fp_type, n_bits=1024):
    """
    Pycaret 노트북 기준 FP 컬럼명 반환
    - 비-MACCS: X1~X1024
    - MACCS: X2~X167 (X1은 항상 0이므로 제외)
    """
    fp_type = fp_type.lower()
    n_fp_bits = 167 if fp_type == 'maccs' else n_bits
    fp_cols_all = [f'X{i+1}' for i in range(n_fp_bits)]
    return fp_cols_all[1:] if fp_type == 'maccs' else fp_cols_all


class FeatureTransformService:

    def __init__(self):
        self.descriptor_selection = self._load_descriptor_selection()

    def _load_descriptor_selection(self):
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

    def transform_to_fingerprint(self, smiles_list: list, fp_type: str = "ecfp4",
                                  fp_size: int = 1024) -> pd.DataFrame:
        """
        SMILES → Fingerprint DataFrame (Pycaret 노트북 기준 컬럼명)
        Returns DataFrame with columns X1~X1024 (or X2~X167 for MACCS)
        """
        fp_type_lower = fp_type.lower()
        get_fp = make_get_fingerprint(fp_type_lower, fp_size)
        fp_data = pd.DataFrame(
            [get_fp(smi) for smi in smiles_list],
            columns=[f'X{i+1}' for i in range(167 if fp_type_lower == 'maccs' else fp_size)]
        )
        fp_cols = get_fp_cols(fp_type_lower, fp_size)
        return fp_data[fp_cols]

    def get_selected_descriptors(self, dataset_key: str):
        return self.descriptor_selection.get(dataset_key, [])

    def transform_to_fingerprint_with_descriptors(self, smiles_list: list,
                                                   fp_type: str = "ecfp4",
                                                   dataset_ratio: str = "5x",
                                                   ignore3D: bool = True) -> pd.DataFrame:
        """
        SMILES → Fingerprint + Descriptor DataFrame (Pycaret 노트북 동일)
        """
        fp_type_lower = fp_type.lower()

        # 1. Fingerprint
        fp_df = self.transform_to_fingerprint(smiles_list, fp_type=fp_type_lower)

        # 2. Descriptor selection
        ignore3D_str = "True" if ignore3D else "False"
        dataset_key = f"descriptors_filtered_FTO_training_{dataset_ratio}_ignore3D_{ignore3D_str}.csv"
        selected_descriptors = self.get_selected_descriptors(dataset_key)

        if not selected_descriptors:
            print(f"Warning: No selected descriptors found for {dataset_key}")
            return fp_df

        # 3. Calculate selected descriptors
        try:
            desc_df = calc_specific_descriptors(
                smiles_list,
                descriptor_list=selected_descriptors,
                ignore3D=ignore3D
            )
        except Exception as e:
            print(f"Error in calc_specific_descriptors: {e}")
            return fp_df

        if 'canonical_SMILES' in desc_df.columns:
            desc_df = desc_df.drop('canonical_SMILES', axis=1)

        # 4. Combine
        result_df = pd.concat([fp_df, desc_df], axis=1)
        n_fp = fp_df.shape[1]
        n_md = desc_df.shape[1]
        print(f"Total features: {result_df.shape[1]} ({n_fp} fingerprint + {n_md} descriptors)")
        return result_df

    def transform_to_descriptors(self, smiles_list: list, descriptor_type: str = "MORDRED_2D",
                                  descriptor_list: list = None):
        ignore3D = descriptor_type != "MORDRED_3D"
        if descriptor_list:
            desc_df = calc_specific_descriptors(smiles_list, descriptor_list=descriptor_list,
                                                ignore3D=ignore3D)
        else:
            desc_df = calc_descriptors(smiles_list, ignore3D=ignore3D)

        desc_df, _ = remove_invalid_descriptors(desc_df)
        desc_df_numeric = desc_df.select_dtypes(include=[np.number])
        desc_df_cleaned, _ = remove_low_variance_descriptors(desc_df_numeric)
        desc_df_final, _ = remove_correlated_descriptors(desc_df_cleaned)
        return desc_df_final

    def prepare_training_data(self, chembl_df: pd.DataFrame, bindingdb_df: pd.DataFrame,
                               protein_name: str, dataset_type: str = "fewshot",
                               pos_threshold: float = 10000, neg_threshold: float = 20000):
        return convert_protein_data(
            chembl_df=chembl_df,
            BDB_df=bindingdb_df,
            protein_name=protein_name,
            dataset_type=dataset_type,
            pos_threshold=pos_threshold,
            neg_threshold=neg_threshold
        )
