"""
분자 특성 변환 모듈 - Fingerprint 및 Molecular Descriptor 계산
"""

import numpy as np
import pandas as pd
from typing import List, Optional, Union
from rdkit import Chem
from rdkit.Chem import (
    AllChem, 
    Descriptors, 
    DataStructs, 
    rdMolDescriptors,
    MACCSkeys
)
from rdkit.ML.Descriptors import MoleculeDescriptors


class MolecularFeatureExtractor:
    """분자 특성 추출 클래스"""
    
    FINGERPRINT_TYPES = {
        'ECFP4': 'Morgan Fingerprint (ECFP4)',
        'ECFP6': 'Morgan Fingerprint (ECFP6)',
        'MACCS': 'MACCS Keys',
        'AtomPair': 'Atom Pair Fingerprint',
        'TopologicalTorsion': 'Topological Torsion',
        'RDKit': 'RDKit Fingerprint'
    }
    
    def __init__(self, fp_type: str = 'ECFP4', fp_size: int = 1024):
        """
        Args:
            fp_type: Fingerprint 타입
            fp_size: Fingerprint 크기 (MACCS는 고정 167)
        """
        self.fp_type = fp_type
        self.fp_size = fp_size if fp_type != 'MACCS' else 167
        
    def smiles_to_fingerprint(self, smiles: str) -> np.ndarray:
        """SMILES를 Fingerprint로 변환"""
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return np.zeros(self.fp_size)
        
        try:
            if self.fp_type == 'ECFP4':
                fp = rdMolDescriptors.GetMorganFingerprintAsBitVect(
                    mol, radius=2, nBits=self.fp_size
                )
            elif self.fp_type == 'ECFP6':
                fp = rdMolDescriptors.GetMorganFingerprintAsBitVect(
                    mol, radius=3, nBits=self.fp_size
                )
            elif self.fp_type == 'MACCS':
                fp = MACCSkeys.GenMACCSKeys(mol)
            elif self.fp_type == 'AtomPair':
                fp = rdMolDescriptors.GetHashedAtomPairFingerprintAsBitVect(
                    mol, nBits=self.fp_size
                )
            elif self.fp_type == 'TopologicalTorsion':
                fp = rdMolDescriptors.GetHashedTopologicalTorsionFingerprintAsBitVect(
                    mol, nBits=self.fp_size
                )
            elif self.fp_type == 'RDKit':
                fp = Chem.RDKFingerprint(mol, fpSize=self.fp_size)
            else:
                return np.zeros(self.fp_size)
            
            arr = np.zeros(self.fp_size)
            DataStructs.ConvertToNumpyArray(fp, arr)
            return arr
        except:
            return np.zeros(self.fp_size)
    
    def calculate_descriptors(
        self, 
        smiles: str, 
        descriptor_names: Optional[List[str]] = None
    ) -> dict:
        """분자 descriptor 계산"""
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return {}
        
        try:
            # 기본 descriptor 리스트
            if descriptor_names is None:
                descriptor_names = [
                    'MolWt', 'LogP', 'NumHDonors', 'NumHAcceptors',
                    'TPSA', 'NumRotatableBonds', 'NumAromaticRings',
                    'FractionCSP3', 'MolMR', 'BalabanJ',
                    'BertzCT', 'Chi0', 'Chi1', 'HallKierAlpha',
                    'Kappa1', 'Kappa2'
                ]
            
            # Descriptor 계산
            calc = MoleculeDescriptors.MolecularDescriptorCalculator(descriptor_names)
            desc_values = calc.CalcDescriptors(mol)
            
            return dict(zip(descriptor_names, desc_values))
        except:
            return {}


def prepare_training_data(
    chembl_df: pd.DataFrame,
    protein_name: str,
    fp_type: str = 'ECFP4',
    fp_size: int = 1024,
    include_descriptors: bool = True,
    pos_threshold: float = 10000,
    neg_threshold: float = 20000,
    dataset_type: str = 'binary'
) -> pd.DataFrame:
    """
    학습용 데이터 준비
    
    Args:
        chembl_df: ChEMBL 데이터프레임
        protein_name: 타겟 단백질 이름
        fp_type: Fingerprint 타입
        fp_size: Fingerprint 크기
        include_descriptors: Molecular descriptor 포함 여부
        pos_threshold: 활성 임계값 (nM)
        neg_threshold: 비활성 임계값 (nM)
        dataset_type: 'binary' 또는 'active_only'
    
    Returns:
        처리된 데이터프레임
    """
    extractor = MolecularFeatureExtractor(fp_type, fp_size)
    
    all_compounds = []
    
    # ChEMBL 데이터 처리
    for _, row in chembl_df.iterrows():
        if pd.notna(row['canonical_smiles']):
            ic50_value = row.get('standard_value', None)
            relation = row.get('standard_relation', '')
            
            if pd.notna(ic50_value):
                try:
                    ic50_num = float(ic50_value)
                    
                    if ic50_num <= pos_threshold:
                        if relation in ['=', '<']:
                            potency = 1 
                        else:
                            potency = '-' 
                    elif ic50_num >= neg_threshold:
                        if relation in ['=', '>']:
                            potency = 0 
                        else:
                            potency = '-' 
                    else:
                        potency = '-'
                except (ValueError, TypeError):
                    potency = '-'
            else:
                potency = '-'
            
            all_compounds.append({
                'smiles': row['canonical_smiles'],
                'potency': potency,
                'ic50': ic50_value
            })
    
    df = pd.DataFrame(all_compounds)
    
    # 활성/비활성 화합물 필터링
    if dataset_type == 'active_only':
        filtered_data = df[df['potency'] == 1].copy()
        print(f"Active compounds: {len(filtered_data)}")
    else:
        active_data = df[df['potency'] == 1].copy()
        inactive_data = df[df['potency'] == 0].copy()
        filtered_data = pd.concat([active_data, inactive_data]).reset_index(drop=True)
        print(f"Active: {len(active_data)}, Inactive: {len(inactive_data)}")
    
    if len(filtered_data) == 0:
        raise ValueError("No valid compounds found")
    
    # Fingerprint 계산
    print("Calculating fingerprints...")
    fingerprints = []
    valid_indices = []
    
    for idx, smiles in enumerate(filtered_data['smiles']):
        fp = extractor.smiles_to_fingerprint(smiles)
        if fp.sum() > 0:
            fingerprints.append(fp)
            valid_indices.append(idx)
    
    if len(valid_indices) == 0:
        raise ValueError("No valid fingerprints generated")
    
    valid_data = filtered_data.iloc[valid_indices].reset_index(drop=True)
    fingerprints = np.array(fingerprints)
    
    # DataFrame 생성
    meta_data = pd.DataFrame({
        'SMILES': valid_data['smiles'],
        'Y': valid_data['potency'],
        'IC50': valid_data['ic50']
    })
    
    fp_columns = [f'FP_{i+1}' for i in range(fingerprints.shape[1])]
    fp_data = pd.DataFrame(fingerprints, columns=fp_columns)
    
    result_df = pd.concat([meta_data, fp_data], axis=1)
    
    # Molecular descriptor 추가
    if include_descriptors:
        print("Calculating molecular descriptors...")
        descriptors_list = []
        for smiles in valid_data['smiles']:
            desc = extractor.calculate_descriptors(smiles)
            descriptors_list.append(desc)
        
        desc_df = pd.DataFrame(descriptors_list)
        if not desc_df.empty:
            desc_df.columns = [f'DESC_{col}' for col in desc_df.columns]
            result_df = pd.concat([result_df, desc_df], axis=1)
    
    print(f"Final dataset: {len(result_df)} compounds, {len(result_df.columns)} features")
    
    return result_df
