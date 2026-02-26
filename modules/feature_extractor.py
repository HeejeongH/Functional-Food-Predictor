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
    dataset_type: str = 'binary',
    use_decoys: bool = False,
    decoy_ratio: float = 50.0,
    decoy_method: str = 'dude',
    decoy_source: Optional[pd.DataFrame] = None
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
        use_decoys: DUDE-style decoy 사용 여부
        decoy_ratio: 활성:비활성 비율 (예: 50 = 1:50)
        decoy_method: 'dude' (DUDE-style) 또는 'random'
        decoy_source: Decoy 소스 데이터 (None이면 chembl_df 사용)
    
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
    
    # 활성 화합물 추출
    active_data = df[df['potency'] == 1].copy()
    print(f"Active compounds: {len(active_data)}")
    
    if len(active_data) == 0:
        raise ValueError("No active compounds found (활성 화합물이 없습니다. IC50 임계값을 조정해보세요.)")
    
    # Decoy 생성 또는 기존 비활성 화합물 사용
    if use_decoys and dataset_type == 'binary':
        from modules.decoy_generator import add_decoys_to_dataset
        
        print(f"\n{'='*60}")
        print(f"Generating DECOY compounds (method: {decoy_method})")
        print(f"{'='*60}")
        
        # Decoy 소스 설정
        if decoy_source is None:
            decoy_source = chembl_df  # 전체 ChEMBL 데이터 사용
        
        # Decoy 추가
        filtered_data = add_decoys_to_dataset(
            active_df=active_data.rename(columns={'smiles': 'SMILES', 'potency': 'Y', 'ic50': 'IC50'}),
            decoy_source=decoy_source,
            decoy_ratio=decoy_ratio,
            decoy_method=decoy_method,
            similarity_threshold=0.75,
            property_tolerance=0.2
        )
        
        # 컬럼명 복원
        filtered_data = filtered_data.rename(columns={'SMILES': 'smiles', 'IC50': 'ic50'})
        print(f"\nFinal dataset: Active {(filtered_data['Y']==1).sum()}, Decoy {(filtered_data['Y']==0).sum()}")
        
    elif dataset_type == 'active_only':
        filtered_data = active_data
        print(f"Using only active compounds (no negatives)")
    else:
        # 기존 방식: IC50 기반 비활성 화합물
        inactive_data = df[df['potency'] == 0].copy()
        filtered_data = pd.concat([active_data, inactive_data]).reset_index(drop=True)
        print(f"Using IC50-based inactive compounds: Active {len(active_data)}, Inactive {len(inactive_data)}")
    
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
    
    # DataFrame 생성 (Y를 명시적으로 정수형으로 변환)
    meta_data = pd.DataFrame({
        'SMILES': valid_data['smiles'],
        'Y': valid_data['potency'].astype(int),  # 명시적 정수 변환
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
    
    # 데이터 검증: Y 컬럼이 0 또는 1만 포함하는지 확인
    unique_y = result_df['Y'].unique()
    if not all(y in [0, 1] for y in unique_y):
        raise ValueError(f"Y 컬럼에 잘못된 값이 있습니다: {unique_y}. 0과 1만 허용됩니다.")
    
    print(f"Data validation passed - Y values: {sorted(unique_y)}")
    
    return result_df
