"""
Decoy 생성 모듈 - DUDE 스타일 네거티브 샘플 생성
물리화학적 특성은 유사하지만 구조적으로 다른 decoy 화합물 생성
"""

import numpy as np
import pandas as pd
from typing import List, Optional, Tuple
from rdkit import Chem
from rdkit.Chem import Descriptors, AllChem, DataStructs
from rdkit import RDLogger
RDLogger.DisableLog('rdApp.*')


class DecoyGenerator:
    """
    DUDE-style Decoy Generator
    활성 화합물과 물리화학적 특성은 유사하지만 구조적으로 다른 decoy 생성
    """
    
    def __init__(
        self,
        similarity_threshold: float = 0.75,
        property_tolerance: float = 0.2,
        random_state: int = 42
    ):
        """
        Args:
            similarity_threshold: Tanimoto 유사도 상한선 (이보다 낮아야 함)
            property_tolerance: 물리화학적 특성 허용 오차 (±20%)
            random_state: 랜덤 시드
        """
        self.similarity_threshold = similarity_threshold
        self.property_tolerance = property_tolerance
        self.random_state = random_state
        np.random.seed(random_state)
    
    def calculate_properties(self, smiles: str) -> dict:
        """분자의 물리화학적 특성 계산"""
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return {}
        
        try:
            return {
                'MW': Descriptors.MolWt(mol),
                'LogP': Descriptors.MolLogP(mol),
                'HBD': Descriptors.NumHDonors(mol),
                'HBA': Descriptors.NumHAcceptors(mol),
                'TPSA': Descriptors.TPSA(mol),
                'RotBonds': Descriptors.NumRotatableBonds(mol),
                'Rings': Descriptors.RingCount(mol)
            }
        except:
            return {}
    
    def calculate_similarity(self, smiles1: str, smiles2: str) -> float:
        """두 화합물 간 Tanimoto 유사도 계산"""
        mol1 = Chem.MolFromSmiles(smiles1)
        mol2 = Chem.MolFromSmiles(smiles2)
        
        if mol1 is None or mol2 is None:
            return 0.0
        
        try:
            fp1 = AllChem.GetMorganFingerprintAsBitVect(mol1, 2, nBits=2048)
            fp2 = AllChem.GetMorganFingerprintAsBitVect(mol2, 2, nBits=2048)
            return DataStructs.TanimotoSimilarity(fp1, fp2)
        except:
            return 0.0
    
    def is_property_similar(self, prop1: dict, prop2: dict) -> bool:
        """두 화합물의 물리화학적 특성이 유사한지 확인"""
        if not prop1 or not prop2:
            return False
        
        tolerance = self.property_tolerance
        
        for key in ['MW', 'LogP', 'TPSA']:
            if key not in prop1 or key not in prop2:
                continue
            
            val1, val2 = prop1[key], prop2[key]
            if val1 == 0:
                continue
            
            diff = abs(val1 - val2) / val1
            if diff > tolerance:
                return False
        
        # 정수형 특성은 정확히 일치하거나 ±1
        for key in ['HBD', 'HBA', 'RotBonds', 'Rings']:
            if key not in prop1 or key not in prop2:
                continue
            
            if abs(prop1[key] - prop2[key]) > 1:
                return False
        
        return True
    
    def generate_decoys_from_database(
        self,
        active_smiles_list: List[str],
        compound_database: pd.DataFrame,
        n_decoys_per_active: int = 50,
        smiles_column: str = 'canonical_smiles'
    ) -> pd.DataFrame:
        """
        데이터베이스에서 DUDE-style decoy 선택
        
        Args:
            active_smiles_list: 활성 화합물 SMILES 리스트
            compound_database: 후보 화합물 데이터베이스 (ChEMBL 등)
            n_decoys_per_active: 활성 화합물당 decoy 개수
            smiles_column: SMILES 컬럼 이름
        
        Returns:
            Decoy 화합물 DataFrame
        """
        print(f"Generating decoys from database with {len(compound_database)} compounds...")
        
        # 활성 화합물 특성 계산
        active_properties = []
        for smiles in active_smiles_list:
            props = self.calculate_properties(smiles)
            if props:
                active_properties.append({'smiles': smiles, 'props': props})
        
        if not active_properties:
            raise ValueError("No valid active compounds")
        
        # Decoy 후보 준비
        decoy_candidates = []
        for _, row in compound_database.iterrows():
            candidate_smiles = row[smiles_column]
            if pd.isna(candidate_smiles):
                continue
            
            # 활성 화합물과 동일하면 제외
            if candidate_smiles in active_smiles_list:
                continue
            
            props = self.calculate_properties(candidate_smiles)
            if not props:
                continue
            
            decoy_candidates.append({
                'smiles': candidate_smiles,
                'props': props
            })
        
        print(f"Found {len(decoy_candidates)} decoy candidates")
        
        # 각 활성 화합물에 대해 decoy 선택
        selected_decoys = []
        
        for active_data in active_properties:
            active_smiles = active_data['smiles']
            active_props = active_data['props']
            
            # 후보 중에서 조건에 맞는 decoy 찾기
            valid_decoys = []
            
            for decoy_data in decoy_candidates:
                # 1. 물리화학적 특성 유사성 확인
                if not self.is_property_similar(active_props, decoy_data['props']):
                    continue
                
                # 2. 구조적 유사도 확인 (너무 유사하면 안됨)
                similarity = self.calculate_similarity(active_smiles, decoy_data['smiles'])
                if similarity >= self.similarity_threshold:
                    continue
                
                valid_decoys.append({
                    'smiles': decoy_data['smiles'],
                    'similarity': similarity,
                    'active_reference': active_smiles
                })
            
            # 유사도가 낮은 순서로 정렬하여 선택
            valid_decoys.sort(key=lambda x: x['similarity'])
            
            # 필요한 만큼 선택
            n_select = min(n_decoys_per_active, len(valid_decoys))
            selected = valid_decoys[:n_select]
            selected_decoys.extend(selected)
            
            print(f"Active {active_smiles[:50]}... -> {n_select} decoys selected")
        
        if not selected_decoys:
            raise ValueError("No valid decoys found. Try adjusting similarity_threshold or property_tolerance")
        
        decoy_df = pd.DataFrame(selected_decoys)
        print(f"Total {len(decoy_df)} decoys generated")
        
        return decoy_df
    
    def generate_random_decoys(
        self,
        active_smiles_list: List[str],
        compound_pool: List[str],
        n_decoys_per_active: int = 50
    ) -> pd.DataFrame:
        """
        랜덤 풀에서 decoy 선택 (빠른 방법)
        
        Args:
            active_smiles_list: 활성 화합물 SMILES 리스트
            compound_pool: 후보 화합물 SMILES 풀
            n_decoys_per_active: 활성 화합물당 decoy 개수
        
        Returns:
            Decoy 화합물 DataFrame
        """
        print(f"Generating random decoys from pool of {len(compound_pool)} compounds...")
        
        # 활성 화합물 제외
        decoy_pool = [s for s in compound_pool if s not in active_smiles_list]
        
        # 랜덤 샘플링
        n_total_decoys = len(active_smiles_list) * n_decoys_per_active
        n_total_decoys = min(n_total_decoys, len(decoy_pool))
        
        selected_smiles = np.random.choice(decoy_pool, n_total_decoys, replace=False)
        
        decoy_df = pd.DataFrame({
            'smiles': selected_smiles,
            'similarity': 0.0,  # 계산 생략
            'active_reference': 'random'
        })
        
        print(f"Total {len(decoy_df)} random decoys generated")
        
        return decoy_df


def add_decoys_to_dataset(
    active_df: pd.DataFrame,
    decoy_source: pd.DataFrame,
    decoy_ratio: float = 50.0,
    decoy_method: str = 'dude',
    similarity_threshold: float = 0.75,
    property_tolerance: float = 0.2
) -> pd.DataFrame:
    """
    활성 화합물 데이터셋에 decoy 추가
    
    Args:
        active_df: 활성 화합물 DataFrame (SMILES, Y=1)
        decoy_source: Decoy 소스 데이터 (ChEMBL 등)
        decoy_ratio: 활성:비활성 비율 (예: 50 = 1:50)
        decoy_method: 'dude' (DUDE-style) 또는 'random'
        similarity_threshold: DUDE 유사도 임계값
        property_tolerance: DUDE 특성 허용오차
    
    Returns:
        활성 + Decoy 통합 DataFrame
    """
    print("="*60)
    print(f"Adding decoys with method: {decoy_method.upper()}")
    print(f"Active compounds: {len(active_df)}")
    print(f"Decoy ratio: 1:{int(decoy_ratio)}")
    print("="*60)
    
    # 활성 화합물 SMILES 추출
    if 'SMILES' in active_df.columns:
        active_smiles = active_df['SMILES'].tolist()
    elif 'canonical_smiles' in active_df.columns:
        active_smiles = active_df['canonical_smiles'].tolist()
    else:
        raise ValueError("No SMILES column found in active_df")
    
    # Decoy 생성기 초기화
    generator = DecoyGenerator(
        similarity_threshold=similarity_threshold,
        property_tolerance=property_tolerance
    )
    
    # Decoy 생성
    n_decoys = int(len(active_smiles) * decoy_ratio)
    n_decoys_per_active = int(decoy_ratio)
    
    if decoy_method == 'dude':
        decoy_df = generator.generate_decoys_from_database(
            active_smiles_list=active_smiles,
            compound_database=decoy_source,
            n_decoys_per_active=n_decoys_per_active,
            smiles_column='canonical_smiles'
        )
    elif decoy_method == 'random':
        if 'canonical_smiles' in decoy_source.columns:
            compound_pool = decoy_source['canonical_smiles'].dropna().tolist()
        else:
            compound_pool = decoy_source['SMILES'].dropna().tolist()
        
        decoy_df = generator.generate_random_decoys(
            active_smiles_list=active_smiles,
            compound_pool=compound_pool,
            n_decoys_per_active=n_decoys_per_active
        )
    else:
        raise ValueError(f"Unknown decoy_method: {decoy_method}")
    
    # Decoy 데이터프레임 생성 (활성 데이터와 동일한 구조)
    decoy_data = pd.DataFrame({
        'SMILES': decoy_df['smiles'],
        'Y': 0,  # Decoy는 비활성
        'IC50': np.nan,
        'source': 'decoy',
        'decoy_method': decoy_method
    })
    
    # 활성 데이터에 source 추가
    active_data = active_df.copy()
    active_data['source'] = 'active'
    active_data['decoy_method'] = 'N/A'
    
    # 통합
    combined_df = pd.concat([active_data, decoy_data], ignore_index=True)
    
    print("\n" + "="*60)
    print("Dataset Summary:")
    print(f"  Active compounds: {(combined_df['Y'] == 1).sum()}")
    print(f"  Decoy compounds:  {(combined_df['Y'] == 0).sum()}")
    print(f"  Total:            {len(combined_df)}")
    print(f"  Ratio:            1:{(combined_df['Y'] == 0).sum() / (combined_df['Y'] == 1).sum():.1f}")
    print("="*60)
    
    return combined_df
