"""
FooDB 예측 모듈 - 실제 식품 화합물 데이터에 대한 예측
"""

import json
import numpy as np
import pandas as pd
from typing import Optional
from rdkit import Chem
from .feature_extractor import MolecularFeatureExtractor


def load_foodb_from_json(json_file: str) -> pd.DataFrame:
    """
    FooDB JSON 파일에서 화합물 데이터 로드
    
    Args:
        json_file: FooDB JSON 파일 경로
    
    Returns:
        화합물 데이터프레임
    """
    compounds = []
    
    with open(json_file, 'r') as f:
        for line in f:
            line = line.strip()
            if line:
                try:
                    compound = json.loads(line)
                    smiles = compound.get('moldb_smiles')
                    if smiles:
                        compounds.append({
                            'id': compound.get('id'),
                            'public_id': compound.get('public_id'),
                            'name': compound.get('name'),
                            'raw_SMILES': smiles
                        })
                except json.JSONDecodeError:
                    continue
    
    df = pd.DataFrame(compounds)
    df['canonical_SMILES'] = df['raw_SMILES'].apply(_to_canonical_smiles)
    
    print(f"Loaded {len(df)} compounds from FooDB")
    
    return df


def load_foodb_from_csv(csv_file: str) -> pd.DataFrame:
    """
    FooDB CSV 파일에서 화합물 데이터 로드
    
    Args:
        csv_file: FooDB CSV 파일 경로
    
    Returns:
        화합물 데이터프레임
    """
    df = pd.read_csv(csv_file)
    print(f"Loaded {len(df)} compounds from FooDB")
    return df


def _to_canonical_smiles(smiles: str) -> Optional[str]:
    """SMILES를 canonical 형태로 변환"""
    if pd.isna(smiles):
        return None
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return None
        return Chem.MolToSmiles(mol, canonical=True)
    except:
        return None


def prepare_foodb_features(
    foodb_df: pd.DataFrame,
    fp_type: str = 'ECFP4',
    fp_size: int = 1024,
    include_descriptors: bool = True
) -> pd.DataFrame:
    """
    FooDB 데이터에 대한 피처 추출
    
    Args:
        foodb_df: FooDB 데이터프레임
        fp_type: Fingerprint 타입
        fp_size: Fingerprint 크기
        include_descriptors: Molecular descriptor 포함 여부
    
    Returns:
        피처가 추가된 데이터프레임
    """
    extractor = MolecularFeatureExtractor(fp_type, fp_size)
    
    print("Extracting features from FooDB compounds...")
    
    # Fingerprint 계산
    fingerprints = []
    valid_indices = []
    
    for idx, row in foodb_df.iterrows():
        smiles = row.get('canonical_SMILES', row.get('raw_SMILES'))
        if pd.notna(smiles):
            fp = extractor.smiles_to_fingerprint(smiles)
            if fp.sum() > 0:
                fingerprints.append(fp)
                valid_indices.append(idx)
    
    if len(valid_indices) == 0:
        raise ValueError("No valid fingerprints generated from FooDB")
    
    valid_data = foodb_df.iloc[valid_indices].reset_index(drop=True)
    fingerprints = np.array(fingerprints)
    
    # DataFrame 생성
    fp_columns = [f'FP_{i+1}' for i in range(fingerprints.shape[1])]
    fp_df = pd.DataFrame(fingerprints, columns=fp_columns)
    
    result_df = pd.concat([valid_data, fp_df], axis=1)
    
    # Molecular descriptor 추가
    if include_descriptors:
        print("Calculating molecular descriptors...")
        descriptors_list = []
        for _, row in valid_data.iterrows():
            smiles = row.get('canonical_SMILES', row.get('raw_SMILES'))
            desc = extractor.calculate_descriptors(smiles)
            descriptors_list.append(desc)
        
        desc_df = pd.DataFrame(descriptors_list)
        if not desc_df.empty:
            desc_df.columns = [f'DESC_{col}' for col in desc_df.columns]
            result_df = pd.concat([result_df, desc_df], axis=1)
    
    print(f"Features extracted for {len(result_df)} compounds")
    
    return result_df


def predict_foodb(
    foodb_df: pd.DataFrame,
    trainer,
    feature_columns: list,
    threshold: float = 0.5
) -> pd.DataFrame:
    """
    FooDB 화합물에 대한 예측 수행
    
    Args:
        foodb_df: 피처가 포함된 FooDB 데이터프레임
        trainer: 학습된 ModelTrainer 객체
        feature_columns: 사용할 피처 컬럼 리스트
        threshold: 예측 임계값
    
    Returns:
        예측 결과가 포함된 데이터프레임
    """
    print("Predicting FooDB compounds...")
    
    # 피처 추출
    X = foodb_df[feature_columns].values
    
    # 예측
    predictions, probabilities = trainer.predict(X)
    
    # 결과 DataFrame 생성
    result_df = foodb_df.copy()
    result_df['prediction'] = predictions
    result_df['probability_class_0'] = probabilities[:, 0]
    result_df['probability_class_1'] = probabilities[:, 1]
    result_df['is_active'] = (probabilities[:, 1] >= threshold).astype(int)
    
    # 활성 화합물만 필터링
    active_df = result_df[result_df['is_active'] == 1].copy()
    active_df = active_df.sort_values('probability_class_1', ascending=False)
    
    print(f"Found {len(active_df)} active compounds (threshold: {threshold})")
    
    return active_df


def batch_predict_foodb(
    foodb_df: pd.DataFrame,
    trainer,
    feature_columns: list,
    batch_size: int = 500,
    threshold: float = 0.5
) -> pd.DataFrame:
    """
    대용량 FooDB 화합물에 대한 배치 예측
    
    Args:
        foodb_df: 피처가 포함된 FooDB 데이터프레임
        trainer: 학습된 ModelTrainer 객체
        feature_columns: 사용할 피처 컬럼 리스트
        batch_size: 배치 크기
        threshold: 예측 임계값
    
    Returns:
        예측 결과가 포함된 데이터프레임
    """
    print(f"Batch predicting {len(foodb_df)} FooDB compounds...")
    
    n_samples = len(foodb_df)
    n_batches = (n_samples + batch_size - 1) // batch_size
    
    all_predictions = []
    all_probabilities = []
    
    for i in range(n_batches):
        start_idx = i * batch_size
        end_idx = min((i + 1) * batch_size, n_samples)
        
        batch_df = foodb_df.iloc[start_idx:end_idx]
        X_batch = batch_df[feature_columns].values
        
        pred, prob = trainer.predict(X_batch)
        
        all_predictions.append(pred)
        all_probabilities.append(prob)
        
        if (i + 1) % 10 == 0:
            print(f"Processed {end_idx}/{n_samples} compounds")
    
    # 결과 합치기
    predictions = np.concatenate(all_predictions)
    probabilities = np.vstack(all_probabilities)
    
    result_df = foodb_df.copy()
    result_df['prediction'] = predictions
    result_df['probability_class_0'] = probabilities[:, 0]
    result_df['probability_class_1'] = probabilities[:, 1]
    result_df['is_active'] = (probabilities[:, 1] >= threshold).astype(int)
    
    # 활성 화합물만 필터링
    active_df = result_df[result_df['is_active'] == 1].copy()
    active_df = active_df.sort_values('probability_class_1', ascending=False)
    
    print(f"Found {len(active_df)} active compounds (threshold: {threshold})")
    
    return active_df
