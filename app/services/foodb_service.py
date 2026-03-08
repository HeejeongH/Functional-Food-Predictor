import pandas as pd
import numpy as np
from typing import List, Dict
from app.services.feature_transform import FeatureTransformService

class FooDBService:
    """
    FooDB 식품 화합물 예측 서비스
    """
    
    def __init__(self):
        self.feature_service = FeatureTransformService()
    
    def predict_batch(
        self,
        model,
        smiles_list: List[str],
        fingerprint_type: str,
        dataset_ratio: str,
        ignore3D: bool,
        expected_features: int,
        batch_size: int = 100
    ) -> Dict:
        """
        FooDB 화합물 배치 예측
        
        Args:
            model: 학습된 모델
            smiles_list: SMILES 리스트
            fingerprint_type: Fingerprint 타입 (MACCS, ECFP4, MORGAN)
            dataset_ratio: 데이터셋 비율 (5x, 10x, 20x)
            ignore3D: 3D descriptor 무시 여부
            expected_features: 예상 특성 개수
            batch_size: 배치 크기
        
        Returns:
            예측 결과 딕셔너리
        """
        
        total_compounds = len(smiles_list)
        n_batches = (total_compounds + batch_size - 1) // batch_size
        
        all_predictions = []
        all_probabilities = []
        
        print(f"Starting batch prediction: {total_compounds} compounds, {n_batches} batches")
        
        for i in range(n_batches):
            start_idx = i * batch_size
            end_idx = min((i + 1) * batch_size, total_compounds)
            batch_smiles = smiles_list[start_idx:end_idx]
            
            print(f"  Batch {i+1}/{n_batches}: {start_idx}-{end_idx}")
            
            try:
                # 특성 변환
                features_df = self.feature_service.transform_to_fingerprint_with_descriptors(
                    smiles_list=batch_smiles,
                    fp_type=fingerprint_type,
                    dataset_ratio=dataset_ratio,
                    ignore3D=ignore3D
                )
                
                X = features_df.values
                
                # Feature shape 검증
                if X.shape[1] != expected_features:
                    print(f"  Warning: Feature mismatch in batch {i+1}, expected {expected_features}, got {X.shape[1]}")
                    # 특성 개수 맞추기 (패딩 또는 자르기)
                    if X.shape[1] < expected_features:
                        # 패딩
                        padding = np.zeros((X.shape[0], expected_features - X.shape[1]))
                        X = np.hstack([X, padding])
                    else:
                        # 자르기
                        X = X[:, :expected_features]
                
                # 예측
                predictions = model.predict(X)
                
                if hasattr(model, 'predict_proba'):
                    probabilities = model.predict_proba(X)
                else:
                    probabilities = np.column_stack([1 - predictions, predictions])
                
                all_predictions.extend(predictions.tolist())
                all_probabilities.extend(probabilities.tolist())
                
            except Exception as e:
                print(f"  Error in batch {i+1}: {str(e)}")
                # 오류 발생 시 기본값으로 채우기
                all_predictions.extend([0] * len(batch_smiles))
                all_probabilities.extend([[1.0, 0.0]] * len(batch_smiles))
        
        # 결과 정리
        probabilities_array = np.array(all_probabilities)
        
        return {
            'predictions': all_predictions,
            'probabilities_inactive': probabilities_array[:, 0].tolist(),
            'probabilities_active': probabilities_array[:, 1].tolist()
        }
