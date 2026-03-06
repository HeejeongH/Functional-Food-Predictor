from fastapi import APIRouter, HTTPException
from app.models.schemas import SHAPAnalysisRequest, SHAPAnalysisResponse
from app.services.model_training import ModelTrainingService
from app.services.shap_analysis import SHAPAnalysisService
from app.services.feature_transform import FeatureTransformService
import pandas as pd
import os

router = APIRouter(prefix="/api/shap", tags=["SHAP Analysis"])

model_service = ModelTrainingService()
shap_service = SHAPAnalysisService()
feature_service = FeatureTransformService()

@router.post("/analyze", response_model=SHAPAnalysisResponse)
async def analyze_shap(request: SHAPAnalysisRequest):
    """
    SHAP 분석을 통한 특성 중요도 분석
    
    - **model_id**: 학습된 모델 ID
    - **feature_type**: fingerprint 또는 descriptor
    - **top_n**: 상위 N개 중요 특성 (기본값: 20)
    """
    try:
        # 모델 로드
        model_info = model_service.load_model(request.model_id)
        model = model_info['model']
        protein_name = model_info['protein_name']
        
        # 학습 데이터 로드 (단일 Dataset 폴더 사용)
        if request.feature_type == "fingerprint":
            data_path = f"raw/Dataset/{protein_name}.csv"
        else:
            data_path = f"raw/descriptors/{protein_name}_descriptors.csv"
        
        if not os.path.exists(data_path):
            raise HTTPException(
                status_code=404,
                detail=f"Training data not found: {data_path}"
            )
        
        df = pd.read_csv(data_path)
        X = df.drop(['SMILES', 'Y', 'potency'], axis=1, errors='ignore').values
        feature_names = df.drop(['SMILES', 'Y', 'potency'], axis=1, errors='ignore').columns.tolist()
        
        # SHAP 분석
        shap_result = shap_service.analyze_features(
            model=model,
            X=X[:500],  # 샘플링 (계산 속도를 위해)
            feature_names=feature_names,
            top_n=request.top_n
        )
        
        return SHAPAnalysisResponse(
            top_features=shap_result['top_features'],
            shap_values_summary={
                'mean_abs_shap': float(abs(shap_result['shap_values']).mean()),
                'max_abs_shap': float(abs(shap_result['shap_values']).max()),
                'samples_analyzed': len(X[:500])
            },
            plot_path=shap_result['plot_path']
        )
    
    except HTTPException:
        raise
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))

@router.get("/feature-importance/{model_id}")
async def get_feature_importance(model_id: str, top_n: int = 20):
    """
    모델의 특성 중요도 조회 (간단 버전)
    """
    try:
        model_info = model_service.load_model(model_id)
        model = model_info['model']
        
        # Tree-based 모델의 경우 feature_importances_ 사용
        if hasattr(model, 'feature_importances_'):
            importances = model.feature_importances_
            
            # Top N 특성
            top_indices = importances.argsort()[-top_n:][::-1]
            
            results = []
            for idx in top_indices:
                results.append({
                    'feature_index': int(idx),
                    'feature_name': f"Feature_{idx}",
                    'importance': float(importances[idx])
                })
            
            return {
                'model_id': model_id,
                'top_features': results,
                'method': 'built-in feature_importances_'
            }
        else:
            return {
                'model_id': model_id,
                'message': 'Model does not support feature importance. Use SHAP analysis instead.',
                'method': 'none'
            }
    
    except FileNotFoundError:
        raise HTTPException(status_code=404, detail=f"Model not found: {model_id}")
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))
