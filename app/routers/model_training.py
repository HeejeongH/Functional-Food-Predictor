from fastapi import APIRouter, HTTPException
from app.models.schemas import (
    ModelTrainRequest, 
    ModelTrainResponse,
    PredictionRequest,
    PredictionResponse
)
from app.services.model_training import ModelTrainingService
from app.services.feature_transform import FeatureTransformService
import os
import glob

router = APIRouter(prefix="/api/models", tags=["Model Training & Prediction"])

model_service = ModelTrainingService()
feature_service = FeatureTransformService()

@router.post("/train", response_model=ModelTrainResponse)
async def train_model(request: ModelTrainRequest):
    """
    머신러닝 모델 학습
    
    - **protein_name**: 단백질 이름
    - **model_type**: TabPFN, XGBoost, LightGBM, CatBoost, RandomForest
    - **feature_type**: fingerprint 또는 descriptor
    - **test_size**: 테스트 세트 비율 (0.0 ~ 1.0)
    - **random_state**: 랜덤 시드
    """
    try:
        # 데이터 경로 찾기
        if request.feature_type == "fingerprint":
            data_path = f"raw/TransferSet/{request.protein_name}.csv"
            if not os.path.exists(data_path):
                data_path = f"raw/FewshotSet/{request.protein_name}.csv"
        else:
            data_path = f"raw/descriptors/{request.protein_name}_descriptors.csv"
        
        if not os.path.exists(data_path):
            raise HTTPException(
                status_code=404,
                detail=f"Training data not found: {data_path}. Please transform features first."
            )
        
        # 모델 학습
        model_id, metrics = model_service.train_model(
            data_path=data_path,
            model_type=request.model_type,
            protein_name=request.protein_name,
            test_size=request.test_size,
            random_state=request.random_state
        )
        
        return ModelTrainResponse(
            model_id=model_id,
            accuracy=metrics['accuracy'],
            precision=metrics['precision'],
            recall=metrics['recall'],
            f1_score=metrics['f1_score'],
            roc_auc=metrics['roc_auc'],
            message=f"Model trained successfully: {model_id}"
        )
    
    except HTTPException:
        raise
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))

@router.post("/predict", response_model=PredictionResponse)
async def predict(request: PredictionRequest):
    """
    학습된 모델로 예측 수행
    
    - **smiles_list**: 예측할 SMILES 리스트
    - **model_id**: 학습된 모델 ID
    - **feature_type**: fingerprint 또는 descriptor
    """
    try:
        # 모델 로드
        model_info = model_service.load_model(request.model_id)
        
        # 특성 변환
        if request.feature_type == "fingerprint":
            X = feature_service.transform_to_fingerprint(
                smiles_list=request.smiles_list,
                fp_type="ECFP4"
            )
        else:
            desc_df = feature_service.transform_to_descriptors(
                smiles_list=request.smiles_list,
                descriptor_type="MORDRED_2D"
            )
            X = desc_df.drop(['canonical_SMILES'], axis=1, errors='ignore').values
        
        # 예측
        predictions, probabilities = model_service.predict(request.model_id, X)
        
        # 결과 구성
        results = []
        for i, smiles in enumerate(request.smiles_list):
            result = {
                'smiles': smiles,
                'prediction': int(predictions[i]),
                'prediction_label': 'Active' if predictions[i] == 1 else 'Inactive'
            }
            
            if probabilities is not None:
                result['probability_inactive'] = float(probabilities[i][0])
                result['probability_active'] = float(probabilities[i][1])
            
            results.append(result)
        
        return PredictionResponse(
            predictions=results,
            model_info={
                'model_id': request.model_id,
                'model_type': model_info['model_type'],
                'protein_name': model_info['protein_name'],
                'metrics': model_info['metrics']
            }
        )
    
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))

@router.get("/list")
async def list_models():
    """
    학습된 모델 리스트 조회
    """
    try:
        model_files = glob.glob("models_trained/*.pkl")
        models = []
        
        for model_path in model_files:
            model_id = os.path.basename(model_path).replace('.pkl', '')
            try:
                model_info = model_service.load_model(model_id)
                models.append({
                    'model_id': model_id,
                    'model_type': model_info['model_type'],
                    'protein_name': model_info['protein_name'],
                    'metrics': model_info['metrics'],
                    'feature_count': model_info['feature_count'],
                    'train_size': model_info['train_size'],
                    'test_size': model_info['test_size']
                })
            except:
                continue
        
        return {
            'total_models': len(models),
            'models': models
        }
    
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))

@router.get("/info/{model_id}")
async def get_model_info(model_id: str):
    """
    특정 모델 정보 조회
    """
    try:
        model_info = model_service.load_model(model_id)
        
        return {
            'model_id': model_id,
            'model_type': model_info['model_type'],
            'protein_name': model_info['protein_name'],
            'metrics': model_info['metrics'],
            'feature_count': model_info['feature_count'],
            'train_size': model_info['train_size'],
            'test_size': model_info['test_size']
        }
    
    except FileNotFoundError:
        raise HTTPException(status_code=404, detail=f"Model not found: {model_id}")
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))
