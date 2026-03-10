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
        # 데이터 경로 찾기 (단일 Dataset 폴더)
        if request.feature_type == "fingerprint":
            data_path = f"raw/Dataset/{request.protein_name}.csv"
        else:
            data_path = f"raw/descriptors/{request.protein_name}_descriptors.csv"
        
        if not os.path.exists(data_path):
            raise HTTPException(
                status_code=404,
                detail=f"Training data not found: {data_path}. Please transform features first."
            )
        
        # CSV에서 학습 메타데이터 추출 (fingerprint_type, dataset_ratio, ignore3D)
        import pandas as pd
        df = pd.read_csv(data_path)
        
        # Feature count로 fingerprint 타입 추정
        feature_count = len(df.columns) - 3  # SMILES, Y, potency 제외
        if feature_count <= 200:  # MACCS (167) + descriptors
            fingerprint_type = "MACCS"
        elif feature_count <= 1100:  # ECFP4/MORGAN (1024) + descriptors
            fingerprint_type = "ECFP4"
        else:
            fingerprint_type = "ECFP4"
        
        # 모델 학습 시 메타데이터 전달
        model_id, metrics = model_service.train_model(
            data_path=data_path,
            model_type=request.model_type,
            protein_name=request.protein_name,
            test_size=request.test_size,
            random_state=request.random_state,
            fingerprint_type=fingerprint_type,
            dataset_ratio="20x",  # 기본값
            ignore3D=True  # 기본값
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
        
        # 모델에 저장된 메타데이터로 동일한 특성 생성
        if request.feature_type == "fingerprint":
            # 모델 메타데이터 사용
            fp_type = model_info.get('fingerprint_type', 'ECFP4')
            dataset_ratio = model_info.get('dataset_ratio', '20x')
            ignore3D = model_info.get('ignore3D', True)
            
            # Fingerprint + Descriptor 생성 (학습 시와 동일)
            features_df = feature_service.transform_to_fingerprint_with_descriptors(
                smiles_list=request.smiles_list,
                fp_type=fp_type,
                dataset_ratio=dataset_ratio,
                ignore3D=ignore3D
            )
            X = features_df.values
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
