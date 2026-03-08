from fastapi import APIRouter, HTTPException, UploadFile, File
from app.models.schemas import FooDBPredictRequest, PredictionResponse
from app.services.model_training import ModelTrainingService
from app.services.feature_transform import FeatureTransformService
from app.services.foodb_service import FooDBService
from app.services.foodb_preprocessing import FooDBPreprocessingService
import pandas as pd
import os
from typing import List

router = APIRouter(prefix="/api/foodb", tags=["FooDB Prediction"])

model_service = ModelTrainingService()
feature_service = FeatureTransformService()
foodb_service = FooDBService()
foodb_preprocessing = FooDBPreprocessingService()

@router.post("/upload")
async def upload_foodb_csv(file: UploadFile = File(...)):
    """
    FooDB CSV 파일 업로드 및 저장
    
    - **file**: FooDB에서 다운로드한 CSV 파일
    """
    try:
        # FooDB 디렉토리 생성
        foodb_dir = "saved_data/FooDB"
        os.makedirs(foodb_dir, exist_ok=True)
        
        # CSV 파일 저장
        file_path = os.path.join(foodb_dir, file.filename)
        
        contents = await file.read()
        with open(file_path, 'wb') as f:
            f.write(contents)
        
        # CSV 파일 검증
        df = pd.read_csv(file_path)
        
        # SMILES 컬럼 확인
        smiles_cols = [col for col in df.columns if 'smiles' in col.lower()]
        if not smiles_cols:
            raise HTTPException(
                status_code=400,
                detail="CSV file must contain a SMILES column (e.g., 'canonical_smiles', 'smiles')"
            )
        
        return {
            "message": "FooDB CSV uploaded successfully",
            "file_path": file_path,
            "total_compounds": len(df),
            "columns": df.columns.tolist(),
            "smiles_column": smiles_cols[0]
        }
    
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"Failed to upload FooDB CSV: {str(e)}")


@router.post("/preprocess")
async def preprocess_foodb(
    foodb_file: str,
    protein_name: str,
    fingerprint_type: str = "MACCS",
    dataset_ratio: str = "20x",
    ignore3D: bool = True
):
    """
    FooDB 전처리 (3D Conformer + Mordred Descriptor)
    
    - **foodb_file**: FooDB CSV 파일 경로
    - **protein_name**: 단백질 이름 (descriptor 선택용)
    - **fingerprint_type**: MACCS, ECFP4, MORGAN
    - **dataset_ratio**: 5x, 10x, 20x
    - **ignore3D**: 3D descriptor 무시 여부 (False = 3D 계산)
    
    ⚠️ 주의: 3D conformer 생성은 시간이 오래 걸립니다 (1000개 화합물 ≈ 10-30분)
    """
    try:
        # Descriptor selection 파일 로드
        descriptor_selection_path = "descriptor_selection.csv"
        if not os.path.exists(descriptor_selection_path):
            raise HTTPException(
                status_code=404,
                detail=f"descriptor_selection.csv not found. Please upload this file first."
            )
        
        selected_descriptor = pd.read_csv(descriptor_selection_path)
        
        # 해당 protein/ratio/3D 조합의 descriptor 찾기
        ignore3D_str = "True" if ignore3D else "False"
        file_name = f'descriptors_filtered_{protein_name}_training_{dataset_ratio}_ignore3D_{ignore3D_str}.csv'
        
        if file_name not in selected_descriptor.columns:
            raise HTTPException(
                status_code=404,
                detail=f"Descriptor column '{file_name}' not found in descriptor_selection.csv"
            )
        
        descriptor_names = selected_descriptor[file_name].dropna().tolist()
        
        if not descriptor_names:
            raise HTTPException(
                status_code=400,
                detail=f"No descriptors found for {file_name}"
            )
        
        print(f"Selected {len(descriptor_names)} descriptors for {protein_name} {dataset_ratio} ignore3D={ignore3D}")
        
        # 출력 경로 설정
        output_path = f'saved_data/FooDB/preprocessed/{protein_name}_{fingerprint_type}_{dataset_ratio}_ignore3D_{ignore3D}.csv'
        cache_3d_path = f'saved_data/FooDB/3d_conformers/{protein_name}_{dataset_ratio}.sdf'
        
        # 전처리 실행
        result_df = foodb_preprocessing.preprocess_foodb(
            foodb_csv_path=foodb_file,
            protein_name=protein_name,
            fingerprint_type=fingerprint_type,
            dataset_ratio=dataset_ratio,
            ignore3D=ignore3D,
            descriptor_names=descriptor_names,
            output_path=output_path,
            cache_3d_path=cache_3d_path
        )
        
        return {
            "message": "FooDB preprocessing completed successfully",
            "input_file": foodb_file,
            "output_file": output_path,
            "total_compounds": len(result_df),
            "fingerprint_type": fingerprint_type,
            "fingerprint_size": 167 if fingerprint_type == 'MACCS' else 1024,
            "descriptor_count": len(descriptor_names),
            "total_features": (167 if fingerprint_type == 'MACCS' else 1024) + len(descriptor_names),
            "3d_cache": cache_3d_path if not ignore3D else None
        }
    
    except HTTPException:
        raise
    except Exception as e:
        import traceback
        traceback.print_exc()
        raise HTTPException(status_code=500, detail=f"FooDB preprocessing failed: {str(e)}")


@router.post("/predict")
async def predict_foodb_compounds(
    model_id: str,
    foodb_file: str = "saved_data/FooDB/foodb_compounds.csv",
    batch_size: int = 100,
    top_n: int = 100
):
    """
    FooDB 화합물에 대한 활성 예측 (배치 처리)
    
    - **model_id**: 사용할 모델 ID
    - **foodb_file**: FooDB CSV 파일 경로
    - **batch_size**: 배치 크기 (기본값: 100)
    - **top_n**: 상위 N개 활성 화합물 반환 (기본값: 100)
    """
    try:
        # 모델 로드
        model_info = model_service.load_model(model_id)
        model = model_info['model']
        fingerprint_type = model_info.get('fingerprint_type', 'ECFP4')
        dataset_ratio = model_info.get('dataset_ratio', '20x')
        ignore3D = model_info.get('ignore3D', True)
        feature_count = model_info['feature_count']
        
        print(f"Model loaded: {model_id}")
        print(f"  Fingerprint: {fingerprint_type}, Ratio: {dataset_ratio}, Ignore3D: {ignore3D}")
        print(f"  Expected features: {feature_count}")
        
        # FooDB 파일 로드
        if not os.path.exists(foodb_file):
            raise HTTPException(
                status_code=404,
                detail=f"FooDB file not found: {foodb_file}. Please upload FooDB CSV first."
            )
        
        foodb_df = pd.read_csv(foodb_file)
        
        # SMILES 컬럼 찾기
        smiles_col = None
        for col in ['canonical_smiles', 'smiles', 'SMILES', 'canonical_SMILES']:
            if col in foodb_df.columns:
                smiles_col = col
                break
        
        if not smiles_col:
            raise HTTPException(
                status_code=400,
                detail="FooDB CSV must contain a SMILES column"
            )
        
        print(f"FooDB loaded: {len(foodb_df)} compounds, SMILES column: {smiles_col}")
        
        # 배치 예측
        predictions = foodb_service.predict_batch(
            model=model,
            smiles_list=foodb_df[smiles_col].tolist(),
            fingerprint_type=fingerprint_type,
            dataset_ratio=dataset_ratio,
            ignore3D=ignore3D,
            expected_features=feature_count,
            batch_size=batch_size
        )
        
        # 결과 데이터프레임 생성
        result_df = pd.DataFrame({
            'smiles': foodb_df[smiles_col],
            'prediction': predictions['predictions'],
            'probability_active': predictions['probabilities_active'],
            'probability_inactive': predictions['probabilities_inactive']
        })
        
        # 추가 컬럼이 있으면 포함
        if 'id' in foodb_df.columns:
            result_df['foodb_id'] = foodb_df['id']
        if 'name' in foodb_df.columns:
            result_df['compound_name'] = foodb_df['name']
        
        # 활성 확률 기준으로 정렬
        result_df = result_df.sort_values('probability_active', ascending=False)
        
        # 상위 N개 추출
        top_active = result_df.head(top_n)
        
        # 결과 저장
        output_dir = f"food_predictions"
        os.makedirs(output_dir, exist_ok=True)
        output_path = os.path.join(output_dir, f"{model_id}_foodb_predictions.csv")
        result_df.to_csv(output_path, index=False)
        
        print(f"Predictions saved: {output_path}")
        
        # 통계 계산
        active_count = (result_df['prediction'] == 1).sum()
        inactive_count = (result_df['prediction'] == 0).sum()
        
        return {
            "message": "FooDB prediction completed successfully",
            "model_id": model_id,
            "total_compounds": len(result_df),
            "active_compounds": int(active_count),
            "inactive_compounds": int(inactive_count),
            "top_active_compounds": top_active.to_dict('records'),
            "output_path": output_path,
            "statistics": {
                "mean_probability_active": float(result_df['probability_active'].mean()),
                "max_probability_active": float(result_df['probability_active'].max()),
                "min_probability_active": float(result_df['probability_active'].min())
            }
        }
    
    except HTTPException:
        raise
    except Exception as e:
        import traceback
        traceback.print_exc()
        raise HTTPException(status_code=500, detail=f"FooDB prediction failed: {str(e)}")


@router.post("/predict-smiles")
async def predict_food_smiles(request: FooDBPredictRequest):
    """
    식품 화합물 SMILES 직접 입력 예측
    
    - **smiles_list**: 예측할 SMILES 리스트
    - **model_id**: 사용할 모델 ID
    """
    smiles_list = request.smiles_list
    model_id = request.model_id
    try:
        # 모델 로드
        model_info = model_service.load_model(model_id)
        model = model_info['model']
        fingerprint_type = model_info.get('fingerprint_type', 'ECFP4')
        dataset_ratio = model_info.get('dataset_ratio', '20x')
        ignore3D = model_info.get('ignore3D', True)
        feature_count = model_info['feature_count']
        
        # 특성 변환
        features_df = feature_service.transform_to_fingerprint_with_descriptors(
            smiles_list=smiles_list,
            fp_type=fingerprint_type,
            dataset_ratio=dataset_ratio,
            ignore3D=ignore3D
        )
        
        X = features_df.values
        
        # Feature shape 검증
        if X.shape[1] != feature_count:
            raise ValueError(
                f"Feature shape mismatch, expected: {feature_count}, got {X.shape[1]}. "
                f"Please use the same fingerprint_type and settings as training."
            )
        
        # 예측
        predictions, probabilities = model_service.predict(model_id, X)
        
        # 결과 구성
        results = []
        for i, smiles in enumerate(smiles_list):
            result = {
                'smiles': smiles,
                'prediction': int(predictions[i]),
                'prediction_label': 'Active' if predictions[i] == 1 else 'Inactive'
            }
            
            if probabilities is not None:
                result['probability_inactive'] = float(probabilities[i][0])
                result['probability_active'] = float(probabilities[i][1])
            
            results.append(result)
        
        return {
            "predictions": results,
            "model_info": {
                "model_id": model_id,
                "model_type": model_info['model_type'],
                "protein_name": model_info['protein_name'],
                "fingerprint_type": fingerprint_type,
                "dataset_ratio": dataset_ratio,
                "ignore3D": ignore3D
            }
        }
    
    except FileNotFoundError:
        raise HTTPException(status_code=404, detail=f"Model not found: {model_id}")
    except ValueError as e:
        raise HTTPException(status_code=400, detail=str(e))
    except Exception as e:
        import traceback
        traceback.print_exc()
        raise HTTPException(status_code=500, detail=f"Prediction failed: {str(e)}")
