from fastapi import APIRouter, HTTPException
from app.models.schemas import (
    ModelTrainRequest,
    ModelTrainResponse,
    PredictionRequest,
    PredictionResponse,
)
from app.services.model_training import ModelTrainingService
from app.services.feature_transform import FeatureTransformService
from app.services.foodb_service import FooDBService
import pandas as pd
import numpy as np
import os
import glob

router = APIRouter(prefix="/api/models", tags=["Model Training & Prediction"])

model_service   = ModelTrainingService()
feature_service = FeatureTransformService()
foodb_service   = FooDBService()


# ── 학습 ─────────────────────────────────────────────────────────────────

@router.post("/train", response_model=ModelTrainResponse)
async def train_model(request: ModelTrainRequest):
    """
    머신러닝 모델 학습 (노트북 final-model 셀과 동일한 파이프라인)

    - **protein_name**    : 단백질 이름
    - **model_type**      : TabPFN, XGBoost, LightGBM, CatBoost, RandomForest
    - **feature_type**    : fingerprint 또는 descriptor
    - **fingerprint_type**: ecfp4 | ecfp6 | fcfp4 | fcfp6 | maccs | rdkit
    - **dataset_ratio**   : 1x | 3x | 5x | 10x | 20x | 50x
    - **fs_method**       : mi | rfe | shap | pca | random | none
    - **n_features**      : 선택 feature 수 (기본 100)
    - **test_size**       : 테스트 비율 (노트북 기본 0.1)
    """
    try:
        # 데이터 경로 결정
        if request.feature_type == "fingerprint":
            data_path = f"raw/Dataset/{request.protein_name}.csv"
        else:
            data_path = f"raw/descriptors/{request.protein_name}_descriptors.csv"

        if not os.path.exists(data_path):
            raise HTTPException(
                status_code=404,
                detail=(
                    f"Training data not found: {data_path}. "
                    "Please run feature transform first."
                ),
            )

        model_id, metrics = model_service.train_model(
            data_path        = data_path,
            model_type       = request.model_type,
            protein_name     = request.protein_name,
            test_size        = request.test_size,
            random_state     = request.random_state,
            fingerprint_type = request.fingerprint_type,
            dataset_ratio    = request.dataset_ratio,
            ignore3D         = request.ignore3D,
            fs_method        = request.fs_method,
            n_features       = request.n_features,
        )

        return ModelTrainResponse(
            model_id   = model_id,
            f1         = metrics["f1"],
            auc        = metrics["auc"],
            mcc        = metrics["mcc"],
            accuracy   = metrics["accuracy"],
            precision  = metrics["precision"],
            recall     = metrics["recall"],
            fs_method  = request.fs_method,
            n_features = request.n_features,
            message    = f"Model trained successfully: {model_id}",
        )

    except HTTPException:
        raise
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))


# ── 예측 ─────────────────────────────────────────────────────────────────

@router.post("/predict", response_model=PredictionResponse)
async def predict(request: PredictionRequest):
    """
    학습된 모델로 예측 수행.
    모델에 저장된 fs_method / transformer / selected_features 를 자동으로 적용.

    - **smiles_list**: 예측할 SMILES 리스트
    - **model_id**   : 학습된 모델 ID
    - **feature_type**: fingerprint 또는 descriptor
    """
    try:
        model_info = model_service.load_model(request.model_id)

        if request.feature_type == "fingerprint":
            fp_type      = model_info.get("fingerprint_type", "ecfp4")
            dataset_ratio = model_info.get("dataset_ratio", "20x")
            ignore3D     = model_info.get("ignore3D", True)

            # 학습 시와 동일한 raw feature 생성
            features_df = feature_service.transform_to_fingerprint_with_descriptors(
                smiles_list   = request.smiles_list,
                fp_type       = fp_type,
                dataset_ratio = dataset_ratio,
                ignore3D      = ignore3D,
            )

            # 학습 컬럼 순서 정렬
            train_cols = model_info.get("X_train_columns", features_df.columns.tolist())
            X_raw = features_df.reindex(columns=train_cols, fill_value=0).fillna(0)

            # FS 변환 적용
            X_pred = foodb_service._apply_fs_transform(
                X_raw,
                fs_method   = model_info.get("fs_method", "none"),
                transformer = model_info.get("transformer"),
                features    = model_info.get("selected_features"),
            )
        else:
            desc_df = feature_service.transform_to_descriptors(
                smiles_list     = request.smiles_list,
                descriptor_type = "MORDRED_2D",
            )
            X_pred = desc_df.drop(["canonical_SMILES"], axis=1, errors="ignore").values

        predictions, probabilities = model_service.predict(request.model_id, X_pred)

        results = []
        for i, smiles in enumerate(request.smiles_list):
            item = {
                "smiles"           : smiles,
                "prediction"       : int(predictions[i]),
                "prediction_label" : "Active" if predictions[i] == 1 else "Inactive",
            }
            if probabilities is not None:
                item["probability_inactive"] = float(probabilities[i][0])
                item["probability_active"]   = float(probabilities[i][1])
            results.append(item)

        return PredictionResponse(
            predictions = results,
            model_info  = {
                "model_id"    : request.model_id,
                "model_type"  : model_info["model_type"],
                "protein_name": model_info["protein_name"],
                "metrics"     : model_info["metrics"],
                "fs_method"   : model_info.get("fs_method", "none"),
                "n_features"  : model_info.get("n_features"),
            },
        )

    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))


# ── 모델 목록 ─────────────────────────────────────────────────────────────

@router.get("/list")
async def list_models():
    """학습된 모델 목록 조회 (.joblib 우선, .pkl 폴백)"""
    try:
        # joblib 파일 우선, 구버전 pkl 도 포함
        joblib_files = glob.glob("models_trained/*.joblib")
        pkl_files    = glob.glob("models_trained/*.pkl")

        # 이미 joblib 이 있는 모델은 pkl 중복 제거
        joblib_ids = {os.path.basename(f).replace(".joblib", "") for f in joblib_files}
        all_files  = joblib_files + [
            f for f in pkl_files
            if os.path.basename(f).replace(".pkl", "") not in joblib_ids
        ]

        models = []
        for model_path in all_files:
            ext      = ".joblib" if model_path.endswith(".joblib") else ".pkl"
            model_id = os.path.basename(model_path).replace(ext, "")
            try:
                info = model_service.load_model(model_id)
                models.append({
                    "model_id"    : model_id,
                    "model_type"  : info["model_type"],
                    "protein_name": info["protein_name"],
                    "metrics"     : info["metrics"],
                    "feature_count": info["feature_count"],
                    "train_size"  : info["train_size"],
                    "test_size_n" : info.get("test_size_n"),
                    "fs_method"   : info.get("fs_method", "unknown"),
                    "n_features"  : info.get("n_features"),
                    "fingerprint_type": info.get("fingerprint_type"),
                    "dataset_ratio"   : info.get("dataset_ratio"),
                })
            except Exception:
                continue

        return {"total_models": len(models), "models": models}

    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))


# ── 모델 상세 ─────────────────────────────────────────────────────────────

@router.get("/info/{model_id}")
async def get_model_info(model_id: str):
    """특정 모델 상세 정보 조회"""
    try:
        info = model_service.load_model(model_id)
        return {
            "model_id"        : model_id,
            "model_type"      : info["model_type"],
            "protein_name"    : info["protein_name"],
            "metrics"         : info["metrics"],
            "feature_count"   : info["feature_count"],
            "train_size"      : info["train_size"],
            "test_size_n"     : info.get("test_size_n"),
            "fingerprint_type": info.get("fingerprint_type"),
            "dataset_ratio"   : info.get("dataset_ratio"),
            "ignore3D"        : info.get("ignore3D"),
            "fs_method"       : info.get("fs_method", "unknown"),
            "n_features"      : info.get("n_features"),
            "has_transformer" : info.get("transformer") is not None,
            "has_selected_features": info.get("selected_features") is not None,
        }
    except FileNotFoundError:
        raise HTTPException(status_code=404, detail=f"Model not found: {model_id}")
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))
