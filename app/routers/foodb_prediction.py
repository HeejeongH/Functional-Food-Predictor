from fastapi import APIRouter, HTTPException, UploadFile, File, Query
from app.models.schemas import FooDBPredictRequest, PredictionResponse
from app.services.model_training import ModelTrainingService
from app.services.feature_transform import FeatureTransformService
from app.services.foodb_service import FooDBService
from app.services.foodb_preprocessing import FooDBPreprocessingService
import pandas as pd
import numpy as np
import os
from typing import List, Optional

router = APIRouter(prefix="/api/foodb", tags=["FooDB Prediction"])

model_service       = ModelTrainingService()
feature_service     = FeatureTransformService()
foodb_service       = FooDBService()
foodb_preprocessing = FooDBPreprocessingService()


# ── /fetch (비활성화) ─────────────────────────────────────────────────────

@router.get("/fetch")
async def fetch_compounds_from_foodb(
    query:     Optional[str] = Query(None, description="검색어"),
    food_name: Optional[str] = Query(None, description="식품 이름"),
    limit:     int           = Query(1000, ge=1, le=5000),
):
    """
    ⚠️ FooDB API 가져오기 (비활성화)

    FooDB는 공식 공개 REST API를 제공하지 않습니다.
    https://foodb.ca/downloads 에서 CSV를 직접 다운로드하세요.
    """
    raise HTTPException(
        status_code=501,
        detail={
            "error"   : "FooDB public API is not available",
            "solution": "Download CSV files from https://foodb.ca/downloads",
            "recommended_files": [
                {"name": "Compound.csv", "description": "70,000+ compounds", "size": "~50 MB"},
                {"name": "Food.csv",     "description": "900+ foods",         "size": "~5 MB"},
                {"name": "Content.csv",  "description": "compound-food relationships", "size": "~200 MB"},
            ],
            "alternative_endpoints": {
                "upload_csv"         : "POST /api/foodb/upload",
                "predict_from_csv"   : "POST /api/foodb/predict",
                "predict_from_smiles": "POST /api/foodb/predict-smiles",
            },
        },
    )


# ── /search-food (비활성화) ───────────────────────────────────────────────

@router.get("/search-food")
async def search_food_compounds(
    food_name: str = Query(..., description="식품 이름"),
    limit:     int = Query(100, ge=1, le=1000),
):
    """⚠️ 특정 식품의 화합물 검색 (비활성화)"""
    raise HTTPException(
        status_code=501,
        detail={
            "error"   : "FooDB search API is not available",
            "solution": "Download CSV from https://foodb.ca/downloads and filter locally",
        },
    )


# ── /upload ───────────────────────────────────────────────────────────────

@router.post("/upload")
async def upload_foodb_csv(file: UploadFile = File(...)):
    """FooDB CSV 파일 업로드"""
    try:
        foodb_dir = "saved_data/FooDB"
        os.makedirs(foodb_dir, exist_ok=True)

        file_path = os.path.join(foodb_dir, file.filename)
        contents  = await file.read()
        with open(file_path, "wb") as f:
            f.write(contents)

        df = pd.read_csv(file_path)

        smiles_cols = [c for c in df.columns if "smiles" in c.lower()]
        if not smiles_cols:
            raise HTTPException(
                status_code=400,
                detail="CSV must contain a SMILES column (e.g. 'canonical_smiles', 'smiles')",
            )

        return {
            "message"       : "FooDB CSV uploaded successfully",
            "file_path"     : file_path,
            "total_compounds": len(df),
            "columns"       : df.columns.tolist(),
            "smiles_column" : smiles_cols[0],
        }

    except HTTPException:
        raise
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"Upload failed: {e}")


# ── /preprocess ───────────────────────────────────────────────────────────

@router.post("/preprocess")
async def preprocess_foodb(
    foodb_file:       str,
    protein_name:     str,
    fingerprint_type: str  = "ecfp4",
    dataset_ratio:    str  = "20x",
    ignore3D:         bool = True,
):
    """
    FooDB 전처리 (Fingerprint + Mordred Descriptor, 청크+체크포인트 지원)

    ⚠️ 3D descriptor 계산은 오래 걸립니다 (1,000개 ≈ 10~30분).
    """
    try:
        descriptor_selection_path = "descriptor_selection.csv"
        if not os.path.exists(descriptor_selection_path):
            raise HTTPException(
                status_code=404,
                detail="descriptor_selection.csv not found. Please upload this file first.",
            )

        sel_df    = pd.read_csv(descriptor_selection_path)
        ignore_str = "True" if ignore3D else "False"
        col_key   = f"descriptors_filtered_{protein_name}_training_{dataset_ratio}_ignore3D_{ignore_str}.csv"

        if col_key not in sel_df.columns:
            raise HTTPException(
                status_code=404,
                detail=f"Descriptor column '{col_key}' not found in descriptor_selection.csv",
            )

        descriptor_names = sel_df[col_key].dropna().tolist()
        if not descriptor_names:
            raise HTTPException(
                status_code=400,
                detail=f"No descriptors found for {col_key}",
            )

        output_path = (
            f"saved_data/FooDB/preprocessed/"
            f"{protein_name}_{fingerprint_type}_{dataset_ratio}_ignore3D_{ignore3D}.csv"
        )
        checkpoint_path = (
            f"saved_data/FooDB/checkpoints/"
            f"{protein_name}_{dataset_ratio}_ignore3D_{ignore3D}_ckpt.csv"
        )
        os.makedirs(os.path.dirname(output_path),    exist_ok=True)
        os.makedirs(os.path.dirname(checkpoint_path), exist_ok=True)

        result_df = foodb_preprocessing.preprocess_foodb(
            foodb_csv_path             = foodb_file,
            protein_name               = protein_name,
            fingerprint_type           = fingerprint_type,
            dataset_ratio              = dataset_ratio,
            ignore3D                   = ignore3D,
            descriptor_names           = descriptor_names,
            output_path                = output_path,
            descriptor_checkpoint_path = checkpoint_path,
        )

        fp_size = 166 if fingerprint_type.lower() == "maccs" else 1024   # FP0 제외
        return {
            "message"         : "FooDB preprocessing completed successfully",
            "input_file"      : foodb_file,
            "output_file"     : output_path,
            "total_compounds" : len(result_df),
            "fingerprint_type": fingerprint_type,
            "fingerprint_size": fp_size,
            "descriptor_count": len(descriptor_names),
            "total_features"  : fp_size + len(descriptor_names),
            "checkpoint"      : checkpoint_path,
        }

    except HTTPException:
        raise
    except Exception as e:
        import traceback; traceback.print_exc()
        raise HTTPException(status_code=500, detail=f"Preprocessing failed: {e}")


# ── /predict ─────────────────────────────────────────────────────────────

@router.post("/predict")
async def predict_foodb_compounds(
    model_id:   str = Query(..., description="사용할 모델 ID"),
    foodb_file: str = Query("saved_data/FooDB/foodb_compounds.csv",
                            description="FooDB CSV 경로"),
    batch_size: int = Query(500,  ge=1,  le=5000, description="배치 크기"),
    top_n:      int = Query(100,  ge=1,  le=10000, description="반환할 상위 N개"),
):
    """
    FooDB CSV 전체 배치 예측.

    모델에 저장된 fs_method / transformer / selected_features 를 자동으로 적용.
    (노트북 FooDB 예측 셀과 동일한 pca / mi / shap / rfe / random 분기)
    """
    try:
        model_info = model_service.load_model(model_id)

        if not os.path.exists(foodb_file):
            raise HTTPException(
                status_code=404,
                detail=f"FooDB file not found: {foodb_file}. Please upload first.",
            )

        foodb_df  = pd.read_csv(foodb_file)
        smiles_col = None
        for col in ["moldb_smiles", "canonical_SMILES", "canonical_smiles",
                    "smiles", "SMILES"]:
            if col in foodb_df.columns:
                smiles_col = col
                break

        if not smiles_col:
            raise HTTPException(
                status_code=400,
                detail=f"No SMILES column found. Columns: {foodb_df.columns.tolist()}",
            )

        print(f"[FooDB /predict] model={model_id}  compounds={len(foodb_df)}  "
              f"smiles_col={smiles_col}")

        # FooDBService.predict_batch → FS 분기 포함
        result = foodb_service.predict_batch(
            model_info  = model_info,
            smiles_list = foodb_df[smiles_col].tolist(),
            batch_size  = batch_size,
        )

        result_df = foodb_df.copy()
        result_df["prediction"]        = result["predictions"]
        result_df["probability_active"] = result["probabilities_active"]
        result_df["probability_inactive"] = result["probabilities_inactive"]
        result_df = result_df.sort_values("probability_active", ascending=False)

        output_dir  = "food_predictions"
        os.makedirs(output_dir, exist_ok=True)
        output_path = os.path.join(output_dir, f"{model_id}_foodb_predictions.csv")
        result_df.to_csv(output_path, index=False)

        active_count   = int((result_df["prediction"] == 1).sum())
        inactive_count = int((result_df["prediction"] == 0).sum())

        return {
            "message"            : "FooDB prediction completed",
            "model_id"           : model_id,
            "fs_method"          : model_info.get("fs_method", "unknown"),
            "n_features"         : model_info.get("n_features"),
            "total_compounds"    : len(result_df),
            "active_compounds"   : active_count,
            "inactive_compounds" : inactive_count,
            "top_active_compounds": result_df.head(top_n).to_dict("records"),
            "output_path"        : output_path,
            "statistics": {
                "mean_probability_active": float(result_df["probability_active"].mean()),
                "max_probability_active" : float(result_df["probability_active"].max()),
                "min_probability_active" : float(result_df["probability_active"].min()),
            },
        }

    except HTTPException:
        raise
    except Exception as e:
        import traceback; traceback.print_exc()
        raise HTTPException(status_code=500, detail=f"FooDB prediction failed: {e}")


# ── /predict-smiles ───────────────────────────────────────────────────────

@router.post("/predict-smiles")
async def predict_food_smiles(request: FooDBPredictRequest):
    """
    식품 화합물 SMILES 직접 입력 예측.

    모델에 저장된 fs_method / transformer / selected_features 를 자동으로 적용.
    (노트북 FooDB/Phyto 예측 셀과 동일한 FS 분기)
    """
    try:
        model_info = model_service.load_model(request.model_id)

        # FooDBService.predict_batch — FS 분기 포함
        result = foodb_service.predict_batch(
            model_info  = model_info,
            smiles_list = request.smiles_list,
            batch_size  = len(request.smiles_list),  # 소규모면 한 배치로
        )

        results = []
        for i, smiles in enumerate(request.smiles_list):
            results.append({
                "smiles"             : smiles,
                "prediction"         : int(result["predictions"][i]),
                "prediction_label"   : "Active" if result["predictions"][i] == 1 else "Inactive",
                "probability_active" : float(result["probabilities_active"][i]),
                "probability_inactive": float(result["probabilities_inactive"][i]),
            })

        return {
            "predictions": results,
            "model_info" : {
                "model_id"       : request.model_id,
                "model_type"     : model_info["model_type"],
                "protein_name"   : model_info["protein_name"],
                "fingerprint_type": model_info.get("fingerprint_type"),
                "dataset_ratio"  : model_info.get("dataset_ratio"),
                "ignore3D"       : model_info.get("ignore3D"),
                "fs_method"      : model_info.get("fs_method", "unknown"),
                "n_features"     : model_info.get("n_features"),
            },
        }

    except FileNotFoundError:
        raise HTTPException(
            status_code=404,
            detail=f"Model not found: {request.model_id}",
        )
    except ValueError as e:
        raise HTTPException(status_code=400, detail=str(e))
    except Exception as e:
        import traceback; traceback.print_exc()
        raise HTTPException(status_code=500, detail=f"Prediction failed: {e}")
