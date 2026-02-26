from pydantic import BaseModel, Field
from typing import List, Optional, Dict, Any
from enum import Enum

class FingerprintType(str, Enum):
    ECFP4 = "ECFP4"
    MACCS = "MACCS"
    MORGAN = "MORGAN"

class DescriptorType(str, Enum):
    MORDRED_2D = "MORDRED_2D"
    MORDRED_3D = "MORDRED_3D"
    CUSTOM = "CUSTOM"

class MLModelType(str, Enum):
    TABPFN = "TabPFN"
    XGBOOST = "XGBoost"
    LIGHTGBM = "LightGBM"
    CATBOOST = "CatBoost"
    RANDOM_FOREST = "RandomForest"

class TargetCollectionRequest(BaseModel):
    target_list: List[str] = Field(..., description="타겟 유전자 리스트")
    standard_type: str = Field(default="IC50", description="표준 타입 (IC50, Ki 등)")
    binding_db_folder: Optional[str] = Field(None, description="BindingDB 폴더 경로")

class DataTransformRequest(BaseModel):
    protein_name: str = Field(..., description="단백질 이름")
    fingerprint_type: Optional[FingerprintType] = Field(None, description="Fingerprint 타입")
    descriptor_type: Optional[DescriptorType] = Field(None, description="Descriptor 타입")
    dataset_type: str = Field(default="fewshot", description="fewshot 또는 transfer")
    pos_threshold: float = Field(default=10000, description="활성 임계값 (nM)")
    neg_threshold: float = Field(default=20000, description="비활성 임계값 (nM)")

class ModelTrainRequest(BaseModel):
    protein_name: str = Field(..., description="단백질 이름")
    model_type: MLModelType = Field(..., description="ML 모델 타입")
    feature_type: str = Field(..., description="특성 타입 (fingerprint/descriptor)")
    test_size: float = Field(default=0.2, description="테스트 세트 비율")
    random_state: int = Field(default=42, description="랜덤 시드")

class PredictionRequest(BaseModel):
    smiles_list: List[str] = Field(..., description="SMILES 리스트")
    model_id: str = Field(..., description="학습된 모델 ID")
    feature_type: str = Field(..., description="특성 타입")

class SHAPAnalysisRequest(BaseModel):
    model_id: str = Field(..., description="학습된 모델 ID")
    feature_type: str = Field(..., description="특성 타입")
    top_n: int = Field(default=20, description="상위 N개 특성")

class DataCollectionResponse(BaseModel):
    chembl_count: int
    bindingdb_count: int
    total_compounds: int
    message: str

class ModelTrainResponse(BaseModel):
    model_id: str
    accuracy: float
    precision: float
    recall: float
    f1_score: float
    roc_auc: float
    message: str

class PredictionResponse(BaseModel):
    predictions: List[Dict[str, Any]]
    model_info: Dict[str, Any]

class SHAPAnalysisResponse(BaseModel):
    top_features: List[Dict[str, Any]]
    shap_values_summary: Dict[str, Any]
    plot_path: Optional[str]
