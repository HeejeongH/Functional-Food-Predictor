from pydantic import BaseModel, Field
from typing import List, Optional, Dict, Any
from enum import Enum


class FingerprintType(str, Enum):
    MACCS = "maccs"
    ECFP4 = "ecfp4"
    ECFP6 = "ecfp6"
    FCFP4 = "fcfp4"
    FCFP6 = "fcfp6"
    RDKIT = "rdkit"


class FeatureSelectionMethod(str, Enum):
    """노트북 Step3 와 동일한 FS 방법 목록"""
    MI     = "mi"
    RFE    = "rfe"
    SHAP   = "shap"
    PCA    = "pca"
    RANDOM = "random"
    NONE   = "none"


class DescriptorType(str, Enum):
    MORDRED_2D = "MORDRED_2D"
    MORDRED_3D = "MORDRED_3D"
    CUSTOM = "CUSTOM"


class MLModelType(str, Enum):
    TABPFN       = "TabPFN"
    XGBOOST      = "XGBoost"
    LIGHTGBM     = "LightGBM"
    CATBOOST     = "CatBoost"
    RANDOM_FOREST = "RandomForest"


# ── 요청 모델 ────────────────────────────────────────────────────────────

class FingerprintRequest(BaseModel):
    smiles_list:   List[str]       = Field(..., description="SMILES 문자열 리스트")
    fp_type:       FingerprintType = Field(default=FingerprintType.ECFP4, description="Fingerprint 타입")
    dataset_ratio: str             = Field(default="5x",  description="데이터셋 비율 (1x~50x)")
    ignore3D:      bool            = Field(default=True,  description="3D descriptor 무시 여부")


class DescriptorRequest(BaseModel):
    smiles_list:     List[str]     = Field(..., description="SMILES 문자열 리스트")
    descriptor_type: DescriptorType = Field(default=DescriptorType.MORDRED_2D, description="Descriptor 타입")


class TargetCollectionRequest(BaseModel):
    target_list:       List[str]    = Field(..., description="타겟 유전자 리스트")
    standard_type:     str          = Field(default="IC50", description="표준 타입 (IC50, Ki 등)")
    binding_db_folder: Optional[str] = Field(None, description="BindingDB 폴더 경로")


class DataTransformRequest(BaseModel):
    protein_name:    str                    = Field(..., description="단백질 이름")
    fingerprint_type: Optional[FingerprintType] = Field(FingerprintType.ECFP4, description="Fingerprint 타입")
    descriptor_type: Optional[DescriptorType]   = Field(None, description="Descriptor 타입")
    dataset_type:    str                    = Field(default="dataset", description="데이터셋 타입")
    dataset_ratio:   str                    = Field(default="5x",  description="데이터셋 비율")
    ignore3D:        bool                   = Field(default=True,  description="3D descriptor 무시 여부")
    pos_threshold:   float                  = Field(default=10000, description="활성 임계값 (nM)")
    neg_threshold:   float                  = Field(default=20000, description="비활성 임계값 (nM)")


class ModelTrainRequest(BaseModel):
    """
    노트북 final-model 셀과 동일한 학습 파라미터.
    기본값은 노트북 Step3/Step4 best 값에 맞춤.
    """
    protein_name:     str                      = Field(..., description="단백질 이름")
    model_type:       MLModelType              = Field(..., description="ML 모델 타입")
    feature_type:     str                      = Field(..., description="특성 타입 (fingerprint/descriptor)")
    fingerprint_type: FingerprintType          = Field(default=FingerprintType.ECFP4, description="Fingerprint 타입")
    dataset_ratio:    str                      = Field(default="20x", description="데이터셋 비율")
    ignore3D:         bool                     = Field(default=True,  description="3D descriptor 무시 여부")
    # ── Feature Selection (노트북 Step3/Step4) ──────────────────
    fs_method:        FeatureSelectionMethod   = Field(default=FeatureSelectionMethod.PCA,
                                                       description="Feature-selection 방법 (mi/rfe/shap/pca/random/none)")
    n_features:       int                      = Field(default=100, ge=1, le=2000,
                                                       description="선택할 feature 수 (노트북 기본 100)")
    # ── Train/Test split ────────────────────────────────────────
    test_size:        float                    = Field(default=0.1, gt=0.0, lt=1.0,
                                                       description="테스트 세트 비율 (노트북 기본 0.1)")
    random_state:     int                      = Field(default=42, description="랜덤 시드")


class PredictionRequest(BaseModel):
    smiles_list:  List[str] = Field(..., description="SMILES 리스트")
    model_id:     str       = Field(..., description="학습된 모델 ID")
    feature_type: str       = Field(..., description="특성 타입")


class SHAPAnalysisRequest(BaseModel):
    model_id:     str = Field(..., description="학습된 모델 ID")
    feature_type: str = Field(..., description="특성 타입")
    top_n:        int = Field(default=20, description="상위 N개 특성")


class FooDBPredictRequest(BaseModel):
    smiles_list: List[str] = Field(..., description="식품 화합물 SMILES 리스트")
    model_id:    str        = Field(..., description="학습된 모델 ID")


class DecoyGenerationRequest(BaseModel):
    protein_name:  str            = Field(..., description="타겟 단백질 이름")
    decoy_ratio:   float          = Field(default=50.0, description="활성 화합물 대비 비활성 화합물 비율")
    pos_threshold: float          = Field(default=10000, description="활성 임계값 (nM)")
    neg_threshold: float          = Field(default=20000, description="비활성 임계값 (nM)")
    zinc_db_path:  Optional[str]  = Field(None, description="ZINC 데이터베이스 파일 경로 (선택적)")


# ── 응답 모델 ────────────────────────────────────────────────────────────

class DataCollectionResponse(BaseModel):
    chembl_count:    int
    bindingdb_count: int
    total_compounds: int
    message:         str


class ModelTrainResponse(BaseModel):
    """
    노트북과 동일한 핵심 지표(F1, AUC, MCC) + 보조 지표 반환.
    """
    model_id:  str
    # 핵심 지표 (노트북 기준)
    f1:        float
    auc:       float
    mcc:       float
    # 보조 지표
    accuracy:  float
    precision: float
    recall:    float
    # 학습 설정 요약
    fs_method: str
    n_features: int
    message:   str


class PredictionResponse(BaseModel):
    predictions: List[Dict[str, Any]]
    model_info:  Dict[str, Any]


class SHAPAnalysisResponse(BaseModel):
    top_features:        List[Dict[str, Any]]
    shap_values_summary: Dict[str, Any]
    plot_path:           Optional[str]


class DecoyGenerationResponse(BaseModel):
    protein_name:         str
    num_actives:          int
    num_inactives_real:   int
    num_decoys_generated: int
    total_compounds:      int
    final_ratio:          str
    output_path:          str
    message:              str
