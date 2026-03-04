from fastapi import APIRouter, HTTPException
from app.models.schemas import DataTransformRequest, FingerprintRequest, DescriptorRequest
from app.services.feature_transform import FeatureTransformService
import pandas as pd
import os

router = APIRouter(prefix="/api/features", tags=["Feature Transformation"])

feature_service = FeatureTransformService()

@router.post("/transform")
async def transform_features(request: DataTransformRequest):
    """
    수집된 데이터를 Fingerprint 또는 Molecular Descriptor로 변환
    
    - **protein_name**: 단백질 이름
    - **fingerprint_type**: ECFP4, MACCS, MORGAN (선택)
    - **descriptor_type**: MORDRED_2D, MORDRED_3D, CUSTOM (선택)
    - **dataset_type**: fewshot 또는 transfer
    - **pos_threshold**: 활성 화합물 IC50 임계값 (nM)
    - **neg_threshold**: 비활성 화합물 IC50 임계값 (nM)
    """
    try:
        # 수집된 데이터 로드
        data_path = f"saved_data/IC50/{request.protein_name}_summary.xlsx"
        
        if not os.path.exists(data_path):
            raise HTTPException(
                status_code=404, 
                detail=f"Data not found for protein: {request.protein_name}. Please collect data first."
            )
        
        chembl_df = pd.read_excel(data_path, sheet_name='ChemBL')
        
        try:
            bindingdb_df = pd.read_excel(data_path, sheet_name='BindingDB')
        except:
            bindingdb_df = pd.DataFrame()
        
        # 데이터 변환
        result_df = feature_service.prepare_training_data(
            chembl_df=chembl_df,
            bindingdb_df=bindingdb_df,
            protein_name=request.protein_name,
            dataset_type=request.dataset_type,
            pos_threshold=request.pos_threshold,
            neg_threshold=request.neg_threshold
        )
        
        if result_df is None or len(result_df) == 0:
            raise HTTPException(
                status_code=400,
                detail="Failed to transform data. Check input data quality."
            )
        
        return {
            "protein_name": request.protein_name,
            "dataset_type": request.dataset_type,
            "total_compounds": len(result_df),
            "active_compounds": int((result_df['Y'] == 1).sum()),
            "inactive_compounds": int((result_df['Y'] == 0).sum()),
            "feature_count": len(result_df.columns) - 3,  # SMILES, Y, potency 제외
            "message": "Data transformation completed successfully"
        }
    
    except HTTPException:
        raise
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))

@router.post("/fingerprint")
async def generate_fingerprint(request: FingerprintRequest):
    """
    SMILES 리스트를 Fingerprint + 선택된 Descriptor로 변환
    (원본 연구 코드 방식: MACCS 167개 + 선택된 Descriptor 10~16개)
    
    - **smiles_list**: SMILES 문자열 리스트
    - **fp_type**: ECFP4 (1024비트), MACCS (167비트), MORGAN (1024비트)
    - **dataset_ratio**: 5x, 10x, 20x (descriptor 선택에 영향)
    - **ignore3D**: True (2D descriptor만), False (3D descriptor 포함)
    
    예시:
    - MACCS + 5x + 2D: 167 (MACCS) + 16개 (descriptor) = 183개
    - ECFP4 + 10x + 2D: 1024 (ECFP4) + 15개 (descriptor) = 1039개
    """
    try:
        result_df = feature_service.transform_to_fingerprint_with_descriptors(
            smiles_list=request.smiles_list,
            fp_type=request.fp_type.value,
            dataset_ratio=request.dataset_ratio,
            ignore3D=request.ignore3D
        )
        
        fp_size = 167 if request.fp_type.value == "MACCS" else 1024
        descriptor_count = len(result_df.columns) - fp_size
        
        return {
            "fingerprint_type": request.fp_type.value,
            "dataset_ratio": request.dataset_ratio,
            "ignore3D": request.ignore3D,
            "count": len(result_df),
            "fingerprint_size": fp_size,
            "descriptor_count": descriptor_count,
            "total_feature_count": len(result_df.columns),
            "features_preview": result_df.head(3).to_dict('records') if len(result_df) > 0 else []
        }
    
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))

@router.post("/descriptors")
async def generate_descriptors(request: DescriptorRequest):
    """
    SMILES 리스트를 Molecular Descriptors로 변환
    
    - **smiles_list**: SMILES 문자열 리스트
    - **descriptor_type**: MORDRED_2D, MORDRED_3D, CUSTOM
    """
    try:
        desc_df = feature_service.transform_to_descriptors(
            smiles_list=request.smiles_list,
            descriptor_type=request.descriptor_type.value
        )
        
        return {
            "descriptor_type": request.descriptor_type.value,
            "count": len(desc_df),
            "descriptor_count": len(desc_df.columns) - 1,  # canonical_SMILES 제외
            "descriptors": desc_df.to_dict('records')
        }
    
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))
