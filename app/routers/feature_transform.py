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
    수집된 데이터를 Fingerprint + Descriptor로 변환 (원본 연구 방식)
    
    - **protein_name**: 단백질 이름
    - **fingerprint_type**: ECFP4, MACCS, MORGAN
    - **dataset_type**: 데이터셋 타입 (단일 dataset 사용)
    - **dataset_ratio**: 5x, 10x, 20x
    - **ignore3D**: True (2D만), False (3D 포함)
    - **pos_threshold**: 활성 화합물 IC50 임계값 (nM)
    - **neg_threshold**: 비활성 화합물 IC50 임계값 (nM)
    """
    try:
        # 캐시 확인: 이미 변환된 데이터가 있으면 재사용
        output_folder = "raw/Dataset"
        os.makedirs(output_folder, exist_ok=True)
        output_path = f"{output_folder}/{request.protein_name}.csv"
        
        # 캐시가 있고, 설정이 동일하면 재사용
        cache_valid = False
        if os.path.exists(output_path):
            try:
                cached_df = pd.read_csv(output_path)
                # Feature count로 설정 확인
                feature_count = len(cached_df.columns) - 3  # SMILES, Y, potency 제외
                expected_fp = 167 if request.fingerprint_type.value == "MACCS" else 1024
                
                # 동일한 설정이면 캐시 사용
                if abs(feature_count - expected_fp) < 100:  # descriptor 범위 고려
                    cache_valid = True
                    features_df = cached_df
                    print(f"✅ Using cached data: {output_path} ({feature_count} features)")
            except:
                pass
        
        if not cache_valid:
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
            
            # 1. 기본 데이터 변환 (SMILES + potency)
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
            
            # 데이터 품질 경고
            if len(result_df) < 100:
                print(f"⚠️ WARNING: Small dataset ({len(result_df)} compounds). Model may overfit!")
                print(f"   Recommendation: Collect more data for reliable predictions.")
            
            # 샘플링 제거: 모든 데이터 사용
            # (이전 코드는 불필요하게 데이터를 줄였습니다)
            
            # 2. Fingerprint + Descriptor 변환 (원본 연구 방식)
            # 컬럼명 확인 (SMILES 또는 canonical_SMILES)
            smiles_col = 'canonical_SMILES' if 'canonical_SMILES' in result_df.columns else 'SMILES'
            smiles_list = result_df[smiles_col].tolist()
            
            print(f"⏳ Computing features for {len(smiles_list)} compounds (this may take 1-2 minutes)...")
            
            # Fingerprint + Descriptor 생성
            features_df = feature_service.transform_to_fingerprint_with_descriptors(
                smiles_list=smiles_list,
                fp_type=request.fingerprint_type.value,
                dataset_ratio=request.dataset_ratio,
                ignore3D=request.ignore3D
            )
            
            # SMILES, Y(label), potency 컬럼 추가
            features_df['SMILES'] = result_df[smiles_col].values
            features_df['Y'] = result_df['Y'].values
            features_df['potency'] = result_df['potency'].values
            
            # 저장
            features_df.to_csv(output_path, index=False)
            print(f"✅ Saved to cache: {output_path}")
        
        # Fingerprint 크기 계산
        fp_size = 167 if request.fingerprint_type == "MACCS" else 1024
        descriptor_count = len(features_df.columns) - 3 - fp_size  # SMILES, Y, potency 제외
        
        return {
            "protein_name": request.protein_name,
            "fingerprint_type": request.fingerprint_type.value,
            "dataset_type": request.dataset_type,
            "dataset_ratio": request.dataset_ratio,
            "ignore3D": request.ignore3D,
            "total_compounds": len(features_df),
            "active_compounds": int((features_df['Y'] == 1).sum()),
            "inactive_compounds": int((features_df['Y'] == 0).sum()),
            "fingerprint_size": fp_size,
            "descriptor_count": descriptor_count,
            "feature_count": len(features_df.columns) - 3,  # SMILES, Y, potency 제외
            "output_path": output_path,
            "message": "Data transformation completed successfully"
        }
    
    except HTTPException:
        raise
    except Exception as e:
        import traceback
        traceback.print_exc()
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
