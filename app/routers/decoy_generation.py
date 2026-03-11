"""
DUD-E Decoy Generation API Router
"""

from fastapi import APIRouter, HTTPException
from app.models.schemas import DecoyGenerationRequest, DecoyGenerationResponse
from app.services.decoy_generation import DecoyGenerationService
import pandas as pd
import os

router = APIRouter(
    prefix="/api/decoy",
    tags=["decoy"],
    responses={404: {"description": "Not found"}}
)


@router.post("/generate", response_model=DecoyGenerationResponse)
async def generate_decoys(request: DecoyGenerationRequest):
    """
    DUD-E 방식으로 Decoy 생성
    
    **전략:**
    1. 실제 데이터베이스의 활성(Active) + 비활성(Inactive) 화합물 우선 사용
    2. 비활성 화합물이 부족한 경우, 활성 화합물과 유사한 물리화학적 특성을 가진 decoy 생성
    3. DUD-E 기준: 분자량, LogP, 회전 가능한 결합, 수소결합 donor/acceptor 등을 매칭
    
    **Parameters:**
    - **protein_name**: 타겟 단백질 이름 (예: "PDE4", "FTO")
    - **decoy_ratio**: 활성 화합물 대비 비활성 화합물 비율 (기본 50:1)
    - **pos_threshold**: 활성 임계값 (nM, 기본 10000)
    - **neg_threshold**: 비활성 임계값 (nM, 기본 20000)
    - **zinc_db_path**: ZINC 데이터베이스 파일 경로 (선택적)
    
    **Returns:**
    - 활성 화합물 개수
    - 실제 비활성 화합물 개수
    - 생성된 decoy 개수
    - 최종 데이터셋 정보
    """
    try:
        # 1. ChEMBL 데이터 로드
        chembl_path = f'saved_data/IC50/{request.protein_name}_summary.xlsx'
        
        if not os.path.exists(chembl_path):
            raise HTTPException(
                status_code=404, 
                detail=f"ChEMBL data not found for {request.protein_name}. "
                       f"Please run data collection first."
            )
        
        # Excel 파일 읽기
        chembl_df = pd.read_excel(chembl_path, sheet_name='ChemBL')
        
        if 'canonical_smiles' not in chembl_df.columns or 'standard_value' not in chembl_df.columns:
            raise HTTPException(
                status_code=400,
                detail="ChEMBL data must contain 'canonical_smiles' and 'standard_value' columns"
            )
        
        # 2. Decoy Generation Service 초기화
        decoy_service = DecoyGenerationService(decoy_ratio=request.decoy_ratio)
        
        # 3. ZINC 데이터베이스 로드 (선택적)
        if request.zinc_db_path and os.path.exists(request.zinc_db_path):
            decoy_service.load_zinc_database(request.zinc_db_path)
        
        # 4. Decoy 생성
        balanced_df = decoy_service.generate_balanced_dataset(
            chembl_df=chembl_df,
            target_protein=request.protein_name,
            pos_threshold=request.pos_threshold,
            neg_threshold=request.neg_threshold
        )
        
        # 5. 결과 저장
        output_dir = f'saved_data/IC50_with_decoys'
        os.makedirs(output_dir, exist_ok=True)
        output_path = f'{output_dir}/{request.protein_name}_balanced.csv'
        
        balanced_df.to_csv(output_path, index=False)
        
        # 6. 통계 계산
        num_actives = len(balanced_df[balanced_df['source'] == 'active'])
        num_inactives_real = len(balanced_df[balanced_df['source'] == 'inactive'])
        num_decoys = len(balanced_df[balanced_df['source'] == 'decoy'])
        total_compounds = len(balanced_df)
        final_ratio = f"1:{(num_inactives_real + num_decoys) / num_actives:.1f}" if num_actives > 0 else "N/A"
        
        return DecoyGenerationResponse(
            protein_name=request.protein_name,
            num_actives=num_actives,
            num_inactives_real=num_inactives_real,
            num_decoys_generated=num_decoys,
            total_compounds=total_compounds,
            final_ratio=final_ratio,
            output_path=output_path,
            message=f"Successfully generated {num_decoys} decoys for {request.protein_name}. "
                   f"Final dataset: {num_actives} actives + {num_inactives_real} inactives (real) + "
                   f"{num_decoys} decoys = {total_compounds} total compounds"
        )
    
    except HTTPException as he:
        raise he
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"Error generating decoys: {str(e)}")


@router.get("/list")
async def list_decoy_datasets():
    """
    생성된 decoy 데이터셋 목록 조회
    
    **Returns:**
    - 생성된 모든 decoy 데이터셋 정보
    """
    try:
        decoy_dir = 'saved_data/IC50_with_decoys'
        
        if not os.path.exists(decoy_dir):
            return {
                "datasets": [],
                "message": "No decoy datasets found"
            }
        
        datasets = []
        for filename in os.listdir(decoy_dir):
            if filename.endswith('_balanced.csv'):
                protein_name = filename.replace('_balanced.csv', '')
                file_path = os.path.join(decoy_dir, filename)
                
                # 데이터셋 정보 읽기
                df = pd.read_csv(file_path)
                
                num_actives = len(df[df['source'] == 'active'])
                num_inactives = len(df[df['source'] == 'inactive'])
                num_decoys = len(df[df['source'] == 'decoy'])
                
                datasets.append({
                    'protein_name': protein_name,
                    'num_actives': num_actives,
                    'num_inactives_real': num_inactives,
                    'num_decoys': num_decoys,
                    'total_compounds': len(df),
                    'file_path': file_path
                })
        
        return {
            "datasets": datasets,
            "count": len(datasets),
            "message": f"Found {len(datasets)} decoy datasets"
        }
    
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"Error listing decoy datasets: {str(e)}")


@router.get("/info/{protein_name}")
async def get_decoy_info(protein_name: str):
    """
    특정 단백질의 decoy 데이터셋 정보 조회
    
    **Parameters:**
    - **protein_name**: 타겟 단백질 이름
    
    **Returns:**
    - 상세 데이터셋 정보 및 통계
    """
    try:
        decoy_path = f'saved_data/IC50_with_decoys/{protein_name}_balanced.csv'
        
        if not os.path.exists(decoy_path):
            raise HTTPException(
                status_code=404,
                detail=f"Decoy dataset not found for {protein_name}"
            )
        
        df = pd.read_csv(decoy_path)
        
        # 상세 통계
        actives = df[df['source'] == 'active']
        inactives = df[df['source'] == 'inactive']
        decoys = df[df['source'] == 'decoy']
        
        info = {
            'protein_name': protein_name,
            'total_compounds': len(df),
            'actives': {
                'count': len(actives),
                'avg_standard_value': actives['standard_value'].mean() if len(actives) > 0 else 0,
                'min_standard_value': actives['standard_value'].min() if len(actives) > 0 else 0,
                'max_standard_value': actives['standard_value'].max() if len(actives) > 0 else 0
            },
            'inactives_real': {
                'count': len(inactives),
                'avg_standard_value': inactives['standard_value'].mean() if len(inactives) > 0 else 0
            },
            'decoys': {
                'count': len(decoys),
                'avg_standard_value': decoys['standard_value'].mean() if len(decoys) > 0 else 0
            },
            'ratio': f"1:{(len(inactives) + len(decoys)) / len(actives):.1f}" if len(actives) > 0 else "N/A",
            'file_path': decoy_path
        }
        
        return info
    
    except HTTPException as he:
        raise he
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"Error getting decoy info: {str(e)}")
