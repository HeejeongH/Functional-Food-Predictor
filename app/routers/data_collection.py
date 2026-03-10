from fastapi import APIRouter, HTTPException
from app.models.schemas import (
    TargetCollectionRequest, 
    DataCollectionResponse
)
from app.services.data_collection import DataCollectionService

router = APIRouter(prefix="/api/data", tags=["Data Collection"])

data_service = DataCollectionService()

@router.post("/collect", response_model=DataCollectionResponse)
async def collect_target_data(request: TargetCollectionRequest):
    """
    타겟 유전자에 대한 ChEMBL 및 BindingDB 데이터 수집
    
    - **target_list**: 타겟 유전자 리스트 (예: ["PDE4", "PDE5"])
    - **standard_type**: IC50, Ki 등의 표준 타입
    - **binding_db_folder**: BindingDB TSV 파일이 있는 폴더 경로 (선택)
    """
    try:
        result = data_service.collect_target_data(
            target_list=request.target_list,
            standard_type=request.standard_type,
            binding_db_folder=request.binding_db_folder
        )
        
        return DataCollectionResponse(
            chembl_count=result['chembl_count'],
            bindingdb_count=result['bindingdb_count'],
            total_compounds=result['total_compounds'],
            message=f"Data collected successfully. Saved to {result['output_path']}"
        )
    
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))

@router.get("/status")
async def get_data_status():
    """
    수집된 데이터 상태 확인
    """
    import os
    import glob
    
    ic50_files = glob.glob("saved_data/IC50/*.xlsx")
    bindingdb_files = glob.glob("saved_data/BindingDB/*/*.tsv")
    
    return {
        "ic50_data_files": len(ic50_files),
        "bindingdb_files": len(bindingdb_files),
        "ic50_files": [os.path.basename(f) for f in ic50_files],
        "status": "ready"
    }
