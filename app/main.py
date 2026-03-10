from fastapi import FastAPI
from fastapi.middleware.cors import CORSMiddleware
from fastapi.responses import JSONResponse, FileResponse
from fastapi.staticfiles import StaticFiles
import uvicorn
import os

# 라우터 임포트
from app.routers import (
    data_collection,
    feature_transform,
    model_training,
    shap_analysis,
    foodb_prediction
)

# FastAPI 앱 생성
app = FastAPI(
    title="PCI Prediction API",
    description="""
    ## Protein-Compound Interaction (PCI) Prediction API
    
    이 API는 푸드테크 연구를 위한 단백질-화합물 상호작용 예측 시스템입니다.
    
    ### 주요 기능:
    1. **데이터 수집**: ChEMBL 및 BindingDB에서 타겟 유전자 데이터 수집
    2. **특성 변환**: Fingerprint 및 Molecular Descriptor 계산
    3. **모델 학습**: TabPFN, XGBoost, LightGBM, CatBoost 등 AutoML 모델 학습
    4. **예측**: 학습된 모델로 화합물 활성 예측
    5. **SHAP 분석**: 중요 화학 특성 추출 및 해석
    6. **FooDB 예측**: 식품 화합물 데이터베이스 활성 예측
    
    ### 워크플로우:
    1. `/api/data/collect` - 타겟 유전자 데이터 수집
    2. `/api/features/transform` - 데이터 변환 및 정제
    3. `/api/models/train` - 모델 학습
    4. `/api/models/predict` - 예측 수행
    5. `/api/shap/analyze` - SHAP 분석
    """,
    version="1.0.0",
    docs_url="/docs",
    redoc_url="/redoc"
)

# CORS 설정
app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],  # 프로덕션에서는 특정 도메인으로 제한
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

# 정적 파일 서빙
app.mount("/static", StaticFiles(directory="static"), name="static")

# 라우터 등록
app.include_router(data_collection.router)
app.include_router(feature_transform.router)
app.include_router(model_training.router)
app.include_router(shap_analysis.router)
app.include_router(foodb_prediction.router)

# 루트 엔드포인트 - 웹사이트 메인 페이지
@app.get("/")
async def root():
    """
    웹사이트 메인 페이지
    """
    return FileResponse('static/index.html')

# 헬스 체크
@app.get("/health")
async def health_check():
    """
    서버 상태 확인
    """
    import glob
    
    return {
        "status": "healthy",
        "services": {
            "data_collection": "available",
            "feature_transform": "available",
            "model_training": "available",
            "shap_analysis": "available",
            "foodb_prediction": "available"
        },
        "data": {
            "collected_datasets": len(glob.glob("saved_data/IC50/*.xlsx")),
            "trained_models": len(glob.glob("models_trained/*.pkl"))
        }
    }

# 에러 핸들러
@app.exception_handler(Exception)
async def global_exception_handler(request, exc):
    return JSONResponse(
        status_code=500,
        content={
            "message": "Internal server error",
            "detail": str(exc)
        }
    )

if __name__ == "__main__":
    # 필요한 디렉토리 생성
    os.makedirs("saved_data/IC50", exist_ok=True)
    os.makedirs("saved_data/BindingDB", exist_ok=True)
    os.makedirs("saved_data/FooDB", exist_ok=True)
    os.makedirs("raw/FewshotSet", exist_ok=True)
    os.makedirs("raw/TransferSet", exist_ok=True)
    os.makedirs("raw/descriptors", exist_ok=True)
    os.makedirs("models_trained", exist_ok=True)
    os.makedirs("shap_outputs", exist_ok=True)
    os.makedirs("food_predictions", exist_ok=True)
    
    # 서버 실행
    uvicorn.run(
        "main:app",
        host="0.0.0.0",
        port=3000,
        reload=True,
        log_level="info"
    )
