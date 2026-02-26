# PCI Prediction API

## 프로젝트 개요

**Protein-Compound Interaction (PCI) Prediction API**는 푸드테크 연구를 위한 단백질-화합물 상호작용 예측 시스템입니다. ChEMBL과 BindingDB 데이터베이스에서 타겟 유전자 데이터를 수집하고, 분자 지문(Fingerprint) 및 분자 기술자(Molecular Descriptor)로 변환하여 머신러닝 모델을 학습시킵니다. 학습된 모델로 화합물의 활성을 예측하고 SHAP 분석을 통해 중요한 화학 특성을 추출합니다.

## 주요 기능

### 현재 구현된 기능

1. ✅ **데이터 수집** - ChEMBL 및 BindingDB에서 IC50 기준 PCI 데이터 수집
2. ✅ **특성 변환** - Fingerprint (ECFP4, MACCS, MORGAN) 및 Molecular Descriptor (Mordred 2D/3D) 계산
3. ✅ **모델 학습** - TabPFN, XGBoost, LightGBM, CatBoost, RandomForest 모델 학습
4. ✅ **예측** - 학습된 모델로 SMILES 리스트에 대한 활성 예측
5. ✅ **SHAP 분석** - 중요 화학 특성 추출 및 시각화
6. ✅ **RESTful API** - FastAPI 기반 표준 REST API
7. ✅ **자동 문서화** - Swagger/OpenAPI 자동 생성
8. ✅ **CORS 지원** - 프론트엔드 연동 가능

### 아직 구현되지 않은 기능

- ⏳ FooDB 데이터 통합 및 실제 식품 화합물 예측
- ⏳ 배치 예측 및 대용량 데이터 처리
- ⏳ 모델 앙상블 및 하이퍼파라미터 튜닝
- ⏳ 사용자 인증 및 API 키 관리

## API 엔드포인트

### 공개 URL

- **Production URL**: https://3000-itvs5ceiquncckj2qhd7g-2b54fc91.sandbox.novita.ai
- **API 문서**: https://3000-itvs5ceiquncckj2qhd7g-2b54fc91.sandbox.novita.ai/docs
- **ReDoc 문서**: https://3000-itvs5ceiquncckj2qhd7g-2b54fc91.sandbox.novita.ai/redoc
- **Health Check**: https://3000-itvs5ceiquncckj2qhd7g-2b54fc91.sandbox.novita.ai/health

### 주요 엔드포인트

#### 1. 데이터 수집
```bash
POST /api/data/collect
{
  "target_list": ["PDE4", "PDE5"],
  "standard_type": "IC50",
  "binding_db_folder": "optional"
}
```

#### 2. 특성 변환
```bash
POST /api/features/transform
{
  "protein_name": "PDE4",
  "fingerprint_type": "ECFP4",
  "dataset_type": "fewshot",
  "pos_threshold": 10000,
  "neg_threshold": 20000
}
```

#### 3. 모델 학습
```bash
POST /api/models/train
{
  "protein_name": "PDE4",
  "model_type": "XGBoost",
  "feature_type": "fingerprint",
  "test_size": 0.2,
  "random_state": 42
}
```

#### 4. 예측
```bash
POST /api/models/predict
{
  "smiles_list": ["CCO", "CC(=O)O"],
  "model_id": "PDE4_XGBoost_20260226_120000",
  "feature_type": "fingerprint"
}
```

#### 5. SHAP 분석
```bash
POST /api/shap/analyze
{
  "model_id": "PDE4_XGBoost_20260226_120000",
  "feature_type": "fingerprint",
  "top_n": 20
}
```

#### 6. 모델 리스트 조회
```bash
GET /api/models/list
```

## 데이터 구조

### 데이터베이스 및 저장소

- **ChEMBL**: 공개 화합물 활성 데이터베이스
- **BindingDB**: 단백질-화합물 결합 데이터베이스
- **저장 위치**: 
  - `saved_data/IC50/`: 수집된 원본 데이터
  - `raw/FewshotSet/`: Few-shot learning용 데이터
  - `raw/TransferSet/`: Transfer learning용 데이터
  - `models_trained/`: 학습된 모델 (.pkl 파일)
  - `shap_outputs/`: SHAP 분석 결과 플롯

### 데이터 플로우

```
1. 타겟 유전자 선택
   ↓
2. ChEMBL/BindingDB 데이터 수집 (IC50 기준)
   ↓
3. SMILES → Fingerprint/Descriptor 변환
   ↓
4. 데이터 정제 (결측값, 낮은 분산, 높은 상관관계 제거)
   ↓
5. ML 모델 학습 (TabPFN/AutoML)
   ↓
6. 예측 및 SHAP 분석
   ↓
7. 중요 화학 특성 추출
```

## 기술 스택

- **Framework**: FastAPI 0.133.1
- **Server**: Uvicorn (ASGI)
- **Chemistry**: RDKit 2025.9.5, Mordred 1.2.0, ChEMBL Web Client 0.10.9
- **Machine Learning**: 
  - XGBoost 3.2.0
  - LightGBM 4.6.0
  - CatBoost 1.2.10
  - Scikit-learn 1.6.1
  - TabPFN 0.1.10 (선택적)
- **Interpretability**: SHAP 0.50.0
- **Data Processing**: Pandas 2.2.3, NumPy 1.26.4
- **Process Manager**: PM2

## 설치 및 실행

### 1. 의존성 설치

```bash
cd /home/user/webapp
pip install -r requirements.txt
```

### 2. 서버 실행

#### PM2로 실행 (권장)
```bash
# 포트 정리
fuser -k 3000/tcp 2>/dev/null || true

# PM2로 시작
pm2 start ecosystem.config.cjs

# 로그 확인
pm2 logs pci-api --nostream

# 서버 재시작
pm2 restart pci-api

# 서버 중지
pm2 stop pci-api

# 서버 삭제
pm2 delete pci-api
```

#### 직접 실행 (개발용)
```bash
cd /home/user/webapp
uvicorn app.main:app --host 0.0.0.0 --port 3000 --reload
```

### 3. API 문서 확인

브라우저에서 다음 URL을 열어 자동 생성된 API 문서를 확인할 수 있습니다:

- Swagger UI: http://localhost:3000/docs
- ReDoc: http://localhost:3000/redoc

## 사용 가이드

### 워크플로우 예시

#### 1단계: 데이터 수집
```bash
curl -X POST "http://localhost:3000/api/data/collect" \
  -H "Content-Type: application/json" \
  -d '{
    "target_list": ["PDE4"],
    "standard_type": "IC50"
  }'
```

#### 2단계: 특성 변환
```bash
curl -X POST "http://localhost:3000/api/features/transform" \
  -H "Content-Type: application/json" \
  -d '{
    "protein_name": "PDE4",
    "fingerprint_type": "ECFP4",
    "dataset_type": "transfer",
    "pos_threshold": 10000,
    "neg_threshold": 20000
  }'
```

#### 3단계: 모델 학습
```bash
curl -X POST "http://localhost:3000/api/models/train" \
  -H "Content-Type: application/json" \
  -d '{
    "protein_name": "PDE4",
    "model_type": "XGBoost",
    "feature_type": "fingerprint",
    "test_size": 0.2,
    "random_state": 42
  }'
```

#### 4단계: 예측
```bash
curl -X POST "http://localhost:3000/api/models/predict" \
  -H "Content-Type: application/json" \
  -d '{
    "smiles_list": ["CCO", "CC(=O)Oc1ccccc1C(=O)O"],
    "model_id": "PDE4_XGBoost_20260226_120000",
    "feature_type": "fingerprint"
  }'
```

#### 5단계: SHAP 분석
```bash
curl -X POST "http://localhost:3000/api/shap/analyze" \
  -H "Content-Type: application/json" \
  -d '{
    "model_id": "PDE4_XGBoost_20260226_120000",
    "feature_type": "fingerprint",
    "top_n": 20
  }'
```

## 프로젝트 구조

```
webapp/
├── app/
│   ├── main.py                 # FastAPI 메인 애플리케이션
│   ├── models/
│   │   └── schemas.py          # Pydantic 스키마 정의
│   ├── routers/
│   │   ├── data_collection.py  # 데이터 수집 API
│   │   ├── feature_transform.py # 특성 변환 API
│   │   ├── model_training.py   # 모델 학습 및 예측 API
│   │   └── shap_analysis.py    # SHAP 분석 API
│   ├── services/
│   │   ├── data_collection.py  # 데이터 수집 비즈니스 로직
│   │   ├── feature_transform.py # 특성 변환 비즈니스 로직
│   │   ├── model_training.py   # 모델 학습 비즈니스 로직
│   │   └── shap_analysis.py    # SHAP 분석 비즈니스 로직
│   └── utils/
│       ├── molecular_descriptor.py # 분자 기술자 계산
│       └── utils.py            # ChEMBL/BindingDB 유틸리티
├── saved_data/                 # 수집된 데이터
│   ├── IC50/
│   └── BindingDB/
├── raw/                        # 변환된 학습 데이터
│   ├── FewshotSet/
│   └── TransferSet/
├── models_trained/             # 학습된 모델
├── shap_outputs/               # SHAP 분석 결과
├── ecosystem.config.cjs        # PM2 설정
├── requirements.txt            # Python 의존성
├── .gitignore                  # Git 제외 파일
└── README.md                   # 프로젝트 문서
```

## 배포 상태

- **플랫폼**: Sandbox Environment (개발/테스트용)
- **상태**: ✅ 활성 (Active)
- **마지막 업데이트**: 2026-02-26

## 다음 개발 단계

### 우선순위: 높음
1. **FooDB 통합**: 실제 식품 데이터베이스와 연동하여 식품 화합물 예측
2. **배치 처리**: 대량 SMILES 리스트 처리 최적화
3. **모델 성능 개선**: 하이퍼파라미터 튜닝 및 앙상블 모델

### 우선순위: 중간
4. **API 인증**: JWT 기반 사용자 인증 시스템
5. **데이터베이스**: PostgreSQL 또는 MongoDB로 데이터 영구 저장
6. **모니터링**: Prometheus + Grafana 모니터링 시스템

### 우선순위: 낮음
7. **프론트엔드**: React 기반 웹 인터페이스 구축
8. **도커화**: Docker 컨테이너화 및 Kubernetes 배포
9. **CI/CD**: GitHub Actions 자동 배포 파이프라인

## 참고 자료

- **ChEMBL**: https://www.ebi.ac.uk/chembl/
- **BindingDB**: https://www.bindingdb.org/
- **FooDB**: https://foodb.ca/
- **RDKit**: https://www.rdkit.org/
- **SHAP**: https://github.com/slundberg/shap
- **FastAPI**: https://fastapi.tiangolo.com/

## 라이선스

이 프로젝트는 연구 목적으로만 사용됩니다.

## 연락처

푸드테크 박사과정 연구 프로젝트
- 전문 분야: 푸드테크, 개인맞춤형식품, AI, 자동화
