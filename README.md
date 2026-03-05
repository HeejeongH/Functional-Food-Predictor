# PCI Prediction API

## 프로젝트 개요

**Protein-Compound Interaction (PCI) Prediction API**는 푸드테크 연구를 위한 단백질-화합물 상호작용 예측 시스템입니다. ChEMBL과 BindingDB 데이터베이스에서 타겟 유전자 데이터를 수집하고, 분자 지문(Fingerprint) 및 분자 기술자(Molecular Descriptor)로 변환하여 머신러닝 모델을 학습시킵니다. 학습된 모델로 화합물의 활성을 예측하고 SHAP 분석을 통해 중요한 화학 특성을 추출합니다.

## 주요 기능

### 현재 구현된 기능

1. ✅ **데이터 수집** - ChEMBL 및 BindingDB에서 IC50 기준 PCI 데이터 수집
2. ✅ **특성 변환** - Fingerprint (ECFP4, MACCS, MORGAN) 및 Molecular Descriptor (Mordred 2D/3D) 계산
3. ✅ **Fingerprint + Descriptor 조합** - descriptor_selection.csv 기반 선택된 descriptor와 fingerprint 결합
4. ✅ **모델 학습** - TabPFN, XGBoost, LightGBM, CatBoost, RandomForest 모델 학습
5. ✅ **예측** - 학습된 모델로 SMILES 리스트에 대한 활성 예측
6. ✅ **SHAP 분석** - 중요 화학 특성 추출 및 시각화
7. ✅ **RESTful API** - FastAPI 기반 표준 REST API
8. ✅ **자동 문서화** - Swagger/OpenAPI 자동 생성
9. ✅ **CORS 지원** - 프론트엔드 연동 가능
10. ✅ **통합 워크플로우 UI** - 단계별 체크박스 선택 및 한번에 실행
11. ✅ **단일 데이터셋 구조** - Transfer/Fewshot 구분 제거, 모델 선택으로 학습 방식 결정

### 최근 업데이트 (2026-03-05)

**🎯 워크플로우 개선**
- Step 1: 데이터 수집 → 단백질 목록 확인 버튼
- Step 2-5: 수집된 단백질만 드롭다운 선택 가능
- Step 3: 모델 학습 후 자동 저장 및 "AUTO" 선택
- Step 4-5: 방금 학습한 모델 자동 사용 (또는 기존 모델 선택)

**📁 데이터셋 구조 단순화**
- `raw/Dataset/` - 모든 데이터를 단일 폴더에 저장
- Transfer Learning vs Few-shot Learning 구분 제거
- 학습 방식은 모델 선택으로 결정:
  - **TabPFN** → Few-shot Learning (소량 데이터 특화)
  - **XGBoost/LightGBM/CatBoost** → 일반 학습

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

**MACCS Fingerprint + Descriptor (권장)**
```bash
POST /api/features/fingerprint
{
  "smiles_list": ["CCO", "CC(=O)O"],
  "fp_type": "MACCS",
  "dataset_ratio": "5x",
  "ignore3D": true
}

# 응답 예시:
{
  "fingerprint_type": "MACCS",
  "dataset_ratio": "5x",
  "ignore3D": true,
  "count": 2,
  "fingerprint_size": 167,
  "descriptor_count": 16,
  "total_feature_count": 183  # 167 (MACCS) + 16 (선택된 descriptor)
}
```

**ECFP4 Fingerprint + Descriptor**
```bash
POST /api/features/fingerprint
{
  "smiles_list": ["CCO", "c1ccccc1"],
  "fp_type": "ECFP4",
  "dataset_ratio": "10x",
  "ignore3D": false
}

# 응답 예시:
{
  "fingerprint_type": "ECFP4",
  "dataset_ratio": "10x",
  "ignore3D": false,
  "count": 2,
  "fingerprint_size": 1024,
  "descriptor_count": 16,
  "total_feature_count": 1040  # 1024 (ECFP4) + 16 (선택된 descriptor)
}
```

**Transfer/Fewshot용 데이터 변환**
```bash
POST /api/features/transform
{
  "protein_name": "PDE4",
  "fingerprint_type": "MACCS",
  "dataset_type": "transfer",
  "dataset_ratio": "5x",
  "ignore3D": true,
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
- **Descriptor Selection**: `descriptor_selection.csv` - 각 데이터셋(5x, 10x, 20x, 2D/3D)별 선택된 descriptor 리스트
- **저장 위치**: 
  - `saved_data/IC50/`: 수집된 원본 데이터
  - `raw/FewshotSet/`: Few-shot learning용 데이터
  - `raw/TransferSet/`: Transfer learning용 데이터
  - `models_trained/`: 학습된 모델 (.pkl 파일)
  - `shap_outputs/`: SHAP 분석 결과 플롯
  - `descriptor_selection.csv`: 데이터셋별 descriptor 선택 정보

### 데이터 플로우

```
1. 타겟 유전자 선택
   ↓
2. ChEMBL/BindingDB 데이터 수집 (IC50 기준)
   ↓
3. SMILES → Fingerprint (MACCS/ECFP4/MORGAN) 변환
   ↓
4. descriptor_selection.csv에서 선택된 Descriptor 추출
   ↓
5. Fingerprint + Descriptor 결합 (예: MACCS 167비트 + 16 descriptor = 183 특성)
   ↓
6. 데이터 정제 (결측값, 이상값 처리)
   ↓
7. ML 모델 학습 (TabPFN/XGBoost/AutoML)
   ↓
8. 예측 및 SHAP 분석
   ↓
9. 중요 화학 특성 추출
```

### Descriptor Selection 방식

`descriptor_selection.csv` 파일은 각 데이터셋 조합별로 최적화된 descriptor 리스트를 포함합니다:

- `descriptors_filtered_FTO_training_5x_ignore3D_True.csv`: Transfer Learning 5x 2D descriptor
- `descriptors_filtered_FTO_training_5x_ignore3D_False.csv`: Transfer Learning 5x 3D descriptor
- `descriptors_filtered_FTO_training_10x_ignore3D_True.csv`: Transfer Learning 10x 2D descriptor
- `descriptors_filtered_FTO_training_10x_ignore3D_False.csv`: Transfer Learning 10x 3D descriptor
- (각 조합별 20x 데이터셋도 동일하게 존재)

**API 사용 시 자동으로 적절한 descriptor 세트를 선택하여 적용합니다.**

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
- **마지막 업데이트**: 2026-03-04
- **주요 개선사항**:
  - ✅ descriptor_selection.csv 기반 특성 선택 시스템 구축
  - ✅ Fingerprint + Descriptor 조합 API 완성
  - ✅ MACCS/ECFP4/MORGAN 모든 fingerprint 타입 지원
  - ✅ 2D/3D descriptor 및 5x/10x/20x 데이터셋 비율 지원
  - ✅ 웹 UI를 통한 사용자 친화적 인터페이스 제공

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
