# PCI Prediction API

## 프로젝트 개요

**Protein-Compound Interaction (PCI) Prediction API**는 푸드테크 연구를 위한 단백질-화합물 상호작용 예측 시스템입니다. ChEMBL과 BindingDB 데이터베이스에서 타겟 유전자 데이터를 수집하고, 분자 지문(Fingerprint) 및 분자 기술자(Molecular Descriptor)로 변환하여 머신러닝 모델을 학습시킵니다. 학습된 모델로 화합물의 활성을 예측하고 SHAP 분석을 통해 중요한 화학 특성을 추출합니다.

## 🐳 빠른 시작 (Docker 권장)

### Docker로 실행 (추천 ⭐)

```bash
# 1. 프로젝트 다운로드 및 압축 해제
wget https://www.genspark.ai/api/files/s/BWcmopJe -O pci-api.tar.gz
tar -xzf pci-api.tar.gz
cd webapp

# 2. Docker Compose로 실행 (원클릭 배포)
docker-compose up -d

# 3. 서비스 확인
curl http://localhost:3000/health
```

**Docker 실행 후 접속:**
- Web UI: http://localhost:3000
- API 문서: http://localhost:3000/docs
- Health Check: http://localhost:3000/health

**Docker 관리 명령어:**
```bash
docker-compose logs -f          # 로그 확인
docker-compose restart          # 재시작
docker-compose down             # 중지 및 삭제
docker-compose up -d --build    # 재빌드 후 실행
```

### PM2로 실행 (레거시 옵션)

Docker를 사용할 수 없는 환경에서만 사용하세요.

```bash
# 1. Python 가상환경 생성
python3 -m venv venv
source venv/bin/activate  # Windows: venv\Scripts\activate

# 2. 의존성 설치
pip install -r requirements.txt

# 3. PM2로 실행
pm2 start legacy/ecosystem.config.cjs

# 4. 로그 확인
pm2 logs pci-api --nostream
```

**주의**: PM2 방식은 RDKit/Mordred 설치 시 환경별 에러가 발생할 수 있습니다. Docker 사용을 강력히 권장합니다.

## 주요 기능

### 현재 구현된 기능

1. ✅ **데이터 수집** - ChEMBL 및 BindingDB에서 IC50 기준 PCI 데이터 수집
2. ✅ **특성 변환** - Fingerprint (ECFP4, MACCS, MORGAN) 및 Molecular Descriptor (Mordred 2D/3D) 계산
3. ✅ **Fingerprint + Descriptor 조합** - descriptor_selection.csv 기반 선택된 descriptor와 fingerprint 결합
4. ✅ **모델 학습** - TabPFN, XGBoost, LightGBM, CatBoost, RandomForest 모델 학습
5. ✅ **예측** - 학습된 모델로 SMILES 리스트에 대한 활성 예측
6. ✅ **SHAP 분석** - 중요 화학 특성 추출 및 시각화
7. ✅ **FooDB 예측** - 식품 화합물 CSV 업로드 및 활성 예측
8. ✅ **배치 예측** - 대량 화합물 처리 (100개 단위 배치)
9. ✅ **RESTful API** - FastAPI 기반 표준 REST API
10. ✅ **자동 문서화** - Swagger/OpenAPI 자동 생성
11. ✅ **CORS 지원** - 프론트엔드 연동 가능
12. ✅ **통합 워크플로우 UI** - 단계별 체크박스 선택 및 한번에 실행
13. ✅ **단일 데이터셋 구조** - Transfer/Fewshot 구분 제거, 모델 선택으로 학습 방식 결정

### 최근 업데이트 (2026-03-06)

**🧠 TabPFN 모델 지원**
- TabPFN 6.4.1 Few-shot Learning 모델 통합
- Hugging Face 인증 가이드 추가 (`DEPLOYMENT.md` 참고)
- 소량 데이터(50-1000개)에서 우수한 성능 제공
- 프론트엔드 UI에 TabPFN 옵션 추가

**🎯 워크플로우 개선** (2026-03-05)
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

### 최근 추가된 기능 (2026-03-08)

**🍎 FooDB 식품 화합물 예측**
- FooDB CSV 업로드 기능
- 식품 화합물 SMILES 직접 입력 예측
- 배치 예측 (대량 화합물 처리)
- 활성 확률 기준 정렬 및 상위 N개 추출

### 아직 구현되지 않은 기능

- ⏳ FooDB API 자동 수집 (현재는 CSV 업로드 방식)
- ⏳ 모델 앙상블 및 하이퍼파라미터 튜닝
- ⏳ 사용자 인증 및 API 키 관리
- ⏳ 식품-단백질 상호작용 네트워크 시각화

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
- **FooDB**: 식품 화합물 데이터베이스 (CSV 업로드 방식)
- **Descriptor Selection**: `descriptor_selection.csv` - 각 데이터셋(5x, 10x, 20x, 2D/3D)별 선택된 descriptor 리스트

### 폴더 구조 (Docker 볼륨 마운트)

```
webapp/
├── Dockerfile                    # Docker 빌드 설정
├── docker-compose.yml            # 원클릭 배포 설정
├── .dockerignore                 # Docker 빌드 최적화
├── descriptor_selection.csv      # 필수 파일 (RO 마운트)
├── saved_data/                   # 수집된 원본 데이터 (볼륨)
│   ├── IC50/                     # ChEMBL/BindingDB IC50 데이터
│   └── FooDB/                    # FooDB CSV 업로드 위치
│       ├── 3d_conformers/        # 3D conformer 캐시 (SDF)
│       └── preprocessed/         # 전처리된 데이터
├── raw/                          # 변환된 특성 데이터 (볼륨)
│   └── Dataset/                  # 통합 데이터셋
├── models_trained/               # 학습된 모델 (볼륨)
├── shap_outputs/                 # SHAP 분석 결과 (볼륨)
├── food_predictions/             # FooDB 예측 결과 (볼륨)
├── app/                          # FastAPI 애플리케이션
└── legacy/                       # PM2 레거시 설정
    ├── ecosystem.config.cjs
    └── README.md
```

**중요**: Docker 볼륨 마운트 방식으로 데이터가 컨테이너 외부에 보존됩니다.

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

- **Deployment**: Docker + Docker Compose (권장) / PM2 (레거시)
- **Framework**: FastAPI 0.109.0
- **Server**: Uvicorn (ASGI)
- **Chemistry**: RDKit 2023.9.4, Mordred 1.2.0, ChEMBL Web Client 0.10.8
- **Machine Learning**: 
  - XGBoost 2.0.3
  - LightGBM 4.3.0
  - CatBoost 1.2.2
  - Scikit-learn 1.4.0
  - **TabPFN 0.1.10** (소량 데이터 Few-shot 학습 전문, Hugging Face 인증 필요)
- **Interpretability**: SHAP 0.44.1
- **Data Processing**: Pandas 2.2.0, NumPy 1.26.3
- **Containerization**: Docker 20.10+, Docker Compose 2.0+

## 설치 및 실행

### ⭐ Docker 방식 (권장)

```bash
# 1. 프로젝트 다운로드
wget https://www.genspark.ai/api/files/s/BWcmopJe -O pci-api.tar.gz
tar -xzf pci-api.tar.gz
cd webapp

# 2. TabPFN 사용 시 환경변수 설정 (선택)
echo "HF_TOKEN=hf_your_token_here" > .env

# 3. Docker Compose 실행
docker-compose up -d

# 4. 로그 확인
docker-compose logs -f

# 5. 서비스 테스트
curl http://localhost:3000/health
```

**Docker 리소스 조정** (docker-compose.yml 수정):
```yaml
deploy:
  resources:
    limits:
      cpus: '16.0'     # CPU 코어 수 조정
      memory: 32G      # 메모리 조정
```

### 🔧 PM2 방식 (레거시)

Docker를 사용할 수 없는 환경에서만 사용하세요. 자세한 내용은 `legacy/README.md` 참조.

```bash
# 1. Python 가상환경 생성
python3 -m venv venv
source venv/bin/activate

# 2. 의존성 설치
pip install -r requirements.txt

# 3. PM2 실행
pm2 start legacy/ecosystem.config.cjs
```

### 0. TabPFN 사용 설정 (선택사항)

TabPFN은 Few-shot Learning에 특화된 모델로, 소량 데이터(50-1000개)에서 우수한 성능을 제공합니다.

**⚠️ 중요**: TabPFN을 사용하려면 Hugging Face 인증이 필요합니다.

#### Docker 환경에서 TabPFN 활성화

```bash
# 1. Hugging Face 토큰 발급 (https://huggingface.co/settings/tokens)
# 2. .env 파일 생성
echo "HF_TOKEN=hf_xxxxxxxxxxxx" > .env

# 3. Docker Compose로 실행 (자동으로 환경변수 주입)
docker-compose up -d
```

#### PM2 환경에서 TabPFN 활성화

```bash
# 1. TabPFN 패키지 설치 (아나콘다 가상환경 권장)
pip install tabpfn

# 2. Hugging Face 로그인
huggingface-cli login
# 또는
hf auth login

# 프롬프트가 나오면 Hugging Face 토큰 입력
# 토큰 생성: https://huggingface.co/settings/tokens

# 3. TabPFN 모델 접근 권한 요청
# 브라우저에서 https://huggingface.co/Prior-Labs/tabpfn_2_5 방문
# "Accept terms" 버튼 클릭하여 약관 동의

# 4. TabPFN 동작 확인
python3 -c "from tabpfn import TabPFNClassifier; print('TabPFN available!')"
```

**TabPFN vs 다른 모델 비교**

| 모델 | 데이터 크기 | 학습 속도 | Few-shot 성능 | 권장 사용 시나리오 |
|------|-----------|---------|--------------|-----------------|
| **TabPFN** | 50-1000개 | 매우 빠름 | ⭐⭐⭐⭐⭐ | 신규 단백질, 프로토타입, 소량 데이터 |
| **XGBoost** | 100+ 개 | 빠름 | ⭐⭐⭐ | 범용 예측, 충분한 데이터 |
| **LightGBM** | 1000+ 개 | 매우 빠름 | ⭐⭐ | 대용량 데이터, 속도 중시 |
| **CatBoost** | 100+ 개 | 중간 | ⭐⭐⭐ | 범주형 특성 많을 때 |
| **RandomForest** | 50+ 개 | 중간 | ⭐⭐ | 기본 baseline 모델 |

✅ **TabPFN 권장 상황**:
- 새로운 단백질에 대한 초기 탐색
- 데이터가 50-1000개 정도로 적을 때
- 빠른 프로토타이핑이 필요할 때
- Transfer Learning 효과를 극대화하고 싶을 때

❌ **TabPFN 비권장 상황**:
- 데이터가 10,000개 이상일 때
- GPU가 없는 환경 (CPU로도 작동하지만 느림)
- Hugging Face 인증이 어려운 환경

**📖 자세한 TabPFN 설정 가이드**: `DEPLOYMENT.md` 참고

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

#### 3단계: 모델 학습 (TabPFN 예시)
```bash
curl -X POST "http://localhost:3000/api/models/train" \
  -H "Content-Type: application/json" \
  -d '{
    "protein_name": "FTO",
    "model_type": "TabPFN",
    "feature_type": "fingerprint",
    "test_size": 0.2,
    "random_state": 42
  }'
```

**💡 Tip**: TabPFN 사용 시 먼저 Hugging Face 인증을 완료해야 합니다. 자세한 내용은 위의 "0. TabPFN 사용 설정" 섹션을 참고하세요.

#### 3단계: 모델 학습 (XGBoost 예시)
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
