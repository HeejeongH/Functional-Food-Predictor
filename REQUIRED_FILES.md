# 📋 프로젝트 실행을 위해 필요한 파일 목록

## 🔴 필수 파일 (반드시 필요)

### 1. **descriptor_selection.csv**
**위치**: `/home/user/webapp/descriptor_selection.csv` (프로젝트 루트)

**내용**: 각 단백질/데이터셋 비율/2D-3D 조합별로 선택된 Mordred descriptor 리스트

**형식**:
```csv
descriptors_filtered_FTO_training_5x_ignore3D_False.csv,descriptors_filtered_FTO_training_5x_ignore3D_True.csv,...
ATSC4c,ATSC7dv,...
ATSC6Z,ATSC6Z,...
ATSC8i,ATSC8i,...
...
```

**컬럼 명명 규칙**:
- `descriptors_filtered_{PROTEIN}_training_{RATIO}_ignore3D_{True/False}.csv`
- 예: `descriptors_filtered_FTO_training_20x_ignore3D_True.csv`

**역할**: FooDB 전처리 시 어떤 descriptor를 계산할지 결정

**✅ 현재 상태**: 이미 존재함 (17줄, 6개 컬럼)

---

## 🟡 선택적 파일 (기능에 따라 필요)

### 2. **FooDB CSV 파일**
**위치**: `saved_data/FooDB/` (사용자가 업로드)

**형식**: CSV 파일, 필수 컬럼: `canonical_smiles` 또는 `smiles`

**예시**:
```csv
id,canonical_smiles,name,moldb_mono_mass
FDB000001,CCO,Ethanol,46.0419
FDB000002,CC(=O)O,Acetic acid,60.0211
FDB000003,c1ccccc1,Benzene,78.0469
...
```

**다운로드**: https://foodb.ca/downloads → `Compounds.csv`

**역할**: 식품 화합물 데이터베이스, 예측할 SMILES 리스트

**📝 참고**: 
- 원본 FooDB CSV는 70,000+ 화합물 포함
- 필요한 화합물만 필터링하여 사용 권장
- 예시 코드의 `filtered_foodb_{ratio}.csv`는 이미 필터링된 버전

---

### 3. **학습 데이터 (ChEMBL/BindingDB)**
**위치**: `saved_data/IC50/{PROTEIN}_summary.xlsx`

**형식**: Excel 파일, 시트: `ChemBL`, `BindingDB` (선택)

**필수 컬럼**:
- `canonical_smiles`: SMILES 구조
- `standard_value`: IC50 값 (nM)
- `standard_relation`: 관계 (=, <, >, <=, >=)
- `target_pref_name`: 단백질 이름
- `molecule_chembl_id`: ChEMBL ID

**생성 방법**: API Step 1 "데이터 수집" 실행 시 자동 생성

**✅ 현재 상태**: `saved_data/IC50/FTO_summary.xlsx` 존재 (46개 화합물)

---

### 4. **변환된 특성 데이터**
**위치**: `raw/Dataset/{PROTEIN}.csv`

**형식**: CSV 파일

**컬럼 구조**:
- `X0, X1, X2, ..., X{N}`: Fingerprint bits + Descriptor values
- `SMILES`: 화합물 구조
- `Y`: 레이블 (0=비활성, 1=활성)
- `potency`: IC50 값

**생성 방법**: API Step 2 "특성 변환" 실행 시 자동 생성

**✅ 현재 상태**: 
- `raw/Dataset/FTO.csv` 존재 (40개 화합물, 1041 특성)
- `raw/TransferSet/FTO.csv` 존재 (39개 화합물, 181 특성)

---

### 5. **학습된 모델**
**위치**: `models_trained/{PROTEIN}_{MODEL}_{TIMESTAMP}.pkl`

**형식**: Python pickle 파일

**내용**:
```python
{
    'model': <XGBoost/LightGBM/CatBoost 모델 객체>,
    'model_type': 'XGBoost',
    'protein_name': 'FTO',
    'metrics': {...},
    'feature_count': 178,
    'fingerprint_type': 'MACCS',
    'dataset_ratio': '20x',
    'ignore3D': True
}
```

**생성 방법**: API Step 3 "모델 학습" 실행 시 자동 생성

**✅ 현재 상태**: 6개 FTO 모델 존재

---

## 📂 디렉토리 구조 (자동 생성)

프로젝트 실행 시 다음 디렉토리가 자동 생성됩니다:

```
webapp/
├── descriptor_selection.csv          # ✅ 필수! 수동으로 추가
│
├── saved_data/                        # 수집된 원본 데이터
│   ├── IC50/                          # ChEMBL/BindingDB 데이터
│   │   └── FTO_summary.xlsx           # Step 1에서 자동 생성
│   ├── BindingDB/                     # BindingDB 전용
│   └── FooDB/                         # FooDB 업로드 및 전처리
│       ├── foodb_compounds.csv        # 업로드한 원본
│       ├── 3d_conformers/             # 3D 구조 캐시 (SDF)
│       │   ├── FTO_5x.sdf
│       │   ├── FTO_10x.sdf
│       │   └── FTO_20x.sdf
│       └── preprocessed/              # 전처리 결과
│           ├── FTO_MACCS_20x_ignore3D_True.csv
│           └── FTO_MACCS_20x_ignore3D_False.csv
│
├── raw/                               # 변환된 특성 데이터
│   ├── Dataset/                       # 통합 데이터셋
│   │   └── FTO.csv                    # Step 2에서 자동 생성
│   ├── TransferSet/                   # (구버전, 호환용)
│   └── descriptors/                   # Descriptor 전용
│
├── models_trained/                    # 학습된 모델
│   └── FTO_XGBoost_20260306_053810.pkl  # Step 3에서 자동 생성
│
├── shap_outputs/                      # SHAP 분석 결과
│   └── shap_summary_20260306_054629.png  # Step 5에서 자동 생성
│
└── food_predictions/                  # FooDB 예측 결과
    └── FTO_XGBoost_..._foodb_predictions.csv
```

---

## 🚀 시작 방법

### Step 1: 프로젝트 다운로드 및 압축 해제

```bash
# 최신 백업 다운로드
wget https://www.genspark.ai/api/files/s/AMOsEmIW -O pci-api.tar.gz

# 압축 해제
tar -xzf pci-api.tar.gz
cd webapp
```

### Step 2: 필수 파일 확인

```bash
# descriptor_selection.csv 확인
ls -lh descriptor_selection.csv

# 없으면 에러! 반드시 추가해야 함
```

**✅ 현재**: 이미 포함되어 있음!

### Step 3: 의존성 설치

```bash
# Python 가상환경 (권장)
python3 -m venv venv
source venv/bin/activate  # Windows: venv\Scripts\activate

# 패키지 설치
pip install -r requirements.txt
```

### Step 4: 서버 실행

```bash
# 개발 모드
uvicorn app.main:app --host 0.0.0.0 --port 3000 --reload

# 또는 PM2 (프로덕션)
pm2 start ecosystem.config.cjs
```

### Step 5: FooDB 데이터 준비 (선택)

```bash
# FooDB 다운로드 (https://foodb.ca/downloads)
# Compounds.csv 다운로드 후 업로드

curl -X POST "http://localhost:3000/api/foodb/upload" \
  -F "file=@Compounds.csv"
```

---

## 🎯 일반적인 워크플로우

### 시나리오 A: 새로운 단백질 학습 및 예측

```bash
# 1. 단백질 데이터 수집 (ChEMBL)
POST /api/data/collect
→ saved_data/IC50/{PROTEIN}_summary.xlsx

# 2. 특성 변환
POST /api/features/transform
→ raw/Dataset/{PROTEIN}.csv

# 3. 모델 학습
POST /api/models/train
→ models_trained/{PROTEIN}_{MODEL}_{TIME}.pkl

# 4. 예측
POST /api/models/predict

# 5. SHAP 분석
POST /api/shap/analyze
→ shap_outputs/shap_summary_{TIME}.png
```

**필요한 파일**: `descriptor_selection.csv` ✅

---

### 시나리오 B: FooDB 식품 화합물 예측

```bash
# 1. FooDB CSV 업로드
POST /api/foodb/upload
→ saved_data/FooDB/foodb_compounds.csv

# 2. FooDB 전처리 (3D + Mordred)
POST /api/foodb/preprocess
→ saved_data/FooDB/preprocessed/{PROTEIN}_{FP}_{RATIO}_ignore3D_{BOOL}.csv
→ saved_data/FooDB/3d_conformers/{PROTEIN}_{RATIO}.sdf (캐시)

# 3. 모델로 예측
POST /api/foodb/predict
→ food_predictions/{MODEL}_foodb_predictions.csv
```

**필요한 파일**: 
- `descriptor_selection.csv` ✅
- FooDB CSV (사용자가 업로드)
- 학습된 모델 (Step 3에서 생성)

---

### 시나리오 C: 간단한 SMILES 예측

```bash
# 학습된 모델로 바로 예측
POST /api/foodb/predict-smiles
{
  "smiles_list": ["CCO", "CC(=O)O"],
  "model_id": "FTO_XGBoost_20260306_053810"
}
```

**필요한 파일**: 
- 학습된 모델 ✅ (이미 6개 존재)
- `descriptor_selection.csv` (내부적으로 사용)

---

## ⚠️ 자주 묻는 질문

### Q1: `descriptor_selection.csv`가 없으면?

**증상**: 
```
HTTPException: descriptor_selection.csv not found
```

**해결**: 
1. 파일을 프로젝트 루트에 추가
2. 형식 확인 (컬럼명 규칙 준수)
3. 서버 재시작

---

### Q2: FooDB CSV 형식이 다르면?

**증상**:
```
HTTPException: CSV file must contain a SMILES column
```

**해결**:
1. `canonical_smiles` 또는 `smiles` 컬럼 추가
2. SMILES 형식 검증 (RDKit으로 파싱 가능한지)

---

### Q3: 3D conformer 생성이 너무 느리면?

**해결**:
1. `ignore3D=True` 사용 (2D만 계산)
2. 화합물 수 줄이기 (필터링)
3. SDF 캐시 재사용

---

## 📝 요약

### ✅ 반드시 필요
- `descriptor_selection.csv` (프로젝트 루트)

### 🟢 API로 자동 생성
- `saved_data/IC50/*.xlsx` (Step 1)
- `raw/Dataset/*.csv` (Step 2)
- `models_trained/*.pkl` (Step 3)
- `shap_outputs/*.png` (Step 5)

### 🔵 사용자가 업로드
- `saved_data/FooDB/*.csv` (FooDB 원본)

### 🟡 자동 캐시
- `saved_data/FooDB/3d_conformers/*.sdf` (3D 구조)
- `saved_data/FooDB/preprocessed/*.csv` (전처리 결과)

---

**프로젝트 백업에 포함된 파일**:
- ✅ `descriptor_selection.csv`
- ✅ `saved_data/IC50/FTO_summary.xlsx`
- ✅ `raw/Dataset/FTO.csv`
- ✅ `models_trained/FTO_*.pkl` (6개)

**추가로 필요한 파일**:
- 🔵 FooDB CSV (https://foodb.ca/downloads에서 다운로드)

모든 준비 완료! 🎉
