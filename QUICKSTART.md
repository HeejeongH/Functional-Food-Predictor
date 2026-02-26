# PCI Prediction Platform - 빠른 시작 가이드

## 🚀 빠른 시작

### 1. 환경 설정

```bash
# 1. 저장소 클론
git clone <repository-url>
cd webapp

# 2. Python 가상환경 생성 (권장)
python -m venv venv
source venv/bin/activate  # Linux/Mac
# 또는
venv\Scripts\activate  # Windows

# 3. 의존성 설치
pip install -r requirements.txt

# 4. RDKit 설치 (conda 권장)
conda install -c conda-forge rdkit
```

### 2. 웹 애플리케이션 실행

```bash
# Linux/Mac
./run.sh

# Windows
streamlit run app/streamlit_app.py
```

브라우저에서 `http://localhost:8501` 접속

### 3. 데모 스크립트 실행

```bash
python demo.py
```

## 📝 사용 예제

### 예제 1: FTO 단백질 PCI 예측

#### Step 1: 데이터 수집
```
타겟 입력: FTO, Alpha-ketoglutarate-dependent dioxygenase
활성 데이터 타입: IC50
```

#### Step 2: 특성 변환
```
Fingerprint: ECFP4
크기: 1024
Molecular Descriptors: 포함
활성 임계값: 10000 nM
비활성 임계값: 20000 nM
```

#### Step 3: 모델 학습
```
모델: TabPFN
피처 선택: PCA
선택 피처 수: 200
테스트 비율: 0.2
```

#### Step 4: SHAP 분석
```
표시 피처 수: 20
SHAP 샘플: 100
```

#### Step 5: FooDB 예측
```
FooDB 파일: foodb_compounds.csv
예측 임계값: 0.5
배치 크기: 500
```

### 예제 2: Python 코드로 직접 사용

```python
from modules.data_collector import collect_pci_data
from modules.feature_extractor import prepare_training_data
from modules.model_trainer import ModelTrainer
from modules.shap_analyzer import analyze_model_with_shap

# 1. 데이터 수집
chembl_df = collect_pci_data(
    chembl_search_list=['FTO', 'Alpha-ketoglutarate-dependent dioxygenase'],
    standard_type='IC50'
)

# 2. 특성 변환
prepared_df = prepare_training_data(
    chembl_df=chembl_df,
    protein_name='FTO',
    fp_type='ECFP4',
    fp_size=1024,
    include_descriptors=True,
    pos_threshold=10000,
    neg_threshold=20000
)

# 3. 모델 학습
trainer = ModelTrainer(
    model_type='TabPFN',
    n_features=200,
    feature_selection_method='pca'
)

metrics = trainer.train(prepared_df, test_size=0.2)
print(f"F1 Score: {metrics['f1']:.4f}")
print(f"AUC: {metrics['auc']:.4f}")

# 4. 모델 저장
trainer.save_model('models/fto_model.joblib')

# 5. SHAP 분석
from sklearn.model_selection import train_test_split
feature_columns = [c for c in prepared_df.columns if c.startswith(('FP_', 'DESC_'))]
X = prepared_df[feature_columns].values
y = prepared_df['Y'].values

X_train, X_test, y_train, y_test = train_test_split(
    X, y, test_size=0.2, random_state=42
)

shap_results = analyze_model_with_shap(
    model=trainer.model,
    X_train=X_train,
    X_test=X_test,
    feature_names=feature_columns,
    save_dir='results/shap'
)

print("Top 10 important features:")
print(shap_results['importance_df'].head(10))
```

## 🔧 고급 설정

### 커스텀 Fingerprint 크기

```python
from modules.feature_extractor import MolecularFeatureExtractor

extractor = MolecularFeatureExtractor(fp_type='ECFP4', fp_size=2048)
fp = extractor.smiles_to_fingerprint('CC(=O)OC1=CC=CC=C1C(=O)O')
```

### 배치 예측 최적화

```python
from modules.foodb_predictor import batch_predict_foodb

# 대용량 데이터는 배치 크기 조정
predictions = batch_predict_foodb(
    foodb_df=foodb_features,
    trainer=trainer,
    feature_columns=feature_columns,
    batch_size=1000,  # 메모리에 따라 조정
    threshold=0.5
)
```

### 다중 모델 앙상블

```python
from sklearn.ensemble import VotingClassifier

# 여러 모델 학습
models = []
for model_type in ['RandomForest', 'GradientBoosting', 'SVM']:
    trainer = ModelTrainer(model_type=model_type)
    trainer.train(prepared_df)
    models.append((model_type, trainer.model))

# 앙상블
ensemble = VotingClassifier(estimators=models, voting='soft')
```

## 📊 결과 해석

### SHAP 값 이해하기

- **양수 SHAP 값**: 해당 특성이 활성 예측에 기여
- **음수 SHAP 값**: 해당 특성이 비활성 예측에 기여
- **절대값 크기**: 특성의 중요도

### 피처 중요도 활용

```python
# 상위 중요 피처 추출
top_features = shap_results['top_features'][:20]

# 중요한 구조적 특성 확인
structural_features = [f for f in top_features if f.startswith('DESC_')]
fingerprint_features = [f for f in top_features if f.startswith('FP_')]

print(f"중요한 구조적 특성: {len(structural_features)}개")
print(f"중요한 Fingerprint 비트: {len(fingerprint_features)}개")
```

## 🐛 문제 해결

### ChEMBL 연결 오류
```bash
pip install --upgrade chembl-webresource-client
```

### RDKit 설치 문제
```bash
# conda를 사용하세요
conda install -c conda-forge rdkit
```

### TabPFN 메모리 오류
- 피처 수를 줄이세요 (예: 200개)
- 샘플 수를 제한하세요
- 다른 모델을 사용하세요 (RandomForest 등)

### SHAP 계산 시간 오류
- `max_background_samples`를 줄이세요 (50-100)
- `nsamples`를 줄이세요 (100)
- 테스트 데이터 샘플 수를 제한하세요

## 📚 추가 자료

### API 문서
- ChEMBL API: https://chembl.gitbook.io/chembl-interface-documentation/
- RDKit: https://www.rdkit.org/docs/
- TabPFN: https://github.com/automl/TabPFN

### 데이터베이스
- ChEMBL: https://www.ebi.ac.uk/chembl/
- BindingDB: https://www.bindingdb.org/
- FooDB: https://foodb.ca/

### 논문
- TabPFN: "TabPFN: A Transformer That Solves Small Tabular Classification Problems in a Second" (2022)
- SHAP: "A Unified Approach to Interpreting Model Predictions" (2017)
- Morgan Fingerprints: "Extended-Connectivity Fingerprints" (2010)

## 💡 팁

1. **데이터 품질이 중요합니다**: IC50 값이 명확한 데이터를 사용하세요
2. **임계값 조정**: 활성/비활성 임계값은 타겟에 따라 조정하세요
3. **피처 선택**: PCA는 빠르고, RFE는 정확하지만 느립니다
4. **모델 선택**: 소규모 데이터는 TabPFN, 대규모는 RandomForest
5. **교차 검증**: 중요한 프로젝트는 별도 검증 세트를 사용하세요

## 🤝 기여하기

버그 리포트, 기능 제안, 코드 기여를 환영합니다!

1. Fork the repository
2. Create your feature branch
3. Commit your changes
4. Push to the branch
5. Create a Pull Request

## 📄 라이센스

MIT License - 자유롭게 사용하고 수정할 수 있습니다.
