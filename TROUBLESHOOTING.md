# 문제 해결 가이드

## ❌ "Unknown label type" 오류

### 문제 설명
```
Unknown label type: unknown. Maybe you are trying to fit a classifier, 
which expects discrete classes on a regression target with continuous values.
```

### 원인
데이터의 Y 컬럼(레이블)에 0, 1이 아닌 다른 값이 포함되어 있을 때 발생합니다.

### ✅ 해결 방법 (이미 수정됨!)

**v1.0.1 업데이트로 자동 수정됨**:
- Y 컬럼이 정수형(0 또는 1)으로만 저장되도록 수정
- 중간값('-')이 자동으로 필터링됨
- 데이터 검증 로직 추가

### 수동 해결 방법 (필요시)

#### 1. IC50 임계값 조정

수집된 데이터가 임계값 범위에 충분히 없을 수 있습니다.

```python
# 웹앱에서 조정:
# 활성 임계값: 10000 → 50000 nM
# 비활성 임계값: 20000 → 10000 nM
```

#### 2. 데이터 확인

데이터 수집 탭에서 IC50 분포를 확인하세요:
- 대부분의 데이터가 10000-20000 범위 밖에 있으면 임계값 조정 필요
- 충분한 활성/비활성 샘플이 필요 (각 최소 10개 이상)

#### 3. 다른 타겟 시도

현재 타겟에서 데이터가 부족하면:
- 더 잘 연구된 타겟 선택 (예: EGFR, VEGFR)
- IC50 데이터가 풍부한 타겟 우선

### 기술적 설명

**수정 전**:
```python
# 문제가 있던 코드
meta_data = pd.DataFrame({
    'Y': valid_data['potency'],  # '-' 문자열이 포함될 수 있음
})
```

**수정 후**:
```python
# 수정된 코드
# 1. '-' 값을 필터링
active_data = df[df['potency'] == 1].copy()
inactive_data = df[df['potency'] == 0].copy()

# 2. 정수형으로 명시적 변환
meta_data = pd.DataFrame({
    'Y': valid_data['potency'].astype(int),  # 0 또는 1만
})

# 3. 검증
if not all(y in [0, 1] for y in result_df['Y'].unique()):
    raise ValueError("Y 컬럼 오류")
```

## 🔍 기타 일반적인 문제들

### 1. ChEMBL 연결 오류

**증상**: "Failed to connect to ChEMBL API"

**해결**:
```bash
# 1. 인터넷 연결 확인
ping www.ebi.ac.uk

# 2. 패키지 업데이트
pip install --upgrade chembl-webresource-client

# 3. SSL 인증서 문제
pip install --upgrade certifi
```

### 2. RDKit 설치 오류

**증상**: "ModuleNotFoundError: No module named 'rdkit'"

**해결**:
```bash
# Conda 사용 (권장)
conda install -c conda-forge rdkit

# Pip 사용 (대안)
pip install rdkit-pypi
```

### 3. TabPFN 메모리 오류

**증상**: "CUDA out of memory" 또는 프로그램 멈춤

**해결**:
```python
# 웹앱에서:
# 1. 피처 선택 방법: PCA 사용
# 2. 선택 피처 수: 200 이하로 감소
# 3. 다른 모델 사용: RandomForest
```

### 4. Streamlit 실행 오류

**증상**: "streamlit: command not found"

**해결**:
```bash
# 1. 설치 확인
pip install streamlit

# 2. PATH 확인
which streamlit

# 3. 직접 실행
python -m streamlit run app/streamlit_app.py
```

### 5. 포트 충돌

**증상**: "Port 8501 is already in use"

**해결**:
```bash
# 1. 다른 포트 사용
streamlit run app/streamlit_app.py --server.port 8502

# 2. 기존 프로세스 종료 (macOS/Linux)
lsof -ti:8501 | xargs kill -9

# 3. 기존 프로세스 종료 (Windows)
netstat -ano | findstr :8501
taskkill /PID <PID번호> /F
```

### 6. 데이터가 수집되지 않음

**증상**: "Found 0 targets" 또는 "ChEMBL data: 0 records"

**해결**:
```python
# 1. 타겟 이름 확인
# 올바른 예: FTO, EGFR, PDE4B
# 잘못된 예: fto (소문자), FTO protein (너무 구체적)

# 2. 여러 이름 시도
# "FTO, Alpha-ketoglutarate-dependent dioxygenase"

# 3. ChEMBL에서 직접 검색
# https://www.ebi.ac.uk/chembl/
```

### 7. SHAP 계산 시간이 너무 오래 걸림

**증상**: SHAP 분석이 몇 분 이상 걸림

**해결**:
```python
# 웹앱에서:
# 1. SHAP 샘플 수: 100 → 50으로 감소
# 2. 표시 피처 수: 20 → 10으로 감소

# 또는 Python 코드:
analyzer.calculate_shap_values(
    X_test, 
    max_background_samples=50,
    nsamples=50
)
```

### 8. FooDB 예측 오류

**증상**: "No valid fingerprints generated from FooDB"

**해결**:
```python
# 1. CSV 파일 형식 확인
# 필수 컬럼: canonical_SMILES 또는 raw_SMILES

# 2. SMILES 유효성 확인
from rdkit import Chem
mol = Chem.MolFromSmiles(smiles)
if mol is None:
    print("Invalid SMILES")

# 3. 샘플 데이터로 테스트
# demo.py 실행
```

## 💡 디버깅 팁

### 1. 로그 확인

```bash
# 로그 파일 확인
tail -f logs/app.log

# 또는 Streamlit 터미널 출력 확인
```

### 2. Python 대화형 모드로 테스트

```python
# 문제가 있는 부분만 테스트
from modules.feature_extractor import MolecularFeatureExtractor

extractor = MolecularFeatureExtractor(fp_type='ECFP4')
fp = extractor.smiles_to_fingerprint('CC(=O)OC1=CC=CC=C1C(=O)O')
print(fp.sum())  # 0이 아닌 값이어야 함
```

### 3. 단계별 실행

웹앱 대신 Python 코드로 각 단계 확인:

```python
# 1. 데이터 수집
from modules.data_collector import collect_pci_data
chembl_df = collect_pci_data(['FTO'])
print(len(chembl_df))

# 2. 특성 변환
from modules.feature_extractor import prepare_training_data
prepared_df = prepare_training_data(chembl_df, 'FTO')
print(prepared_df['Y'].value_counts())

# 3. 모델 학습
from modules.model_trainer import ModelTrainer
trainer = ModelTrainer(model_type='RandomForest')
metrics = trainer.train(prepared_df)
print(metrics)
```

## 🆘 여전히 문제가 있나요?

### GitHub Issues

버그 리포트나 질문:
https://github.com/HeejeongH/Functional-Food-Predictor/issues

다음 정보를 포함해주세요:
1. 오류 메시지 전체
2. 사용한 타겟 이름
3. 운영체제 (macOS/Windows/Linux)
4. Python 버전
5. 실행한 명령어

### 커뮤니티

- Streamlit 포럼: https://discuss.streamlit.io/
- RDKit 포럼: https://github.com/rdkit/rdkit/discussions
- ChEMBL 지원: https://www.ebi.ac.uk/chembl/

## 📋 체크리스트

문제 해결 전 확인:

- [ ] Python 3.8 이상
- [ ] RDKit 설치됨
- [ ] 모든 패키지 최신 버전
- [ ] 인터넷 연결 정상
- [ ] 충분한 메모리 (8GB+)
- [ ] 최신 코드 (git pull)
- [ ] demo.py 실행 성공

---

**마지막 업데이트**: 2024-02-26 (v1.0.1)

**주요 버그 수정**: "Unknown label type" 오류 완전 해결
