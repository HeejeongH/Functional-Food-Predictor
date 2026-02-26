# PCI Prediction Platform - 사용 안내서

## 🎯 프로젝트 개요

이 프로젝트는 **단백질-화합물 상호작용(Protein-Compound Interaction, PCI)** 예측을 위한 통합 웹 플랫폼입니다.

### 주요 기능

1. **데이터 수집**: ChEMBL/BindingDB에서 IC50 기준 PCI 데이터 자동 수집
2. **특성 변환**: 6가지 Fingerprint 및 Molecular Descriptor 계산
3. **모델 학습**: TabPFN, RandomForest 등 5가지 ML 모델 지원
4. **SHAP 분석**: 화합물의 어떤 특성이 중요한지 해석
5. **FooDB 예측**: 실제 식품 화합물 데이터로 활성 예측

## 📁 프로젝트 구조

```
webapp/
├── app/
│   └── streamlit_app.py          # Streamlit 웹 애플리케이션
├── modules/
│   ├── data_collector.py         # ChEMBL 데이터 수집
│   ├── feature_extractor.py      # Fingerprint & Descriptor 계산
│   ├── model_trainer.py          # ML 모델 학습 (TabPFN, RF, etc.)
│   ├── shap_analyzer.py          # SHAP 분석
│   └── foodb_predictor.py        # FooDB 예측
├── DB/
│   └── utils.py                  # 기존 유틸리티 (호환성)
├── data/
│   ├── temp/                     # 임시 데이터
│   └── results/                  # 결과 저장
├── models/                       # 학습된 모델
├── logs/                         # 로그
├── requirements.txt              # Python 패키지
├── run.sh                        # 실행 스크립트
├── demo.py                       # 데모 스크립트
├── README.md                     # 프로젝트 설명
└── QUICKSTART.md                 # 빠른 시작 가이드
```

## 🚀 실행 방법

### 방법 1: 웹 애플리케이션 (추천)

```bash
# 실행
./run.sh

# 또는
streamlit run app/streamlit_app.py

# 브라우저에서 http://localhost:8501 접속
```

### 방법 2: Python 코드로 직접 사용

```python
# demo.py 참고
python demo.py
```

## 📖 사용 가이드

### 1단계: 데이터 수집 🔍

웹앱에서:
1. 타겟 유전자 입력 (예: `FTO, Alpha-ketoglutarate-dependent dioxygenase`)
2. IC50/EC50 선택
3. "데이터 수집 시작" 클릭

결과: ChEMBL에서 수집된 화합물 데이터

### 2단계: 특성 변환 ⚙️

1. **Fingerprint 선택**:
   - ECFP4 (추천): Morgan Fingerprint, radius=2
   - ECFP6: Morgan Fingerprint, radius=3
   - MACCS: 167개 고정 키
   - AtomPair: 원자 쌍 기반
   - TopologicalTorsion: 위상학적 특성
   - RDKit: RDKit 기본 FP

2. **설정**:
   - Fingerprint 크기: 1024 (추천)
   - Molecular Descriptors 포함: 체크
   - 활성 임계값: 10000 nM
   - 비활성 임계값: 20000 nM

3. "특성 변환 시작" 클릭

결과: Fingerprint + Descriptor가 포함된 학습 데이터

### 3단계: 모델 학습 🤖

1. **모델 선택**:
   - **TabPFN** (추천): 최신 트랜스포머 기반, 빠르고 정확
   - RandomForest: 전통적 앙상블, 안정적
   - GradientBoosting: 높은 성능
   - SVM: 소규모 데이터에 적합
   - LogisticRegression: 빠른 기본 모델

2. **피처 선택**:
   - PCA (추천): 빠르고 효과적
   - RFE: 정확하지만 느림
   - Mutual Information: 중간 수준
   - Univariate: 빠른 기본 방법
   
3. **선택 피처 수**: 200 (추천)

4. "모델 학습 시작" 클릭

결과: F1 Score, AUC 등 성능 메트릭

### 4단계: SHAP 분석 📊

1. 표시할 피처 수: 20
2. SHAP 샘플 수: 100
3. "SHAP 분석 시작" 클릭

결과:
- Summary Plot: 각 피처의 SHAP 값 분포
- Bar Plot: 피처 중요도 순위
- Feature Importance CSV: 다운로드 가능

**해석 방법**:
- 높은 SHAP 값 = 해당 특성이 활성 예측에 중요
- FP_XXX: 특정 화학 구조 패턴
- DESC_XXX: 물리화학적 특성 (MW, LogP 등)

### 5단계: FooDB 예측 🍎

1. FooDB CSV 파일 업로드
2. 예측 임계값: 0.5
3. 배치 크기: 500
4. "예측 시작" 클릭

결과:
- 활성 예측된 식품 화합물 리스트
- 확률 기반 정렬
- CSV 다운로드 가능

## 💡 사용 팁

### 데이터 수집
- 타겟 이름은 여러 개 입력 가능 (쉼표로 구분)
- IC50 데이터가 가장 풍부함
- 최소 100개 이상의 화합물 필요

### 특성 변환
- ECFP4가 가장 범용적
- Molecular Descriptors는 해석에 유용
- 임계값은 문헌 참고하여 조정

### 모델 학습
- TabPFN은 1000개 미만 샘플에 최적
- 피처 수가 많으면 PCA 추천
- 클래스 불균형이 심하면 SMOTE 고려

### SHAP 분석
- 계산 시간이 오래 걸리면 샘플 수 줄이기
- DESC로 시작하는 피처가 해석하기 쉬움
- 상위 20개 피처로 충분

### FooDB 예측
- 대용량 데이터는 배치 크기 조정
- 임계값 0.5는 일반적 기준
- 확률 0.7 이상이 높은 신뢰도

## 🔧 문제 해결

### Q1: ChEMBL 데이터가 수집되지 않음
```bash
pip install --upgrade chembl-webresource-client
```

### Q2: RDKit 설치 오류
```bash
conda install -c conda-forge rdkit
```

### Q3: TabPFN 메모리 부족
- 피처 수를 200 이하로 줄이기
- 다른 모델 사용 (RandomForest)

### Q4: SHAP 계산이 너무 느림
- 샘플 수를 50으로 줄이기
- 배경 샘플을 50으로 제한

### Q5: Streamlit이 실행되지 않음
```bash
pip install --upgrade streamlit
streamlit run app/streamlit_app.py --server.port 8502
```

## 📊 예상 결과

### 성능 지표 (FTO 예제)
- **데이터 수**: 200-300개 화합물
- **F1 Score**: 0.95-0.98
- **AUC**: 0.95-0.99
- **학습 시간**: 1-3분

### 중요 피처 (일반적)
1. LogP (소수성)
2. TPSA (극성 표면적)
3. MolWt (분자량)
4. NumHDonors (수소 공여체)
5. ECFP4 특정 비트 (구조 패턴)

## 🎓 추가 학습

### 권장 논문
1. TabPFN: "TabPFN: A Transformer That Solves Small Tabular Classification Problems in a Second"
2. SHAP: "A Unified Approach to Interpreting Model Predictions"
3. Morgan FP: "Extended-Connectivity Fingerprints"

### 권장 자료
- [ChEMBL Tutorial](https://chembl.gitbook.io/)
- [RDKit Cookbook](https://www.rdkit.org/docs/Cookbook.html)
- [SHAP Documentation](https://shap.readthedocs.io/)

## 📞 지원

- 버그 리포트: GitHub Issues
- 기능 제안: GitHub Discussions
- 이메일: pci-research@example.com

## 📝 라이센스

MIT License - 자유롭게 사용 및 수정 가능

---

**마지막 업데이트**: 2024-02-26
**버전**: 1.0.0
