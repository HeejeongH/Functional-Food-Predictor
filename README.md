# PCI Prediction Platform

단백질-화합물 상호작용(Protein-Compound Interaction) 예측을 위한 통합 웹 플랫폼

## 주요 기능

### 1. 데이터 수집 🔍
- **ChEMBL API**: 타겟 유전자/단백질에 대한 IC50/EC50/Ki 데이터 자동 수집
- **BindingDB**: 추가 PCI 데이터 통합 (수동 업로드)
- **데이터 정제**: Canonical SMILES 변환 및 중복 제거

### 2. 분자 특성 변환 ⚙️
- **Fingerprints**:
  - ECFP4 (Morgan, radius=2)
  - ECFP6 (Morgan, radius=3)
  - MACCS Keys
  - Atom Pair Fingerprint
  - Topological Torsion
  - RDKit Fingerprint
  
- **Molecular Descriptors**:
  - 물리화학적 특성 (MW, LogP, TPSA 등)
  - 구조적 특성 (회전 가능 결합, 방향족 고리 등)
  - 분자 복잡도 지표

### 3. 머신러닝 모델 학습 🤖
- **지원 모델**:
  - TabPFN (v2.5) - 최신 트랜스포머 기반 모델
  - Random Forest
  - Gradient Boosting
  - SVM
  - Logistic Regression

- **피처 선택**:
  - Mutual Information
  - Recursive Feature Elimination (RFE)
  - Principal Component Analysis (PCA)
  - Univariate Feature Selection

- **성능 메트릭**:
  - Accuracy, Precision, Recall
  - F1 Score, AUC-ROC
  - Confusion Matrix

### 4. SHAP 분석 📊
- **모델 설명가능성**: SHAP (SHapley Additive exPlanations)
- **피처 중요도**: 상위 중요 화학적 특성 추출
- **시각화**:
  - Summary Plot
  - Bar Plot
  - Feature Importance Table

### 5. FooDB 예측 🍎
- **실제 식품 화합물**: FooDB 데이터베이스 (70,000+ 화합물)
- **배치 예측**: 대용량 데이터 효율적 처리
- **결과 필터링**: 활성 확률 기반 화합물 선별
- **CSV 출력**: 예측 결과 다운로드

## 설치 방법

### 1. 필수 요구사항
- Python 3.8+
- pip

### 2. 의존성 설치
```bash
pip install -r requirements.txt
```

### 3. RDKit 설치 (별도)
RDKit은 conda를 통한 설치를 권장합니다:
```bash
conda install -c conda-forge rdkit
```

또는 pip으로:
```bash
pip install rdkit-pypi
```

## 실행 방법

### Linux/Mac
```bash
chmod +x run.sh
./run.sh
```

### Windows
```bash
streamlit run app/streamlit_app.py
```

### Python으로 직접 실행
```python
import streamlit.cli as stcli
import sys

sys.argv = ["streamlit", "run", "app/streamlit_app.py"]
sys.exit(stcli.main())
```

## 프로젝트 구조

```
webapp/
├── app/
│   └── streamlit_app.py      # Streamlit 웹 애플리케이션
├── modules/
│   ├── __init__.py
│   ├── data_collector.py     # ChEMBL/BindingDB 데이터 수집
│   ├── feature_extractor.py  # 분자 특성 변환
│   ├── model_trainer.py      # ML 모델 학습
│   ├── shap_analyzer.py      # SHAP 분석
│   └── foodb_predictor.py    # FooDB 예측
├── DB/
│   └── utils.py              # 기존 유틸리티 (호환성)
├── data/
│   ├── temp/                 # 임시 데이터
│   └── results/              # 결과 저장
├── models/                   # 학습된 모델 저장
├── logs/                     # 로그 파일
├── requirements.txt          # Python 의존성
├── run.sh                    # 실행 스크립트 (Linux/Mac)
└── README.md                 # 이 파일
```

## 사용 가이드

### 1. 데이터 수집
1. 타겟 유전자/단백질 이름 입력 (예: FTO, PDE4B)
2. 활성 데이터 타입 선택 (IC50, EC50, Ki 등)
3. "데이터 수집 시작" 클릭
4. 수집된 데이터 확인 및 다운로드

### 2. 특성 변환
1. Fingerprint 타입 선택
2. Molecular Descriptor 포함 여부 선택
3. 활성/비활성 임계값 설정
4. "특성 변환 시작" 클릭
5. 준비된 학습 데이터 확인

### 3. 모델 학습
1. 모델 타입 선택 (TabPFN 권장)
2. 피처 선택 방법 및 개수 설정
3. 테스트 데이터 비율 설정
4. "모델 학습 시작" 클릭
5. 성능 메트릭 및 Confusion Matrix 확인

### 4. SHAP 분석
1. 표시할 피처 수 설정
2. SHAP 계산 샘플 수 설정
3. "SHAP 분석 시작" 클릭
4. 피처 중요도 시각화 확인
5. 상위 중요 특성 CSV 다운로드

### 5. FooDB 예측
1. FooDB CSV 파일 업로드
2. 예측 임계값 설정
3. "예측 시작" 클릭
4. 활성 화합물 결과 확인
5. 전체 예측 결과 CSV 다운로드

## 기술 스택

- **데이터 수집**: ChEMBL Web Resource Client
- **화학 정보**: RDKit
- **머신러닝**: Scikit-learn, TabPFN
- **설명가능 AI**: SHAP
- **시각화**: Plotly, Matplotlib
- **웹 프레임워크**: Streamlit
- **데이터 처리**: Pandas, NumPy

## 참고 자료

- [ChEMBL Database](https://www.ebi.ac.uk/chembl/)
- [BindingDB](https://www.bindingdb.org/)
- [FooDB](https://foodb.ca/)
- [RDKit Documentation](https://www.rdkit.org/docs/)
- [TabPFN Paper](https://arxiv.org/abs/2207.01848)
- [SHAP Documentation](https://shap.readthedocs.io/)

## 라이센스

MIT License

## 문의

PCI Research Team
- Email: pci-research@example.com
- GitHub: https://github.com/pci-research/pci-platform

## 버전

**v1.0.0** (2024-02-26)
- 초기 릴리스
- ChEMBL 데이터 수집
- 6가지 Fingerprint 타입 지원
- 5가지 ML 모델 지원
- SHAP 분석 통합
- FooDB 예측 기능
