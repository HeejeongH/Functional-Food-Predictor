# 🚀 API 사용 가이드

## 실행 방법

### 📌 방법 1: 웹 브라우저로 API 문서 사용 (가장 쉬움!)

1. **브라우저에서 API 문서 열기**
   ```
   https://3000-itvs5ceiquncckj2qhd7g-2b54fc91.sandbox.novita.ai/docs
   ```

2. **각 API 테스트하기**
   - 원하는 엔드포인트 클릭
   - "Try it out" 버튼 클릭
   - 파라미터 입력
   - "Execute" 버튼 클릭
   - 결과 확인!

---

### 📌 방법 2: curl 명령어로 테스트 (터미널)

#### 1️⃣ 서버 상태 확인
```bash
curl https://3000-itvs5ceiquncckj2qhd7g-2b54fc91.sandbox.novita.ai/health
```

#### 2️⃣ SMILES를 Fingerprint로 변환
```bash
curl -X POST "https://3000-itvs5ceiquncckj2qhd7g-2b54fc91.sandbox.novita.ai/api/features/fingerprint?fp_type=ECFP4" \
  -H "Content-Type: application/json" \
  -d '["CCO", "CC(=O)O", "c1ccccc1"]'
```

#### 3️⃣ 데이터 상태 확인
```bash
curl https://3000-itvs5ceiquncckj2qhd7g-2b54fc91.sandbox.novita.ai/api/data/status
```

#### 4️⃣ 모델 리스트 확인
```bash
curl https://3000-itvs5ceiquncckj2qhd7g-2b54fc91.sandbox.novita.ai/api/models/list
```

---

### 📌 방법 3: Python 스크립트로 테스트

샌드박스에서 제공된 `test_example.py`를 실행:

```bash
cd /home/user/webapp
python3 test_example.py
```

또는 직접 Python 코드 작성:

```python
import requests

BASE_URL = "https://3000-itvs5ceiquncckj2qhd7g-2b54fc91.sandbox.novita.ai"

# Health Check
response = requests.get(f"{BASE_URL}/health")
print(response.json())

# Fingerprint 생성
response = requests.post(
    f"{BASE_URL}/api/features/fingerprint?fp_type=ECFP4",
    json=["CCO", "CC(=O)O"]
)
print(response.json())
```

---

## 📋 전체 워크플로우 예제

### 실제 연구 시나리오: PDE4 타겟 화합물 예측

```bash
# 1단계: 데이터 수집 (ChEMBL/BindingDB)
curl -X POST "https://3000-itvs5ceiquncckj2qhd7g-2b54fc91.sandbox.novita.ai/api/data/collect" \
  -H "Content-Type: application/json" \
  -d '{
    "target_list": ["PDE4"],
    "standard_type": "IC50"
  }'

# 2단계: 특성 변환 (Fingerprint)
curl -X POST "https://3000-itvs5ceiquncckj2qhd7g-2b54fc91.sandbox.novita.ai/api/features/transform" \
  -H "Content-Type: application/json" \
  -d '{
    "protein_name": "PDE4",
    "fingerprint_type": "ECFP4",
    "dataset_type": "transfer",
    "pos_threshold": 10000,
    "neg_threshold": 20000
  }'

# 3단계: 모델 학습
curl -X POST "https://3000-itvs5ceiquncckj2qhd7g-2b54fc91.sandbox.novita.ai/api/models/train" \
  -H "Content-Type: application/json" \
  -d '{
    "protein_name": "PDE4",
    "model_type": "XGBoost",
    "feature_type": "fingerprint",
    "test_size": 0.2,
    "random_state": 42
  }'

# 4단계: 새로운 화합물 예측
# (model_id는 3단계에서 받은 값 사용)
curl -X POST "https://3000-itvs5ceiquncckj2qhd7g-2b54fc91.sandbox.novita.ai/api/models/predict" \
  -H "Content-Type: application/json" \
  -d '{
    "smiles_list": ["CCO", "CC(=O)Oc1ccccc1C(=O)O"],
    "model_id": "PDE4_XGBoost_20260226_120000",
    "feature_type": "fingerprint"
  }'

# 5단계: SHAP 분석
curl -X POST "https://3000-itvs5ceiquncckj2qhd7g-2b54fc91.sandbox.novita.ai/api/shap/analyze" \
  -H "Content-Type: application/json" \
  -d '{
    "model_id": "PDE4_XGBoost_20260226_120000",
    "feature_type": "fingerprint",
    "top_n": 20
  }'
```

---

## 🧪 간단한 테스트 예제

### 지금 바로 테스트할 수 있는 API

#### ✅ Health Check
```bash
curl https://3000-itvs5ceiquncckj2qhd7g-2b54fc91.sandbox.novita.ai/health
```

**결과:**
```json
{
  "status": "healthy",
  "services": {
    "data_collection": "available",
    "feature_transform": "available",
    "model_training": "available",
    "shap_analysis": "available"
  },
  "data": {
    "collected_datasets": 0,
    "trained_models": 0
  }
}
```

#### ✅ SMILES → Fingerprint 변환
```bash
curl -X POST "https://3000-itvs5ceiquncckj2qhd7g-2b54fc91.sandbox.novita.ai/api/features/fingerprint?fp_type=ECFP4" \
  -H "Content-Type: application/json" \
  -d '["CCO"]'
```

---

## 📚 추가 리소스

- **Swagger UI**: https://3000-itvs5ceiquncckj2qhd7g-2b54fc91.sandbox.novita.ai/docs
- **ReDoc**: https://3000-itvs5ceiquncckj2qhd7g-2b54fc91.sandbox.novita.ai/redoc
- **README**: `/home/user/webapp/README.md`

---

## ⚠️ 주의사항

1. **데이터 수집**: ChEMBL API는 외부 네트워크 연결이 필요합니다.
2. **모델 학습**: 데이터가 충분해야 합니다 (최소 100개 이상 권장).
3. **예측**: 학습된 모델이 있어야 예측 가능합니다.
4. **SHAP 분석**: 계산 시간이 걸릴 수 있습니다.

---

## 💡 팁

- **API 문서**: 모든 파라미터와 응답 형식을 확인할 수 있습니다.
- **에러 메시지**: 자세한 에러 정보가 응답에 포함됩니다.
- **로그 확인**: `pm2 logs pci-api --nostream`로 서버 로그 확인.
