# 🎉 API 실행 방법 - 간단 가이드

## 가장 쉬운 방법: 웹 브라우저 사용 🌐

### 1단계: API 문서 열기
브라우저에서 이 URL을 여세요:
```
https://3000-itvs5ceiquncckj2qhd7g-2b54fc91.sandbox.novita.ai/docs
```

### 2단계: API 테스트하기
1. **원하는 API 선택** (예: `/api/features/fingerprint`)
2. **"Try it out" 버튼 클릭**
3. **데이터 입력**:
   ```json
   {
     "smiles_list": ["CCO", "CC(=O)O"],
     "fp_type": "ECFP4"
   }
   ```
4. **"Execute" 버튼 클릭**
5. **결과 확인!** ✅

---

## 명령어로 테스트하기 💻

### 간단한 테스트 (지금 바로 해보세요!)

```bash
# 1. 서버 상태 확인
curl https://3000-itvs5ceiquncckj2qhd7g-2b54fc91.sandbox.novita.ai/health

# 2. SMILES → Fingerprint 변환
curl -X POST "https://3000-itvs5ceiquncckj2qhd7g-2b54fc91.sandbox.novita.ai/api/features/fingerprint" \
  -H "Content-Type: application/json" \
  -d '{"smiles_list": ["CCO", "CC(=O)O"], "fp_type": "ECFP4"}'
```

---

## Python으로 테스트하기 🐍

### 방법 1: 제공된 스크립트 실행

```bash
cd /home/user/webapp
python3 test_example.py
```

### 방법 2: 직접 코드 작성

```python
import requests

# API URL
url = "https://3000-itvs5ceiquncckj2qhd7g-2b54fc91.sandbox.novita.ai"

# 1. Health Check
response = requests.get(f"{url}/health")
print(response.json())

# 2. SMILES를 Fingerprint로 변환
data = {
    "smiles_list": ["CCO", "CC(=O)O"],
    "fp_type": "ECFP4"
}
response = requests.post(f"{url}/api/features/fingerprint", json=data)
print(response.json())
```

---

## 🔍 어떤 API를 사용할 수 있나요?

| API | 설명 | 예제 |
|-----|------|------|
| `GET /health` | 서버 상태 확인 | 지금 바로 테스트 가능! |
| `GET /api/data/status` | 데이터 상태 | 지금 바로 테스트 가능! |
| `POST /api/features/fingerprint` | SMILES → Fingerprint | ✅ 지금 바로 테스트 가능! |
| `POST /api/data/collect` | ChEMBL 데이터 수집 | 실제 데이터 필요 |
| `POST /api/models/train` | 모델 학습 | 데이터 수집 후 가능 |
| `POST /api/models/predict` | 화합물 예측 | 모델 학습 후 가능 |

---

## 📚 더 자세한 정보

- **전체 가이드**: `/home/user/webapp/USAGE_GUIDE.md`
- **README**: `/home/user/webapp/README.md`
- **API 문서**: https://3000-itvs5ceiquncckj2qhd7g-2b54fc91.sandbox.novita.ai/docs

---

## ✨ 바로 시작하기 (추천!)

가장 쉬운 방법은 **웹 브라우저로 API 문서**를 여는 것입니다:

👉 https://3000-itvs5ceiquncckj2qhd7g-2b54fc91.sandbox.novita.ai/docs

거기서 모든 API를 클릭 한 번으로 테스트할 수 있습니다!
