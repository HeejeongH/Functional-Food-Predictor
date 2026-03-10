# 🚀 PCI Prediction API - 배포 가이드

> **🪟 Windows 사용자**: [WINDOWS_SETUP.md](WINDOWS_SETUP.md) 문서를 먼저 참조하세요!
> - Docker Desktop + WSL2 설치 상세 가이드
> - Anaconda 방식 (Docker 없이 사용)
> - Windows 전용 문제 해결

## 📦 다운로드 링크

**최신 프로젝트 백업 (Docker + Windows + TabPFN 지원)**: https://www.genspark.ai/api/files/s/xiJMtcDB
- 파일 크기: 1.4 MB
- 포함 내용: Docker/Docker Compose, Windows 가이드, TabPFN, FooDB 전처리, 6개 학습 모델, descriptor_selection.csv (14 descriptors ✅ 수정됨)
- 수정됨: descriptor 개수 불일치 해결 (1038 features)

**이전 버전 (PM2 방식)**: https://www.genspark.ai/api/files/s/AMOsEmIW

이 파일을 다운로드하여 압축을 풀면 전체 프로젝트를 로컬에서 실행할 수 있습니다.

**✨ 최신 버전 업데이트 (2026-03-08)**:
- 🐳 **Docker + Docker Compose 지원** (권장 배포 방식)
- 🔄 PM2 설정은 레거시 옵션으로 이동 (`legacy/` 폴더)
- 🧠 TabPFN Few-shot Learning 모델 완전 지원
- 🍎 FooDB 식품 화합물 예측 기능 추가
- 📖 Hugging Face 인증 단계별 가이드
- 🎯 소량 데이터(50-1000개)에 최적화된 학습

---

## 💻 시스템 요구사항

### Docker 방식 (권장 ⭐)

#### 최소 사양
- **OS**: Linux (Ubuntu 20.04+), macOS (Intel/M1), Windows 10/11 with WSL2
- **Docker**: 20.10 이상
- **Docker Compose**: 2.0 이상
- **RAM**: 4GB 이상
- **저장공간**: 10GB 이상

#### 권장 사양
- **OS**: Ubuntu 22.04 LTS
- **Docker**: 최신 버전
- **Docker Compose**: 최신 버전
- **RAM**: 16GB 이상 (3D conformer 생성 시)
- **CPU**: 8코어 이상
- **저장공간**: 20GB 이상

### PM2 방식 (레거시)

#### 최소 사양
- **OS**: Linux (Ubuntu 20.04+), macOS, Windows 10/11
- **Python**: 3.10 이상
- **RAM**: 4GB 이상
- **저장공간**: 5GB 이상

#### 권장 사양
- **OS**: Ubuntu 22.04 LTS
- **Python**: 3.11 or 3.12
- **RAM**: 8GB 이상
- **저장공간**: 10GB 이상

---

## 🐳 Docker 배포 방법 (권장)

### 1. Docker 설치

#### Ubuntu/Debian
```bash
# Docker 설치
curl -fsSL https://get.docker.com -o get-docker.sh
sudo sh get-docker.sh

# 현재 사용자를 docker 그룹에 추가
sudo usermod -aG docker $USER
newgrp docker

# Docker Compose 설치 (이미 포함됨)
docker compose version
```

#### macOS
```bash
# Homebrew로 설치
brew install --cask docker

# 또는 공식 사이트에서 다운로드
# https://www.docker.com/products/docker-desktop
```

#### Windows
```bash
# WSL2 활성화 후 Docker Desktop 설치
# https://www.docker.com/products/docker-desktop
```

### 2. 프로젝트 다운로드 및 압축 해제

```bash
# 최신 버전 다운로드 (wget 사용)
wget https://www.genspark.ai/api/files/s/xiJMtcDB -O pci-api-docker.tar.gz

# 또는 curl 사용
curl -L https://www.genspark.ai/api/files/s/xiJMtcDB -o pci-api-docker.tar.gz

# 압축 해제
tar -xzf pci-api-docker.tar.gz

# 프로젝트 디렉토리로 이동
cd webapp
```

### 3. TabPFN 환경변수 설정 (선택사항)

TabPFN 모델을 사용하려면 Hugging Face 토큰이 필요합니다.

```bash
# .env 파일 생성
echo "HF_TOKEN=hf_your_token_here" > .env

# Hugging Face 토큰 발급 방법:
# 1. https://huggingface.co/settings/tokens 접속
# 2. "New token" 클릭 → Read 권한으로 토큰 생성
# 3. 토큰을 복사하여 위 명령어에 붙여넣기
```

### 4. Docker Compose로 실행

```bash
# 컨테이너 빌드 및 실행 (첫 실행 시 5-10분 소요)
docker-compose up -d

# 빌드 로그 확인
docker-compose logs -f

# 서비스 상태 확인
docker-compose ps
```

### 5. 서비스 테스트

```bash
# Health Check
curl http://localhost:3000/health

# 웹 UI 접속
# http://localhost:3000

# API 문서 접속
# http://localhost:3000/docs
```

### 6. Docker 관리 명령어

```bash
# 로그 확인 (실시간)
docker-compose logs -f

# 컨테이너 재시작
docker-compose restart

# 컨테이너 중지
docker-compose down

# 컨테이너 중지 및 볼륨 삭제 (데이터 삭제 주의!)
docker-compose down -v

# 이미지 재빌드 후 실행
docker-compose up -d --build

# 리소스 사용량 확인
docker stats pci-prediction-api
```

### 7. 리소스 제한 조정

`docker-compose.yml` 파일을 수정하여 리소스 제한을 조정할 수 있습니다:

```yaml
deploy:
  resources:
    limits:
      cpus: '16.0'     # CPU 코어 수 조정 (워크스테이션 사양에 맞게)
      memory: 32G      # 메모리 조정
    reservations:
      cpus: '4.0'      # 최소 예약 코어
      memory: 8G       # 최소 예약 메모리
```

수정 후 재시작:
```bash
docker-compose down
docker-compose up -d
```

---

## 🔧 PM2 배포 방법 (레거시)

**⚠️ 주의**: Docker를 사용할 수 없는 환경에서만 사용하세요. RDKit/Mordred 설치 시 환경별 에러가 발생할 수 있습니다.

자세한 내용은 `legacy/README.md` 파일을 참조하세요.

### 1. Python 설치 확인

```bash
python3 --version
# Python 3.10 이상이어야 함
```

### 2. 프로젝트 다운로드 및 압축 해제

```bash
# 최신 버전 다운로드 (wget 사용)
wget https://www.genspark.ai/api/files/s/xiJMtcDB -O pci-api.tar.gz

# 또는 curl 사용
curl -L https://www.genspark.ai/api/files/s/xiJMtcDB -o pci-api.tar.gz

# 압축 해제
tar -xzf pci-api.tar.gz

# 프로젝트 디렉토리로 이동
cd webapp
```

### 3. Python 가상환경 생성 (권장)

```bash
# venv 생성
python3 -m venv venv

# 활성화 (Linux/macOS)
source venv/bin/activate

# 활성화 (Windows)
venv\Scripts\activate
```

### 4. 의존성 패키지 설치

```bash
# requirements.txt 설치
pip install -r requirements.txt

# 설치 시간: 약 5-10분 소요
```

### 5. 서버 실행

#### Option A: 직접 실행 (개발용)

```bash
# Uvicorn으로 직접 실행
uvicorn app.main:app --host 0.0.0.0 --port 3000 --reload

# 브라우저에서 접속
# http://localhost:3000
```

#### Option B: PM2 사용 (프로덕션 권장)

```bash
# PM2 설치 (Node.js 필요)
npm install -g pm2

# PM2로 서버 시작
pm2 start ecosystem.config.cjs

# 서버 상태 확인
pm2 status

# 로그 확인
pm2 logs pci-api

# 서버 재시작
pm2 restart pci-api

# 서버 중지
pm2 stop pci-api
```

---

## 🌐 접속 URL

서버가 정상적으로 실행되면 다음 URL로 접속 가능합니다:

- **웹 UI**: http://localhost:3000
- **API 문서 (Swagger)**: http://localhost:3000/docs
- **API 문서 (ReDoc)**: http://localhost:3000/redoc
- **Health Check**: http://localhost:3000/health

---

## 🧠 TabPFN 모델 사용 설정 (선택사항)

TabPFN은 Few-shot Learning에 특화된 모델로, 소량의 데이터(50개)에서 우수한 성능을 보입니다.

### TabPFN 활성화 방법

#### 1. Hugging Face 계정 생성 및 인증

```bash
# Hugging Face CLI 설치 (이미 requirements.txt에 포함됨)
pip install huggingface-hub

# Hugging Face 로그인
huggingface-cli login
# 또는
hf auth login

# 프롬프트가 나오면 Hugging Face 토큰 입력
# 토큰 생성: https://huggingface.co/settings/tokens
```

#### 2. TabPFN 모델 접근 권한 요청

1. **모델 페이지 방문**: https://huggingface.co/Prior-Labs/tabpfn_2_5
2. **"Accept terms"** 버튼 클릭하여 약관 동의
3. **접근 권한 승인 대기** (보통 즉시 승인됨)

#### 3. 환경변수 설정 (선택)

```bash
# ~/.bashrc 또는 ~/.zshrc에 추가
export HF_TOKEN="your_huggingface_token"

# 또는 .env 파일 생성
echo "HF_TOKEN=your_huggingface_token" > .env
```

#### 4. TabPFN 모델 테스트

```bash
# Python으로 테스트
python3 -c "from tabpfn import TabPFNClassifier; print('TabPFN available!')"

# 또는 API로 테스트
curl -X POST http://localhost:3000/api/models/train \
  -H "Content-Type: application/json" \
  -d '{
    "protein_name": "FTO",
    "model_type": "TabPFN",
    "feature_type": "fingerprint"
  }'
```

### TabPFN vs 다른 모델 비교

| 모델 | 데이터 크기 | 학습 속도 | 정확도 | Few-shot 성능 |
|------|-----------|---------|--------|--------------|
| **TabPFN** | 소량 (50개) | 매우 빠름 | 높음 | ⭐⭐⭐⭐⭐ |
| **XGBoost** | 중간~대량 | 빠름 | 높음 | ⭐⭐⭐ |
| **LightGBM** | 대량 | 매우 빠름 | 높음 | ⭐⭐ |
| **CatBoost** | 중간~대량 | 중간 | 매우 높음 | ⭐⭐⭐ |
| **RandomForest** | 소량~중간 | 중간 | 중간 | ⭐⭐ |

### TabPFN 권장 사용 시나리오

✅ **사용 권장**:
- 데이터가 50-1000개 정도로 적을 때
- 빠른 프로토타이핑이 필요할 때
- Transfer Learning 효과를 극대화하고 싶을 때
- 새로운 단백질에 대한 초기 탐색

❌ **사용 비권장**:
- 데이터가 10,000개 이상일 때
- 매우 복잡한 feature engineering이 필요할 때
- GPU가 없는 환경 (CPU로도 작동하지만 느림)

### TabPFN 문제 해결

#### "Failed to download TabPFN model" 오류

```bash
# 1. 인증 확인
huggingface-cli whoami

# 2. 모델 접근 권한 확인
# https://huggingface.co/Prior-Labs/tabpfn_2_5 에서 Accept terms 클릭

# 3. 토큰 재설정
huggingface-cli login --token YOUR_NEW_TOKEN
```

#### "GPU out of memory" 오류

```python
# CPU 모드로 강제 실행 (app/services/model_training.py)
return TabPFNClassifier(device='cpu')
```

---

## 🔒 방화벽 설정 (외부 접속 허용 시)

### Linux (ufw)

```bash
# 포트 3000 열기
sudo ufw allow 3000/tcp

# 방화벽 상태 확인
sudo ufw status
```

### Windows (PowerShell 관리자 권한)

```powershell
# 포트 3000 인바운드 규칙 추가
New-NetFirewallRule -DisplayName "PCI Prediction API" -Direction Inbound -Protocol TCP -LocalPort 3000 -Action Allow
```

---

## 📂 프로젝트 구조

```
webapp/
├── app/                      # FastAPI 백엔드
│   ├── main.py              # 메인 애플리케이션
│   ├── routers/             # API 라우터
│   │   ├── data_collection.py
│   │   ├── feature_transform.py
│   │   ├── model_training.py
│   │   └── shap_analysis.py
│   ├── services/            # 비즈니스 로직
│   ├── models/              # Pydantic 스키마
│   └── utils/               # 유틸리티 함수
├── static/                  # 프론트엔드
│   ├── index.html
│   ├── css/
│   └── js/
├── saved_data/              # 수집된 데이터
│   └── IC50/
├── raw/                     # 변환된 특성 데이터
│   └── Dataset/
├── models_trained/          # 학습된 모델
├── shap_outputs/            # SHAP 분석 결과
├── requirements.txt         # Python 패키지
├── ecosystem.config.cjs     # PM2 설정
└── README.md               # 프로젝트 문서
```

---

## 🎯 사용 방법

### 1. 데이터 수집

```bash
# API로 직접 호출
curl -X POST http://localhost:3000/api/data/collect \
  -H "Content-Type: application/json" \
  -d '{
    "target_list": ["PDE4", "PDE5", "FTO"],
    "standard_type": "IC50"
  }'
```

### 2. 특성 변환

```bash
curl -X POST http://localhost:3000/api/features/transform \
  -H "Content-Type: application/json" \
  -d '{
    "protein_name": "FTO",
    "fingerprint_type": "MACCS",
    "dataset_ratio": "20x",
    "ignore3D": true
  }'
```

### 3. 모델 학습

```bash
curl -X POST http://localhost:3000/api/models/train \
  -H "Content-Type: application/json" \
  -d '{
    "protein_name": "FTO",
    "model_type": "XGBoost",
    "feature_type": "fingerprint"
  }'
```

### 4. 예측

```bash
curl -X POST http://localhost:3000/api/models/predict \
  -H "Content-Type: application/json" \
  -d '{
    "smiles_list": ["CCO", "CC(=O)O"],
    "model_id": "FTO_MLModelType.XGBOOST_20260306_053810",
    "feature_type": "fingerprint"
  }'
```

---

## 🐳 Docker 배포 (선택사항)

Docker를 사용하면 더 쉽게 배포할 수 있습니다.

### Dockerfile 생성

```dockerfile
FROM python:3.11-slim

WORKDIR /app

# 시스템 패키지 설치
RUN apt-get update && apt-get install -y \
    build-essential \
    && rm -rf /var/lib/apt/lists/*

# Python 패키지 설치
COPY requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt

# 프로젝트 복사
COPY . .

# 포트 노출
EXPOSE 3000

# 서버 실행
CMD ["uvicorn", "app.main:app", "--host", "0.0.0.0", "--port", "3000"]
```

### Docker 빌드 및 실행

```bash
# 이미지 빌드
docker build -t pci-prediction-api .

# 컨테이너 실행
docker run -d -p 3000:3000 --name pci-api pci-prediction-api

# 로그 확인
docker logs -f pci-api
```

---

## 🌍 외부 서버 배포 (AWS, GCP, Azure 등)

### 1. 클라우드 VM 생성
- AWS EC2, GCP Compute Engine, Azure VM 등
- Ubuntu 22.04 LTS 권장
- 최소 t2.medium (2 vCPU, 4GB RAM)

### 2. 보안 그룹 설정
- 포트 3000 인바운드 허용
- SSH 포트 22 허용

### 3. 서버 설정

```bash
# SSH 접속
ssh user@your-server-ip

# 최신 버전 다운로드 (TabPFN 지원)
wget https://www.genspark.ai/api/files/s/AMOsEmIW -O pci-prediction-api-tabpfn.tar.gz
tar -xzf pci-prediction-api-tabpfn.tar.gz
cd webapp

# Python 3.11 설치
sudo apt update
sudo apt install -y python3.11 python3.11-venv python3-pip

# 가상환경 생성 및 패키지 설치
python3.11 -m venv venv
source venv/bin/activate
pip install -r requirements.txt

# (선택) TabPFN 사용 시 Hugging Face 인증
huggingface-cli login
# 토큰 입력 후 https://huggingface.co/Prior-Labs/tabpfn_2_5 에서 Accept terms

# PM2로 서버 시작
npm install -g pm2
pm2 start ecosystem.config.cjs

# PM2 자동 재시작 설정
pm2 startup
pm2 save
```

### 4. Nginx 리버스 프록시 (선택)

```bash
# Nginx 설치
sudo apt install -y nginx

# 설정 파일 생성
sudo nano /etc/nginx/sites-available/pci-api
```

```nginx
server {
    listen 80;
    server_name your-domain.com;

    location / {
        proxy_pass http://localhost:3000;
        proxy_http_version 1.1;
        proxy_set_header Upgrade $http_upgrade;
        proxy_set_header Connection 'upgrade';
        proxy_set_header Host $host;
        proxy_cache_bypass $http_upgrade;
    }
}
```

```bash
# 설정 활성화
sudo ln -s /etc/nginx/sites-available/pci-api /etc/nginx/sites-enabled/
sudo nginx -t
sudo systemctl restart nginx
```

---

## 🔧 문제 해결

### 포트가 이미 사용 중

```bash
# 포트 3000 사용 프로세스 확인
lsof -i :3000

# 프로세스 종료
kill -9 <PID>
```

### 패키지 설치 오류

```bash
# pip 업그레이드
pip install --upgrade pip

# 캐시 클리어 후 재설치
pip cache purge
pip install -r requirements.txt
```

### RDKit 설치 오류

```bash
# conda 사용 (권장)
conda install -c conda-forge rdkit

# 또는 시스템 패키지 설치 후
sudo apt install -y librdkit-dev python3-rdkit
```

---

## 📊 성능 최적화

### 1. Worker 수 증가

```bash
# ecosystem.config.cjs 수정
instances: 4,  # CPU 코어 수에 맞게 조정
exec_mode: 'cluster'
```

### 2. Gunicorn 사용

```bash
pip install gunicorn

# Gunicorn으로 실행 (4 workers)
gunicorn app.main:app -w 4 -k uvicorn.workers.UvicornWorker --bind 0.0.0.0:3000
```

### 3. 캐싱 설정

- Redis 설치하여 descriptor 계산 결과 캐싱
- 대용량 데이터셋은 사전 변환 후 저장

---

## 🔐 보안 권장사항

1. **환경변수 사용**: API 키, DB 비밀번호 등은 `.env` 파일로 관리
2. **HTTPS 설정**: Let's Encrypt로 SSL 인증서 발급
3. **인증 추가**: JWT 토큰 기반 API 인증 구현
4. **Rate Limiting**: 요청 제한으로 DDoS 방어

---

## 📞 지원

문제가 발생하면:
1. 로그 확인: `pm2 logs pci-api --lines 50`
2. Health Check: `curl http://localhost:3000/health`
3. API 문서: http://localhost:3000/docs

---

## 📄 라이선스

이 프로젝트는 연구 목적으로 사용 가능합니다.

---

**배포 성공을 기원합니다!** 🎉
