# 🪟 Windows 환경 설치 가이드

## 목차
1. [Docker Desktop 방식 (권장)](#option-1-docker-desktop-권장-)
2. [Anaconda 방식 (Docker 없이)](#option-2-anaconda-방식-docker-없이)
3. [문제 해결](#문제-해결)

---

## Option 1: Docker Desktop (권장 ⭐)

### 시스템 요구사항
- **Windows 10 64-bit**: Pro, Enterprise, Education (Build 19041 이상)
- **Windows 11 64-bit**: Home, Pro, Enterprise, Education
- **BIOS에서 가상화 활성화** (CPU virtualization/VT-x)
- **최소 4GB RAM** (권장 8GB+)

### 1. WSL2 설치

**PowerShell을 관리자 권한으로 실행**:

```powershell
# WSL2 설치
wsl --install

# 재부팅
Restart-Computer
```

**재부팅 후 확인**:
```powershell
wsl --list --verbose
```

**예상 출력**:
```
  NAME      STATE           VERSION
* Ubuntu    Running         2
```

**만약 VERSION이 1이면**:
```powershell
wsl --set-default-version 2
wsl --set-version Ubuntu 2
```

### 2. Docker Desktop 설치

1. **다운로드**: https://www.docker.com/products/docker-desktop
2. **실행**: `Docker Desktop Installer.exe`
3. **설치 옵션**:
   - ✅ "Use WSL 2 instead of Hyper-V" 체크
   - ✅ "Add shortcut to desktop" 체크
4. **재부팅**
5. **Docker Desktop 실행**

**설정 확인**:
- Docker Desktop → Settings (⚙️) → General
- ✅ "Use the WSL 2 based engine" 활성화
- Docker Desktop → Settings → Resources → WSL Integration
- ✅ "Ubuntu" 체크

### 3. 프로젝트 다운로드 및 실행

**PowerShell 또는 Windows Terminal에서**:

```powershell
# WSL2 Ubuntu로 진입
wsl

# 홈 디렉토리로 이동
cd ~

# 프로젝트 다운로드
wget https://www.genspark.ai/api/files/s/zcdmXoey -O pci-api.tar.gz

# 압축 해제
tar -xzf pci-api.tar.gz
cd webapp

# TabPFN 사용 시 (선택)
echo "HF_TOKEN=hf_your_token_here" > .env

# Docker Compose 실행
docker-compose up -d

# 로그 확인
docker-compose logs -f
```

### 4. 접속 확인

**Windows 브라우저에서**:
- Web UI: http://localhost:3000
- API 문서: http://localhost:3000/docs
- Health Check: http://localhost:3000/health

### 5. 관리 명령어

```bash
# WSL2 내부에서 실행
docker-compose logs -f         # 로그 확인
docker-compose restart         # 재시작
docker-compose down            # 중지
docker-compose up -d --build   # 재빌드
```

---

## Option 2: Anaconda 방식 (Docker 없이)

Docker Desktop을 사용할 수 없는 경우 (Windows Home Edition Build < 19041, 회사 정책 등)

### 1. Anaconda 설치

**다운로드**: https://www.anaconda.com/download

**설치 옵션**:
- ✅ "Add Anaconda to my PATH environment variable" (권장하지 않지만 편의상 체크)
- ✅ "Register Anaconda as my default Python 3.x"

### 2. Anaconda Prompt 실행

**시작 메뉴 → Anaconda3 → Anaconda Prompt**

### 3. 프로젝트 다운로드

```bash
# 작업 디렉토리로 이동 (예: C:\Users\YourName\Documents)
cd %USERPROFILE%\Documents

# 프로젝트 다운로드 (curl 사용)
curl -L https://www.genspark.ai/api/files/s/zcdmXoey -o pci-api.tar.gz

# 압축 해제 (7-Zip 또는 Windows 기본 도구 사용)
tar -xzf pci-api.tar.gz
cd webapp
```

**또는 Windows 탐색기에서**:
1. 브라우저에서 https://www.genspark.ai/api/files/s/zcdmXoey 다운로드
2. 파일 우클릭 → 압축 풀기
3. `webapp` 폴더로 이동

### 4. Conda 가상환경 생성

```bash
# Python 3.11 가상환경 생성
conda create -n pci-api python=3.11 -y

# 가상환경 활성화
conda activate pci-api

# 기본 과학 패키지 설치 (권장)
conda install -c conda-forge rdkit -y
```

### 5. Python 패키지 설치

```bash
# requirements.txt 설치
pip install -r requirements.txt

# 설치 시간: 약 10-15분 소요
```

**자주 발생하는 에러**:
- **RDKit 설치 실패**: `conda install -c conda-forge rdkit` 먼저 실행
- **Microsoft Visual C++ 14.0 required**: https://visualstudio.microsoft.com/downloads/ 에서 "Build Tools for Visual Studio" 설치

### 6. 서버 실행

**Option A: Uvicorn 직접 실행 (개발용)**

```bash
# 가상환경이 활성화된 상태에서
uvicorn app.main:app --host 0.0.0.0 --port 3000 --reload
```

**Option B: PM2 사용 (프로덕션용)**

```bash
# Node.js 설치 (https://nodejs.org/)
# 설치 후 Anaconda Prompt 재시작

# PM2 설치
npm install -g pm2

# PM2로 서버 실행
pm2 start legacy/ecosystem.config.cjs

# 로그 확인
pm2 logs pci-api --nostream

# 재시작
pm2 restart pci-api

# 중지
pm2 delete pci-api
```

### 7. TabPFN 사용 설정 (선택)

```bash
# Hugging Face CLI 설치
pip install huggingface-hub

# 로그인
huggingface-cli login
# 프롬프트에 토큰 입력 (https://huggingface.co/settings/tokens)

# TabPFN 접근 권한 요청
# https://huggingface.co/Prior-Labs/tabpfn_2_5 접속 → "Agree and access repository"

# 테스트
python -c "from tabpfn import TabPFNClassifier; print('TabPFN available!')"
```

---

## 문제 해결

### Docker Desktop 관련

#### 문제: "WSL 2 installation is incomplete"
```powershell
# PowerShell 관리자 권한
wsl --update
wsl --set-default-version 2
```

#### 문제: "Docker daemon is not running"
- Docker Desktop 재시작
- Windows 재부팅
- BIOS에서 가상화(VT-x/AMD-V) 활성화 확인

#### 문제: "error during connect: This error may indicate that the docker daemon is not running"
```powershell
# Docker Desktop 완전 재시작
taskkill /IM "Docker Desktop.exe" /F
Start-Process "C:\Program Files\Docker\Docker\Docker Desktop.exe"
```

#### 문제: 포트 3000이 이미 사용 중
```powershell
# 포트 사용 프로세스 확인
netstat -ano | findstr :3000

# PID로 프로세스 종료 (예: PID 1234)
taskkill /PID 1234 /F
```

### Anaconda 방식 관련

#### 문제: "error: Microsoft Visual C++ 14.0 or greater is required"
**해결**:
1. https://visualstudio.microsoft.com/downloads/
2. "Build Tools for Visual Studio 2022" 다운로드
3. "Desktop development with C++" 워크로드 설치

#### 문제: RDKit 설치 실패
```bash
# Conda로 먼저 설치
conda install -c conda-forge rdkit -y

# 그 후 pip install
pip install -r requirements.txt
```

#### 문제: "ModuleNotFoundError: No module named 'rdkit'"
```bash
# 가상환경 재생성
conda deactivate
conda remove -n pci-api --all -y
conda create -n pci-api python=3.11 -y
conda activate pci-api
conda install -c conda-forge rdkit -y
pip install -r requirements.txt
```

#### 문제: Mordred 계산 시 "RuntimeError: 3D coordinates are required"
- `ignore3D=True` 옵션 사용
- API 호출 시 `"ignore3D": true` 파라미터 추가

### 일반 문제

#### 문제: "Address already in use" (포트 3000 충돌)
**Windows에서 포트 확인**:
```powershell
netstat -ano | findstr :3000
```

**프로세스 종료**:
```powershell
taskkill /PID <PID> /F
```

#### 문제: 메모리 부족 (SHAP 분석 또는 3D conformer 생성 시)
**Docker 방식**:
```yaml
# docker-compose.yml 수정
deploy:
  resources:
    limits:
      memory: 8G  # 메모리 증가
```

**Anaconda 방식**:
- 배치 크기 줄이기 (100 → 50개)
- `ignore3D=True` 사용 (2D descriptor만)

---

## 권장 방식 요약

| 상황 | 권장 방식 |
|------|----------|
| **Windows 10/11 Pro/Enterprise** | Docker Desktop ⭐ |
| **Windows 10/11 Home** | Docker Desktop ⭐ |
| **회사 PC (Docker 설치 불가)** | Anaconda + PM2 |
| **저사양 PC (RAM < 8GB)** | Anaconda + Uvicorn |
| **개발/디버깅** | Anaconda + Uvicorn |

---

## 추가 리소스

- **Docker Desktop 공식 문서**: https://docs.docker.com/desktop/windows/
- **WSL2 공식 문서**: https://learn.microsoft.com/en-us/windows/wsl/
- **Anaconda 공식 문서**: https://docs.anaconda.com/
- **RDKit 설치 가이드**: https://www.rdkit.org/docs/Install.html
- **Hugging Face 토큰 발급**: https://huggingface.co/settings/tokens

---

## 연락처

문제가 계속 발생하면:
1. GitHub Issues: (프로젝트 GitHub URL)
2. Docker Desktop 로그: `%LOCALAPPDATA%\Docker\log.txt`
3. PM2 로그: `pm2 logs pci-api --lines 50`
