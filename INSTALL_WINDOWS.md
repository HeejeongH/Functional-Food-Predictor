# Windows 설치 가이드

## 🪟 Windows에서 PCI Prediction Platform 설치하기

### 필수 요구사항
- Windows 10/11 (64-bit)
- Python 3.8 이상
- 8GB RAM 이상 (16GB 권장)

⚠️ **주의**: Windows에서는 몇 가지 추가 단계가 필요합니다.

## 📦 1단계: Anaconda 설치 (필수!)

Windows에서는 RDKit 설치를 위해 **Anaconda가 필수**입니다.

```bash
# 1. Anaconda 다운로드
# https://www.anaconda.com/download 방문
# Windows 64-bit 버전 다운로드

# 2. 설치 시 옵션
☑ Add Anaconda to PATH (권장)
☑ Register Anaconda as default Python

# 3. 설치 완료 후 Anaconda Prompt 실행
# 시작 메뉴 > Anaconda3 > Anaconda Prompt
```

## 🚀 2단계: 프로젝트 설치

**Anaconda Prompt**를 관리자 권한으로 실행하세요!

```bash
# 1. 저장소 클론
git clone https://github.com/HeejeongH/Functional-Food-Predictor.git
cd Functional-Food-Predictor

# 2. Conda 환경 생성
conda create -n pci_platform python=3.9
conda activate pci_platform

# 3. RDKit 설치 (중요!)
conda install -c conda-forge rdkit

# 4. 나머지 패키지 설치
pip install -r requirements.txt
```

## 🎯 3단계: 실행

### 방법 1: Anaconda Prompt에서 실행 (권장)

```bash
# 1. Anaconda Prompt 실행
# 2. 환경 활성화
conda activate pci_platform

# 3. 프로젝트 폴더로 이동
cd C:\Users\YourName\Functional-Food-Predictor

# 4. Streamlit 실행
streamlit run app/streamlit_app.py
```

### 방법 2: PowerShell 스크립트 생성

프로젝트 폴더에 `run.ps1` 파일 생성:

```powershell
# run.ps1
Write-Host "PCI Prediction Platform Starting..." -ForegroundColor Green

# Conda 환경 활성화
& "C:\Users\YourName\anaconda3\Scripts\activate.bat" pci_platform

# Streamlit 실행
streamlit run app/streamlit_app.py --server.port 8501
```

실행:
```powershell
# PowerShell에서
.\run.ps1
```

### 방법 3: 배치 파일 생성 (가장 쉬움!)

프로젝트 폴더에 `run.bat` 파일 생성:

```batch
@echo off
echo =========================================
echo   PCI Prediction Platform
echo   Starting Streamlit Application...
echo =========================================

call C:\Users\YourName\anaconda3\Scripts\activate.bat pci_platform

streamlit run app\streamlit_app.py --server.port 8501

pause
```

더블클릭으로 실행!

## 🔧 4단계: Git 설치 (선택사항)

Git이 없다면:

```bash
# 1. Git 다운로드
# https://git-scm.com/download/win

# 2. 설치 시 옵션
☑ Git Bash
☑ Git GUI
☑ Add to PATH
```

## 🐛 문제 해결

### Q1: "conda: command not found"

**해결책**:
```bash
# 1. Anaconda Prompt를 사용하세요 (일반 CMD가 아닌)
# 2. 또는 PATH에 Anaconda 추가:
# 시스템 환경 변수에 추가:
C:\Users\YourName\anaconda3
C:\Users\YourName\anaconda3\Scripts
C:\Users\YourName\anaconda3\Library\bin
```

### Q2: RDKit 설치 오류

**해결책**:
```bash
# Anaconda Prompt에서:
conda config --add channels conda-forge
conda config --set channel_priority strict
conda install rdkit
```

### Q3: "Microsoft Visual C++ 14.0 is required"

**해결책**:
```bash
# 1. Visual Studio Build Tools 설치
# https://visualstudio.microsoft.com/downloads/
# "Build Tools for Visual Studio" 다운로드

# 2. 설치 시 "C++ build tools" 선택
```

### Q4: ChEMBL API 연결 오류

**해결책**:
```bash
# 1. 방화벽 확인
# Windows Defender 방화벽 > Python 허용

# 2. SSL 인증서 업데이트
pip install --upgrade certifi
```

### Q5: Streamlit 실행 시 포트 오류

**해결책**:
```bash
# 다른 포트 사용
streamlit run app/streamlit_app.py --server.port 8502
```

### Q6: run.sh 실행 안됨

Windows에서는 `.sh` 파일이 기본적으로 실행되지 않습니다.

**해결책**:
- **방법 1**: Git Bash 사용
  ```bash
  # Git Bash 실행
  ./run.sh
  ```

- **방법 2**: `run.bat` 사용 (위의 배치 파일 생성)

- **방법 3**: 직접 명령어 입력
  ```bash
  streamlit run app/streamlit_app.py
  ```

## 💡 Windows 전용 팁

### 1. 바탕화면 바로가기 만들기

`PCI Platform.bat` 파일 생성:
```batch
@echo off
cd C:\Users\YourName\Functional-Food-Predictor
call C:\Users\YourName\anaconda3\Scripts\activate.bat pci_platform
streamlit run app\streamlit_app.py
```

바탕화면에 바로가기 생성!

### 2. Windows Terminal 사용 (권장)

Microsoft Store에서 "Windows Terminal" 설치
- 탭 지원
- 더 나은 색상
- 복사/붙여넣기 편리

### 3. WSL2 사용 (고급)

Linux 환경을 원한다면:
```powershell
# PowerShell (관리자)
wsl --install
wsl --set-default-version 2
```

그 다음 Ubuntu에서 Linux처럼 설치!

### 4. 경로 주의사항

Windows에서는 백슬래시(`\`) 사용:
```bash
# ❌ 잘못된 예
cd /home/user/project

# ✅ 올바른 예
cd C:\Users\YourName\Functional-Food-Predictor
```

## 📊 성능 비교

| 항목 | Windows | macOS |
|------|---------|-------|
| 설치 난이도 | ⭐⭐⭐ | ⭐⭐ |
| RDKit 설치 | Anaconda 필수 | Conda 권장 |
| 실행 속도 | 보통 | 빠름 |
| 안정성 | 보통 | 높음 |
| 추천도 | ⭐⭐⭐ | ⭐⭐⭐⭐⭐ |

## ✅ 설치 완료 체크리스트

- [ ] Anaconda 설치됨
- [ ] Anaconda Prompt 실행 가능
- [ ] Python 3.8+ 설치됨
- [ ] RDKit 설치됨
- [ ] 모든 패키지 설치됨
- [ ] 배치 파일 생성됨 (선택)
- [ ] 데모 실행 성공
- [ ] 웹앱 실행 성공

## 🎓 권장 개발 환경

### VSCode 설치
```bash
# 1. VSCode 다운로드
# https://code.visualstudio.com/

# 2. Extensions 설치
- Python (Microsoft)
- Pylance
- Jupyter
```

### Anaconda Navigator 사용

GUI로 쉽게 관리:
```bash
# 시작 메뉴 > Anaconda Navigator
1. Environments > pci_platform 선택
2. 패키지 관리 UI
3. Jupyter Lab 실행 가능
```

## 🆘 추가 도움말

### 일반적인 Windows 문제

1. **경로에 공백 있음**: 
   ```bash
   # ❌
   cd C:\Program Files\project
   
   # ✅
   cd "C:\Program Files\project"
   ```

2. **관리자 권한 필요**:
   - Anaconda Prompt를 우클릭 > 관리자 권한으로 실행

3. **방화벽/백신**:
   - Python.exe 허용
   - Streamlit 포트 허용

4. **인코딩 문제**:
   ```bash
   # 파일 저장 시 UTF-8 사용
   # VSCode에서: 하단 인코딩 > UTF-8 선택
   ```

## 🔄 업데이트

```bash
# Anaconda Prompt에서
cd C:\Users\YourName\Functional-Food-Predictor
git pull origin main
conda activate pci_platform
pip install -r requirements.txt --upgrade
```

---

**설치 시간**: 약 20-30분 (인터넷 속도에 따라 다름)

**Windows보다 macOS가 더 권장됩니다!** 

하지만 Windows에서도 충분히 사용 가능합니다. 🪟✨
