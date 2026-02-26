# macOS 설치 가이드

## 🍎 macOS에서 PCI Prediction Platform 설치하기

### 필수 요구사항
- macOS 10.15 (Catalina) 이상
- Python 3.8 이상
- 8GB RAM 이상 (16GB 권장)

## 📦 1단계: Homebrew 설치 (선택사항)

```bash
# Homebrew가 없다면 설치
/bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"
```

## 🐍 2단계: Conda 설치 (권장)

### Miniforge 설치 (Apple Silicon/Intel 모두 지원)

```bash
# Miniforge 다운로드 및 설치
curl -L -O "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-MacOSX-$(uname -m).sh"
bash Miniforge3-MacOSX-$(uname -m).sh

# 설치 후 터미널 재시작
source ~/.zshrc  # zsh 사용 시
# 또는
source ~/.bash_profile  # bash 사용 시
```

### 또는 Anaconda 설치

```bash
# Anaconda 다운로드 (https://www.anaconda.com/download)
# 설치 후 터미널 재시작
```

## 🚀 3단계: 프로젝트 설치

```bash
# 1. 저장소 클론
git clone https://github.com/HeejeongH/Functional-Food-Predictor.git
cd Functional-Food-Predictor

# 2. Conda 환경 생성 (권장)
conda create -n pci_platform python=3.9
conda activate pci_platform

# 3. RDKit 설치 (가장 중요!)
conda install -c conda-forge rdkit

# 4. 나머지 패키지 설치
pip install -r requirements.txt

# 5. 실행 권한 부여
chmod +x run.sh
```

## 🎯 4단계: 실행

```bash
# Conda 환경 활성화 (매번 실행 시 필요)
conda activate pci_platform

# 웹앱 실행
./run.sh

# 또는 직접 실행
streamlit run app/streamlit_app.py
```

브라우저에서 자동으로 `http://localhost:8501` 열림

## 🔧 Apple Silicon (M1/M2/M3) 특별 가이드

Apple Silicon Mac을 사용하는 경우:

```bash
# 1. Rosetta 없이 네이티브로 설치 (권장)
# Miniforge를 사용하면 자동으로 ARM64 네이티브로 설치됨

# 2. 환경 생성
conda create -n pci_platform python=3.9
conda activate pci_platform

# 3. RDKit 설치
conda install -c conda-forge rdkit

# 4. 나머지 패키지
pip install -r requirements.txt

# 5. TabPFN 설치 (중요!)
pip install tabpfn --no-deps
pip install torch torchvision
```

## ⚡ 5단계: 테스트

```bash
# 데모 스크립트 실행
python demo.py

# 출력 예시:
# ============================================================
# Testing Molecular Feature Extractor
# ============================================================
# ECFP4 Fingerprint:
# Aspirin         - FP Size: 1024, Non-zero: 45
#                   MW: 180.16, LogP: 1.19
# ...
```

## 🐛 문제 해결

### Q1: "command not found: conda"
```bash
# Conda 경로를 수동으로 추가
echo 'export PATH="$HOME/miniforge3/bin:$PATH"' >> ~/.zshrc
source ~/.zshrc
```

### Q2: RDKit 설치 오류
```bash
# Conda-forge 채널 우선순위 설정
conda config --add channels conda-forge
conda config --set channel_priority strict
conda install rdkit
```

### Q3: TabPFN 설치 오류
```bash
# PyTorch를 먼저 설치
conda install pytorch torchvision -c pytorch
pip install tabpfn
```

### Q4: "Permission denied: ./run.sh"
```bash
chmod +x run.sh
```

### Q5: Streamlit 포트 충돌
```bash
# 다른 포트로 실행
streamlit run app/streamlit_app.py --server.port 8502
```

## 💡 성능 최적화 (Apple Silicon)

Apple Silicon Mac에서 최적 성능을 위해:

```bash
# 1. PyTorch Metal 가속 활용
pip install torch torchvision torchaudio

# 2. NumPy 최적화
conda install "libblas=*=*accelerate"

# 3. 메모리 설정
export PYTORCH_MPS_HIGH_WATERMARK_RATIO=0.0
```

## 🎓 추가 팁

### 1. Jupyter Notebook도 사용하고 싶다면
```bash
conda install jupyter
jupyter notebook
```

### 2. VSCode에서 개발하려면
```bash
# VSCode 설치
brew install --cask visual-studio-code

# Python Extension 설치
code --install-extension ms-python.python
```

### 3. 가상환경 자동 활성화
```bash
# .zshrc 또는 .bash_profile에 추가
echo 'conda activate pci_platform' >> ~/.zshrc
```

### 4. 빠른 시작 alias 만들기
```bash
# ~/.zshrc에 추가
alias pci='cd ~/Functional-Food-Predictor && conda activate pci_platform && ./run.sh'

# 이후 터미널에서 'pci' 입력만으로 실행!
```

## 📊 시스템 요구사항 확인

```bash
# Python 버전 확인
python --version  # 3.8 이상 필요

# 메모리 확인
sysctl hw.memsize  # 8GB 이상 권장

# 디스크 공간
df -h  # 최소 5GB 여유 공간
```

## 🔄 업데이트

```bash
# 최신 버전으로 업데이트
cd Functional-Food-Predictor
git pull origin main
pip install -r requirements.txt --upgrade
```

## ✅ 설치 완료 체크리스트

- [ ] Conda 설치됨
- [ ] Python 3.8+ 설치됨
- [ ] RDKit 설치됨 (`conda install -c conda-forge rdkit`)
- [ ] 모든 패키지 설치됨 (`pip install -r requirements.txt`)
- [ ] 실행 권한 부여됨 (`chmod +x run.sh`)
- [ ] 데모 실행 성공 (`python demo.py`)
- [ ] 웹앱 실행 성공 (`./run.sh`)

## 🆘 추가 도움말

문제가 계속되면:
1. [GitHub Issues](https://github.com/HeejeongH/Functional-Food-Predictor/issues)
2. 에러 메시지 전체 복사해서 문의

---

**설치 시간**: 약 10-15분 (인터넷 속도에 따라 다름)

**macOS에서 가장 안정적으로 작동합니다!** 🍎✨
