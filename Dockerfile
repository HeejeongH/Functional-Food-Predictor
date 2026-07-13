# ============================================
# PCI Prediction API - Dockerfile (GPU)
# ============================================
# Base: NVIDIA CUDA 12.9 + Python 3.11
# GPU: RTX 5070 Ti (CUDA 12.9, cuDNN 9)
# Includes: RDKit, Mordred, XGBoost, LightGBM, CatBoost, TabPFN (GPU)

# Stage 1: Builder
FROM nvidia/cuda:12.9.0-cudnn-devel-ubuntu22.04 AS builder

# 시간대 설정 (apt 설치 중 대화창 방지)
ENV DEBIAN_FRONTEND=noninteractive
ENV TZ=Asia/Seoul

RUN apt-get update && apt-get install -y --no-install-recommends \
    python3.11 \
    python3.11-venv \
    python3.11-dev \
    python3-pip \
    build-essential \
    gcc \
    g++ \
    cmake \
    libxrender1 \
    libxext6 \
    libgomp1 \
    && rm -rf /var/lib/apt/lists/*

# python3.11을 기본 python으로 설정
RUN update-alternatives --install /usr/bin/python python /usr/bin/python3.11 1 && \
    update-alternatives --install /usr/bin/python3 python3 /usr/bin/python3.11 1

# 가상환경 생성
RUN python -m venv /opt/venv
ENV PATH="/opt/venv/bin:$PATH"

# pip 업그레이드
RUN pip install --no-cache-dir --upgrade pip setuptools wheel

# PyTorch GPU (CUDA 12.9 호환 — cu128 빌드 사용)
RUN pip install --no-cache-dir \
    torch torchvision torchaudio \
    --index-url https://download.pytorch.org/whl/cu128

# 나머지 패키지 설치
COPY requirements.txt /tmp/requirements.txt
RUN pip install --no-cache-dir -r /tmp/requirements.txt

# TabPFN
# TabPFN v2 (최신, Nature 2025)
RUN pip install --no-cache-dir "tabpfn>=2.0"


# Stage 2: Runtime
FROM nvidia/cuda:12.9.0-cudnn-runtime-ubuntu22.04

ENV DEBIAN_FRONTEND=noninteractive
ENV TZ=Asia/Seoul

RUN apt-get update && apt-get install -y --no-install-recommends \
    python3.11 \
    python3.11-venv \
    libxrender1 \
    libxext6 \
    libgomp1 \
    curl \
    && rm -rf /var/lib/apt/lists/*

RUN update-alternatives --install /usr/bin/python python /usr/bin/python3.11 1 && \
    update-alternatives --install /usr/bin/python3 python3 /usr/bin/python3.11 1

# Builder에서 가상환경 복사
COPY --from=builder /opt/venv /opt/venv
ENV PATH="/opt/venv/bin:$PATH"

# CUDA 환경변수
ENV CUDA_HOME=/usr/local/cuda
ENV LD_LIBRARY_PATH=/usr/local/cuda/lib64:$LD_LIBRARY_PATH

WORKDIR /app
COPY . .

# 필요 디렉토리 생성
RUN mkdir -p \
    saved_data/IC50 \
    saved_data/BindingDB \
    saved_data/FooDB \
    saved_data/FooDB/3d_conformers \
    saved_data/FooDB/preprocessed \
    raw/FewshotSet \
    raw/TransferSet \
    raw/Dataset \
    raw/descriptors \
    models_trained \
    shap_outputs \
    food_predictions

ENV PYTHONUNBUFFERED=1 \
    PYTHONDONTWRITEBYTECODE=1 \
    PORT=3000

HEALTHCHECK --interval=30s --timeout=10s --start-period=30s --retries=3 \
    CMD curl -f http://localhost:3000/health || exit 1

EXPOSE 3000

CMD ["uvicorn", "app.main:app", "--host", "0.0.0.0", "--port", "3000", "--workers", "1"]
