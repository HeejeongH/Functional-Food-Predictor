# ============================================
# PCI Prediction API - Dockerfile
# ============================================
# Multi-stage build for optimized image size
# Includes: Python 3.11, RDKit, Mordred, XGBoost, LightGBM, CatBoost, TabPFN
# Target: Research-grade protein-compound interaction prediction

# Stage 1: Builder (compile dependencies)
FROM python:3.11-slim as builder

# Install system build dependencies
RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential \
    gcc \
    g++ \
    cmake \
    libxrender1 \
    libxext6 \
    libgomp1 \
    && rm -rf /var/lib/apt/lists/*

# Create virtual environment
RUN python -m venv /opt/venv
ENV PATH="/opt/venv/bin:$PATH"

# Copy requirements and install Python packages
COPY requirements.txt /tmp/requirements.txt
RUN pip install --no-cache-dir --upgrade pip setuptools wheel && \
    pip install --no-cache-dir -r /tmp/requirements.txt

# Stage 2: Runtime (minimal final image)
FROM python:3.11-slim

# Install runtime dependencies only
RUN apt-get update && apt-get install -y --no-install-recommends \
    libxrender1 \
    libxext6 \
    libgomp1 \
    curl \
    && rm -rf /var/lib/apt/lists/*

# Copy virtual environment from builder
COPY --from=builder /opt/venv /opt/venv
ENV PATH="/opt/venv/bin:$PATH"

# Set working directory
WORKDIR /app

# Copy application code
COPY . .

# Create necessary directories
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

# Set environment variables
ENV PYTHONUNBUFFERED=1 \
    PYTHONDONTWRITEBYTECODE=1 \
    PORT=3000

# Health check
HEALTHCHECK --interval=30s --timeout=10s --start-period=5s --retries=3 \
    CMD curl -f http://localhost:3000/health || exit 1

# Expose port
EXPOSE 3000

# Run application
CMD ["uvicorn", "app.main:app", "--host", "0.0.0.0", "--port", "3000", "--workers", "1"]
