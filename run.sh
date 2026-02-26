#!/bin/bash

# PCI Prediction Platform - Run Script

echo "========================================="
echo "  PCI Prediction Platform"
echo "  Starting Streamlit Application..."
echo "========================================="

# 현재 디렉토리 확인
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"
cd "$SCRIPT_DIR"

# Python 가상환경 확인 (선택사항)
if [ -d "venv" ]; then
    echo "Activating virtual environment..."
    source venv/bin/activate
fi

# 필요한 디렉토리 생성
mkdir -p data/temp data/results models logs

# Streamlit 실행
echo ""
echo "Starting Streamlit on http://localhost:8501"
echo ""

streamlit run app/streamlit_app.py \
    --server.port 8501 \
    --server.address localhost \
    --server.headless true \
    --browser.serverAddress localhost \
    --browser.gatherUsageStats false

echo ""
echo "Application stopped."
