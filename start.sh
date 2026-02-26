#!/bin/bash

echo "======================================================"
echo "🚀 PCI Prediction Platform - Professional Edition 🚀"
echo "======================================================"
echo ""

# 색상 정의
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

# 프로젝트 루트로 이동
cd "$(dirname "$0")"

echo -e "${BLUE}📦 Step 1: 의존성 확인 중...${NC}"

# Python 의존성 확인
if ! python3 -c "import flask" 2>/dev/null; then
    echo -e "${YELLOW}⚠️  Flask가 설치되지 않았습니다. 설치 중...${NC}"
    pip install -r requirements.txt
fi

# Node.js 확인
if ! command -v node &> /dev/null; then
    echo -e "${RED}❌ Node.js가 설치되지 않았습니다!${NC}"
    echo "Node.js를 먼저 설치해주세요: https://nodejs.org/"
    exit 1
fi

echo -e "${GREEN}✅ Python 의존성 확인 완료${NC}"

# 프론트엔드 의존성 설치
if [ ! -d "frontend/node_modules" ]; then
    echo -e "${BLUE}📦 Step 2: 프론트엔드 의존성 설치 중...${NC}"
    cd frontend
    npm install
    cd ..
    echo -e "${GREEN}✅ 프론트엔드 의존성 설치 완료${NC}"
else
    echo -e "${GREEN}✅ 프론트엔드 의존성 이미 설치됨${NC}"
fi

echo ""
echo -e "${BLUE}🚀 Step 3: 서버 시작 중...${NC}"
echo ""

# 백엔드 서버 시작 (백그라운드)
echo -e "${BLUE}🔧 Flask 백엔드 서버 시작 중... (포트 5000)${NC}"
cd backend
python3 app.py &
BACKEND_PID=$!
cd ..

# 잠시 대기 (백엔드가 시작될 시간을 줌)
sleep 3

# 프론트엔드 서버 시작
echo -e "${BLUE}🎨 React 프론트엔드 서버 시작 중... (포트 3000)${NC}"
cd frontend
npm run dev &
FRONTEND_PID=$!
cd ..

echo ""
echo -e "${GREEN}======================================================"
echo "✅ 서버가 성공적으로 시작되었습니다!"
echo "======================================================"
echo ""
echo -e "🌐 웹 애플리케이션: ${BLUE}http://localhost:3000${NC}"
echo -e "📡 API 백엔드: ${BLUE}http://localhost:5000${NC}"
echo ""
echo -e "종료하려면 ${RED}Ctrl+C${NC}를 누르세요"
echo -e "${NC}"

# 종료 시그널 처리
trap "echo ''; echo '🛑 서버 종료 중...'; kill $BACKEND_PID $FRONTEND_PID 2>/dev/null; exit" INT TERM

# 서버가 실행되는 동안 대기
wait
