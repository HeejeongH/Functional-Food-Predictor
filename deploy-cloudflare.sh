#!/bin/bash
# ============================================
# Cloudflare 배포 자동화 스크립트
# ============================================
# 방법 1: Cloudflare Tunnel (권장 — 서버가 필요함)
# 방법 2: Cloudflare Pages + Worker Proxy (정적 프론트엔드만)
#
# 사용법: ./deploy-cloudflare.sh [tunnel|pages] [프로젝트명]

set -e

METHOD=${1:-tunnel}
PROJECT=${2:-pci-prediction}

echo "======================================"
echo "  PCI Prediction — Cloudflare 배포"
echo "======================================"
echo ""

# ── 공통 사전 체크 ──────────────────────
check_env() {
  if [ ! -f .env ]; then
    echo "⚠  .env 파일이 없습니다. .env.example을 복사하세요:"
    echo "   cp .env.example .env && nano .env"
    exit 1
  fi
  source .env
}

# ──────────────────────────────────────────
# 방법 1: Cloudflare Tunnel
# ──────────────────────────────────────────
deploy_tunnel() {
  echo "▶ 방법 1: Cloudflare Tunnel 배포"
  echo ""

  if [ -z "$CLOUDFLARE_TUNNEL_TOKEN" ]; then
    echo "❌ CLOUDFLARE_TUNNEL_TOKEN 이 .env에 설정되지 않았습니다."
    echo ""
    echo "토큰 발급 방법:"
    echo "  1. https://one.dash.cloudflare.com 접속"
    echo "  2. Networks → Tunnels → Create a tunnel"
    echo "  3. 터널 이름 입력 (예: pci-prediction)"
    echo "  4. Docker 선택 후 토큰 복사"
    echo "  5. .env 파일에 추가: CLOUDFLARE_TUNNEL_TOKEN=eyJ..."
    echo "  6. Cloudflare 대시보드에서 Public Hostname 설정:"
    echo "     - Subdomain: pci-api"
    echo "     - Domain: 본인 도메인"
    echo "     - Service: http://pci-api:3000"
    echo ""
    exit 1
  fi

  echo "✅ CLOUDFLARE_TUNNEL_TOKEN 확인됨"
  echo ""

  # Docker Compose로 실행
  echo "🐳 Docker Compose 시작 (API + Cloudflare Tunnel)..."
  docker-compose -f docker-compose.yml -f docker-compose.cloudflare.yml up -d --build

  echo ""
  echo "✅ 배포 완료!"
  echo ""
  echo "📋 상태 확인:"
  docker-compose -f docker-compose.yml -f docker-compose.cloudflare.yml ps
  echo ""
  echo "📝 다음 단계:"
  echo "  - Cloudflare Zero Trust 대시보드에서 Public Hostname 설정"
  echo "  - 도메인: https://pci-api.your-domain.com → Service: http://pci-api:3000"
}

# ──────────────────────────────────────────
# 방법 2: Cloudflare Pages (프론트엔드만)
# ──────────────────────────────────────────
deploy_pages() {
  echo "▶ 방법 2: Cloudflare Pages 배포 (프론트엔드 전용)"
  echo ""
  echo "⚠  주의: 이 방법은 정적 파일(static/)만 배포합니다."
  echo "   백엔드 API는 별도 서버에서 실행되어야 합니다."
  echo ""

  # wrangler 설치 확인
  if ! command -v npx &> /dev/null; then
    echo "❌ npx가 없습니다. Node.js를 설치하세요."
    exit 1
  fi

  # Cloudflare API 토큰 확인
  if [ -z "$CLOUDFLARE_API_TOKEN" ]; then
    echo "❌ CLOUDFLARE_API_TOKEN 이 .env에 설정되지 않았습니다."
    echo ""
    echo "API 토큰 발급:"
    echo "  1. https://dash.cloudflare.com/profile/api-tokens"
    echo "  2. Create Token → Custom Token"
    echo "  3. Permissions: Cloudflare Pages - Edit"
    echo "  4. .env에 추가: CLOUDFLARE_API_TOKEN=your_token"
    exit 1
  fi

  export CLOUDFLARE_API_TOKEN

  # Pages 배포
  echo "📦 static/ 폴더를 Cloudflare Pages에 배포합니다..."
  
  # API_URL을 실제 백엔드 URL로 교체
  if [ -n "$BACKEND_URL" ]; then
    echo "🔗 API URL을 ${BACKEND_URL}으로 설정합니다..."
    sed -i "s|const API = window.location.origin;|const API = '${BACKEND_URL}';|g" static/js/workflow.js
  fi

  npx wrangler pages deploy static \
    --project-name "$PROJECT" \
    --branch main

  echo ""
  echo "✅ Pages 배포 완료!"
  echo "🌐 URL: https://${PROJECT}.pages.dev"

  # API_URL 원복
  if [ -n "$BACKEND_URL" ]; then
    sed -i "s|const API = '${BACKEND_URL}';|const API = window.location.origin;|g" static/js/workflow.js
  fi
}

# ── 실행 ──────────────────────────────────
check_env

case "$METHOD" in
  tunnel) deploy_tunnel ;;
  pages)  deploy_pages  ;;
  *)
    echo "사용법: $0 [tunnel|pages] [프로젝트명]"
    echo ""
    echo "  tunnel  — Cloudflare Zero Trust Tunnel (백엔드 포함 전체 배포, 권장)"
    echo "  pages   — Cloudflare Pages (정적 프론트엔드만 배포)"
    echo ""
    echo "예시:"
    echo "  $0 tunnel                   # 터널로 전체 배포"
    echo "  $0 pages pci-prediction    # Pages로 프론트엔드 배포"
    exit 1
    ;;
esac
