# Cloudflare 배포 가이드

PCI Prediction API를 Cloudflare를 통해 배포하는 방법은 **두 가지**입니다.

---

## 방법 1: Cloudflare Zero Trust Tunnel ⭐ 권장

**전체 스택(백엔드 + 프론트엔드)을 하나의 명령으로 배포합니다.**
공개 IP가 없는 환경(홈 서버, 연구실 PC 등)에서도 사용 가능합니다.

### 아키텍처
```
사용자 브라우저
    ↕ HTTPS (your-domain.com)
Cloudflare Edge Network
    ↕ Encrypted Tunnel
cloudflared container (로컬/서버)
    ↕ HTTP (localhost)
PCI API container (FastAPI)
```

### 설정 절차

1. **Cloudflare Zero Trust 대시보드에서 터널 생성**
   ```
   https://one.dash.cloudflare.com/
   → Networks → Tunnels → Create a tunnel
   → Cloudflared 선택 → 터널 이름 입력 (예: pci-prediction)
   → Install and run a connector → Docker 선택
   → 토큰(eyJ...) 복사
   ```

2. **도메인 연결 (Public Hostname 설정)**
   ```
   → Public Hostname 탭 → Add a public hostname
   → Subdomain: pci-api
   → Domain: 본인 도메인 선택 (예: yourdomain.com)
   → Service Type: HTTP
   → URL: pci-api:3000
   → Save
   ```
   
   결과: `https://pci-api.yourdomain.com` 으로 접근 가능

3. **.env 파일 설정**
   ```bash
   cp .env.example .env
   # .env 파일 수정:
   CLOUDFLARE_TUNNEL_TOKEN=eyJhIjoixxxxxxxx...   # 복사한 터널 토큰
   ```

4. **배포 실행**
   ```bash
   # 방법 A: 자동화 스크립트
   ./deploy-cloudflare.sh tunnel

   # 방법 B: 직접 실행
   docker-compose -f docker-compose.yml -f docker-compose.cloudflare.yml up -d --build
   ```

5. **상태 확인**
   ```bash
   docker-compose -f docker-compose.yml -f docker-compose.cloudflare.yml ps
   docker-compose logs cloudflared
   ```

6. **접속 확인**
   ```bash
   curl https://pci-api.yourdomain.com/health
   ```

### 관리 명령어
```bash
# 중지
docker-compose -f docker-compose.yml -f docker-compose.cloudflare.yml down

# 로그 확인
docker-compose -f docker-compose.yml -f docker-compose.cloudflare.yml logs -f

# 재시작
docker-compose -f docker-compose.yml -f docker-compose.cloudflare.yml restart
```

---

## 방법 2: Cloudflare Pages (프론트엔드만)

정적 HTML/CSS/JS 파일만 Cloudflare Pages CDN에 배포합니다.
백엔드 API는 별도 서버(VPS, 홈서버 등)에서 실행되어야 합니다.

### 설정 절차

1. **Cloudflare API 토큰 발급**
   ```
   https://dash.cloudflare.com/profile/api-tokens
   → Create Token → Custom Token
   → Permissions: Cloudflare Pages - Edit
   → Create Token → 토큰 복사
   ```

2. **.env 파일 설정**
   ```bash
   CLOUDFLARE_API_TOKEN=your_api_token_here
   BACKEND_URL=https://your-backend-server.com   # 백엔드 서버 URL
   ```

3. **Pages 배포**
   ```bash
   ./deploy-cloudflare.sh pages pci-prediction
   ```
   
   또는 수동으로:
   ```bash
   export CLOUDFLARE_API_TOKEN=your_token
   npx wrangler pages deploy static --project-name pci-prediction
   ```

4. **접속 확인**
   ```
   https://pci-prediction.pages.dev
   ```

### 주의사항
- 프론트엔드만 배포되므로 API 기능은 `BACKEND_URL`에서 제공해야 합니다
- CORS 설정이 올바르게 되어 있어야 합니다

---

## Cloudflare 추가 기능 (선택사항)

### Access Policy (접근 제한)
Zero Trust 대시보드에서 인증된 사용자만 접근하도록 설정:
```
https://one.dash.cloudflare.com/
→ Access → Applications → Add an application
→ Self-hosted 선택
→ 도메인 및 접근 정책 설정
```

### Rate Limiting
```
Cloudflare 대시보드 → Security → WAF → Rate limiting rules
→ API 엔드포인트(/api/*)에 분당 요청 제한 설정
```

### DDoS 보호
Cloudflare를 통해 자동으로 DDoS 보호가 활성화됩니다.

---

## 트러블슈팅

| 문제 | 해결 방법 |
|------|-----------|
| 터널이 연결되지 않음 | `CLOUDFLARE_TUNNEL_TOKEN` 확인, cloudflared 컨테이너 로그 확인 |
| API 응답 없음 | `pci-api` 컨테이너 상태 확인: `docker-compose ps` |
| CORS 오류 | `app/main.py`의 `allow_origins` 설정 확인 |
| Pages 배포 실패 | API 토큰 권한(Pages:Edit) 확인 |
| SSL 인증서 오류 | Cloudflare에서 SSL/TLS 모드를 "Full" 또는 "Flexible"로 설정 |

---

*마지막 업데이트: 2026-07-13*
