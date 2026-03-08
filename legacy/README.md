# Legacy PM2 Configuration
# ========================
# This folder contains PM2-based deployment files for environments where Docker is not available.
# 
# **IMPORTANT**: Docker is the recommended deployment method.
# Only use PM2 if:
# - Docker cannot be installed on your system
# - You need quick local development without Docker overhead
# - You're debugging environment-specific issues
#
# Files:
# - ecosystem.config.cjs: PM2 process manager configuration
#
# Usage:
# ```bash
# pm2 start legacy/ecosystem.config.cjs
# pm2 logs pci-api --nostream
# pm2 restart pci-api
# pm2 delete pci-api
# ```
#
# For production deployment, use Docker instead:
# ```bash
# docker-compose up -d
# ```
