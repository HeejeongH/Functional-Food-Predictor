module.exports = {
  apps: [
    {
      name: 'pci-api',
      script: 'uvicorn',
      args: 'app.main:app --host 0.0.0.0 --port 3000 --reload --timeout-keep-alive 300',
      cwd: '/home/user/webapp',
      interpreter: 'python3',
      env: {
        NODE_ENV: 'development',
        PORT: 3000,
        PYTHONPATH: '/home/user/webapp'
      },
      watch: false,
      instances: 1,
      exec_mode: 'fork',
      autorestart: true,
      max_restarts: 10,
      min_uptime: '10s'
    }
  ]
}
