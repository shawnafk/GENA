#!/usr/bin/env bash
set -euo pipefail
DIGITALQQ_ROOT="$(cd "$(dirname "$0")/.." && pwd)"
cd "$DIGITALQQ_ROOT"
if [ ! -d digitalQQ/frontend/dist ]; then
  (cd digitalQQ/frontend && npm ci && npm run build)
fi
exec python -m uvicorn digitalQQ.backend.app:app --host 127.0.0.1 --port 8765
