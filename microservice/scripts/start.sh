#!/usr/bin/env bash
set -euo pipefail

export APP_ROOT=${APP_ROOT:-/app}
export PLUMBER_PORT=${PLUMBER_PORT:-8002}
export PLUMBER_HOST=${PLUMBER_HOST:-0.0.0.0}

echo "Starting AmPEP30 Final microservice on ${PLUMBER_HOST}:${PLUMBER_PORT}"
cd /app
R -e "pr <- plumber::plumb('microservice/api/ampep30_final_api.R'); pr\$run(host='${PLUMBER_HOST}', port=as.integer('${PLUMBER_PORT}'))"

