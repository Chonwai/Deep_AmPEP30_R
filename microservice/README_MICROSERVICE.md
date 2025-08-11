## Deep-AmPEP30 Microservice

### Quickstart

Local (requires R and packages):

```bash
R -e "pr <- plumber::plumb('microservice/api/plumber.R'); pr$run(host='0.0.0.0', port=8001)"
```

Docker:

```bash
docker build -f microservice/docker/Dockerfile -t deep-ampep30 .
docker run --rm -p 8001:8001 deep-ampep30
```

Compose:

```bash
docker compose -f microservice/docker/docker-compose.yml up --build
```

### API

- GET /health
- POST /api/predict

Example request:

```bash
curl -X POST http://localhost:8001/api/predict \
  -H 'Content-Type: application/json' \
  -d '{"fasta": ">seq1\\nALWKTMLKKLGTMALHAGKAALGAAADTISQGTQ", "method": "rf"}'
```

