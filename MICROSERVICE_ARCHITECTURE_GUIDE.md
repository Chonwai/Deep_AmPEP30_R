# 🏗️ 機器學習模型微服務化架構指南

## 📋 概述

本文檔基於 **AmPEP30 抗菌胜肽預測服務** 的微服務化經驗，提供通用的架構設計指南，適用於任何需要從單體應用轉換為微服務的機器學習項目。

**適用場景**：
- 原本需要在系統上保存文件才能進行預測的應用
- 基於文件輸入/輸出的機器學習模型
- 需要容器化部署的預測服務
- 支持多種預測算法的服務

## 🎯 核心目標

### 從這樣 ❌
```
用戶上傳文件 → 保存到服務器 → 調用預測腳本 → 讀取結果文件 → 返回結果
```

### 變成這樣 ✅
```
HTTP 請求 → API 驗證 → 內存處理 → 直接預測 → JSON 響應
```

---

## 🏛️ 整體架構

### 📦 容器化架構
```
┌─────────────────────────────────────────┐
│                Docker Container          │
│  ┌─────────────────────────────────────┐ │
│  │           Web API Layer             │ │
│  │  ┌─────────────┐ ┌─────────────┐   │ │
│  │  │   REST API  │ │ Health Check│   │ │
│  │  │  (plumber)  │ │   Endpoint  │   │ │
│  │  └─────────────┘ └─────────────┘   │ │
│  └─────────────────────────────────────┘ │
│  ┌─────────────────────────────────────┐ │
│  │        Business Logic Layer        │ │
│  │  ┌─────────────┐ ┌─────────────┐   │ │
│  │  │ Input       │ │ Model       │   │ │
│  │  │ Validation  │ │ Prediction  │   │ │
│  │  └─────────────┘ └─────────────┘   │ │
│  └─────────────────────────────────────┘ │
│  ┌─────────────────────────────────────┐ │
│  │         Model Layer                 │ │
│  │  ┌─────────────┐ ┌─────────────┐   │ │
│  │  │ Random      │ │ Deep        │   │ │
│  │  │ Forest      │ │ Learning    │   │ │
│  │  │ Model       │ │ Model       │   │ │
│  │  └─────────────┘ └─────────────┘   │ │
│  └─────────────────────────────────────┘ │
│  ┌─────────────────────────────────────┐ │
│  │      Configuration Layer            │ │
│  │  ┌─────────────┐ ┌─────────────┐   │ │
│  │  │ Environment │ │ Model       │   │ │
│  │  │ Variables   │ │ Paths       │   │ │
│  │  └─────────────┘ └─────────────┘   │ │
│  └─────────────────────────────────────┘ │
└─────────────────────────────────────────┘
```

### 🔄 數據流程
```
┌─────────────┐    ┌─────────────┐    ┌─────────────┐    ┌─────────────┐
│   Client    │───▶│  API Gateway│───▶│  Validation │───▶│  Prediction │
│  Request    │    │  (plumber)  │    │   Layer     │    │   Engine    │
└─────────────┘    └─────────────┘    └─────────────┘    └─────────────┘
                            │                                     │
                            ▼                                     ▼
                   ┌─────────────┐                      ┌─────────────┐
                   │   Error     │                      │   Success   │
                   │  Response   │                      │  Response   │
                   └─────────────┘                      └─────────────┘
```

---

## 📁 目錄結構設計

### 標準微服務目錄結構
```
project-root/
├── microservice/                    # 微服務相關文件
│   ├── api/                        # API 層
│   │   ├── main_api.R              # 主 API 文件
│   │   └── endpoints/              # 端點定義
│   ├── config/                     # 配置管理
│   │   ├── config.R                # 配置文件
│   │   └── validation.R            # 驗證規則
│   ├── docker/                     # 容器化配置
│   │   ├── Dockerfile              # 容器構建文件
│   │   └── docker-compose.yml      # 編排配置
│   ├── scripts/                    # 啟動腳本
│   │   └── start.sh                # 服務啟動腳本
│   └── docs/                       # 文檔
│       └── API_DOCUMENTATION.md    # API 文檔
├── models/                         # 模型文件
│   ├── model_v1.pkl                # 訓練好的模型
│   └── model_v2.h5                 # 深度學習模型
├── core/                           # 核心業務邏輯
│   ├── prediction.R                # 預測邏輯
│   └── preprocessing.R             # 數據預處理
└── tests/                          # 測試文件
    ├── test_api.R                  # API 測試
    └── test_models.R               # 模型測試
```

---

## 🔧 核心組件設計

### 1. 🌐 API 層 (Web Layer)
**職責**: 處理 HTTP 請求和響應

**關鍵特性**:
- RESTful API 設計
- 統一的響應格式
- 錯誤處理機制
- 請求驗證

**技術選擇**:
- **R**: Plumber 框架
- **Python**: FastAPI, Flask
- **Node.js**: Express
- **Java**: Spring Boot

### 2. 🛡️ 驗證層 (Validation Layer)
**職責**: 輸入數據驗證和清理

**驗證項目**:
- 數據格式檢查
- 參數範圍驗證
- 業務規則檢查
- 安全性驗證

### 3. 🧠 預測層 (Prediction Layer)
**職責**: 執行機器學習預測

**設計原則**:
- 模型抽象化
- 支持多模型切換
- 內存中處理
- 結果標準化

### 4. ⚙️ 配置層 (Configuration Layer)
**職責**: 管理服務配置和環境變量

**配置項目**:
- 服務端口和主機
- 模型路徑配置
- 驗證規則參數
- 日誌級別設置

---

## 🐳 容器化設計

### Dockerfile 最佳實踐

```dockerfile
# 1. 選擇合適的基礎鏡像
FROM python:3.9-slim  # 或 rocker/r-ver:4.3.2

# 2. 設置環境變量
ENV APP_ROOT=/app \
    API_PORT=8000 \
    PYTHONUNBUFFERED=1

# 3. 安裝系統依賴
RUN apt-get update && apt-get install -y \
    build-essential \
    curl \
    && rm -rf /var/lib/apt/lists/*

# 4. 設置工作目錄
WORKDIR /app

# 5. 複製並安裝依賴
COPY requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt

# 6. 複製應用代碼
COPY . .

# 7. 暴露端口
EXPOSE 8000

# 8. 健康檢查
HEALTHCHECK --interval=30s --timeout=10s --start-period=60s --retries=3 \
  CMD curl -f http://localhost:8000/health || exit 1

# 9. 啟動命令
CMD ["python", "app.py"]
```

### Docker Compose 編排

```yaml
version: '3.8'
services:
  ml-prediction-service:
    build:
      context: .
      dockerfile: Dockerfile
    container_name: ml-service
    environment:
      - APP_ROOT=/app
      - API_PORT=8000
      - MODEL_PATH=/app/models
    ports:
      - "8000:8000"
    volumes:
      - ./models:/app/models:ro  # 只讀掛載模型
    healthcheck:
      test: ["CMD", "curl", "-f", "http://localhost:8000/health"]
      interval: 30s
      timeout: 10s
      retries: 3
    restart: unless-stopped
```

---

## 📊 API 設計模式

### 統一響應格式

```json
{
  "status": "success|error",
  "data": {
    "prediction": 1,
    "confidence": 0.95,
    "model_used": "random_forest"
  },
  "error": null,
  "metadata": {
    "timestamp": "2024-08-17T10:30:00Z",
    "request_id": "req_123456"
  }
}
```

### 標準端點設計

| 端點 | 方法 | 功能 | 用途 |
|------|------|------|------|
| `/health` | GET | 健康檢查 | 容器編排監控 |
| `/predict` | POST | 單次預測 | 主要功能 |
| `/predict/batch` | POST | 批量預測 | 效率優化 |
| `/models` | GET | 模型信息 | 運維監控 |
| `/metrics` | GET | 服務指標 | 性能監控 |

---

## 🔐 安全性考慮

### 1. 輸入驗證
- 參數類型檢查
- 數據長度限制
- 特殊字符過濾
- SQL 注入防護

### 2. 訪問控制
- API 密鑰驗證
- 速率限制
- IP 白名單
- CORS 配置

### 3. 數據保護
- 敏感信息脫敏
- 日誌安全
- 傳輸加密
- 存儲加密

---

## 📈 監控和日誌

### 監控指標
- **服務健康**: 響應時間、錯誤率
- **業務指標**: 預測準確率、請求量
- **資源使用**: CPU、內存、磁盤
- **模型性能**: 推理時間、準確度

### 日誌設計
```
[TIMESTAMP] [LEVEL] [COMPONENT] [REQUEST_ID] MESSAGE
2024-08-17 10:30:00 INFO API req_123 Prediction request received
2024-08-17 10:30:01 DEBUG MODEL req_123 Using random_forest model
2024-08-17 10:30:02 INFO API req_123 Prediction completed: confidence=0.95
```

---

## 🚀 部署策略

### 1. 本地開發
```bash
# 開發環境啟動
docker-compose up -d
curl http://localhost:8000/health
```

### 2. 測試環境
```bash
# 自動化測試
docker-compose -f docker-compose.test.yml up --abort-on-container-exit
```

### 3. 生產環境
```bash
# 生產部署
docker-compose -f docker-compose.prod.yml up -d
# 滾動更新
docker-compose -f docker-compose.prod.yml up -d --no-deps service-name
```

---

## 🔄 版本管理

### API 版本控制
```
/v1/predict  # 版本 1
/v2/predict  # 版本 2
```

### 模型版本管理
```
models/
├── v1/
│   └── model.pkl
├── v2/
│   └── model.pkl
└── current -> v2/  # 符號連結指向當前版本
```

---

## 💡 最佳實踐總結

### ✅ 做這些
1. **統一響應格式** - 便於客戶端集成
2. **完整的錯誤處理** - 包含詳細上下文
3. **健康檢查端點** - 支持容器編排
4. **配置外部化** - 使用環境變量
5. **內存處理** - 避免文件 I/O
6. **模型抽象** - 支持多模型切換
7. **詳細日誌** - 便於問題診斷

### ❌ 避免這些
1. **硬編碼配置** - 降低靈活性
2. **文件依賴** - 增加部署複雜性
3. **同步阻塞** - 影響並發性能
4. **單一錯誤格式** - 客戶端處理困難
5. **缺少監控** - 運維困難
6. **不一致的 API** - 增加集成成本

---

## 🎓 技術棧選擇指南

### 語言和框架對比

| 語言 | 框架 | 優勢 | 適用場景 |
|------|------|------|----------|
| **Python** | FastAPI | 高性能、自動文檔 | 深度學習、數據科學 |
| **Python** | Flask | 輕量、靈活 | 簡單 API、快速原型 |
| **R** | Plumber | 統計建模友好 | 統計分析、生物信息 |
| **Node.js** | Express | 高並發、JavaScript 生態 | 實時應用、前後端統一 |
| **Java** | Spring Boot | 企業級、生態完整 | 大型系統、企業應用 |

---

**文檔版本**: 1.0.0  
**更新日期**: 2024年8月17日  
**基於項目**: AmPEP30 抗菌胜肽預測服務