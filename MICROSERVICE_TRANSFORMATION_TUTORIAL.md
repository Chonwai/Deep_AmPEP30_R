# 🔄 機器學習應用微服務化改造教學指南

## 📋 概述

本指南基於 **AmPEP30** 項目的實際改造經驗，提供完整的微服務化改造教學，幫助團隊將基於文件處理的機器學習應用轉換為現代化的微服務架構。

**適用對象**:
- 需要微服務化改造的機器學習團隊
- 從文件處理模式轉向 API 模式的項目
- 希望容器化部署的應用

---

## 🎯 改造目標

### 改造前 ❌ 
```
用戶操作：上傳文件到服務器
系統處理：保存文件 → 調用腳本 → 讀取結果文件
用戶獲取：下載結果文件
```

**問題**:
- 文件 I/O 開銷大
- 部署環境依賴複雜
- 難以水平擴展
- 缺乏標準化接口

### 改造後 ✅
```
用戶操作：發送 HTTP 請求
系統處理：接收數據 → 內存處理 → 直接預測
用戶獲取：實時 JSON 響應
```

**優勢**:
- 高性能內存處理
- 容器化標準部署
- 水平擴展友好
- RESTful API 標準

---

## 🚀 改造步驟

## 步驟 1: 📊 現狀分析

### 1.1 識別核心功能
```python
# 原始應用分析清單
功能分析：
□ 主要預測算法（如：隨機森林、深度學習）
□ 輸入數據格式（如：FASTA、CSV、TXT）
□ 輸出結果格式（如：分類結果、概率值）
□ 預處理步驟（如：序列清理、特徵提取）
□ 後處理邏輯（如：結果解釋、置信度計算）

依賴分析：
□ 機器學習庫（如：scikit-learn、tensorflow、keras）
□ 數據處理庫（如：pandas、numpy、seqinr）
□ 文件處理庫（如：Bio.SeqIO、readr）
□ 系統依賴（如：特定版本的 Python/R）
```

### 1.2 性能基準測試
```bash
# 記錄原始性能指標
echo "測試原始應用性能..."
time python original_predict.py input.fasta output.txt
# 記錄：處理時間、內存使用、文件大小
```

---

## 步驟 2: 🏗️ 架構設計

### 2.1 API 端點設計

```yaml
# API 設計規劃
端點規劃:
  健康檢查: GET /health
  單次預測: POST /predict/single
  批量預測: POST /predict/batch  
  模型信息: GET /model/info
  服務指標: GET /metrics

數據流設計:
  輸入: JSON 格式的結構化數據
  處理: 內存中的數據流處理
  輸出: 統一的 JSON 響應格式
```

### 2.2 目錄結構規劃

```
project/
├── microservice/           # 新增微服務目錄
│   ├── api/               # API 層
│   ├── config/            # 配置管理
│   ├── docker/            # 容器化配置
│   └── scripts/           # 啟動腳本
├── core/                  # 重構後的核心邏輯
│   ├── models/            # 模型加載和預測
│   ├── preprocessing/     # 數據預處理
│   └── validation/        # 輸入驗證
├── original/              # 保留原始代碼（備份）
└── tests/                 # 測試代碼
```

---

## 步驟 3: 🔧 核心邏輯重構

### 3.1 模型加載優化

**原始方式** ❌:
```python
# 每次預測都重新加載模型
def predict(input_file):
    model = load_model('model.pkl')  # 每次都加載
    data = read_file(input_file)
    result = model.predict(data)
    return result
```

**微服務方式** ✅:
```python
# 服務啟動時加載，內存中保持
class PredictionService:
    def __init__(self):
        self.models = {
            'rf': load_model('rf_model.pkl'),
            'cnn': load_model('cnn_model.h5')
        }
    
    def predict(self, data, model_type='rf'):
        model = self.models[model_type]
        return model.predict(data)

# 全局實例，服務啟動時初始化
prediction_service = PredictionService()
```

### 3.2 數據處理流水線

**原始方式** ❌:
```python
def process_fasta_file(file_path):
    # 讀取文件
    sequences = []
    with open(file_path, 'r') as f:
        # 解析 FASTA 格式
        ...
    return sequences
```

**微服務方式** ✅:
```python
def process_fasta_content(fasta_string):
    """直接處理 FASTA 內容，無需文件 I/O"""
    sequences = []
    lines = fasta_string.strip().split('\n')
    current_seq = ""
    current_name = ""
    
    for line in lines:
        if line.startswith('>'):
            if current_seq:
                sequences.append({
                    'name': current_name,
                    'sequence': current_seq
                })
            current_name = line[1:]
            current_seq = ""
        else:
            current_seq += line.strip()
    
    if current_seq:
        sequences.append({
            'name': current_name,
            'sequence': current_seq
        })
    
    return sequences
```

### 3.3 錯誤處理標準化

```python
class APIError(Exception):
    def __init__(self, message, error_code, status_code=400):
        self.message = message
        self.error_code = error_code
        self.status_code = status_code

def validate_sequence(sequence):
    """統一的序列驗證邏輯"""
    if not sequence:
        raise APIError("序列不能為空", "EMPTY_SEQUENCE")
    
    if len(sequence) < 5 or len(sequence) > 30:
        raise APIError(
            f"序列長度必須在 5-30 之間，當前長度: {len(sequence)}", 
            "INVALID_LENGTH"
        )
    
    valid_chars = set('ACDEFGHIKLMNPQRSTVWY')
    invalid_chars = set(sequence.upper()) - valid_chars
    if invalid_chars:
        raise APIError(
            f"序列包含無效字符: {', '.join(invalid_chars)}", 
            "INVALID_CHARACTERS"
        )
```

---

## 步驟 4: 🌐 API 層實現

### 4.1 統一響應格式

```python
def create_response(success=True, data=None, error=None, metadata=None):
    """統一的響應格式生成器"""
    response = {
        "success": success,
        "timestamp": datetime.utcnow().isoformat(),
        "data": data,
        "error": error,
        "metadata": metadata or {}
    }
    return response

# 成功響應示例
def success_response(prediction_result):
    return create_response(
        success=True,
        data={
            "prediction": prediction_result["prediction"],
            "confidence": prediction_result["confidence"],
            "model_used": prediction_result["model"]
        },
        metadata={
            "processing_time": prediction_result["time"],
            "version": "1.0.0"
        }
    )

# 錯誤響應示例
def error_response(error):
    return create_response(
        success=False,
        error={
            "message": error.message,
            "code": error.error_code,
            "details": getattr(error, 'details', None)
        }
    )
```

### 4.2 API 端點實現

```python
from flask import Flask, request, jsonify

app = Flask(__name__)

@app.route('/health', methods=['GET'])
def health_check():
    """健康檢查端點"""
    return jsonify({
        "status": "healthy",
        "service": "ML-Prediction-Service",
        "version": "1.0.0",
        "timestamp": datetime.utcnow().isoformat()
    })

@app.route('/predict/single', methods=['POST'])
def predict_single():
    """單序列預測端點"""
    try:
        # 解析請求
        data = request.get_json()
        sequence = data.get('sequence')
        method = data.get('method', 'rf')
        
        # 驗證輸入
        validate_sequence(sequence)
        
        # 執行預測
        result = prediction_service.predict(sequence, method)
        
        # 返回結果
        return jsonify(success_response(result))
        
    except APIError as e:
        return jsonify(error_response(e)), e.status_code
    except Exception as e:
        return jsonify(error_response(
            APIError("內部服務錯誤", "INTERNAL_ERROR", 500)
        )), 500

@app.route('/predict/batch', methods=['POST'])
def predict_batch():
    """批量預測端點"""
    try:
        data = request.get_json()
        fasta_content = data.get('fasta_content')
        method = data.get('method', 'rf')
        
        # 解析 FASTA 內容
        sequences = process_fasta_content(fasta_content)
        
        # 批量預測
        results = []
        for seq_info in sequences:
            try:
                validate_sequence(seq_info['sequence'])
                result = prediction_service.predict(
                    seq_info['sequence'], 
                    method
                )
                result['sequence_name'] = seq_info['name']
                results.append(result)
            except APIError as e:
                # 記錄錯誤但繼續處理其他序列
                results.append({
                    'sequence_name': seq_info['name'],
                    'error': e.message
                })
        
        return jsonify(success_response({'predictions': results}))
        
    except Exception as e:
        return jsonify(error_response(
            APIError("批量預測失敗", "BATCH_PREDICTION_ERROR", 500)
        )), 500
```

---

## 步驟 5: 🐳 容器化配置

### 5.1 Dockerfile 編寫

```dockerfile
# 基於官方 Python 鏡像
FROM python:3.9-slim

# 設置環境變量
ENV PYTHONUNBUFFERED=1 \
    APP_ROOT=/app \
    API_PORT=8000

# 安裝系統依賴
RUN apt-get update && apt-get install -y \
    build-essential \
    curl \
    && rm -rf /var/lib/apt/lists/*

# 設置工作目錄
WORKDIR /app

# 複製依賴文件並安裝
COPY requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt

# 複製應用代碼
COPY . .

# 創建非 root 用戶
RUN useradd --create-home --shell /bin/bash app && \
    chown -R app:app /app
USER app

# 暴露端口
EXPOSE 8000

# 健康檢查
HEALTHCHECK --interval=30s --timeout=10s --start-period=60s --retries=3 \
  CMD curl -f http://localhost:8000/health || exit 1

# 啟動命令
CMD ["python", "app.py"]
```

### 5.2 Docker Compose 配置

```yaml
version: '3.8'

services:
  ml-service:
    build:
      context: .
      dockerfile: Dockerfile
    container_name: ml-prediction-service
    environment:
      - APP_ROOT=/app
      - API_PORT=8000
      - MODEL_PATH=/app/models
      - LOG_LEVEL=INFO
    ports:
      - "8000:8000"
    volumes:
      - ./models:/app/models:ro
      - ./logs:/app/logs
    healthcheck:
      test: ["CMD", "curl", "-f", "http://localhost:8000/health"]
      interval: 30s
      timeout: 10s
      retries: 3
      start_period: 60s
    restart: unless-stopped
    networks:
      - ml-network

networks:
  ml-network:
    driver: bridge
```

---

## 步驟 6: ⚙️ 配置管理

### 6.1 環境變量配置

```python
import os
from dataclasses import dataclass

@dataclass
class Config:
    # 服務配置
    SERVICE_NAME: str = os.getenv('SERVICE_NAME', 'ML-Prediction-Service')
    SERVICE_VERSION: str = os.getenv('SERVICE_VERSION', '1.0.0')
    API_HOST: str = os.getenv('API_HOST', '0.0.0.0')
    API_PORT: int = int(os.getenv('API_PORT', '8000'))
    
    # 模型配置
    MODEL_PATH: str = os.getenv('MODEL_PATH', '/app/models')
    DEFAULT_MODEL: str = os.getenv('DEFAULT_MODEL', 'random_forest')
    
    # 驗證配置
    MIN_SEQUENCE_LENGTH: int = int(os.getenv('MIN_SEQUENCE_LENGTH', '5'))
    MAX_SEQUENCE_LENGTH: int = int(os.getenv('MAX_SEQUENCE_LENGTH', '30'))
    MAX_BATCH_SIZE: int = int(os.getenv('MAX_BATCH_SIZE', '100'))
    
    # 日誌配置
    LOG_LEVEL: str = os.getenv('LOG_LEVEL', 'INFO')
    LOG_FORMAT: str = os.getenv('LOG_FORMAT', 
        '%(asctime)s - %(name)s - %(levelname)s - %(message)s')

# 全局配置實例
config = Config()
```

### 6.2 多環境配置

```bash
# .env.development
SERVICE_NAME=ML-Service-Dev
API_PORT=8000
LOG_LEVEL=DEBUG
MODEL_PATH=./models

# .env.production  
SERVICE_NAME=ML-Service-Prod
API_PORT=8000
LOG_LEVEL=INFO
MODEL_PATH=/app/models
```

---

## 步驟 7: 🧪 測試策略

### 7.1 單元測試

```python
import unittest
from unittest.mock import patch, MagicMock

class TestPredictionService(unittest.TestCase):
    
    def setUp(self):
        self.service = PredictionService()
    
    def test_valid_sequence_prediction(self):
        """測試有效序列預測"""
        sequence = "GLFDIVKKVVGALGSL"
        result = self.service.predict(sequence, 'rf')
        
        self.assertIn('prediction', result)
        self.assertIn('confidence', result)
        self.assertIsInstance(result['confidence'], float)
        self.assertTrue(0 <= result['confidence'] <= 1)
    
    def test_invalid_sequence_length(self):
        """測試無效序列長度"""
        with self.assertRaises(APIError):
            validate_sequence("ABC")  # 太短
        
        with self.assertRaises(APIError):
            validate_sequence("A" * 50)  # 太長
    
    def test_invalid_characters(self):
        """測試無效字符"""
        with self.assertRaises(APIError):
            validate_sequence("GLFDIVXKVVGALGSL")  # 包含 X
```

### 7.2 集成測試

```python
import requests
import json

class TestAPIEndpoints(unittest.TestCase):
    
    def setUp(self):
        self.base_url = "http://localhost:8000"
    
    def test_health_endpoint(self):
        """測試健康檢查端點"""
        response = requests.get(f"{self.base_url}/health")
        self.assertEqual(response.status_code, 200)
        data = response.json()
        self.assertEqual(data['status'], 'healthy')
    
    def test_single_prediction(self):
        """測試單序列預測"""
        payload = {
            "sequence": "GLFDIVKKVVGALGSL",
            "method": "rf"
        }
        response = requests.post(
            f"{self.base_url}/predict/single",
            json=payload
        )
        self.assertEqual(response.status_code, 200)
        data = response.json()
        self.assertTrue(data['success'])
        self.assertIn('prediction', data['data'])
```

### 7.3 性能測試

```python
import time
import concurrent.futures

def performance_test():
    """性能基準測試"""
    test_sequence = "GLFDIVKKVVGALGSL"
    num_requests = 100
    
    def single_request():
        start_time = time.time()
        response = requests.post(
            "http://localhost:8000/predict/single",
            json={"sequence": test_sequence}
        )
        end_time = time.time()
        return end_time - start_time, response.status_code
    
    # 並發測試
    with concurrent.futures.ThreadPoolExecutor(max_workers=10) as executor:
        futures = [executor.submit(single_request) for _ in range(num_requests)]
        results = [future.result() for future in futures]
    
    times = [r[0] for r in results]
    success_count = sum(1 for r in results if r[1] == 200)
    
    print(f"成功請求: {success_count}/{num_requests}")
    print(f"平均響應時間: {sum(times)/len(times):.3f}s")
    print(f"最大響應時間: {max(times):.3f}s")
    print(f"最小響應時間: {min(times):.3f}s")
```

---

## 步驟 8: 📊 監控和日誌

### 8.1 結構化日誌

```python
import logging
import json
from datetime import datetime

class StructuredLogger:
    def __init__(self, name):
        self.logger = logging.getLogger(name)
        self.logger.setLevel(logging.INFO)
        
        handler = logging.StreamHandler()
        formatter = logging.Formatter(
            '%(asctime)s - %(name)s - %(levelname)s - %(message)s'
        )
        handler.setFormatter(formatter)
        self.logger.addHandler(handler)
    
    def log_request(self, request_id, endpoint, method, data=None):
        self.logger.info(json.dumps({
            "event": "request_received",
            "request_id": request_id,
            "endpoint": endpoint,
            "method": method,
            "timestamp": datetime.utcnow().isoformat(),
            "data_size": len(str(data)) if data else 0
        }))
    
    def log_prediction(self, request_id, model_used, processing_time, success=True):
        self.logger.info(json.dumps({
            "event": "prediction_completed",
            "request_id": request_id,
            "model_used": model_used,
            "processing_time": processing_time,
            "success": success,
            "timestamp": datetime.utcnow().isoformat()
        }))

# 全局日誌實例
logger = StructuredLogger(__name__)
```

### 8.2 性能監控

```python
from functools import wraps
import time

def monitor_performance(func):
    """性能監控裝飾器"""
    @wraps(func)
    def wrapper(*args, **kwargs):
        start_time = time.time()
        try:
            result = func(*args, **kwargs)
            success = True
            return result
        except Exception as e:
            success = False
            raise
        finally:
            end_time = time.time()
            processing_time = end_time - start_time
            
            # 記錄性能指標
            logger.log_prediction(
                request_id=kwargs.get('request_id', 'unknown'),
                model_used=kwargs.get('model_type', 'unknown'),
                processing_time=processing_time,
                success=success
            )
    
    return wrapper

# 使用示例
@monitor_performance
def predict_with_monitoring(sequence, model_type='rf', request_id=None):
    return prediction_service.predict(sequence, model_type)
```

---

## 步驟 9: 🚀 部署和運維

### 9.1 部署腳本

```bash
#!/bin/bash
# deploy.sh

set -e

echo "開始部署 ML 微服務..."

# 1. 構建鏡像
echo "構建 Docker 鏡像..."
docker build -t ml-service:latest .

# 2. 停止舊服務
echo "停止現有服務..."
docker-compose down

# 3. 啟動新服務
echo "啟動新服務..."
docker-compose up -d

# 4. 健康檢查
echo "等待服務啟動..."
sleep 30

# 5. 驗證部署
echo "驗證服務健康..."
curl -f http://localhost:8000/health || {
    echo "健康檢查失敗，回滾部署"
    docker-compose down
    exit 1
}

echo "部署成功！"
```

### 9.2 監控腳本

```bash
#!/bin/bash
# monitor.sh

# 檢查容器狀態
check_container() {
    if docker ps | grep -q ml-service; then
        echo "✅ 容器運行正常"
    else
        echo "❌ 容器未運行"
        return 1
    fi
}

# 檢查 API 健康
check_api() {
    if curl -s -f http://localhost:8000/health > /dev/null; then
        echo "✅ API 健康檢查通過"
    else
        echo "❌ API 健康檢查失敗"
        return 1
    fi
}

# 檢查資源使用
check_resources() {
    echo "📊 資源使用情況："
    docker stats ml-service --no-stream --format "table {{.CPUPerc}}\t{{.MemUsage}}\t{{.MemPerc}}"
}

# 主監控邏輯
main() {
    echo "🔍 開始監控檢查..."
    
    check_container || exit 1
    check_api || exit 1
    check_resources
    
    echo "✅ 所有檢查通過"
}

main
```

---

## 步驟 10: 📈 性能優化

### 10.1 模型預加載優化

```python
class OptimizedPredictionService:
    def __init__(self):
        # 預加載所有模型到內存
        self.models = {}
        self.model_metadata = {}
        self._load_all_models()
    
    def _load_all_models(self):
        """啟動時預加載所有模型"""
        model_configs = [
            {'name': 'rf', 'path': 'models/rf_model.pkl'},
            {'name': 'cnn', 'path': 'models/cnn_model.h5'}
        ]
        
        for config in model_configs:
            start_time = time.time()
            self.models[config['name']] = load_model(config['path'])
            load_time = time.time() - start_time
            
            self.model_metadata[config['name']] = {
                'load_time': load_time,
                'model_size': os.path.getsize(config['path']),
                'loaded_at': datetime.utcnow().isoformat()
            }
            
            logger.info(f"模型 {config['name']} 加載完成，耗時 {load_time:.2f}s")
```

### 10.2 批處理優化

```python
def batch_predict_optimized(sequences, model_type='rf', batch_size=32):
    """優化的批處理預測"""
    model = prediction_service.models[model_type]
    results = []
    
    # 分批處理，避免內存溢出
    for i in range(0, len(sequences), batch_size):
        batch = sequences[i:i + batch_size]
        
        # 批量特徵提取
        features = extract_features_batch(batch)
        
        # 批量預測
        predictions = model.predict(features)
        
        # 處理結果
        for j, (seq_info, pred) in enumerate(zip(batch, predictions)):
            results.append({
                'sequence_name': seq_info['name'],
                'sequence': seq_info['sequence'],
                'prediction': int(pred),
                'confidence': float(max(pred)) if hasattr(pred, '__iter__') else float(pred)
            })
    
    return results
```

### 10.3 緩存策略

```python
from functools import lru_cache
import hashlib

class CachedPredictionService:
    def __init__(self):
        self.models = {}
        self._load_models()
    
    def _get_sequence_hash(self, sequence, model_type):
        """生成序列和模型的哈希值用於緩存"""
        content = f"{sequence}:{model_type}"
        return hashlib.md5(content.encode()).hexdigest()
    
    @lru_cache(maxsize=1000)
    def predict_cached(self, sequence_hash, sequence, model_type):
        """帶緩存的預測方法"""
        return self._predict_internal(sequence, model_type)
    
    def predict(self, sequence, model_type='rf'):
        """對外的預測接口"""
        sequence_hash = self._get_sequence_hash(sequence, model_type)
        return self.predict_cached(sequence_hash, sequence, model_type)
```

---

## 📋 改造檢查清單

### ✅ 技術改造檢查
- [ ] **核心邏輯重構**
  - [ ] 模型加載邏輯改為啟動時預加載
  - [ ] 文件處理改為內存處理
  - [ ] 錯誤處理標準化
  
- [ ] **API 層實現**
  - [ ] RESTful API 設計
  - [ ] 統一響應格式
  - [ ] 輸入驗證機制
  - [ ] 錯誤處理中間件
  
- [ ] **容器化配置**
  - [ ] Dockerfile 編寫
  - [ ] Docker Compose 配置
  - [ ] 環境變量管理
  - [ ] 健康檢查配置

### ✅ 質量保證檢查
- [ ] **測試覆蓋**
  - [ ] 單元測試 (>80% 覆蓋率)
  - [ ] 集成測試
  - [ ] 性能測試
  - [ ] 負載測試
  
- [ ] **監控日誌**
  - [ ] 結構化日誌
  - [ ] 性能監控
  - [ ] 錯誤追蹤
  - [ ] 業務指標監控

### ✅ 運維部署檢查
- [ ] **部署流程**
  - [ ] 自動化部署腳本
  - [ ] 藍綠部署支持
  - [ ] 回滾機制
  - [ ] 健康檢查驗證
  
- [ ] **文檔完整性**
  - [ ] API 文檔
  - [ ] 部署文檔
  - [ ] 運維手冊
  - [ ] 故障排除指南

---

## 🎯 成功指標

### 性能指標
- **響應時間**: < 100ms (單次預測)
- **吞吐量**: > 1000 請求/分鐘
- **可用性**: > 99.9%
- **錯誤率**: < 0.1%

### 開發效率指標
- **部署時間**: < 5 分鐘
- **新功能開發週期**: 縮短 50%
- **Bug 修復時間**: 縮短 60%
- **環境一致性**: 100%

### 運維指標
- **自動化程度**: > 90%
- **監控覆蓋**: 100%
- **故障恢復時間**: < 5 分鐘
- **擴展能力**: 支持水平擴展

---

## 💡 常見問題和解決方案

### Q1: 模型加載時間過長
**問題**: 服務啟動時間過長，影響部署效率
**解決方案**:
```python
# 使用模型預熱和分階段加載
def warm_up_models():
    """模型預熱"""
    dummy_input = "GLFDIVKKVVGALGSL"
    for model_name in prediction_service.models:
        prediction_service.predict(dummy_input, model_name)
```

### Q2: 內存使用過高
**問題**: 多個大模型同時加載導致內存不足
**解決方案**:
```python
# 實現模型懶加載和緩存淘汰
class LazyModelService:
    def __init__(self, max_models_in_memory=2):
        self.max_models = max_models_in_memory
        self.loaded_models = {}
        self.model_usage = {}
    
    def get_model(self, model_type):
        if model_type not in self.loaded_models:
            if len(self.loaded_models) >= self.max_models:
                # 淘汰最少使用的模型
                least_used = min(self.model_usage, key=self.model_usage.get)
                del self.loaded_models[least_used]
                del self.model_usage[least_used]
            
            self.loaded_models[model_type] = load_model(f'{model_type}_model.pkl')
        
        self.model_usage[model_type] = self.model_usage.get(model_type, 0) + 1
        return self.loaded_models[model_type]
```

### Q3: 批量處理性能瓶頸
**問題**: 大批量數據處理時性能下降明顯
**解決方案**:
```python
# 實現流式處理和異步處理
import asyncio
from concurrent.futures import ThreadPoolExecutor

async def batch_predict_async(sequences, model_type='rf'):
    """異步批量預測"""
    loop = asyncio.get_event_loop()
    executor = ThreadPoolExecutor(max_workers=4)
    
    tasks = []
    for seq_batch in chunk_sequences(sequences, batch_size=100):
        task = loop.run_in_executor(
            executor, 
            prediction_service.predict_batch, 
            seq_batch, 
            model_type
        )
        tasks.append(task)
    
    results = await asyncio.gather(*tasks)
    return flatten_results(results)
```

---

**文檔版本**: 1.0.0  
**更新日期**: 2024年8月17日  
**基於項目**: AmPEP30 抗菌胜肽預測服務微服務化改造