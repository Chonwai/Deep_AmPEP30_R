# 🎨 機器學習微服務 API 設計最佳實踐

## 📋 概述

本文檔基於 **AmPEP30** 項目的 API 設計經驗，提供機器學習微服務 API 設計的最佳實踐指南。這些原則適用於任何語言和框架的實現。

**核心理念**: 
- **一致性** - 統一的設計模式
- **可預測性** - 直觀的 API 行為
- **可擴展性** - 支持未來需求變化
- **開發者友好** - 易於理解和使用

---

## 🏗️ API 設計原則

### 1. RESTful 設計原則

#### ✅ 推薦的端點設計
```
GET    /health              # 健康檢查
POST   /predict/single      # 單次預測
POST   /predict/batch       # 批量預測
GET    /models              # 模型列表
GET    /models/{id}         # 特定模型信息
POST   /models/{id}/predict # 指定模型預測
GET    /metrics             # 服務指標
GET    /version             # 版本信息
```

#### ❌ 避免的設計
```
GET    /getPrediction       # 動詞命名
POST   /predict_single      # 下劃線分隔
GET    /model_info          # 不一致的命名
POST   /do-prediction       # 混合分隔符
```

### 2. HTTP 狀態碼使用

```python
# 標準狀態碼映射
STATUS_CODES = {
    200: "請求成功",           # 預測成功
    201: "資源已創建",         # 模型上傳成功
    400: "請求參數錯誤",       # 無效輸入
    401: "未授權訪問",         # 缺少 API 密鑰
    403: "禁止訪問",          # 權限不足
    404: "資源不存在",         # 模型不存在
    422: "實體無法處理",       # 數據驗證失敗
    429: "請求過於頻繁",       # 速率限制
    500: "內部服務器錯誤",     # 預測失敗
    503: "服務不可用"          # 模型加載中
}
```

---

## 📊 統一響應格式

### 基礎響應結構

```json
{
  "success": boolean,
  "timestamp": "ISO8601字符串",
  "data": object|null,
  "error": object|null,
  "metadata": object
}
```

### 成功響應示例

```json
{
  "success": true,
  "timestamp": "2024-08-17T10:30:00Z",
  "data": {
    "prediction": 1,
    "confidence": 0.95,
    "model_used": "random_forest_v2",
    "processing_time_ms": 45
  },
  "error": null,
  "metadata": {
    "request_id": "req_abc123",
    "version": "1.0.0",
    "model_version": "2.1.0"
  }
}
```

### 錯誤響應示例

```json
{
  "success": false,
  "timestamp": "2024-08-17T10:30:00Z",
  "data": null,
  "error": {
    "code": "INVALID_SEQUENCE_LENGTH",
    "message": "序列長度必須在 5-30 個氨基酸之間",
    "details": {
      "provided_length": 3,
      "min_length": 5,
      "max_length": 30,
      "provided_sequence": "ABC"
    }
  },
  "metadata": {
    "request_id": "req_def456",
    "version": "1.0.0"
  }
}
```

---

## 🔍 輸入驗證設計

### 1. 分層驗證策略

```python
class ValidationLayer:
    """分層驗證架構"""
    
    def validate_request(self, data):
        """請求級驗證"""
        self._validate_json_structure(data)
        self._validate_required_fields(data)
        self._validate_data_types(data)
    
    def validate_business_rules(self, data):
        """業務規則驗證"""
        self._validate_sequence_format(data.get('sequence'))
        self._validate_model_availability(data.get('method'))
        self._validate_parameter_ranges(data)
    
    def validate_security(self, request):
        """安全性驗證"""
        self._validate_api_key(request.headers.get('Authorization'))
        self._validate_rate_limit(request.remote_addr)
        self._validate_request_size(request.content_length)
```

### 2. 錯誤碼設計

```python
class ErrorCodes:
    """標準化錯誤碼"""
    
    # 請求格式錯誤 (4000-4099)
    INVALID_JSON = "E4001"
    MISSING_REQUIRED_FIELD = "E4002"
    INVALID_DATA_TYPE = "E4003"
    
    # 業務規則錯誤 (4100-4199)
    INVALID_SEQUENCE_LENGTH = "E4101"
    INVALID_SEQUENCE_CHARACTERS = "E4102"
    UNSUPPORTED_MODEL = "E4103"
    INVALID_PARAMETER_RANGE = "E4104"
    
    # 系統錯誤 (5000-5099)
    MODEL_LOADING_ERROR = "E5001"
    PREDICTION_ERROR = "E5002"
    INTERNAL_ERROR = "E5003"
    
    # 資源錯誤 (5100-5199)
    MODEL_NOT_FOUND = "E5101"
    INSUFFICIENT_MEMORY = "E5102"
    TIMEOUT_ERROR = "E5103"
```

### 3. 驗證規則配置

```yaml
# validation_rules.yaml
sequence_validation:
  min_length: 5
  max_length: 30
  allowed_characters: "ACDEFGHIKLMNPQRSTVWY"
  required: true

method_validation:
  allowed_values: ["rf", "cnn", "svm"]
  default: "rf"
  required: false

precision_validation:
  min_value: 0
  max_value: 6
  default: 3
  required: false

batch_validation:
  max_sequences: 100
  max_total_length: 10000
```

---

## 📈 批量處理設計

### 1. 批量請求格式

```json
{
  "sequences": [
    {
      "id": "seq_001",
      "name": "Peptide A",
      "sequence": "GLFDIVKKVVGALGSL"
    },
    {
      "id": "seq_002", 
      "name": "Peptide B",
      "sequence": "ALWKTMLKKLGTMALH"
    }
  ],
  "method": "rf",
  "precision": 3,
  "options": {
    "include_features": false,
    "return_probabilities": true
  }
}
```

### 2. 批量響應格式

```json
{
  "success": true,
  "timestamp": "2024-08-17T10:30:00Z",
  "data": {
    "total_processed": 2,
    "successful": 2,
    "failed": 0,
    "results": [
      {
        "id": "seq_001",
        "name": "Peptide A",
        "sequence": "GLFDIVKKVVGALGSL",
        "prediction": 1,
        "confidence": 0.95,
        "status": "success"
      },
      {
        "id": "seq_002",
        "name": "Peptide B", 
        "sequence": "ALWKTMLKKLGTMALH",
        "prediction": 0,
        "confidence": 0.78,
        "status": "success"
      }
    ]
  },
  "metadata": {
    "processing_time_ms": 156,
    "model_used": "random_forest_v2",
    "batch_id": "batch_xyz789"
  }
}
```

### 3. 部分失敗處理

```json
{
  "success": true,
  "timestamp": "2024-08-17T10:30:00Z",
  "data": {
    "total_processed": 3,
    "successful": 2,
    "failed": 1,
    "results": [
      {
        "id": "seq_001",
        "status": "success",
        "prediction": 1,
        "confidence": 0.95
      },
      {
        "id": "seq_002",
        "status": "error",
        "error": {
          "code": "INVALID_SEQUENCE_LENGTH",
          "message": "序列過短"
        }
      },
      {
        "id": "seq_003", 
        "status": "success",
        "prediction": 0,
        "confidence": 0.82
      }
    ]
  }
}
```

---

## 🔐 安全性設計

### 1. API 密鑰認證

```python
class APIKeyAuth:
    def __init__(self):
        self.valid_keys = self._load_api_keys()
    
    def authenticate(self, request):
        auth_header = request.headers.get('Authorization')
        if not auth_header:
            raise UnauthorizedError("缺少 Authorization 頭")
        
        if not auth_header.startswith('Bearer '):
            raise UnauthorizedError("無效的 Authorization 格式")
        
        api_key = auth_header[7:]  # 移除 'Bearer ' 前綴
        if api_key not in self.valid_keys:
            raise UnauthorizedError("無效的 API 密鑰")
        
        return self.valid_keys[api_key]  # 返回用戶信息
```

### 2. 速率限制

```python
from collections import defaultdict
import time

class RateLimiter:
    def __init__(self, max_requests=100, time_window=3600):
        self.max_requests = max_requests
        self.time_window = time_window
        self.requests = defaultdict(list)
    
    def is_allowed(self, client_id):
        now = time.time()
        # 清理過期記錄
        self.requests[client_id] = [
            req_time for req_time in self.requests[client_id]
            if now - req_time < self.time_window
        ]
        
        # 檢查是否超過限制
        if len(self.requests[client_id]) >= self.max_requests:
            return False
        
        # 記錄本次請求
        self.requests[client_id].append(now)
        return True
```

### 3. 輸入清理

```python
import re
import html

class InputSanitizer:
    def __init__(self):
        self.sequence_pattern = re.compile(r'^[ACDEFGHIKLMNPQRSTVWY]+$')
        self.name_pattern = re.compile(r'^[\w\s\-_\.]+$')
    
    def sanitize_sequence(self, sequence):
        """清理氨基酸序列"""
        if not sequence:
            raise ValidationError("序列不能為空")
        
        # 移除空格和換行符
        cleaned = re.sub(r'\s+', '', sequence.upper())
        
        # 驗證字符
        if not self.sequence_pattern.match(cleaned):
            raise ValidationError("序列包含無效字符")
        
        return cleaned
    
    def sanitize_name(self, name):
        """清理序列名稱"""
        if not name:
            return "unnamed_sequence"
        
        # HTML 轉義
        cleaned = html.escape(name.strip())
        
        # 長度限制
        if len(cleaned) > 100:
            cleaned = cleaned[:97] + "..."
        
        # 字符驗證
        if not self.name_pattern.match(cleaned):
            raise ValidationError("名稱包含無效字符")
        
        return cleaned
```

---

## 📊 性能優化設計

### 1. 響應壓縮

```python
from flask import Flask
from flask_compress import Compress

app = Flask(__name__)
Compress(app)  # 自動壓縮響應

# 或手動控制
@app.route('/predict/batch', methods=['POST'])
@compress.compressed()
def predict_batch():
    # 大量數據的響應會被自動壓縮
    return jsonify(large_response_data)
```

### 2. 分頁設計

```python
@app.route('/predictions', methods=['GET'])
def get_predictions():
    page = request.args.get('page', 1, type=int)
    per_page = min(request.args.get('per_page', 20, type=int), 100)
    
    predictions = PredictionHistory.query.paginate(
        page=page, 
        per_page=per_page,
        error_out=False
    )
    
    return jsonify({
        'success': True,
        'data': {
            'predictions': [p.to_dict() for p in predictions.items],
            'pagination': {
                'page': page,
                'per_page': per_page,
                'total': predictions.total,
                'pages': predictions.pages,
                'has_next': predictions.has_next,
                'has_prev': predictions.has_prev,
                'next_url': url_for('get_predictions', page=predictions.next_num) if predictions.has_next else None,
                'prev_url': url_for('get_predictions', page=predictions.prev_num) if predictions.has_prev else None
            }
        }
    })
```

### 3. 緩存策略

```python
from functools import wraps
import hashlib
import json

def cache_prediction(expiry=3600):
    """預測結果緩存裝飾器"""
    def decorator(func):
        @wraps(func)
        def wrapper(*args, **kwargs):
            # 生成緩存鍵
            cache_key = hashlib.md5(
                json.dumps(kwargs, sort_keys=True).encode()
            ).hexdigest()
            
            # 檢查緩存
            cached_result = redis_client.get(f"prediction:{cache_key}")
            if cached_result:
                return json.loads(cached_result)
            
            # 執行預測
            result = func(*args, **kwargs)
            
            # 存儲緩存
            redis_client.setex(
                f"prediction:{cache_key}", 
                expiry, 
                json.dumps(result)
            )
            
            return result
        return wrapper
    return decorator
```

---

## 📚 API 文檔設計

### 1. OpenAPI/Swagger 規範

```yaml
openapi: 3.0.0
info:
  title: ML Prediction API
  description: 機器學習預測微服務 API
  version: 1.0.0
  contact:
    name: API Support
    email: api-support@example.com

servers:
  - url: https://api.example.com/v1
    description: 生產環境
  - url: https://staging-api.example.com/v1
    description: 測試環境

paths:
  /predict/single:
    post:
      summary: 單序列預測
      description: 對單個氨基酸序列進行抗菌胜肽預測
      tags:
        - Prediction
      requestBody:
        required: true
        content:
          application/json:
            schema:
              $ref: '#/components/schemas/SinglePredictionRequest'
            example:
              sequence: "GLFDIVKKVVGALGSL"
              method: "rf"
              precision: 3
      responses:
        '200':
          description: 預測成功
          content:
            application/json:
              schema:
                $ref: '#/components/schemas/PredictionResponse'
        '400':
          description: 請求參數錯誤
          content:
            application/json:
              schema:
                $ref: '#/components/schemas/ErrorResponse'

components:
  schemas:
    SinglePredictionRequest:
      type: object
      required:
        - sequence
      properties:
        sequence:
          type: string
          description: 氨基酸序列
          pattern: '^[ACDEFGHIKLMNPQRSTVWY]{5,30}$'
          example: "GLFDIVKKVVGALGSL"
        method:
          type: string
          description: 預測方法
          enum: [rf, cnn]
          default: rf
        precision:
          type: integer
          description: 數值精度
          minimum: 0
          maximum: 6
          default: 3
```

### 2. 代碼示例

```python
# Python 示例
import requests

def predict_sequence(sequence, api_key, method='rf'):
    """調用 API 進行序列預測"""
    url = "https://api.example.com/v1/predict/single"
    headers = {
        'Authorization': f'Bearer {api_key}',
        'Content-Type': 'application/json'
    }
    data = {
        'sequence': sequence,
        'method': method
    }
    
    response = requests.post(url, json=data, headers=headers)
    
    if response.status_code == 200:
        result = response.json()
        if result['success']:
            return result['data']
        else:
            raise Exception(f"預測失敗: {result['error']['message']}")
    else:
        raise Exception(f"HTTP 錯誤: {response.status_code}")

# 使用示例
try:
    result = predict_sequence("GLFDIVKKVVGALGSL", "your_api_key")
    print(f"預測結果: {result['prediction']}")
    print(f"置信度: {result['confidence']}")
except Exception as e:
    print(f"錯誤: {e}")
```

```javascript
// JavaScript 示例
async function predictSequence(sequence, apiKey, method = 'rf') {
    const response = await fetch('https://api.example.com/v1/predict/single', {
        method: 'POST',
        headers: {
            'Authorization': `Bearer ${apiKey}`,
            'Content-Type': 'application/json'
        },
        body: JSON.stringify({
            sequence: sequence,
            method: method
        })
    });

    const result = await response.json();
    
    if (!response.ok) {
        throw new Error(`HTTP error! status: ${response.status}`);
    }
    
    if (!result.success) {
        throw new Error(`Prediction failed: ${result.error.message}`);
    }
    
    return result.data;
}

// 使用示例
try {
    const result = await predictSequence("GLFDIVKKVVGALGSL", "your_api_key");
    console.log(`預測結果: ${result.prediction}`);
    console.log(`置信度: ${result.confidence}`);
} catch (error) {
    console.error(`錯誤: ${error.message}`);
}
```

---

## 🔍 監控和指標設計

### 1. 業務指標

```python
class MetricsCollector:
    def __init__(self):
        self.metrics = {
            'total_predictions': 0,
            'successful_predictions': 0,
            'failed_predictions': 0,
            'average_response_time': 0,
            'model_usage': defaultdict(int),
            'error_types': defaultdict(int)
        }
    
    def record_prediction(self, success, model_used, response_time, error_type=None):
        self.metrics['total_predictions'] += 1
        
        if success:
            self.metrics['successful_predictions'] += 1
        else:
            self.metrics['failed_predictions'] += 1
            if error_type:
                self.metrics['error_types'][error_type] += 1
        
        self.metrics['model_usage'][model_used] += 1
        
        # 更新平均響應時間
        total_time = (self.metrics['average_response_time'] * 
                     (self.metrics['total_predictions'] - 1) + response_time)
        self.metrics['average_response_time'] = total_time / self.metrics['total_predictions']
```

### 2. 健康檢查端點

```python
@app.route('/health', methods=['GET'])
def health_check():
    """詳細健康檢查"""
    health_status = {
        'status': 'healthy',
        'timestamp': datetime.utcnow().isoformat(),
        'service': {
            'name': 'ML-Prediction-Service',
            'version': '1.0.0',
            'uptime': get_uptime()
        },
        'models': {},
        'dependencies': {}
    }
    
    # 檢查模型狀態
    for model_name, model in prediction_service.models.items():
        try:
            # 簡單預測測試
            test_result = model.predict(['GLFDIVKKVVGALGSL'])
            health_status['models'][model_name] = {
                'status': 'healthy',
                'last_prediction': datetime.utcnow().isoformat()
            }
        except Exception as e:
            health_status['models'][model_name] = {
                'status': 'unhealthy',
                'error': str(e)
            }
            health_status['status'] = 'degraded'
    
    # 檢查依賴服務
    health_status['dependencies']['redis'] = check_redis_connection()
    health_status['dependencies']['database'] = check_database_connection()
    
    status_code = 200 if health_status['status'] == 'healthy' else 503
    return jsonify(health_status), status_code
```

---

## 🎯 版本控制設計

### 1. URL 版本控制

```python
# 推薦方式：URL 路徑版本控制
@app.route('/v1/predict/single', methods=['POST'])
def predict_single_v1():
    # 版本 1 的實現
    pass

@app.route('/v2/predict/single', methods=['POST'])  
def predict_single_v2():
    # 版本 2 的實現，支持新功能
    pass

# 向後兼容處理
@app.route('/predict/single', methods=['POST'])
def predict_single():
    # 重定向到最新版本或默認版本
    return predict_single_v1()
```

### 2. Header 版本控制

```python
@app.route('/predict/single', methods=['POST'])
def predict_single():
    api_version = request.headers.get('API-Version', 'v1')
    
    if api_version == 'v1':
        return predict_single_v1()
    elif api_version == 'v2':
        return predict_single_v2()
    else:
        return jsonify({
            'success': False,
            'error': {
                'code': 'UNSUPPORTED_VERSION',
                'message': f'不支持的 API 版本: {api_version}',
                'supported_versions': ['v1', 'v2']
            }
        }), 400
```

### 3. 版本棄用策略

```python
@app.route('/v1/predict/single', methods=['POST'])
def predict_single_v1():
    # 添加棄用警告
    response = make_response(jsonify(prediction_result))
    response.headers['Warning'] = '299 - "API v1 將於 2024-12-31 棄用，請升級到 v2"'
    response.headers['Sunset'] = 'Tue, 31 Dec 2024 23:59:59 GMT'
    return response
```

---

## 💡 最佳實踐總結

### ✅ 應該做的

1. **統一響應格式**
   - 使用一致的 JSON 結構
   - 包含成功/失敗標識
   - 提供詳細的錯誤信息

2. **完善的錯誤處理**
   - 標準化錯誤碼
   - 多語言錯誤消息
   - 包含調試信息

3. **安全性考慮**
   - API 密鑰認證
   - 輸入驗證和清理
   - 速率限制

4. **性能優化**
   - 響應壓縮
   - 結果緩存
   - 批量處理優化

5. **完善的文檔**
   - OpenAPI 規範
   - 代碼示例
   - 錯誤處理指南

### ❌ 應該避免的

1. **不一致的設計**
   - 混合的命名約定
   - 不同的響應格式
   - 不規則的錯誤處理

2. **安全漏洞**
   - 缺少輸入驗證
   - 暴露敏感信息
   - 沒有訪問控制

3. **性能問題**
   - 同步阻塞操作
   - 缺少緩存機制
   - 低效的批量處理

4. **維護困難**
   - 缺少版本控制
   - 不充分的監控
   - 文檔不完整

---

## 🔧 實現框架對比

### Python 框架

| 框架 | 優勢 | 適用場景 | 學習曲線 |
|------|------|----------|----------|
| **FastAPI** | 自動文檔、高性能、類型檢查 | 現代 API、高並發 | 中等 |
| **Flask** | 輕量、靈活、生態豐富 | 快速原型、中小項目 | 簡單 |
| **Django REST** | 功能完整、ORM 集成 | 企業應用、數據庫密集 | 複雜 |

### R 框架

| 框架 | 優勢 | 適用場景 | 學習曲線 |
|------|------|----------|----------|
| **Plumber** | R 原生、統計友好 | 統計分析 API | 簡單 |
| **OpenCPU** | 安全、擴展性好 | 生產環境 R API | 中等 |

### Node.js 框架

| 框架 | 優勢 | 適用場景 | 學習曲線 |
|------|------|----------|----------|
| **Express** | 成熟、生態豐富 | 通用 API 服務 | 簡單 |
| **Fastify** | 高性能、低開銷 | 高並發場景 | 中等 |
| **NestJS** | 企業級、TypeScript | 大型項目 | 複雜 |

---

**文檔版本**: 1.0.0  
**更新日期**: 2024年8月17日  
**基於項目**: AmPEP30 抗菌胜肽預測服務 API 設計經驗