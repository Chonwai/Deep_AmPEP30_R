# Deep-AmPEP30: 基於深度學習的抗菌肽預測系統

[English](README.md) | [繁體中文](README_zh-TW.md)

![R](https://img.shields.io/badge/R-276DC3.svg?style=flat&logo=r&logoColor=white)
![Keras](https://img.shields.io/badge/Keras-FF0000.svg?style=flat&logo=keras&logoColor=white)
![TensorFlow](https://img.shields.io/badge/TensorFlow-FF6F00.svg?style=flat&logo=tensorflow&logoColor=white)

**Deep-AmPEP30** 是一個先進的機器學習系統，專門用於從蛋白質序列中預測抗菌肽（AMPs）。此系統作為 AxPEP 網站伺服器的核心預測引擎，為長度在 5 到 30 個氨基酸之間的肽段提供準確且高效的抗菌肽預測。

## 🌟 核心特色

- **雙模型架構**：實現卷積神經網路（CNN）和隨機森林（RF）模型，提供穩健的預測結果
- **先進特徵工程**：運用 PseKRAAC（偽K元組簡化氨基酸組成）和多種氨基酸分類方案
- **高效能表現**：針對大規模序列處理進行最佳化，每序列預測時間低於1秒
- **生產就緒**：專為實際應用設計，具備完整的錯誤處理和驗證機制
- **批次處理**：支援 FASTA 格式輸入，可同時處理多個序列

## 🛠️ 技術架構

### 機器學習模型

#### CNN 模型（卷積神經網路）
- **框架**：Keras/TensorFlow
- **架構**：128+128 個過濾器，3×3 卷積核，10 個密集單元
- **訓練時間**：約3分鐘（3,297個訓練序列）
- **預測速度**：約22秒（10,000個序列）

#### 隨機森林模型
- **配置**：1,200 棵決策樹
- **訓練時間**：約27秒
- **預測速度**：約23秒（10,000個序列）

### 特徵提取
- **氨基酸組成（AAC）**：20種標準氨基酸的頻率分析
- **PseKRAAC**：偽K元組簡化氨基酸組成
- **多重分類系統**：Type 8-17、Type 3A-19、Type 12-18、Type 7-15、Type 12-17
- **最佳化特徵選擇**：最佳5種特徵組合策略

## 📊 效能指標

### 模型準確度
- **UniProt 非抗菌肽資料集**：75.4% 準確度（24.6% 錯誤率）
- **正樣本錯誤率**：5.9%
- **大規模預測**：38% 序列被預測為抗菌肽

### 計算效能
- **特徵生成**：0.52±0.02 秒（10,000個序列）
- **CNN 訓練**：197.5±0.8 秒
- **RF 訓練**：28.1±0.3 秒
- **預測速度**：兩種模型約1.7秒（10,000個序列）

## 📋 系統需求

### R 套件依賴
```r
# 核心套件
library(seqinr)      # 序列分析
library(keras)       # 深度學習
library(kerasR)      # Keras 介面
library(caret)       # 機器學習
library(randomForest) # 隨機森林
library(ROCR)        # 效能評估
library(protr)       # 蛋白質分析
```

### 系統需求
- **R**：版本 ≥ 3.6.0
- **Python**：版本 ≥ 3.6（用於 TensorFlow 後端）
- **記憶體**：建議 ≥ 4GB RAM
- **儲存空間**：模型和資料集需要 ≥ 1GB 可用空間

## 🚀 安裝指南

1. **複製專案庫**：
```bash
git clone https://github.com/your-username/Deep-AmPEP30.git
cd Deep-AmPEP30
```

2. **安裝 R 套件依賴**：
```r
# 安裝必需套件
install.packages(c("seqinr", "caret", "randomForest", "ROCR", "protr"))

# 安裝 Keras 和 TensorFlow
install.packages("keras")
library(keras)
install_keras()
```

3. **初始化 Keras**：
```r
library(kerasR)
keras_init()
```

## 📖 使用方法

### 基本預測

#### 使用 CNN 模型
```r
source("Deep-AmPEP30.R")

# 使用 CNN 預測抗菌肽
results <- cnn_yan_predict("your_sequences.fasta")
```

#### 使用隨機森林模型
```r
# 使用隨機森林預測
results <- rf_yan_predict("your_sequences.fasta")
```

### 訓練自訂模型

#### 訓練 CNN 模型
```r
# 使用預設參數訓練 CNN 模型
develop_cnn_mdl_yan(rf = FALSE)
```

#### 訓練隨機森林模型
```r
# 訓練隨機森林模型
develop_cnn_mdl_yan(rf = TRUE)
```

### 特徵生成
```r
# 為自訂序列生成特徵
features <- totally_psekraac_best5(
    test_path = "your_sequences.fasta",
    test = TRUE,
    without_fasta_name = FALSE,
    check_len = TRUE
)
```

## 📁 資料集結構

```
dataset/
├── model_po.fasta          # 正樣本訓練資料（3,303個AMPs）
├── model_ne.fasta          # 負樣本訓練資料（3,293個非AMPs）
├── train_po.fasta          # 額外正樣本訓練資料
├── train_ne.fasta          # 額外負樣本訓練資料
├── test_po.fasta           # 正樣本測試資料
└── test_ne.fasta           # 負樣本測試資料
```

### 資料格式
序列應以 FASTA 格式提供：
```
>sequence_id_1
ACSAG
>sequence_id_2
FRWWHR
>sequence_id_3
RKKWFW
```

## 🔬 模型驗證

系統使用 10 重交叉驗證進行模型評估：

```r
# 交叉驗證結果在訓練期間自動產生
# 效能指標包括：
# - 準確度
# - 敏感度（真陽性率）
# - 特異度（真陰性率）
# - AUC（曲線下面積）
```

## 📈 效能基準測試

所有效能基準測試都是基於長度為 5-30 個氨基酸的序列進行：

| 操作 | 時間（秒） | 序列數量 | 模型 |
|------|------------|----------|------|
| 特徵生成 | 0.52±0.02 | 10,000 | - |
| CNN 訓練 | 197.5±0.8 | 3,297 | CNN |
| RF 訓練 | 28.1±0.3 | 3,297 | RF |
| CNN 預測 | 1.73±0.06 | 10,000 | CNN |
| RF 預測 | 1.70±0.04 | 10,000 | RF |

## 🔧 配置設定

### 氨基酸分類系統
系統使用針對抗菌肽預測最佳化的多種分類方案：

```r
# 可用的分類類型
type8.cluster17 <- c('AT','C','DE','F','G','H','IV','K','L','M','N','P','Q','R','S','V','W')
type3a.cluster19 <- c('FA','P','G','S','T','D','E','Q','N','K','R','H','W','Y','M','L','I','V','C')
type12.cluster18 <- c('TVLI','M','F','W','Y','C','A','H','G','N','Q','P','R','K','S','T','D','E')
type7.cluster15 <- c('C','K','R','W','Y','A','FILV','M','D','E','Q','H','TP','GS','N')
type12.cluster17 <- c('TVLI','M','F','W','Y','C','A','H','G','N','Q','P','R','K','S','T','DE')
```

## 🧪 應用領域

此系統在各種研究和工業應用中具有重要價值：

- **藥物發現**：篩選新型抗菌肽候選分子
- **蛋白質組學**：大規模蛋白質功能註釋
- **抗生素替代**：識別新型抗菌療法
- **生物技術**：設計具有特定抗菌活性的肽段

## 📚 檔案結構

```
Deep-AmPEP30/
├── Deep-AmPEP30.R                 # 主要預測腳本
├── dataset/                       # 訓練和測試資料集
├── note                          # 開發筆記和觀察記錄
├── run_time_performance.txt      # 詳細效能基準測試
├── rf_cnn_mdl_use_time.note     # 模型時間比較
└── C_glabarta_protein/          # 應用範例資料
```

## 🤝 貢獻指南

我們歡迎對 Deep-AmPEP30 的改進貢獻：

1. Fork 此專案庫
2. 建立功能分支（`git checkout -b feature/improvement`）
3. 提交您的更改（`git commit -am 'Add new feature'`）
4. 推送到分支（`git push origin feature/improvement`）
5. 建立 Pull Request

## 📄 授權條款

本專案採用 MIT 授權條款 - 詳見 [LICENSE](LICENSE) 檔案。

## 🙏 致謝

- **蛋白質序列分析**：使用 `seqinr` 和 `protr` 套件建構
- **機器學習**：由 `keras`、`caret` 和 `randomForest` 提供技術支援
- **效能最佳化**：經過廣泛的基準測試和生產使用最佳化

---

*Deep-AmPEP30：透過深度學習和全面特徵工程推進抗菌肽預測技術。*