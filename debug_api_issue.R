# 調試API問題 - 測試單獨的特徵提取函數

# 載入必要套件
library(randomForest)
library(seqinr)
library(protr)

# 從我們的API複製常量和函數進行測試
type8.cluster17 <- c("AT", "C", "DE", "F", "G", "H", "IV", "K", "L", "M", "N", "P", "Q", "R", "S", "V", "W")
type3a.cluster19 <- c("FA", "P", "G", "S", "T", "D", "E", "Q", "N", "K", "R", "H", "W", "Y", "M", "L", "I", "V", "C")
type12.cluster18 <- c("TVLI", "M", "F", "W", "Y", "C", "A", "H", "G", "N", "Q", "P", "R", "K", "S", "T", "D", "E")
type7.cluster15 <- c("C", "K", "R", "W", "Y", "A", "FILV", "M", "D", "E", "Q", "H", "TP", "GS", "N")
type12.cluster17 <- c("TVLI", "M", "F", "W", "Y", "C", "A", "H", "G", "N", "Q", "P", "R", "K", "S", "T", "DE")

# 測試gene_aac函數
gene_aac <- function(seq_name, label = NULL) {
  ft <- NULL
  seqs <- seq_name$seq
  name <- seq_name$name
  n <- length(seqs)

  for (seq in seqs) {
    t <- extractAAC(seq)
    t <- t * nchar(seq)
    ft <- rbind(ft, t)
  }

  rownames(ft) <- name
  if (!is.null(label)) {
    class <- rep(label, n)
    ft <- cbind(ft, class)
  }
  ft <- data.frame(ft)
  return(ft)
}

# 測試totally_psekraac_best5函數
totally_psekraac_best5 <- function(sequences, seq_names) {
  cat("輸入序列:", sequences, "\n")
  cat("輸入名稱:", seq_names, "\n")

  # 為單一序列創建proper格式
  seqs <- list(seq = sequences, name = seq_names)
  cat("創建的seqs結構:\n")
  print(seqs)

  data.aac <- gene_aac(seqs)
  cat("AAC數據維度:", dim(data.aac), "\n")
  cat("AAC數據列名:", colnames(data.aac), "\n")

  data.temp <- NULL
  descriptors <- c(type8.cluster17, type3a.cluster19, type12.cluster18, type7.cluster15, type12.cluster17)
  cat("描述符數量:", length(descriptors), "\n")

  for (i in seq_along(descriptors)) {
    x <- descriptors[i]
    cat("處理描述符", i, ":", x, "\n")

    if (nchar(x) > 1) {
      char_multi <- substring(x, seq(1, nchar(x), 1), seq(1, nchar(x), 1))
      cat("  多字符分解:", char_multi, "\n")

      # 檢查列是否存在
      missing_cols <- char_multi[!char_multi %in% colnames(data.aac)]
      if (length(missing_cols) > 0) {
        cat("  警告: 缺少列:", missing_cols, "\n")
        next
      }

      t <- apply(data.aac[, char_multi], 1, sum)
      data.temp <- cbind(data.temp, t)
    } else {
      if (!x %in% colnames(data.aac)) {
        cat("  警告: 缺少列:", x, "\n")
        next
      }
      t <- data.aac[, x]
      data.temp <- cbind(data.temp, t)
    }
  }

  cat("最終數據維度:", dim(data.temp), "\n")

  n <- length(descriptors)
  name <- paste0(descriptors, "_", 1:n)
  data <- data.temp
  whole_name <- name

  colnames(data) <- whole_name
  data <- data.frame(data)
  return(data)
}

# 測試案例
cat("=== 調試特徵提取過程 ===\n")

test_sequence <- "GIGKFLHSAKKFGKAFVGEIMNS"
test_name <- "Magainin-2"

cat("測試序列:", test_sequence, "\n")
cat("序列長度:", nchar(test_sequence), "\n")

# 測試特徵提取
tryCatch(
  {
    feature_vector <- totally_psekraac_best5(c(test_sequence), c(test_name))
    cat("✓ 特徵提取成功\n")
    cat("特徵向量維度:", dim(feature_vector), "\n")
    cat("前10個特徵值:", head(feature_vector, 1)[1:min(10, ncol(feature_vector))], "\n")
  },
  error = function(e) {
    cat("✗ 特徵提取失敗:", e$message, "\n")
  }
)

cat("\n=== 比較與原始RF-AmPEP30.R ===\n")

# 載入和測試原始函數
source("RF-AmPEP30.R")

# 創建測試文件
writeLines(paste0(">", test_name, "\n", test_sequence), "temp_test.fasta")

# 使用原始函數
tryCatch(
  {
    original_data <- totally_psekraac_best5(test_path = "temp_test.fasta", test = TRUE)
    cat("✓ 原始函數成功\n")
    cat("原始數據維度:", dim(original_data), "\n")
    cat("前10個特徵值:", head(original_data, 1)[1:min(10, ncol(original_data))], "\n")

    # 載入模型並測試預測
    model <- readRDS("AmPEP30-RF-1200tree.mdl")
    prediction <- predict(model, original_data, type = "prob")
    cat("✓ 原始預測成功\n")
    cat("預測結果:", prediction, "\n")
  },
  error = function(e) {
    cat("✗ 原始函數失敗:", e$message, "\n")
  }
)

# 清理
file.remove("temp_test.fasta")
