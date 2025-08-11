# AmPEP30 驗證API測試客戶端
# 測試API與本地結果的一致性

library(httr)
library(jsonlite)

# API設定
API_BASE_URL <- "http://localhost:8002"

# 測試函數
test_api_endpoint <- function(endpoint, method = "GET", body = NULL, description = "") {
  cat("\n=== 測試:", description, "===\n")
  cat("端點:", endpoint, "\n")

  url <- paste0(API_BASE_URL, endpoint)

  tryCatch(
    {
      if (method == "GET") {
        response <- GET(url)
      } else if (method == "POST") {
        response <- POST(url, body = body, encode = "json")
      }

      if (status_code(response) == 200) {
        result <- fromJSON(content(response, "text"))
        cat("✓ 狀態: 成功\n")
        return(result)
      } else {
        cat("✗ 狀態碼:", status_code(response), "\n")
        cat("錯誤:", content(response, "text"), "\n")
        return(NULL)
      }
    },
    error = function(e) {
      cat("✗ 連接錯誤:", e$message, "\n")
      return(NULL)
    }
  )
}

# 比較函數 - 比較API結果與本地結果
compare_predictions <- function(api_result, local_result, sequence_name) {
  cat("\n--- 比較結果:", sequence_name, "---\n")

  if (is.null(api_result) || api_result$status == "error") {
    cat("✗ API預測失敗\n")
    return(FALSE)
  }

  # 比較主要字段
  api_pred <- api_result$prediction
  api_prob <- api_result$amp_probability
  local_pred <- local_result$prediction
  local_prob <- local_result$amp_probability

  cat("API預測:   ", api_pred, "(", round(api_prob * 100, 1), "%)\n")
  cat("本地預測:  ", local_pred, "(", round(local_prob * 100, 1), "%)\n")

  # 檢查一致性 (允許小數點後3位的差異)
  pred_match <- api_pred == local_pred
  prob_match <- abs(api_prob - local_prob) < 0.001

  if (pred_match && prob_match) {
    cat("✓ 結果一致\n")
    return(TRUE)
  } else {
    cat("✗ 結果不一致\n")
    cat("  預測一致:", pred_match, "\n")
    cat("  機率一致:", prob_match, "(差異:", abs(api_prob - local_prob), ")\n")
    return(FALSE)
  }
}

# 主要測試流程
main_test <- function() {
  cat("========================================\n")
  cat("  AmPEP30 API與本地結果一致性測試\n")
  cat("========================================\n")

  # 1. 健康檢查
  health_result <- test_api_endpoint("/health", "GET", description = "健康檢查")
  if (is.null(health_result)) {
    cat("API服務未響應，請確保服務正在運行\n")
    return()
  }

  # 2. 模型資訊
  model_info <- test_api_endpoint("/model/info", "GET", description = "模型資訊")

  # 3. 測試已知序列
  cat("\n\n=== 與本地結果比較 ===\n")

  # 載入本地預測函數
  source("AmPEP30_Demo.R")

  # 測試序列
  test_cases <- list(
    list(name = "Magainin-2", sequence = "GIGKFLHSAKKFGKAFVGEIMNS"),
    list(name = "Melittin", sequence = "GIGAVLKVLTTGLPALISWIKRKRQQ"),
    list(name = "Random_peptide", sequence = "ACDEFGHIKL"),
    list(name = "Short_peptide", sequence = "KLMNP")
  )

  all_consistent <- TRUE

  for (test_case in test_cases) {
    # API預測
    api_result <- test_api_endpoint(
      "/predict/single",
      "POST",
      list(sequence = test_case$sequence, name = test_case$name),
      paste("API預測", test_case$name)
    )

    # 本地預測
    local_result <- predict_peptide(test_case$sequence, test_case$name)

    # 比較結果
    is_consistent <- compare_predictions(api_result, local_result, test_case$name)
    all_consistent <- all_consistent && is_consistent
  }

  # 4. 測試FASTA格式
  cat("\n\n=== 測試FASTA格式 ===\n")
  fasta_content <- ">Magainin-2
GIGKFLHSAKKFGKAFVGEIMNS
>Melittin
GIGAVLKVLTTGLPALISWIKRKRQQ"

  fasta_result <- test_api_endpoint(
    "/predict/fasta",
    "POST",
    list(fasta_content = fasta_content),
    "FASTA多序列預測"
  )

  # 5. 測試演示端點
  cat("\n\n=== 測試演示端點 ===\n")
  demo_result <- test_api_endpoint("/test/demo", "GET", description = "演示測試")

  # 總結
  cat("\n\n========================================\n")
  cat("  測試總結\n")
  cat("========================================\n")

  if (all_consistent) {
    cat("✅ 所有測試通過！API結果與本地結果完全一致\n")
    cat("🎉 系統已準備好進行生產部署\n")
  } else {
    cat("❌ 部分測試失敗，需要檢查不一致的原因\n")
  }

  cat("\nAPI服務狀態:", ifelse(is.null(health_result), "離線", "在線"), "\n")
  cat("模型載入狀態:", ifelse(is.null(model_info), "未知", ifelse(model_info$model_loaded, "已載入", "未載入")), "\n")
}

# 執行測試
if (!interactive()) {
  main_test()
} else {
  cat("API測試客戶端已載入\n")
  cat("執行 main_test() 來開始測試\n")
}
