#!/usr/bin/env Rscript

# API 客戶端測試
cat("=== API 客戶端測試 ===\n\n")

library(httr)
library(jsonlite)

# API 基本 URL
base_url <- "http://localhost:8006"

# 測試健康檢查
cat("1. 測試健康檢查...\n")
tryCatch(
  {
    response <- GET(paste0(base_url, "/health"))
    if (status_code(response) == 200) {
      result <- content(response, "text", encoding = "UTF-8")
      cat("✓ 健康檢查成功\n")
      cat("響應:", result, "\n\n")
    } else {
      cat("✗ 健康檢查失敗，狀態碼:", status_code(response), "\n\n")
    }
  },
  error = function(e) {
    cat("✗ 健康檢查錯誤:", e$message, "\n\n")
  }
)

# 測試預測功能
cat("2. 測試預測功能...\n")

test_sequences <- list(
  list(name = "Magainin", seq = "GIGKFLHSAKKFGKAFVGEIMNS", expected = "AMP"),
  list(name = "Short_seq", seq = "KLMNP", expected = "Unknown"),
  list(name = "Random", seq = "ACDEFGHIKL", expected = "Non-AMP")
)

for (i in 1:length(test_sequences)) {
  test_seq <- test_sequences[[i]]

  cat(sprintf("測試 %s (%s)...\n", test_seq$name, test_seq$seq))

  tryCatch(
    {
      # 使用 form data 方式
      response <- POST(
        url = paste0(base_url, "/predict"),
        body = list(sequence = test_seq$seq),
        encode = "form"
      )

      if (status_code(response) == 200) {
        result_text <- content(response, "text", encoding = "UTF-8")
        cat("✓ 請求成功\n")
        cat("響應:", result_text, "\n")

        # 嘗試解析 JSON
        tryCatch(
          {
            result_json <- fromJSON(result_text)
            if (!is.null(result_json$success) && result_json$success) {
              cat(sprintf(
                "  預測: %s (機率: %.3f)\n",
                result_json$prediction, result_json$probability
              ))
            } else {
              cat("  預測失敗:", result_json$error, "\n")
            }
          },
          error = function(e) {
            cat("  JSON 解析錯誤:", e$message, "\n")
          }
        )
      } else {
        cat("✗ 請求失敗，狀態碼:", status_code(response), "\n")
      }
    },
    error = function(e) {
      cat("✗ 請求錯誤:", e$message, "\n")
    }
  )

  cat("\n")
}

cat("測試完成！\n")
