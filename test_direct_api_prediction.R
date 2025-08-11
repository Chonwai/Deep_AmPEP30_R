# 直接測試API風格的預測，使用已驗證的函數

# 載入本地工作的預測函數
source("AmPEP30_Demo.R")

# 測試案例
test_sequence <- "GIGKFLHSAKKFGKAFVGEIMNS"
test_name <- "Magainin-2"

cat("=== 使用本地函數預測 ===\n")
result <- predict_peptide(test_sequence, test_name)
print(result)

cat("\n=== 轉換為API格式 ===\n")
api_result <- list(
  sequence_name = result$sequence_name,
  sequence = result$sequence,
  length = result$length,
  prediction = result$prediction,
  amp_probability = result$amp_probability,
  non_amp_probability = result$non_amp_probability,
  confidence = result$confidence,
  interpretation = result$interpretation,
  status = "success"
)

print(api_result)

cat("\n=== JSON格式 ===\n")
library(jsonlite)
json_result <- toJSON(api_result, pretty = TRUE)
cat(json_result)
